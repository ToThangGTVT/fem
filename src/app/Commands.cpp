#include "Commands.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <map>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <functional>

#include "App.h"
#include "core/Mesher.h"

using namespace fem;

namespace {

constexpr double kPi = 3.14159265358979323846;

std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

bool number(const std::string& s, double& v) {
    const std::string t = trim(s);
    if (t.empty()) return false;
    char* end = nullptr;
    v = std::strtod(t.c_str(), &end);
    return end && *end == 0;
}

// "a,b,c" → các số; trả về số phần tử đọc được (0 nếu lỗi)
int numbers(const std::string& s, double* out, int maxN) {
    int n = 0;
    size_t i = 0;
    while (i <= s.size() && n < maxN) {
        size_t j = s.find(',', i);
        if (j == std::string::npos) j = s.size();
        if (!number(s.substr(i, j - i), out[n])) return 0;
        ++n;
        i = j + 1;
        if (j == s.size()) break;
    }
    return n;
}

std::string fmt(const char* f, ...) __attribute__((format(printf, 1, 2)));
std::string fmt(const char* f, ...) {
    char b[512];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(b, sizeof b, f, ap);
    va_end(ap);
    return b;
}

std::string num(double v) { return fmt("%g", v); }

double tolOf(const Model& m) {
    double x0, y0, x1, y1;
    m.bounds(x0, y0, x1, y1);
    return 1e-6 * std::max({1.0, x1 - x0, y1 - y0});
}

// ---------------------------------------------------------------- chọn đối tượng rồi thực hiện
struct SelectCmd : Command {
    bool selecting = true;
    int filter = 0; // 0: bất kỳ, 1: nút, 2: thanh

    bool usable(App& a) const {
        if (filter == 1) return !a.selNodes.empty();
        if (filter == 2) {
            for (int e : a.selElems)
                if (a.model.elements[e].isBar()) return true;
            return false;
        }
        return !a.selNodes.empty() || !a.selElems.empty();
    }
    std::vector<int> bars(App& a) const {
        std::vector<int> v;
        for (int e : a.selElems)
            if (a.model.elements[e].isBar()) v.push_back(e);
        return v;
    }
    void start(App& a) override {
        if (usable(a)) {
            selecting = false;
            onSelected(a);
        }
    }
    bool wantsSelection() const override { return selecting; }
    std::string prompt(App& a) override {
        if (selecting) {
            const char* what = filter == 1 ? "nút" : filter == 2 ? "thanh" : "đối tượng";
            return fmt("Chọn %s (đã chọn %zu), Enter để kết thúc:", what, a.selNodes.size() + a.selElems.size());
        }
        return afterPrompt(a);
    }
    void onEnter(App& a) override {
        if (selecting) {
            if (!usable(a)) {
                a.echo("Không có đối tượng phù hợp.");
                done = true;
                return;
            }
            selecting = false;
            onSelected(a);
        } else {
            afterEnter(a);
        }
    }
    virtual void onSelected(App& a) = 0;
    virtual std::string afterPrompt(App&) { return ""; }
    virtual void afterEnter(App&) { done = true; }
};

// ---------------------------------------------------------------- vẽ thanh
struct LineCmd : Command {
    bool truss;
    std::vector<SnapPoint> chain;
    explicit LineCmd(bool t) : truss(t) { name = t ? "THANH DÀN" : "THANH"; }
    bool wantsPoint() const override { return true; }
    bool base(double& x, double& y) const override {
        if (chain.empty()) return false;
        x = chain.back().x;
        y = chain.back().y;
        return true;
    }
    std::string prompt(App&) override {
        if (chain.empty()) return "Chỉ định điểm đầu:";
        if (chain.size() >= 3) return "Chỉ định điểm tiếp theo hoặc [Khép kín(C)/Hoàn tác(U)]:";
        return "Chỉ định điểm tiếp theo hoặc [Hoàn tác(U)]:";
    }
    void segment(App& a, const SnapPoint& p0, const SnapPoint& p1) {
        a.snapshot();
        const int n0 = a.nodeForSnap(p0);
        const int n1 = a.nodeForSnap(p1);
        if (n0 >= 0 && n1 >= 0 && n0 != n1) {
            if (a.model.materials.empty()) a.model.materials.push_back(Material{});
            if (a.model.sections.empty()) a.model.sections.push_back(Section{});
            Element el;
            el.type = truss ? ElemType::Truss : ElemType::Frame;
            el.n = {n0, n1, -1, -1};
            el.mat = a.curMat;
            el.sec = a.curSec;
            if (a.model.findElement(el.type, el.n) < 0) a.model.elements.push_back(el);
        }
        a.changed();
    }
    void onPoint(App& a, const SnapPoint& p) override {
        if (!chain.empty()) {
            if (std::hypot(p.x - chain.back().x, p.y - chain.back().y) < tolOf(a.model)) return;
            segment(a, chain.back(), p);
        }
        SnapPoint q = p;
        q.kind = SnapKind::None; // vị trí đã thành nút
        q.elem = -1;
        q.node = -1;
        chain.push_back(q);
    }
    bool onText(App& a, const std::string& s) override {
        const std::string u = upper(s);
        if (u == "U") {
            if (chain.size() > 1) {
                a.undo();
                chain.pop_back();
            } else chain.clear();
            return true;
        }
        if (u == "C" && chain.size() >= 3) {
            segment(a, chain.back(), chain.front());
            done = true;
            return true;
        }
        return false;
    }
    void preview(App& a, ImDrawList* dl, const SnapPoint& c) override {
        if (chain.empty()) return;
        dl->AddLine(a.cam.toScreen(chain.back().x, chain.back().y), a.cam.toScreen(c.x, c.y), a.pal.hover, 1.5f);
    }
};

struct NodeCmd : Command {
    NodeCmd() { name = "NÚT"; }
    bool wantsPoint() const override { return true; }
    std::string prompt(App&) override { return "Chỉ định vị trí nút (Enter để kết thúc):"; }
    void onPoint(App& a, const SnapPoint& p) override {
        a.snapshot();
        a.nodeForSnap(p);
        a.changed();
    }
};

// ---------------------------------------------------------------- xóa
struct EraseCmd : SelectCmd {
    EraseCmd() { name = "XÓA"; }
    void onSelected(App& a) override {
        const size_t n = a.selNodes.size() + a.selElems.size();
        a.deleteSelection();
        a.echo(fmt("Đã xóa %zu đối tượng.", n));
        done = true;
    }
};

// ---------------------------------------------------------------- biến đổi hình học
using Xform = std::function<void(double&, double&)>;

void drawGhost(App& a, ImDrawList* dl, const Xform& f) {
    const Model& m = a.model;
    auto P = [&](int n) {
        double x = m.nodes[n].x, y = m.nodes[n].y;
        f(x, y);
        return a.cam.toScreen(x, y);
    };
    for (int e : a.selElems) {
        const Element& el = m.elements[e];
        const int nc = el.nodeCount();
        for (int k = 0; k < (el.isBar() ? 1 : nc); ++k)
            dl->AddLine(P(el.n[k]), P(el.n[(k + 1) % nc]), a.pal.hover, 1.5f);
    }
    for (int n : a.selNodes) dl->AddCircle(P(n), 4, a.pal.hover, 0, 1.5f);
}

// copy = false: dời các nút; copy = true: sao chép phần tử (kèm gối và tải) sang vị trí mới.
// reflect: phép đối xứng (đổi chiều mômen, phản xạ véc tơ lực nút).
void applyXform(App& a, const Xform& f, bool copy, bool reflect = false, double ux = 0, double uy = 0) {
    Model& m = a.model;
    const double tol = tolOf(m);
    const std::set<int> nodes = a.selectionNodes();
    if (!copy) {
        for (int n : nodes) f(m.nodes[n].x, m.nodes[n].y);
        mergeNodes(m, tol);
        return;
    }
    std::map<int, int> nmap;
    std::set<int> fresh;
    for (int n : nodes) {
        double x = m.nodes[n].x, y = m.nodes[n].y;
        f(x, y);
        const size_t before = m.nodes.size();
        const int nn = m.addNode(x, y, tol);
        if (m.nodes.size() > before) fresh.insert(nn);
        nmap[n] = nn;
    }
    std::map<int, int> emap;
    const std::vector<int> elems(a.selElems.begin(), a.selElems.end());
    for (int e : elems) {
        Element c = m.elements[e];
        for (int k = 0; k < c.nodeCount(); ++k) c.n[k] = nmap[c.n[k]];
        if (m.findElement(c.type, c.n) >= 0) continue;
        m.elements.push_back(c);
        emap[e] = (int)m.elements.size() - 1;
    }
    for (const auto& [o, nn] : nmap) {
        const Support* s = m.supportAt(o);
        if (s && fresh.count(nn)) {
            Support c = *s;
            c.node = nn;
            m.supports.push_back(c);
        }
    }
    // phản xạ véc tơ qua trục có hướng (ux, uy): F' = 2(F·u)u − F
    auto refl = [&](double& fx, double& fy) {
        if (!reflect) return;
        const double d = fx * ux + fy * uy;
        fx = 2 * d * ux - fx;
        fy = 2 * d * uy - fy;
    };
    for (auto& lc : m.cases) {
        const size_t nN = lc.nodal.size(), nE = lc.elem.size();
        for (size_t i = 0; i < nN; ++i) {
            auto it = nmap.find(lc.nodal[i].node);
            if (it == nmap.end() || !fresh.count(it->second)) continue;
            NodalLoad l = lc.nodal[i];
            l.node = it->second;
            refl(l.F[0], l.F[1]);
            if (reflect) l.F[2] = -l.F[2];
            lc.nodal.push_back(l);
        }
        for (size_t i = 0; i < nE; ++i) {
            auto it = emap.find(lc.elem[i].elem);
            if (it == emap.end()) continue;
            ElemLoad l = lc.elem[i];
            l.elem = it->second;
            if (reflect && l.kind == ElemLoadKind::Moment) l.q1 = -l.q1;
            lc.elem.push_back(l);
        }
    }
}

struct MoveCopyCmd : SelectCmd {
    bool copy;
    int step = 0; // 0: điểm gốc, 1: điểm thứ hai
    double bx = 0, by = 0;
    explicit MoveCopyCmd(bool c) : copy(c) { name = c ? "SAO CHÉP" : "DI CHUYỂN"; }
    void onSelected(App&) override { step = 0; }
    bool wantsPoint() const override { return !selecting; }
    bool base(double& x, double& y) const override {
        if (selecting || step == 0) return false;
        x = bx;
        y = by;
        return true;
    }
    std::string afterPrompt(App&) override {
        if (step == 0) return "Chỉ định điểm gốc:";
        return copy ? "Chỉ định điểm thứ hai (Enter để kết thúc):" : "Chỉ định điểm thứ hai:";
    }
    void onPoint(App& a, const SnapPoint& p) override {
        if (step == 0) {
            bx = p.x;
            by = p.y;
            step = 1;
            return;
        }
        const double dx = p.x - bx, dy = p.y - by;
        a.snapshot();
        applyXform(a, [&](double& x, double& y) { x += dx; y += dy; }, copy);
        a.changed();
        if (!copy) {
            a.selNodes.clear();
            a.selElems.clear();
            done = true;
        }
    }
    void preview(App& a, ImDrawList* dl, const SnapPoint& c) override {
        if (selecting || step == 0) return;
        const double dx = c.x - bx, dy = c.y - by;
        dl->AddLine(a.cam.toScreen(bx, by), a.cam.toScreen(c.x, c.y), a.pal.textDim, 1.0f);
        drawGhost(a, dl, [&](double& x, double& y) { x += dx; y += dy; });
    }
    void afterEnter(App& a) override {
        if (copy) {
            a.selNodes.clear();
            a.selElems.clear();
        }
        done = true;
    }
};

struct RotateCmd : SelectCmd {
    int step = 0;
    double bx = 0, by = 0;
    RotateCmd() { name = "XOAY"; }
    void onSelected(App&) override { step = 0; }
    bool wantsPoint() const override { return !selecting; }
    bool base(double& x, double& y) const override {
        if (selecting || step == 0) return false;
        x = bx;
        y = by;
        return true;
    }
    std::string afterPrompt(App&) override {
        return step == 0 ? "Chỉ định tâm xoay:" : "Góc xoay (độ, dương ngược chiều kim đồng hồ) hoặc chỉ điểm:";
    }
    Xform rot(double deg) const {
        const double c = std::cos(deg * kPi / 180), s = std::sin(deg * kPi / 180);
        const double X = bx, Y = by;
        return [=](double& x, double& y) {
            const double dx = x - X, dy = y - Y;
            x = X + c * dx - s * dy;
            y = Y + s * dx + c * dy;
        };
    }
    void apply(App& a, double deg) {
        a.snapshot();
        applyXform(a, rot(deg), false);
        a.changed();
        a.selNodes.clear();
        a.selElems.clear();
        done = true;
    }
    void onPoint(App& a, const SnapPoint& p) override {
        if (step == 0) {
            bx = p.x;
            by = p.y;
            step = 1;
            return;
        }
        apply(a, std::atan2(p.y - by, p.x - bx) * 180 / kPi);
    }
    bool onText(App& a, const std::string& s) override {
        double d;
        if (selecting || step == 0 || !number(s, d)) return false;
        apply(a, d);
        return true;
    }
    void preview(App& a, ImDrawList* dl, const SnapPoint& c) override {
        if (selecting || step == 0) return;
        dl->AddLine(a.cam.toScreen(bx, by), a.cam.toScreen(c.x, c.y), a.pal.textDim, 1.0f);
        drawGhost(a, dl, rot(std::atan2(c.y - by, c.x - bx) * 180 / kPi));
    }
};

struct MirrorCmd : SelectCmd {
    int step = 0;
    double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    MirrorCmd() { name = "ĐỐI XỨNG"; }
    void onSelected(App&) override { step = 0; }
    bool wantsPoint() const override { return !selecting && step < 2; }
    bool base(double& x, double& y) const override {
        if (selecting || step != 1) return false;
        x = x1;
        y = y1;
        return true;
    }
    std::string afterPrompt(App&) override {
        if (step == 0) return "Điểm thứ nhất của trục đối xứng:";
        if (step == 1) return "Điểm thứ hai của trục đối xứng:";
        return "Xóa đối tượng gốc? [Có(C)/Không(K)] <K>:";
    }
    static Xform mir(double x1, double y1, double x2, double y2) {
        const double L = std::hypot(x2 - x1, y2 - y1), ux = (x2 - x1) / L, uy = (y2 - y1) / L;
        return [=](double& x, double& y) {
            const double dx = x - x1, dy = y - y1, d = dx * ux + dy * uy;
            x = x1 + 2 * d * ux - dx;
            y = y1 + 2 * d * uy - dy;
        };
    }
    void onPoint(App&, const SnapPoint& p) override {
        if (step == 0) {
            x1 = p.x;
            y1 = p.y;
            step = 1;
        } else if (std::hypot(p.x - x1, p.y - y1) > 0) {
            x2 = p.x;
            y2 = p.y;
            step = 2;
        }
    }
    void finish(App& a, bool erase) {
        const double L = std::hypot(x2 - x1, y2 - y1);
        a.snapshot();
        applyXform(a, mir(x1, y1, x2, y2), !erase, true, (x2 - x1) / L, (y2 - y1) / L);
        a.changed();
        a.selNodes.clear();
        a.selElems.clear();
        done = true;
    }
    bool onText(App& a, const std::string& s) override {
        if (selecting || step != 2) return false;
        const std::string u = upper(s);
        if (u == "C" || u == "Y") finish(a, true);
        else if (u == "K" || u == "N") finish(a, false);
        else return false;
        return true;
    }
    void afterEnter(App& a) override {
        if (step == 2) finish(a, false);
        else done = true;
    }
    void preview(App& a, ImDrawList* dl, const SnapPoint& c) override {
        if (selecting || step != 1 || std::hypot(c.x - x1, c.y - y1) <= 0) return;
        dl->AddLine(a.cam.toScreen(x1, y1), a.cam.toScreen(c.x, c.y), a.pal.textDim, 1.0f);
        drawGhost(a, dl, mir(x1, y1, c.x, c.y));
    }
};

// ---------------------------------------------------------------- chia thanh
struct DivideCmd : SelectCmd {
    DivideCmd() {
        name = "CHIA THANH";
        filter = 2;
    }
    bool asking = false;
    void onSelected(App&) override { asking = true; }
    std::string afterPrompt(App& a) override { return fmt("Số đoạn <%d>:", a.defDiv); }
    void apply(App& a, int n) {
        a.defDiv = n;
        std::vector<int> b = bars(a);
        std::sort(b.rbegin(), b.rend());
        a.snapshot();
        for (int e : b) subdivideBar(a.model, e, n);
        a.changed();
        a.selElems.clear();
        a.echo(fmt("Đã chia %zu thanh thành %d đoạn.", b.size(), n));
        done = true;
    }
    bool onText(App& a, const std::string& s) override {
        double d;
        if (selecting || !number(s, d) || d < 2 || d > 500) return false;
        apply(a, (int)d);
        return true;
    }
    void afterEnter(App& a) override { apply(a, a.defDiv); }
};

// ---------------------------------------------------------------- gối tựa
struct SupportCmd : SelectCmd {
    std::array<bool, 3> fix;
    bool clear;
    SupportCmd(const char* n, std::array<bool, 3> f, bool c = false) : fix(f), clear(c) {
        name = n;
        filter = 1;
    }
    void onSelected(App& a) override {
        a.snapshot();
        for (int n : a.selNodes) {
            if (clear) {
                if (Support* s = a.model.supportAt(n)) *s = Support{n};
            } else {
                a.model.supportFor(n).fix = fix;
            }
        }
        a.model.cleanSupports();
        a.changed();
        a.echo(fmt("%s: %zu nút.", name.c_str(), a.selNodes.size()));
        a.selNodes.clear();
        done = true;
    }
};

// ---------------------------------------------------------------- tải trọng
struct NodalLoadCmd : SelectCmd {
    NodalLoadCmd() {
        name = "TẢI NÚT";
        filter = 1;
    }
    void onSelected(App&) override {}
    std::string afterPrompt(App& a) override {
        return fmt("Fx,Fy,Mz (kN, kNm) cho TH \"%s\" <%s,%s,%s>:", a.model.cases[a.activeCase].name.c_str(),
                   num(a.defNodal[0]).c_str(), num(a.defNodal[1]).c_str(), num(a.defNodal[2]).c_str());
    }
    void apply(App& a) {
        a.snapshot();
        for (int n : a.selNodes) a.model.cases[a.activeCase].nodal.push_back({n, {a.defNodal[0], a.defNodal[1], a.defNodal[2]}});
        a.changed();
        a.selNodes.clear();
        done = true;
    }
    bool onText(App& a, const std::string& s) override {
        if (selecting) return false;
        double v[3] = {0, 0, 0};
        const int n = numbers(s, v, 3);
        if (n == 0) return false;
        for (int i = 0; i < 3; ++i) a.defNodal[i] = i < n ? v[i] : 0.0;
        if (n == 1) { // một số: hiểu là Fy
            a.defNodal[1] = v[0];
            a.defNodal[0] = a.defNodal[2] = 0;
        }
        apply(a);
        return true;
    }
    void afterEnter(App& a) override { apply(a); }
};

bool parseDir(const std::string& s, int& dir, bool& projected) {
    const std::string u = upper(s);
    projected = false;
    if (u == "GX" || u == "X") dir = (int)LoadDir::GlobalX;
    else if (u == "GY" || u == "Y") dir = (int)LoadDir::GlobalY;
    else if (u == "LX") dir = (int)LoadDir::LocalX;
    else if (u == "LY") dir = (int)LoadDir::LocalY;
    else if (u == "HC") {
        dir = (int)LoadDir::GlobalY;
        projected = true;
    } else return false;
    return true;
}

const char* dirToken(int d, bool projected) {
    if (projected) return "HC";
    switch ((LoadDir)d) {
    case LoadDir::GlobalX: return "GX";
    case LoadDir::GlobalY: return "GY";
    case LoadDir::LocalX: return "LX";
    case LoadDir::LocalY: return "LY";
    }
    return "GY";
}

struct DistLoadCmd : SelectCmd {
    int step = 0;
    DistLoadCmd() {
        name = "TẢI PHÂN BỐ";
        filter = 2;
    }
    void onSelected(App&) override { step = 0; }
    std::string afterPrompt(App& a) override {
        if (step == 0) return fmt("q (kN/m) hoặc q_i,q_j <%s,%s>:", num(a.defQ1).c_str(), num(a.defQ2).c_str());
        return fmt("Hướng [GX/GY/LX/LY/HC(đứng, theo hình chiếu)] <%s>:", dirToken(a.defDir, a.defProj));
    }
    void apply(App& a) {
        a.snapshot();
        for (int e : bars(a)) {
            ElemLoad l;
            l.elem = e;
            l.kind = ElemLoadKind::Distributed;
            l.dir = (LoadDir)a.defDir;
            l.q1 = a.defQ1;
            l.q2 = a.defQ2;
            l.projected = a.defProj;
            a.model.cases[a.activeCase].elem.push_back(l);
        }
        a.changed();
        a.selElems.clear();
        done = true;
    }
    bool onText(App& a, const std::string& s) override {
        if (selecting) return false;
        if (step == 0) {
            double v[2];
            const int n = numbers(s, v, 2);
            if (n == 0) return false;
            a.defQ1 = v[0];
            a.defQ2 = n == 2 ? v[1] : v[0];
            step = 1;
            return true;
        }
        bool pr;
        int d;
        if (!parseDir(s, d, pr)) return false;
        a.defDir = d;
        a.defProj = pr;
        apply(a);
        return true;
    }
    void afterEnter(App& a) override {
        if (step == 0) step = 1;
        else apply(a);
    }
};

struct PointLoadCmd : Command {
    int step = 0;
    bool moment;
    explicit PointLoadCmd(bool m) : moment(m) { name = m ? "MÔMEN TẬP TRUNG" : "LỰC TẬP TRUNG"; }
    bool wantsPoint() const override { return step == 2; }
    std::string prompt(App& a) override {
        if (step == 0)
            return moment ? fmt("M (kNm, dương ngược chiều kim đồng hồ) <%s>:", num(a.defM).c_str())
                          : fmt("P (kN) <%s>:", num(a.defP).c_str());
        if (step == 1) return fmt("Hướng [GX/GY/LX/LY] <%s>:", dirToken(a.defPDir, false));
        return "Chỉ định điểm đặt trên thanh hoặc nút (Enter để kết thúc):";
    }
    void start(App&) override {
        if (moment) step = 0;
    }
    bool onText(App& a, const std::string& s) override {
        double v;
        if (step == 0 && number(s, v)) {
            (moment ? a.defM : a.defP) = v;
            step = moment ? 2 : 1;
            return true;
        }
        bool pr;
        int d;
        if (step == 1 && parseDir(s, d, pr)) {
            a.defPDir = d;
            step = 2;
            return true;
        }
        return false;
    }
    void onEnter(App&) override {
        if (step == 0) step = moment ? 2 : 1;
        else if (step == 1) step = 2;
        else done = true;
    }
    void onPoint(App& a, const SnapPoint& p) override {
        Model& m = a.model;
        LoadCase& lc = m.cases[a.activeCase];
        const double tol = tolOf(m);
        // tại nút → tải nút; trên thanh → tải trên thanh
        const int node = m.findNode(p.x, p.y, tol);
        if (node >= 0) {
            a.snapshot();
            NodalLoad l{node, {0, 0, 0}};
            if (moment) l.F[2] = a.defM;
            else if (a.defPDir == (int)LoadDir::GlobalX) l.F[0] = a.defP;
            else l.F[1] = a.defP;
            lc.nodal.push_back(l);
            a.changed();
            return;
        }
        int best = -1;
        double bt = 0, bd = 1e300;
        for (size_t e = 0; e < m.elements.size(); ++e) {
            const Element& el = m.elements[e];
            if (!el.isBar()) continue;
            const Node& A = m.nodes[el.n[0]];
            const Node& B = m.nodes[el.n[1]];
            const double dx = B.x - A.x, dy = B.y - A.y, L2 = dx * dx + dy * dy;
            const double t = std::clamp(((p.x - A.x) * dx + (p.y - A.y) * dy) / L2, 0.0, 1.0);
            const double d = std::hypot(A.x + t * dx - p.x, A.y + t * dy - p.y);
            if (d < bd) {
                bd = d;
                best = (int)e;
                bt = t;
            }
        }
        if (best < 0 || bd > 1e-3 * std::max(1.0, m.barLength(best))) {
            a.echo("Điểm không nằm trên thanh nào.");
            return;
        }
        a.snapshot();
        ElemLoad l;
        l.elem = best;
        l.kind = moment ? ElemLoadKind::Moment : ElemLoadKind::Point;
        l.dir = (LoadDir)a.defPDir;
        l.q1 = moment ? a.defM : a.defP;
        l.a = bt * m.barLength(best);
        lc.elem.push_back(l);
        a.changed();
    }
};

// ---------------------------------------------------------------- tiện ích
struct DistCmd : Command {
    int step = 0;
    double x1 = 0, y1 = 0;
    DistCmd() { name = "KHOẢNG CÁCH"; }
    bool wantsPoint() const override { return true; }
    bool base(double& x, double& y) const override {
        if (step == 0) return false;
        x = x1;
        y = y1;
        return true;
    }
    std::string prompt(App&) override { return step == 0 ? "Điểm thứ nhất:" : "Điểm thứ hai:"; }
    void onPoint(App& a, const SnapPoint& p) override {
        if (step == 0) {
            x1 = p.x;
            y1 = p.y;
            step = 1;
            return;
        }
        const double dx = p.x - x1, dy = p.y - y1;
        a.echo(fmt("Khoảng cách = %.4f m,  góc = %.2f°,  ΔX = %.4f,  ΔY = %.4f", std::hypot(dx, dy),
                   std::atan2(dy, dx) * 180 / kPi, dx, dy));
        done = true;
    }
    void preview(App& a, ImDrawList* dl, const SnapPoint& c) override {
        if (step == 1) dl->AddLine(a.cam.toScreen(x1, y1), a.cam.toScreen(c.x, c.y), a.pal.hover, 1.0f);
    }
};

struct IdCmd : Command {
    IdCmd() { name = "TỌA ĐỘ"; }
    bool wantsPoint() const override { return true; }
    std::string prompt(App&) override { return "Chỉ định điểm:"; }
    void onPoint(App& a, const SnapPoint& p) override {
        a.echo(fmt("X = %.4f   Y = %.4f%s", p.x, p.y, p.node >= 0 ? fmt("   (nút %d)", p.node + 1).c_str() : ""));
        done = true;
    }
};

struct ListCmd : SelectCmd {
    ListCmd() { name = "LIỆT KÊ"; }
    void onSelected(App& a) override {
        const Model& m = a.model;
        int shown = 0;
        for (int n : a.selNodes) {
            if (++shown > 30) break;
            const Support* s = m.supportAt(n);
            a.echo(fmt("  Nút %d: (%.4f; %.4f)%s", n + 1, m.nodes[n].x, m.nodes[n].y,
                       s ? fmt("  gối [%d%d%d]", s->fix[0], s->fix[1], s->fix[2]).c_str() : ""));
        }
        for (int e : a.selElems) {
            if (++shown > 30) break;
            const Element& el = m.elements[e];
            std::string ns;
            for (int k = 0; k < el.nodeCount(); ++k) ns += (k ? "-" : "") + std::to_string(el.n[k] + 1);
            std::string extra;
            if (el.isBar()) {
                extra = fmt("L = %.4f m, %s, %s", m.barLength(e), m.materials[el.mat].name.c_str(),
                            el.sec < (int)m.sections.size() ? m.sections[el.sec].name.c_str() : "?");
                if (el.hingeI || el.hingeJ) extra += el.hingeI && el.hingeJ ? ", khớp 2 đầu" : el.hingeI ? ", khớp i" : ", khớp j";
            } else {
                extra = fmt("t = %.3g m, %s", el.t, m.materials[el.mat].name.c_str());
            }
            a.echo(fmt("  Phần tử %d (%s): nút %s, %s", e + 1, elemTypeName(el.type), ns.c_str(), extra.c_str()));
        }
        done = true;
    }
};

struct ZoomCmd : Command {
    int step = 0;
    double x1 = 0, y1 = 0;
    ZoomCmd() { name = "ZOOM"; }
    bool wantsPoint() const override { return step > 0; }
    std::string prompt(App&) override {
        if (step == 0) return "[Toàn bộ(E)/Cửa sổ(W)] <E>:";
        return step == 1 ? "Góc thứ nhất:" : "Góc đối diện:";
    }
    bool onText(App& a, const std::string& s) override {
        const std::string u = upper(s);
        if (step != 0) return false;
        if (u == "E" || u == "A") {
            a.fitRequested = true;
            done = true;
        } else if (u == "W") step = 1;
        else return false;
        return true;
    }
    void onEnter(App& a) override {
        if (step == 0) a.fitRequested = true;
        done = true;
    }
    void onPoint(App& a, const SnapPoint& p) override {
        if (step == 1) {
            x1 = p.x;
            y1 = p.y;
            step = 2;
            return;
        }
        const double w = std::abs(p.x - x1), h = std::abs(p.y - y1);
        if (w > 0 && h > 0) {
            a.cam.cx = 0.5 * (p.x + x1);
            a.cam.cy = 0.5 * (p.y + y1);
            a.cam.scale = std::min((a.cam.c1.x - a.cam.c0.x) / w, (a.cam.c1.y - a.cam.c0.y) / h);
        }
        done = true;
    }
    void preview(App& a, ImDrawList* dl, const SnapPoint& c) override {
        if (step == 2) dl->AddRect(a.cam.toScreen(x1, y1), a.cam.toScreen(c.x, c.y), a.pal.hover);
    }
};

} // namespace

// ---------------------------------------------------------------- bảng lệnh
std::string upper(const std::string& s) {
    std::string r = s;
    for (char& c : r)
        if (c >= 'a' && c <= 'z') c = (char)(c - 'a' + 'A');
    return r;
}

const char* snapName(SnapKind k) {
    switch (k) {
    case SnapKind::Node: return "Nút";
    case SnapKind::Mid: return "Trung điểm";
    case SnapKind::Intersection: return "Giao điểm";
    case SnapKind::Perpendicular: return "Vuông góc";
    case SnapKind::Nearest: return "Gần nhất";
    case SnapKind::Grid: return "Lưới";
    default: return "";
    }
}

std::unique_ptr<Command> makeCommand(const std::string& u, std::string& canonical) {
    auto is = [&](std::initializer_list<const char*> names) {
        for (const char* n : names)
            if (u == n) {
                canonical = *names.begin();
                return true;
            }
        return false;
    };
    if (is({"L", "LINE", "THANH", "TK"})) return std::make_unique<LineCmd>(false);
    if (is({"TRUSS", "DAN", "TD"})) return std::make_unique<LineCmd>(true);
    if (is({"PO", "POINT", "NODE", "NUT"})) return std::make_unique<NodeCmd>();
    if (is({"E", "ERASE", "XOA"})) return std::make_unique<EraseCmd>();
    if (is({"M", "MOVE"})) return std::make_unique<MoveCopyCmd>(false);
    if (is({"CO", "CP", "COPY"})) return std::make_unique<MoveCopyCmd>(true);
    if (is({"RO", "ROTATE", "XOAY"})) return std::make_unique<RotateCmd>();
    if (is({"MI", "MIRROR"})) return std::make_unique<MirrorCmd>();
    if (is({"DIV", "DIVIDE", "CHIA"})) return std::make_unique<DivideCmd>();
    if (is({"NGAM", "FIX"})) return std::make_unique<SupportCmd>("NGÀM", std::array<bool, 3>{true, true, true});
    if (is({"KHOP", "PIN"})) return std::make_unique<SupportCmd>("GỐI CỐ ĐỊNH", std::array<bool, 3>{true, true, false});
    if (is({"DIDONG", "ROLLER"})) return std::make_unique<SupportCmd>("GỐI DI ĐỘNG", std::array<bool, 3>{false, true, false});
    if (is({"DIDONGX", "ROLLERX"}))
        return std::make_unique<SupportCmd>("GỐI DI ĐỘNG NGANG", std::array<bool, 3>{true, false, false});
    if (is({"NGAMTRUOT", "GUIDE"})) return std::make_unique<SupportCmd>("NGÀM TRƯỢT", std::array<bool, 3>{true, false, true});
    if (is({"BOGOI", "FREE"})) return std::make_unique<SupportCmd>("BỎ GỐI", std::array<bool, 3>{false, false, false}, true);
    if (is({"TAINUT", "NLOAD"})) return std::make_unique<NodalLoadCmd>();
    if (is({"TAIPB", "DLOAD"})) return std::make_unique<DistLoadCmd>();
    if (is({"TAITT", "PLOAD"})) return std::make_unique<PointLoadCmd>(false);
    if (is({"MOMEN", "MLOAD"})) return std::make_unique<PointLoadCmd>(true);
    if (is({"DI", "DIST", "KC"})) return std::make_unique<DistCmd>();
    if (is({"ID"})) return std::make_unique<IdCmd>();
    if (is({"LI", "LIST"})) return std::make_unique<ListCmd>();
    if (is({"Z", "ZOOM"})) return std::make_unique<ZoomCmd>();
    return nullptr;
}

std::vector<std::pair<std::string, std::string>> commandHelp() {
    return {
        {"L, THANH", "vẽ thanh khung (chuỗi điểm)"},
        {"TRUSS, DAN", "vẽ thanh dàn"},
        {"PO, NUT", "thêm nút"},
        {"E, XOA", "xóa đối tượng"},
        {"M / CO", "di chuyển / sao chép"},
        {"RO / MI", "xoay / đối xứng"},
        {"DIV, CHIA", "chia thanh thành n đoạn"},
        {"GOP", "gộp nút trùng"},
        {"NGAM, KHOP, DIDONG, DIDONGX, NGAMTRUOT, BOGOI", "gối tựa"},
        {"TAINUT / TAIPB / TAITT / MOMEN", "tải nút / phân bố / lực tập trung / mômen"},
        {"TINH (F5)", "tính toán"},
        {"ONDINH / DAODONG", "ổn định / dao động riêng"},
        {"DI, ID, LI", "đo khoảng cách / tọa độ / liệt kê"},
        {"Z", "zoom (E: toàn bộ, W: cửa sổ)"},
        {"U / REDO", "hoàn tác / làm lại"},
        {"LUOI, TAOKHUNG, TAODAN", "tạo lưới phẳng / khung / dàn"},
        {"PR (Ctrl+1)", "bảng thuộc tính"},
        {"NEW, OPEN, SAVE, SAVEAS", "tệp"},
    };
}

bool parsePoint(const std::string& in, bool hasBase, double bx, double by, double dirX, double dirY, double& x, double& y) {
    std::string t = trim(in);
    if (t.empty()) return false;
    bool rel = false;
    if (t[0] == '@') {
        if (!hasBase) return false;
        rel = true;
        t = t.substr(1);
    } else if (t[0] == '#') {
        t = t.substr(1);
    }
    double dx, dy;
    const size_t lt = t.find('<');
    if (lt != std::string::npos) {
        double d, ang;
        if (!number(t.substr(0, lt), d) || !number(t.substr(lt + 1), ang)) return false;
        dx = d * std::cos(ang * kPi / 180);
        dy = d * std::sin(ang * kPi / 180);
    } else if (t.find(',') != std::string::npos) {
        double v[2];
        if (numbers(t, v, 2) != 2) return false;
        dx = v[0];
        dy = v[1];
    } else {
        // nhập khoảng cách trực tiếp theo hướng con trỏ
        double d;
        if (!hasBase || rel || !number(t, d)) return false;
        const double L = std::hypot(dirX, dirY);
        if (L <= 0) return false;
        x = bx + d * dirX / L;
        y = by + d * dirY / L;
        return true;
    }
    x = rel ? bx + dx : dx;
    y = rel ? by + dy : dy;
    return true;
}
