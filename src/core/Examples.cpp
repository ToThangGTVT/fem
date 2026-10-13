#include "core/Examples.h"

#include <cmath>

#include "core/Library.h"
#include "core/Mesher.h"

namespace fem {

namespace {

const MaterialPreset& preset(const char* name) {
    for (const auto& p : materialPresets())
        if (std::string(p.name) == name) return p;
    return materialPresets()[0];
}

const ProfilePreset& profile(const char* name) {
    for (const auto& p : profilePresets())
        if (std::string(p.name) == name) return p;
    return profilePresets()[0];
}

int frame(Model& m, int a, int b, int sec = 0, bool hi = false, bool hj = false) {
    Element e;
    e.type = ElemType::Frame;
    e.n = {a, b, -1, -1};
    e.mat = 0;
    e.sec = sec;
    e.hingeI = hi;
    e.hingeJ = hj;
    m.elements.push_back(e);
    return (int)m.elements.size() - 1;
}

int truss(Model& m, int a, int b, int sec = 0) {
    Element e;
    e.type = ElemType::Truss;
    e.n = {a, b, -1, -1};
    e.sec = sec;
    m.elements.push_back(e);
    return (int)m.elements.size() - 1;
}

void support(Model& m, int node, bool fx, bool fy, bool fr) {
    Support& s = m.supportFor(node);
    s.fix = {fx, fy, fr};
}

ElemLoad dist(int e, double q1, double q2, LoadDir dir = LoadDir::GlobalY) {
    ElemLoad l;
    l.elem = e;
    l.kind = ElemLoadKind::Distributed;
    l.dir = dir;
    l.q1 = q1;
    l.q2 = q2;
    return l;
}

ElemLoad point(int e, double P, double a, LoadDir dir = LoadDir::GlobalY) {
    ElemLoad l;
    l.elem = e;
    l.kind = ElemLoadKind::Point;
    l.dir = dir;
    l.q1 = P;
    l.a = a;
    return l;
}

Model simpleBeam() {
    Model m;
    m.title = "Dầm đơn giản nhịp 6 m, q = 10 kN/m (Mmax = qL²/8 = 45 kNm)";
    m.materials.push_back(toMaterial(preset("Thép CCT34")));
    m.sections.push_back(toSection(profile("I30")));
    int a = m.addNode(0, 0), b = m.addNode(6, 0);
    int e = frame(m, a, b);
    support(m, a, true, true, false);
    support(m, b, false, true, false);
    m.cases.push_back({"Tĩnh tải"});
    m.cases[0].elem.push_back(dist(e, -10, -10));
    return m;
}

Model continuousBeam() {
    Model m;
    m.title = "Dầm liên tục 3 nhịp 5 + 6 + 5 m — hoạt tải cách nhịp, biểu đồ bao";
    m.materials.push_back(toMaterial(preset("Bê tông B25")));
    m.sections.push_back(rectSection(0.3, 0.6));
    const double xs[4] = {0, 5, 11, 16};
    int n[4];
    for (int i = 0; i < 4; ++i) n[i] = m.addNode(xs[i], 0);
    int e[3];
    for (int i = 0; i < 3; ++i) e[i] = frame(m, n[i], n[i + 1]);
    support(m, n[0], true, true, false);
    for (int i = 1; i < 4; ++i) support(m, n[i], false, true, false);
    LoadCase tt{"Tĩnh tải", 1.0};
    for (int i = 0; i < 3; ++i) tt.elem.push_back(dist(e[i], -20, -20));
    // các sơ đồ chất hoạt tải bất lợi: nhịp 1+3 (Mmax nhịp biên), nhịp 2 (Mmax nhịp giữa), nhịp 1+2 (Mmin gối 2)
    const double p = -15;
    LoadCase h1{"Hoạt tải nhịp 1 + 3"};
    h1.elem = {dist(e[0], p, p), dist(e[2], p, p)};
    LoadCase h2{"Hoạt tải nhịp 2"};
    h2.elem = {dist(e[1], p, p), point(e[1], -50, 3)};
    LoadCase h3{"Hoạt tải nhịp 1 + 2"};
    h3.elem = {dist(e[0], p, p), dist(e[1], p, p)};
    LoadCase h4{"Hoạt tải nhịp 2 + 3"};
    h4.elem = {dist(e[1], p, p), dist(e[2], p, p)};
    m.cases = {tt, h1, h2, h3, h4};
    for (int i = 1; i <= 4; ++i) m.combos.push_back({"1,1TT + 1,2·" + m.cases[i].name, {{0, 1.1}, {i, 1.2}}});
    return m;
}

Model portalFrame() {
    Model m;
    m.title = "Khung một tầng một nhịp, chân ngàm";
    m.materials.push_back(toMaterial(preset("Bê tông B20")));
    m.sections.push_back(rectSection(0.3, 0.4)); // cột
    m.sections.push_back(rectSection(0.3, 0.6)); // dầm
    int a = m.addNode(0, 0), b = m.addNode(0, 4), c = m.addNode(6, 4), d = m.addNode(6, 0);
    frame(m, a, b, 0);
    int beam = frame(m, b, c, 1);
    frame(m, d, c, 0);
    support(m, a, true, true, true);
    support(m, d, true, true, true);
    LoadCase tt{"Tĩnh tải"};
    tt.elem.push_back(dist(beam, -15, -15));
    LoadCase gio{"Gió trái"};
    gio.nodal.push_back({b, {10, 0, 0}});
    gio.elem.push_back(dist(0, 2.4, 2.4, LoadDir::GlobalX));
    m.cases = {tt, gio};
    m.combos.push_back({"TT + 0,9Gió", {{0, 1.0}, {1, 0.9}}});
    return m;
}

Model threeHingedFrame() {
    Model m;
    m.title = "Khung ba khớp nhịp 8 m, cao 5 m";
    m.materials.push_back(toMaterial(preset("Thép CCT38")));
    m.sections.push_back(toSection(profile("I36")));
    int a = m.addNode(0, 0), b = m.addNode(0, 5), c = m.addNode(4, 6), d = m.addNode(8, 5), e = m.addNode(8, 0);
    frame(m, a, b);
    int r1 = frame(m, b, c, 0, false, true);
    int r2 = frame(m, c, d);
    frame(m, e, d);
    support(m, a, true, true, false);
    support(m, e, true, true, false);
    LoadCase lc{"Tải mái"};
    ElemLoad q1 = dist(r1, -12, -12);
    q1.projected = true;
    ElemLoad q2 = dist(r2, -12, -12);
    q2.projected = true;
    lc.elem = {q1, q2};
    m.cases = {lc};
    return m;
}

Model twoStoreyFrame() {
    Model m;
    m.title = "Khung 2 tầng 2 nhịp";
    m.materials.push_back(toMaterial(preset("Bê tông B25")));
    m.sections.push_back(rectSection(0.3, 0.3)); // cột
    m.sections.push_back(rectSection(0.25, 0.5)); // dầm
    const double xs[3] = {0, 6, 12}, ys[3] = {0, 3.6, 7.2};
    int n[3][3];
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) n[j][i] = m.addNode(xs[i], ys[j]);
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 3; ++i) frame(m, n[j][i], n[j + 1][i], 0);
    std::vector<int> beams;
    for (int j = 1; j < 3; ++j)
        for (int i = 0; i < 2; ++i) beams.push_back(frame(m, n[j][i], n[j][i + 1], 1));
    for (int i = 0; i < 3; ++i) support(m, n[0][i], true, true, true);
    LoadCase tt{"Tĩnh tải", 1.0};
    for (int b : beams) tt.elem.push_back(dist(b, -18, -18));
    LoadCase ht{"Hoạt tải"};
    for (int b : beams) ht.elem.push_back(dist(b, -9, -9));
    LoadCase gio{"Gió X"};
    gio.nodal.push_back({n[1][0], {12, 0, 0}});
    gio.nodal.push_back({n[2][0], {7, 0, 0}});
    m.cases = {tt, ht, gio};
    m.combos.push_back({"TH cơ bản 1: TT + HT", {{0, 1.0}, {1, 1.0}}});
    m.combos.push_back({"TH cơ bản 2: TT + 0,9(HT + Gió)", {{0, 1.0}, {1, 0.9}, {2, 0.9}}});
    return m;
}

Model bridgeTruss() {
    Model m;
    m.title = "Dàn cầu 6 khoang, nhịp 18 m, cao 3 m";
    m.materials.push_back(toMaterial(preset("Thép CCT38")));
    Section s;
    s.name = "2L75x6";
    s.A = 17.56e-4;
    s.I = 0;
    s.h = 0.075;
    m.sections.push_back(s);
    const int nP = 6;
    const double d = 3, h = 3;
    std::vector<int> bot, top;
    for (int i = 0; i <= nP; ++i) bot.push_back(m.addNode(i * d, 0));
    for (int i = 1; i < nP; ++i) top.push_back(m.addNode(i * d, h));
    for (int i = 0; i < nP; ++i) truss(m, bot[i], bot[i + 1]);
    for (int i = 0; i + 1 < (int)top.size(); ++i) truss(m, top[i], top[i + 1]);
    truss(m, bot[0], top[0]);
    truss(m, bot[nP], top.back());
    for (int i = 1; i < nP; ++i) truss(m, bot[i], top[i - 1]); // thanh đứng
    for (int i = 1; i < nP - 1; ++i) {
        // thanh xiên (dàn Pratt: xiên hướng về giữa nhịp)
        if (i < nP / 2) truss(m, top[i - 1], bot[i + 1]);
        else truss(m, bot[i], top[i]);
    }
    support(m, bot[0], true, true, false);
    support(m, bot[nP], false, true, false);
    LoadCase lc{"Tải nút dàn dưới"};
    for (int i = 1; i < nP; ++i) lc.nodal.push_back({bot[i], {0, -40, 0}});
    m.cases = {lc};
    return m;
}

Model gerberBeam() {
    Model m;
    m.title = "Dầm tĩnh định nhiều nhịp (dầm Gerber)";
    m.materials.push_back(toMaterial(preset("Thép CCT34")));
    m.sections.push_back(toSection(profile("I27")));
    const double xs[6] = {0, 6, 8, 12, 14, 20};
    int n[6];
    for (int i = 0; i < 6; ++i) n[i] = m.addNode(xs[i], 0);
    std::vector<int> e;
    e.push_back(frame(m, n[0], n[1]));
    e.push_back(frame(m, n[1], n[2], 0, false, true)); // khớp tại x = 8
    e.push_back(frame(m, n[2], n[3]));
    e.push_back(frame(m, n[3], n[4], 0, false, true)); // khớp tại x = 14
    e.push_back(frame(m, n[4], n[5]));
    support(m, n[0], true, true, false);
    support(m, n[1], false, true, false);
    support(m, n[3], false, true, false);
    support(m, n[5], false, true, false);
    LoadCase lc{"Tải trọng"};
    for (int k : e) lc.elem.push_back(dist(k, -10, -10));
    lc.elem.push_back(point(e[2], -30, 2));
    m.cases = {lc};
    return m;
}

Model winklerBeam() {
    Model m;
    m.title = "Dầm trên nền đàn hồi Winkler (k = 20000 kN/m³, b = 1 m)";
    m.materials.push_back(toMaterial(preset("Bê tông B20")));
    m.sections.push_back(rectSection(1.0, 0.6));
    const int n = 40;
    const double L = 12;
    std::vector<int> nodes;
    for (int i = 0; i <= n; ++i) nodes.push_back(m.addNode(L * i / n, 0));
    LoadCase lc{"Tải tập trung"};
    for (int i = 0; i < n; ++i) {
        int e = frame(m, nodes[i], nodes[i + 1]);
        m.elements[e].kf = 20000;
    }
    support(m, nodes[n / 2], true, false, false); // chỉ chặn chuyển vị dọc trục
    lc.nodal.push_back({nodes[n / 2], {0, -400, 0}});
    lc.nodal.push_back({nodes[n / 8], {0, -150, 0}});
    m.cases = {lc};
    return m;
}

Model deepBeam() {
    Model m;
    m.title = "Dầm-tường (балка-стенка) 6 × 3 m, ứng suất phẳng";
    m.materials.push_back(toMaterial(preset("Bê tông B25")));
    const int nx = 24, ny = 12;
    const double W = 6, H = 3, t = 0.2;
    std::vector<int> g = meshRectangle(m, 0, 0, W, H, nx, ny, ElemType::Quad4, 0, t);
    // gối có bản đệm rộng 0,5 m ở mỗi đầu (giảm ứng suất tập trung tại điểm gối)
    support(m, g[0], true, true, false);
    support(m, g[1], false, true, false);
    support(m, g[2], false, true, false);
    for (int i = nx - 2; i <= nx; ++i) support(m, g[i], false, true, false);
    LoadCase lc{"Tải phân bố trên mép trên", 1.0};
    const double q = -100; // kN/m
    for (int i = 0; i <= nx; ++i) {
        const double w = (i == 0 || i == nx) ? 0.5 : 1.0;
        lc.nodal.push_back({g[ny * (nx + 1) + i], {0, q * w * W / nx, 0}});
    }
    m.cases = {lc};
    return m;
}

Model gravityDam() {
    Model m;
    m.title = "Đập bê tông trọng lực, biến dạng phẳng (1 m dài)";
    m.planeMode = PlaneMode::Strain;
    m.materials.push_back(toMaterial(preset("Bê tông B15")));
    const int nx = 8, ny = 12;
    const double cx[4] = {0, 8, 2, 0}, cy[4] = {0, 0, 12, 12};
    std::vector<int> g = meshPatch(m, cx, cy, nx, ny, ElemType::Tri3, 0, 1.0);
    for (int i = 0; i <= nx; ++i) support(m, g[i], true, true, false);
    LoadCase tl{"Trọng lượng bản thân", 1.0};
    LoadCase water{"Áp lực nước thượng lưu"};
    const double Hw = 11, gw = 10;
    for (int j = 0; j < ny; ++j) {
        const int na = g[j * (nx + 1)], nb = g[(j + 1) * (nx + 1)];
        const double ya = m.nodes[na].y, yb = m.nodes[nb].y, L = yb - ya;
        const double pa = gw * std::max(0.0, Hw - ya), pb = gw * std::max(0.0, Hw - yb);
        water.nodal.push_back({na, {L * (2 * pa + pb) / 6, 0, 0}});
        water.nodal.push_back({nb, {L * (pa + 2 * pb) / 6, 0, 0}});
    }
    m.cases = {tl, water};
    m.combos.push_back({"TLBT + Nước", {{0, 1.0}, {1, 1.0}}});
    return m;
}

Model settlementFrame() {
    Model m;
    m.title = "Khung chịu lún gối và nhiệt độ";
    m.materials.push_back(toMaterial(preset("Bê tông B25")));
    m.sections.push_back(rectSection(0.3, 0.5));
    int a = m.addNode(0, 0), b = m.addNode(0, 4), c = m.addNode(8, 4), d = m.addNode(8, 0);
    frame(m, a, b);
    int beam = frame(m, b, c);
    frame(m, d, c);
    support(m, a, true, true, true);
    support(m, d, true, true, false);
    LoadCase lun{"Lún gối phải 2 cm"};
    lun.disp.push_back({d, UY, -0.02});
    LoadCase nhiet{"Nhiệt độ dầm: t = +20°C, Δt = 15°C"};
    ElemLoad tl;
    tl.elem = beam;
    tl.kind = ElemLoadKind::Temperature;
    tl.q1 = 20;
    tl.q2 = 15;
    nhiet.elem.push_back(tl);
    m.cases = {lun, nhiet};
    return m;
}

} // namespace

const std::vector<Example>& examples() {
    static const std::vector<Example> v = {
        {"Dầm đơn giản", "Kiểm tra Mmax = qL²/8", simpleBeam},
        {"Dầm liên tục 3 nhịp", "Hoạt tải cách nhịp, biểu đồ bao", continuousBeam},
        {"Dầm tĩnh định nhiều nhịp", "Dầm Gerber có khớp", gerberBeam},
        {"Khung 1 tầng 1 nhịp", "Tĩnh tải + gió", portalFrame},
        {"Khung ba khớp", "Tải phân bố theo hình chiếu", threeHingedFrame},
        {"Khung 2 tầng 2 nhịp", "3 trường hợp tải, 2 tổ hợp", twoStoreyFrame},
        {"Khung chịu lún và nhiệt", "Chuyển vị cưỡng bức, tải nhiệt độ", settlementFrame},
        {"Dàn cầu", "Dàn Pratt 6 khoang", bridgeTruss},
        {"Dầm trên nền đàn hồi", "Mô hình nền Winkler", winklerBeam},
        {"Dầm-tường (ứng suất phẳng)", "Lưới Q4 24×12", deepBeam},
        {"Đập trọng lực (biến dạng phẳng)", "Lưới tam giác CST", gravityDam},
    };
    return v;
}

} // namespace fem
