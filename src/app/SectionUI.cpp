// Sức bền vật liệu: thiết kế tiết diện (đặc trưng hình học, lõi tiết diện) và ứng suất trên tiết diện (σ, τ, Mohr).
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <sstream>

#include "App.h"
#include "core/Library.h"
#include "core/PostProcess.h"
#include "core/SectionGeom.h"
#include "imgui_stdlib.h"

using namespace fem;

namespace {

constexpr double kPi = 3.14159265358979323846;

ImU32 withA(ImU32 c, int a) { return (c & 0x00FFFFFF) | ((ImU32)a << 24); }

// Khung vẽ: ánh xạ tọa độ tiết diện (m) sang màn hình
struct Frame2 {
    ImVec2 o;
    double s = 1, x0 = 0, y0 = 0;
    ImVec2 P(double x, double y) const { return ImVec2((float)(o.x + (x - x0) * s), (float)(o.y - (y - y0) * s)); }
};

void drawContours(ImDrawList* dl, const std::vector<Contour>& cs, const Frame2& f, ImU32 fill, ImU32 /*bg*/, ImU32 line) {
    // tô bằng quét hàng (chẵn–lẻ trên mọi chu tuyến: lỗ được trừ tự nhiên)
    std::vector<std::vector<ImVec2>> scr;
    float y0 = FLT_MAX, y1 = -FLT_MAX;
    for (const auto& c : cs) {
        std::vector<ImVec2> pts;
        for (const auto& v : c) {
            pts.push_back(f.P(v[0], v[1]));
            y0 = std::min(y0, pts.back().y);
            y1 = std::max(y1, pts.back().y);
        }
        scr.push_back(pts);
    }
    std::vector<float> xs;
    for (float y = std::floor(y0) + 0.5f; y < y1; y += 1.0f) {
        xs.clear();
        for (const auto& pts : scr)
            for (size_t i = 0; i < pts.size(); ++i) {
                const ImVec2 a = pts[i], b = pts[(i + 1) % pts.size()];
                if ((a.y <= y && b.y > y) || (b.y <= y && a.y > y)) xs.push_back(a.x + (y - a.y) / (b.y - a.y) * (b.x - a.x));
            }
        std::sort(xs.begin(), xs.end());
        for (size_t i = 0; i + 1 < xs.size(); i += 2) dl->AddRectFilled(ImVec2(xs[i], y - 0.5f), ImVec2(xs[i + 1], y + 0.5f), fill);
    }
    for (const auto& pts : scr) dl->AddPolyline(pts.data(), (int)pts.size(), line, 1.5f, ImDrawFlags_Closed);
}

void dashed(ImDrawList* dl, ImVec2 a, ImVec2 b, ImU32 c, float th = 1.0f) {
    const float dx = b.x - a.x, dy = b.y - a.y, L = std::sqrt(dx * dx + dy * dy);
    if (L <= 0) return;
    for (float t = 0; t < L; t += 9) {
        const float t1 = std::min(t + 5, L);
        dl->AddLine(ImVec2(a.x + dx * t / L, a.y + dy * t / L), ImVec2(a.x + dx * t1 / L, a.y + dy * t1 / L), c, th);
    }
}

std::string fmtCm(double v, int pw, double ref = 0) { // giá trị m^pw → cm^pw; |v| ≪ ref coi là 0
    char b[48];
    if (ref > 0 && std::abs(v) < 1e-9 * ref) v = 0;
    std::snprintf(b, sizeof b, "%.4g", v * std::pow(100.0, pw));
    return b;
}

// Phân tích chuỗi chu tuyến: mỗi dòng "x,y x,y ..." (cm)
std::vector<Contour> parseContours(const std::string& text) {
    std::vector<Contour> out;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        Contour c;
        std::istringstream ls(line);
        std::string tok;
        while (ls >> tok) {
            const size_t k = tok.find(',');
            if (k == std::string::npos) continue;
            c.push_back({std::strtod(tok.substr(0, k).c_str(), nullptr) / 100, std::strtod(tok.substr(k + 1).c_str(), nullptr) / 100});
        }
        if (c.size() >= 3) out.push_back(c);
    }
    return out;
}

std::string contoursText(const std::vector<Contour>& cs) {
    std::string s;
    char b[64];
    for (const auto& c : cs) {
        for (const auto& v : c) {
            std::snprintf(b, sizeof b, "%g,%g ", v[0] * 100, v[1] * 100);
            s += b;
        }
        s += "\n";
    }
    return s;
}

} // namespace

// ============================================================================ thiết kế tiết diện
void App::sectionDesigner() {
    ImGui::SetNextWindowSize(ImVec2(820, 560), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Thiết kế tiết diện", &showSectionDesigner)) {
        ImGui::End();
        return;
    }
    static SectionShape sh = [] {
        SectionShape s;
        s.kind = ShapeKind::Tee;
        s.p = {0.4, 0.1, 0.2, 0.6, 0, 0};
        return s;
    }();
    static std::string polyText = "0,0 30,0 30,50 0,50\n5,5 5,45 25,45 25,5\n";
    static std::string name;
    static int target = 0;

    ImGui::BeginChild("##left", ImVec2(330, 0), ImGuiChildFlags_Borders);
    ImGui::TextDisabled("Kích thước theo cm");
    const ShapeKind kinds[] = {ShapeKind::Rect, ShapeKind::Tee, ShapeKind::IBeam, ShapeKind::Box, ShapeKind::Pipe,
                               ShapeKind::Circle, ShapeKind::Channel, ShapeKind::Angle, ShapeKind::Polygon};
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::BeginCombo("##kind", shapeName(sh.kind))) {
        for (auto k : kinds)
            if (ImGui::Selectable(shapeName(k), k == sh.kind)) {
                sh.kind = k;
                // giá trị mặc định hợp lý
                switch (k) {
                case ShapeKind::Rect: sh.p = {0.3, 0.6, 0, 0, 0, 0}; break;
                case ShapeKind::Tee: sh.p = {0.4, 0.1, 0.2, 0.6, 0, 0}; break;
                case ShapeKind::IBeam: sh.p = {0.2, 0.012, 0.008, 0.4, 0.3, 0.016}; break;
                case ShapeKind::Box: sh.p = {0.2, 0.3, 0.01, 0, 0, 0}; break;
                case ShapeKind::Pipe: sh.p = {0.219, 0.008, 0, 0, 0, 0}; break;
                case ShapeKind::Circle: sh.p = {0.3, 0, 0, 0, 0, 0}; break;
                case ShapeKind::Channel: sh.p = {0.08, 0.2, 0.0052, 0.009, 0, 0}; break;
                case ShapeKind::Angle: sh.p = {0.1, 0.1, 0.01, 0, 0, 0}; break;
                default: break;
                }
            }
        ImGui::EndCombo();
    }
    for (int i = 0; i < shapeParamCount(sh.kind); ++i) {
        double cm = sh.p[i] * 100;
        ImGui::SetNextItemWidth(110);
        if (ImGui::InputDouble(shapeParamName(sh.kind, i), &cm, 0, 0, "%.4g")) sh.p[i] = std::max(0.0, cm / 100);
    }
    if (sh.kind == ShapeKind::Polygon) {
        ImGui::TextWrapped("Mỗi dòng một chu tuyến \"x,y x,y ...\" (cm). Ngược chiều kim đồng hồ: phần đặc; thuận chiều: lỗ.");
        ImGui::InputTextMultiline("##poly", &polyText, ImVec2(-FLT_MIN, 140));
        sh.contours = parseContours(polyText);
    }
    const std::vector<Contour> cs = shapeContours(sh);
    const SectionProps p = sectionProps(cs);
    ImGui::SeparatorText("Đưa vào mô hình");
    ImGui::SetNextItemWidth(-60);
    ImGui::InputTextWithHint("Tên", "tự đặt theo kích thước", &name);
    ImGui::BeginDisabled(!p.ok);
    if (ImGui::Button("Thêm tiết diện mới")) {
        snapshot();
        Section s = sh.kind == ShapeKind::Polygon ? sectionFromShape("Đa giác", sh) : shapeSection(sh.kind, sh.p);
        if (sh.kind == ShapeKind::Rect) s = rectSection(sh.p[0], sh.p[1]);
        if (!name.empty()) s.name = name;
        model.sections.push_back(s);
        curSec = (int)model.sections.size() - 1;
        echo("Đã thêm tiết diện \"" + s.name + "\".");
        changed();
    }
    if (!model.sections.empty()) {
        target = std::clamp(target, 0, (int)model.sections.size() - 1);
        ImGui::SetNextItemWidth(170);
        if (ImGui::BeginCombo("##target", model.sections[target].name.c_str())) {
            for (int i = 0; i < (int)model.sections.size(); ++i) {
                ImGui::PushID(i);
                if (ImGui::Selectable(model.sections[i].name.c_str(), i == target)) target = i;
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        if (ImGui::Button("Ghi đè")) {
            snapshot();
            Section s = sectionFromShape(name.empty() ? model.sections[target].name : name, sh);
            if (sh.kind == ShapeKind::Rect) {
                s.A = sh.p[0] * sh.p[1];
                s.I = sh.p[0] * std::pow(sh.p[1], 3) / 12;
            }
            model.sections[target] = s;
            changed();
        }
    }
    ImGui::EndDisabled();
    if (!model.sections.empty() && model.sections[target].shape.kind != ShapeKind::None) {
        if (ImGui::Button("Nạp hình dạng từ tiết diện này")) {
            sh = model.sections[target].shape;
            if (sh.kind == ShapeKind::Polygon) polyText = contoursText(sh.contours);
            name = model.sections[target].name;
        }
    }
    ImGui::EndChild();
    ImGui::SameLine();

    // ---- hình vẽ + bảng đặc trưng
    ImGui::BeginChild("##right", ImVec2(0, 0));
    const ImVec2 c0 = ImGui::GetCursorScreenPos();
    const float W = ImGui::GetContentRegionAvail().x, H = std::max(200.0f, ImGui::GetContentRegionAvail().y * 0.55f);
    ImGui::InvisibleButton("##canvas", ImVec2(W, H));
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->AddRectFilled(c0, ImVec2(c0.x + W, c0.y + H), pal.bg);
    if (p.ok) {
        const double bw = p.xmax - p.xmin, bh = p.ymax - p.ymin;
        Frame2 f;
        f.s = std::min((W - 80) / std::max(bw, 1e-9), (H - 60) / std::max(bh, 1e-9));
        f.x0 = 0.5 * (p.xmin + p.xmax);
        f.y0 = 0.5 * (p.ymin + p.ymax);
        f.o = ImVec2(c0.x + W * 0.5f, c0.y + H * 0.5f);
        drawContours(dl, cs, f, withA(pal.plane, 255), pal.bg, pal.bar);
        // lõi tiết diện
        if (p.kern.size() >= 3) {
            std::vector<ImVec2> k;
            for (const auto& v : p.kern) k.push_back(f.P(v[0], v[1]));
            dl->AddConvexPolyFilled(k.data(), (int)k.size(), withA(pal.load, 90));
            dl->AddPolyline(k.data(), (int)k.size(), pal.load, 1.5f, ImDrawFlags_Closed);
        }
        // trục trung tâm và trục chính
        const ImVec2 C = f.P(p.xc, p.yc);
        const float R = (float)(0.6 * std::max(bw, bh) * f.s);
        dashed(dl, ImVec2(C.x - R, C.y), ImVec2(C.x + R, C.y), pal.textDim);
        dashed(dl, ImVec2(C.x, C.y - R), ImVec2(C.x, C.y + R), pal.textDim);
        dl->AddText(ImVec2(C.x + R + 4, C.y - 8), pal.textDim, "x");
        dl->AddText(ImVec2(C.x + 4, C.y - R - 16), pal.textDim, "y");
        if (std::abs(p.Jxy) > 1e-9 * (p.Jx + p.Jy)) {
            const float a = (float)(p.alpha * kPi / 180);
            const ImVec2 d(std::cos(a) * R, -std::sin(a) * R);
            dl->AddLine(ImVec2(C.x - d.x, C.y - d.y), ImVec2(C.x + d.x, C.y + d.y), pal.select, 1.5f);
            dl->AddText(ImVec2(C.x + d.x + 4, C.y + d.y - 8), pal.select, "1");
        }
        dl->AddCircleFilled(C, 4, pal.load);
        char b[64];
        std::snprintf(b, sizeof b, "h = %.4g cm", bh * 100);
        dl->AddText(ImVec2(c0.x + 8, c0.y + 6), pal.textDim, b);
        std::snprintf(b, sizeof b, "b = %.4g cm", bw * 100);
        dl->AddText(ImVec2(c0.x + 8, c0.y + 24), pal.textDim, b);
        dl->AddText(ImVec2(c0.x + 8, c0.y + H - 20), pal.load, "Lõi tiết diện (đỏ)");
    } else {
        dl->AddText(ImVec2(c0.x + 10, c0.y + 10), pal.load, "Hình dạng không hợp lệ (diện tích ≤ 0)");
    }
    if (p.ok && ImGui::BeginTable("##props", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        auto row = [&](const char* n1, const std::string& v1, const char* n2, const std::string& v2) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(n1);
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(v1.c_str());
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(n2);
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(v2.c_str());
        };
        char a[48];
        std::snprintf(a, sizeof a, "%.2f°", p.alpha);
        row("A (cm²)", fmtCm(p.A, 2), "Trọng tâm y_c từ đáy (cm)", fmtCm(p.yc - p.ymin, 1));
        row("Jx (cm⁴)", fmtCm(p.Jx, 4), "Jy (cm⁴)", fmtCm(p.Jy, 4));
        row("Jxy (cm⁴)", fmtCm(p.Jxy, 4, p.Jx + p.Jy), "Góc trục chính α", a);
        row("J1 = Jmax (cm⁴)", fmtCm(p.J1, 4), "J2 = Jmin (cm⁴)", fmtCm(p.J2, 4));
        row("W thớ trên (cm³)", fmtCm(p.Wtop, 3), "W thớ dưới (cm³)", fmtCm(p.Wbot, 3));
        row("ix (cm)", fmtCm(p.ix, 1), "iy (cm)", fmtCm(p.iy, 1));
        row("y thớ trên (cm)", fmtCm(p.yTop, 1), "y thớ dưới (cm)", fmtCm(p.yBot, 1));
        ImGui::EndTable();
    }
    ImGui::EndChild();
    ImGui::End();
}

// ============================================================================ ứng suất trên tiết diện
void App::sectionStressDialog() {
    ImGui::SetNextWindowSize(ImVec2(900, 600), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Ứng suất trên tiết diện", &showSectionStress)) {
        ImGui::End();
        return;
    }
    static int elem = 0;
    static float xr = 0.0f, yr = 0.5f;
    static double Rd = 210; // cường độ tính toán, MPa
    const LoadSetResult* rs = currentSet();
    for (int e : selElems)
        if (e < (int)model.elements.size() && model.elements[e].isBar()) {
            elem = e;
            break;
        }
    if (!rs || elem < 0 || elem >= (int)model.elements.size() || !model.elements[elem].isBar()) {
        ImGui::TextWrapped("Chọn một thanh trên sơ đồ (và tính kết quả) để xem ứng suất trên tiết diện.");
        ImGui::End();
        return;
    }
    const Element& el = model.elements[elem];
    const Section& sec = model.sections[el.sec];
    const BarResult& br = rs->bars[elem];
    const double L = model.barLength(elem);
    ImGui::Text("Thanh %d — %s — %s", elem + 1, sec.name.c_str(), rs->name.c_str());
    ImGui::SetNextItemWidth(260);
    ImGui::SliderFloat("Vị trí tiết diện x/L", &xr, 0, 1, "%.3f");
    // nhảy tới tiết diện có |M| lớn nhất
    ImGui::SameLine();
    if (ImGui::Button("|M|max")) {
        size_t k = 0;
        for (size_t i = 0; i < br.M.size(); ++i)
            if (std::abs(br.M[i]) > std::abs(br.M[k])) k = i;
        xr = (float)(br.x[k] / L);
    }
    ImGui::SameLine();
    if (ImGui::Button("|Q|max")) {
        size_t k = 0;
        for (size_t i = 0; i < br.Q.size(); ++i)
            if (std::abs(br.Q[i]) > std::abs(br.Q[k])) k = i;
        xr = (float)(br.x[k] / L);
    }
    const double x = xr * L;
    const double N = barValueAt(br, BarQty::N, x), Q = barValueAt(br, BarQty::Q, x), M = barValueAt(br, BarQty::M, x);
    ImGui::Text("x = %.3f m:   N = %.3f kN   Q = %.3f kN   M = %.3f kNm", x, N, Q, M);

    std::vector<Contour> cs = shapeContours(sec.shape);
    SectionProps p = sectionProps(cs);
    const bool haveShape = p.ok;
    double A = sec.A, J = sec.I, yTop = sec.top(), yBot = sec.bot();
    if (!haveShape) {
        // không có hình dạng: chỉ ứng suất pháp tuyến tính, τ trung bình
        p = SectionProps{};
        p.ymin = -yBot;
        p.ymax = yTop;
        p.yc = 0;
    } else {
        // dùng A, J của tiết diện trong mô hình (có thể là giá trị nhập tay) nhưng phân bố theo hình dạng
        yTop = p.yTop;
        yBot = p.yBot;
    }
    auto sigma = [&](double y) { return (N / A - M * (y - p.yc) / J) * 1e-3; }; // MPa, y hướng lên
    auto tau = [&](double y) {
        if (!haveShape) return Q / A * 1e-3;
        const double b = widthAt(cs, y);
        return b > 0 ? Q * staticMomentAbove(cs, p.yc, y) / (J * b) * 1e-3 : 0.0;
    };
    ImGui::SetNextItemWidth(260);
    ImGui::SliderFloat("Thớ khảo sát (từ dưới lên)", &yr, 0, 1, "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(90);
    ImGui::InputDouble("R (MPa)", &Rd, 0, 0, "%.4g");
    if (!haveShape) ImGui::TextDisabled("Tiết diện không có hình dạng: τ lấy trung bình Q/A. Dùng \"Thiết kế tiết diện\" để có phân bố τ.");

    // ---- vẽ: tiết diện | σ(y) | τ(y) | vòng Mohr
    const ImVec2 c0 = ImGui::GetCursorScreenPos();
    const float W = ImGui::GetContentRegionAvail().x, H = std::max(220.0f, ImGui::GetContentRegionAvail().y - 110);
    ImGui::InvisibleButton("##ss", ImVec2(W, H));
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->AddRectFilled(c0, ImVec2(c0.x + W, c0.y + H), pal.bg);
    const float pw = W / 4;
    const double ymin = p.ymin, ymax = p.ymax;
    const float top = c0.y + 30, bot = c0.y + H - 30;
    auto Y = [&](double y) { return (float)(bot - (y - ymin) / (ymax - ymin) * (bot - top)); };
    const double yf = ymin + yr * (ymax - ymin);
    // tiết diện
    if (haveShape) {
        Frame2 f;
        f.s = std::min((pw - 30) / std::max(p.xmax - p.xmin, 1e-9), (bot - top) / std::max(ymax - ymin, 1e-9));
        f.x0 = 0.5 * (p.xmin + p.xmax);
        f.y0 = 0.5 * (ymin + ymax);
        f.o = ImVec2(c0.x + pw * 0.5f, 0.5f * (top + bot));
        drawContours(dl, cs, f, withA(pal.plane, 255), pal.bg, pal.bar);
    } else {
        dl->AddRect(ImVec2(c0.x + pw * 0.3f, top), ImVec2(c0.x + pw * 0.7f, bot), pal.bar, 0, 0, 1.5f);
    }
    dashed(dl, ImVec2(c0.x + 4, Y(p.yc)), ImVec2(c0.x + 3 * pw, Y(p.yc)), pal.textDim);
    dl->AddLine(ImVec2(c0.x + 4, Y(yf)), ImVec2(c0.x + 3 * pw, Y(yf)), pal.snap, 1.5f);
    dl->AddText(ImVec2(c0.x + 6, c0.y + 6), pal.textDim, "Tiết diện");
    // biểu đồ σ và τ
    const int ns = 80;
    auto diagram = [&](float x0, const char* title, auto fn, ImU32 col) {
        double vmax = 1e-12;
        std::vector<double> v(ns + 1);
        for (int i = 0; i <= ns; ++i) {
            v[i] = fn(ymin + (ymax - ymin) * i / ns);
            vmax = std::max(vmax, std::abs(v[i]));
        }
        const float cx = x0 + pw * 0.5f, k = (float)((pw * 0.4f) / vmax);
        dl->AddLine(ImVec2(cx, top), ImVec2(cx, bot), pal.textDim);
        std::vector<ImVec2> pts;
        for (int i = 0; i <= ns; ++i) {
            const double y = ymin + (ymax - ymin) * i / ns;
            const ImVec2 a(cx, Y(y)), b(cx + (float)v[i] * k, Y(y));
            if (i % 4 == 0) dl->AddLine(a, b, withA(v[i] >= 0 ? pal.posLine : pal.negLine, 120));
            pts.push_back(b);
        }
        dl->AddPolyline(pts.data(), (int)pts.size(), col, 2.0f);
        char b[64];
        std::snprintf(b, sizeof b, "%.3f", v[ns]);
        dl->AddText(ImVec2(cx + 4, top - 18), pal.text, b);
        std::snprintf(b, sizeof b, "%.3f", v[0]);
        dl->AddText(ImVec2(cx + 4, bot + 2), pal.text, b);
        dl->AddText(ImVec2(x0 + 6, c0.y + 6), pal.textDim, title);
        return vmax;
    };
    diagram(c0.x + pw, "σ (MPa)", sigma, pal.posLine);
    const double tmax = diagram(c0.x + 2 * pw, "τ (MPa)", tau, pal.negLine);
    // vòng tròn Mohr cho thớ khảo sát
    const double s = sigma(yf), t = tau(yf);
    const double cM = s / 2, rM = std::hypot(s / 2, t), s1 = cM + rM, s3 = cM - rM;
    {
        const float mx = c0.x + 3.5f * pw, my = 0.5f * (top + bot);
        const double span = std::max({std::abs(s1), std::abs(s3), rM, 1e-9});
        const float k = (float)((pw * 0.42f) / span);
        dl->AddLine(ImVec2(c0.x + 3 * pw + 6, my), ImVec2(c0.x + 4 * pw - 6, my), pal.textDim);
        dl->AddLine(ImVec2(mx, top), ImVec2(mx, bot), pal.textDim);
        dl->AddCircle(ImVec2(mx + (float)cM * k, my), (float)rM * k, pal.select, 64, 2.0f);
        const ImVec2 A1(mx + (float)s * k, my - (float)t * k), A2(mx, my + (float)t * k);
        dl->AddLine(A1, A2, pal.textDim);
        dl->AddCircleFilled(A1, 4, pal.load);
        dl->AddCircleFilled(ImVec2(mx + (float)s1 * k, my), 4, pal.snap);
        dl->AddCircleFilled(ImVec2(mx + (float)s3 * k, my), 4, pal.snap);
        dl->AddText(ImVec2(c0.x + 3 * pw + 6, c0.y + 6), pal.textDim, "Vòng tròn Mohr");
        dl->AddText(ImVec2(mx + (float)s1 * k + 4, my + 4), pal.snap, "σ1");
        dl->AddText(ImVec2(mx + (float)s3 * k - 22, my + 4), pal.snap, "σ3");
    }
    // ---- kết quả số
    const double a0 = 0.5 * std::atan2(-2 * t, s) * 180 / kPi;
    const double td3 = std::sqrt(s * s + 4 * t * t), td4 = std::sqrt(s * s + 3 * t * t);
    ImGui::Text("Thớ y = %.2f cm so với trọng tâm:  σ = %.3f MPa,  τ = %.3f MPa", (yf - p.yc) * 100, s, t);
    ImGui::Text("Ứng suất chính: σ1 = %.3f, σ3 = %.3f MPa, τmax = %.3f MPa, góc phương chính α0 = %.2f°", s1, s3, rM, a0);
    auto ratio = [&](double v) { return Rd > 0 ? v / Rd : 0.0; };
    ImGui::Text("Thuyết bền 3 (Tresca): σtđ = √(σ² + 4τ²) = %.3f MPa  (%.1f%% R)", td3, 100 * ratio(td3));
    ImGui::Text("Thuyết bền 4 (von Mises): σtđ = √(σ² + 3τ²) = %.3f MPa  (%.1f%% R)", td4, 100 * ratio(td4));
    ImGui::Text("Toàn tiết diện: σ thớ trên %.3f, thớ dưới %.3f MPa;  |τ|max = %.3f MPa", sigma(ymax), sigma(ymin), tmax);
    ImGui::End();
}
