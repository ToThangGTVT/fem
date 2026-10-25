// Vẽ sơ đồ kết cấu, kết quả và xử lý tương tác chuột trên canvas.
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>

#include "App.h"
#include "core/Elements.h"
#include "core/Mesher.h"
#include "core/PlaneElements.h"

using namespace fem;

namespace {

constexpr float kPi = 3.14159265358979f;

Palette P;

ImVec2 operator+(ImVec2 a, ImVec2 b) { return ImVec2(a.x + b.x, a.y + b.y); }
ImVec2 operator-(ImVec2 a, ImVec2 b) { return ImVec2(a.x - b.x, a.y - b.y); }
ImVec2 operator*(ImVec2 a, float s) { return ImVec2(a.x * s, a.y * s); }
float len(ImVec2 a) { return std::sqrt(a.x * a.x + a.y * a.y); }
ImVec2 norm(ImVec2 a) {
    const float l = len(a);
    return l > 1e-6f ? a * (1.0f / l) : ImVec2(0, 0);
}

ImU32 withAlpha(ImU32 c, int a) { return (c & 0x00FFFFFF) | ((ImU32)a << 24); }

std::string fmtVal(double v) {
    char b[32];
    const double a = std::abs(v);
    if (a >= 1000) std::snprintf(b, sizeof b, "%.0f", v);
    else if (a >= 100) std::snprintf(b, sizeof b, "%.1f", v);
    else if (a >= 1) std::snprintf(b, sizeof b, "%.2f", v);
    else if (a >= 0.01) std::snprintf(b, sizeof b, "%.3f", v);
    else if (a == 0) std::snprintf(b, sizeof b, "0");
    else std::snprintf(b, sizeof b, "%.2e", v);
    return b;
}

// Nhãn chữ có nền mờ; align: (0..1, 0..1) điểm neo trong hộp chữ
void label(ImDrawList* dl, ImVec2 pos, const std::string& s, ImU32 col, ImVec2 align = ImVec2(0.5f, 0.5f),
           float scale = 0.9f) {
    ImFont* f = ImGui::GetFont();
    const float fs = ImGui::GetFontSize() * scale;
    const ImVec2 ts = f->CalcTextSizeA(fs, FLT_MAX, 0, s.c_str());
    const ImVec2 p(std::floor(pos.x - ts.x * align.x), std::floor(pos.y - ts.y * align.y));
    dl->AddRectFilled(p - ImVec2(2, 0), p + ts + ImVec2(2, 0), withAlpha(P.bg, 190), 2);
    dl->AddText(f, fs, p, col, s.c_str());
}

void arrow(ImDrawList* dl, ImVec2 from, ImVec2 to, ImU32 col, float th = 1.6f, float head = 8) {
    const ImVec2 d = norm(to - from);
    if (d.x == 0 && d.y == 0) return;
    const ImVec2 n(-d.y, d.x);
    dl->AddLine(from, to - d * (head * 0.6f), col, th);
    dl->AddTriangleFilled(to, to - d * head + n * (head * 0.38f), to - d * head - n * (head * 0.38f), col);
}

// Cung mômen quanh điểm c. ccw: chiều ngược kim đồng hồ (trong hệ tọa độ thực).
void momentArc(ImDrawList* dl, ImVec2 c, float r, bool ccw, ImU32 col, float th = 1.6f) {
    const int n = 24;
    const float a0 = -0.75f * kPi, a1 = 0.75f * kPi;
    ImVec2 pts[25];
    for (int i = 0; i <= n; ++i) {
        const float a = a0 + (a1 - a0) * i / n;
        pts[i] = c + ImVec2(std::cos(a) * r, -std::sin(a) * r);
    }
    dl->AddPolyline(pts, n + 1, col, th);
    const ImVec2 tip = ccw ? pts[n] : pts[0];
    const ImVec2 prev = ccw ? pts[n - 3] : pts[3];
    arrow(dl, prev, tip, col, th, 8);
}

void colormap(double t, float rgb[3]) {
    static const float stops[9][3] = {{0.00f, 0.00f, 0.55f}, {0.00f, 0.15f, 1.00f}, {0.00f, 0.60f, 1.00f},
                                      {0.00f, 0.90f, 0.85f}, {0.30f, 0.90f, 0.30f}, {0.85f, 0.95f, 0.10f},
                                      {1.00f, 0.70f, 0.00f}, {1.00f, 0.30f, 0.00f}, {0.75f, 0.00f, 0.00f}};
    t = std::clamp(t, 0.0, 1.0) * 8;
    const int i = std::min(7, (int)t);
    const float f = (float)(t - i);
    for (int k = 0; k < 3; ++k) rgb[k] = stops[i][k] + (stops[i + 1][k] - stops[i][k]) * f;
}

ImU32 colormapU32(double t) {
    float c[3];
    colormap(t, c);
    return IM_COL32((int)(c[0] * 255), (int)(c[1] * 255), (int)(c[2] * 255), 255);
}

// Đường biên phần tử phẳng trên màn hình: p[] là vị trí các nút; cạnh bậc hai được chia nhỏ để vẽ cong.
int planeOutline(const Element& el, const ImVec2* p, ImVec2 out[24]) {
    const int nn = el.nodeCount(), nc = el.cornerCount();
    const int seg = nn > nc ? 6 : 1;
    int n = 0;
    for (int k = 0; k < nc; ++k) {
        double x0, y0, x1, y1;
        planeNodeNatural(el.type, k, x0, y0);
        planeNodeNatural(el.type, (k + 1) % nc, x1, y1);
        for (int s = 0; s < seg; ++s) {
            if (seg == 1) {
                out[n++] = p[k];
                continue;
            }
            double N[8], a[8], b[8];
            planeShape(el.type, x0 + (x1 - x0) * s / seg, y0 + (y1 - y0) * s / seg, N, a, b);
            ImVec2 q(0, 0);
            for (int i = 0; i < nn; ++i) q = q + p[i] * (float)N[i];
            out[n++] = q;
        }
    }
    return n;
}

bool g_axisym = false; // tên đại lượng theo bài toán đối xứng trục (r, z, θ)

const char* fieldName(ContourField f) {
    switch (f) {
    case ContourField::Ux: return g_axisym ? "ur (mm)" : "ux (mm)";
    case ContourField::Uy: return g_axisym ? "uz (mm)" : "uy (mm)";
    case ContourField::Umag: return "|u| (mm)";
    case ContourField::Sx: return g_axisym ? "σr (MPa)" : "σx (MPa)";
    case ContourField::Sy: return g_axisym ? "σz (MPa)" : "σy (MPa)";
    case ContourField::Txy: return g_axisym ? "τrz (MPa)" : "τxy (MPa)";
    case ContourField::Sz: return g_axisym ? "σθ (MPa)" : "σz ngoài mặt phẳng (MPa)";
    case ContourField::S1: return "σ1 (MPa)";
    case ContourField::S2: return "σ2 (MPa)";
    case ContourField::VonMises: return "σ von Mises (MPa)";
    }
    return "";
}

bool isDispField(ContourField f) { return f == ContourField::Ux || f == ContourField::Uy || f == ContourField::Umag; }

double stressField(const std::array<double, 4>& s, ContourField f) {
    const StressInfo si = stressInfo(s[0], s[1], s[2], s[3]);
    switch (f) {
    case ContourField::Sx: return si.sx * 1e-3;
    case ContourField::Sy: return si.sy * 1e-3;
    case ContourField::Txy: return si.txy * 1e-3;
    case ContourField::Sz: return si.sz * 1e-3;
    case ContourField::S1: return si.s1 * 1e-3;
    case ContourField::S2: return si.s2 * 1e-3;
    case ContourField::VonMises: return si.vonMises * 1e-3;
    default: return 0;
    }
}

double dispField(const LoadSetResult& s, int node, ContourField f) {
    const double ux = s.U[3 * node] * 1e3, uy = s.U[3 * node + 1] * 1e3;
    return f == ContourField::Ux ? ux : f == ContourField::Uy ? uy : std::hypot(ux, uy);
}

// Trạng thái tương tác chuột của canvas
struct Interaction {
    bool leftDown = false, boxing = false, rightDown = false;
    ImVec2 start{0, 0};
} g_in;

} // namespace

const char* contourFieldName(ContourField f, bool axisym) {
    g_axisym = axisym;
    return fieldName(f);
}

void App::viewport() {
    P = pal;
    {
        ImGuiWindowClass wc;
        wc.DockNodeFlagsOverrideSet = ImGuiDockNodeFlags_AutoHideTabBar;
        ImGui::SetNextWindowClass(&wc);
    }
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("Sơ đồ kết cấu", nullptr, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    ImGui::PopStyleVar();

    // ---------------------------------------------------------------- canvas
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    ImVec2 sz = ImGui::GetContentRegionAvail();
    sz.x = std::max(sz.x, 50.0f);
    sz.y = std::max(sz.y, 50.0f);
    ImGui::InvisibleButton("canvas", sz,
                           ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight |
                               ImGuiButtonFlags_MouseButtonMiddle);
    const bool hovered = ImGui::IsItemHovered();
    cam.c0 = p0;
    cam.c1 = p0 + sz;
    ImGuiIO& io = ImGui::GetIO();
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->PushClipRect(cam.c0, cam.c1, true);
    dl->AddRectFilled(cam.c0, cam.c1, P.bg);

    const Model& m = model;
    double bx0, by0, bx1, by1;
    m.bounds(bx0, by0, bx1, by1);
    const double modelSize = std::max({bx1 - bx0, by1 - by0, 1.0});
    if (fitRequested) {
        fitRequested = false;
        fitFrames = 4;
    }
    if (fitFrames > 0) {
        --fitFrames;
        const double w = std::max(bx1 - bx0, modelSize * 0.05), h = std::max(by1 - by0, modelSize * 0.05);
        cam.cx = 0.5 * (bx0 + bx1);
        cam.cy = 0.5 * (by0 + by1);
        cam.scale = std::min(sz.x / (w * 1.35), sz.y / (h * 1.5));
    }

    // ---- lưới
    if (showGrid) {
        double step = gridStep;
        while (step * cam.scale < 10) step *= 5;
        double wx0, wy0, wx1, wy1;
        cam.toWorld(cam.c0, wx0, wy1);
        cam.toWorld(cam.c1, wx1, wy0);
        const long i0 = (long)std::floor(wx0 / step), i1 = (long)std::ceil(wx1 / step);
        const long j0 = (long)std::floor(wy0 / step), j1 = (long)std::ceil(wy1 / step);
        if (i1 - i0 < 600 && j1 - j0 < 600) {
            for (long i = i0; i <= i1; ++i) {
                const float x = cam.toScreen(i * step, 0).x;
                dl->AddLine(ImVec2(x, cam.c0.y), ImVec2(x, cam.c1.y), i == 0 ? P.axis : (i % 5 ? P.grid : P.gridMajor));
            }
            for (long j = j0; j <= j1; ++j) {
                const float y = cam.toScreen(0, j * step).y;
                dl->AddLine(ImVec2(cam.c0.x, y), ImVec2(cam.c1.x, y), j == 0 ? P.axis : (j % 5 ? P.grid : P.gridMajor));
            }
        }
    }

    const LoadSetResult* rs = currentSet();
    std::string fmTitle;
    if (const LoadSetResult* o = overrideShown(&fmTitle)) rs = o; // biểu đồ phương pháp lực / động lực học
    const bool showRes = rs != nullptr && view != ResultView::None;

    // ---- hệ số phóng đại biến dạng
    double maxU = 0;
    if (rs) {
        for (size_t i = 0; i < m.nodes.size(); ++i) maxU = std::max(maxU, std::hypot(rs->U[3 * i], rs->U[3 * i + 1]));
        for (size_t e = 0; e < m.elements.size(); ++e)
            for (size_t k = 0; k < rs->bars[e].u.size(); ++k)
                maxU = std::max(maxU, std::hypot(rs->bars[e].u[k], rs->bars[e].v[k]));
    }
    deformAuto = maxU > 0 ? 0.08 * modelSize / maxU : 1.0;
    const double S = deformAuto * deformFactor;
    const bool drawDeformed = rs && (view == ResultView::Deformed || (view == ResultView::Contour && contourOnDeformed) ||
                                     (showDeformedWithDiagram && view != ResultView::None));
    auto nodePos = [&](int i, bool deformed) {
        double x = m.nodes[i].x, y = m.nodes[i].y;
        if (deformed && rs) {
            x += rs->U[3 * i] * S;
            y += rs->U[3 * i + 1] * S;
        }
        return cam.toScreen(x, y);
    };

    // ---- phần tử phẳng: tô màu hoặc nền
    const bool contourOn = showRes && view == ResultView::Contour && m.hasPlane();
    double cmin = 0, cmax = 0;
    if (contourOn) {
        cmin = DBL_MAX;
        cmax = -DBL_MAX;
        for (size_t e = 0; e < m.elements.size(); ++e) {
            const Element& el = m.elements[e];
            if (!el.isPlane()) continue;
            for (int k = 0; k < el.nodeCount(); ++k) {
                const double v = isDispField(contour) ? dispField(*rs, el.n[k], contour)
                                 : contourAveraged        ? stressField(rs->nodeStress[el.n[k]], contour)
                                                          : stressField(rs->elemStress[e][k], contour);
                cmin = std::min(cmin, v);
                cmax = std::max(cmax, v);
            }
        }
        if (cmax - cmin < 1e-12 * std::max(1.0, std::abs(cmax))) {
            cmax += 1e-9;
            cmin -= 1e-9;
        }
    }
    const ImVec2 uvw = ImGui::GetFontTexUvWhitePixel();
    for (size_t e = 0; e < m.elements.size(); ++e) {
        const Element& el = m.elements[e];
        if (!el.isPlane()) continue;
        const int nc = el.nodeCount();
        const bool def = contourOn ? contourOnDeformed : false;
        ImVec2 pts[8], ol[24];
        for (int k = 0; k < nc; ++k) pts[k] = nodePos(el.n[k], def);
        const int no = planeOutline(el, pts, ol);
        if (contourOn) {
            double val[8];
            for (int k = 0; k < nc; ++k)
                val[k] = isDispField(contour) ? dispField(*rs, el.n[k], contour)
                         : contourAveraged     ? stressField(rs->nodeStress[el.n[k]], contour)
                                               : stressField(rs->elemStress[e][k], contour);
            const bool tri = isTriangle(el.type);
            const int n = nc > el.cornerCount() ? 8 : 6;
            auto vtx = [&](double u, double v, ImVec2& p, ImU32& c) {
                double N[8], a[8], b[8], val2 = 0;
                if (tri) planeShape(el.type, u, v, N, a, b);
                else planeShape(el.type, 2 * u - 1, 2 * v - 1, N, a, b);
                p = ImVec2(0, 0);
                for (int k = 0; k < nc; ++k) {
                    p = p + pts[k] * (float)N[k];
                    val2 += N[k] * val[k];
                }
                c = colormapU32((val2 - cmin) / (cmax - cmin));
            };
            auto tri3 = [&](double u0, double v0, double u1, double v1, double u2, double v2) {
                ImVec2 a, b, c;
                ImU32 ca, cb, cc;
                vtx(u0, v0, a, ca);
                vtx(u1, v1, b, cb);
                vtx(u2, v2, c, cc);
                dl->PrimReserve(3, 3);
                dl->PrimVtx(a, uvw, ca);
                dl->PrimVtx(b, uvw, cb);
                dl->PrimVtx(c, uvw, cc);
            };
            const double h = 1.0 / n;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) {
                    if (!tri) {
                        tri3(i * h, j * h, (i + 1) * h, j * h, (i + 1) * h, (j + 1) * h);
                        tri3(i * h, j * h, (i + 1) * h, (j + 1) * h, i * h, (j + 1) * h);
                    } else if (i + j < n) {
                        tri3(i * h, j * h, (i + 1) * h, j * h, i * h, (j + 1) * h);
                        if (i + j + 1 < n) tri3((i + 1) * h, j * h, (i + 1) * h, (j + 1) * h, i * h, (j + 1) * h);
                    }
                }
            dl->AddPolyline(ol, no, withAlpha(P.planeEdge, 90), 1.0f, ImDrawFlags_Closed);
        } else {
            if (no == nc) dl->AddConvexPolyFilled(ol, no, withAlpha(P.plane, drawDeformed ? 70 : 150));
            else dl->AddConcavePolyFilled(ol, no, withAlpha(P.plane, drawDeformed ? 70 : 150));
            dl->AddPolyline(ol, no, withAlpha(P.planeEdge, drawDeformed ? 110 : 255), 1.0f, ImDrawFlags_Closed);
        }
    }
    // trục đối xứng (bài toán đối xứng trục: r = x, trục z ≡ x = 0)
    if (m.planeMode == PlaneMode::Axisymmetric && m.hasPlane()) {
        const ImVec2 a = cam.toScreen(0, 0);
        const float y0 = cam.c0.y, y1 = cam.c1.y;
        for (float y = y0; y < y1; y += 18) dl->AddLine(ImVec2(a.x, y), ImVec2(a.x, std::min(y + 10, y1)), withAlpha(P.textDim, 200), 1.2f);
        for (float y = y0 + 13; y < y1; y += 18) dl->AddLine(ImVec2(a.x, y), ImVec2(a.x, std::min(y + 2, y1)), withAlpha(P.textDim, 200), 1.2f);
        label(dl, ImVec2(a.x + 4, y0 + 6), "z (trục đối xứng)", P.textDim, ImVec2(0, 0), 0.8f);
    }

    // ---- thanh
    for (size_t e = 0; e < m.elements.size(); ++e) {
        const Element& el = m.elements[e];
        if (!el.isBar()) continue;
        const ImVec2 a = nodePos(el.n[0], false), b = nodePos(el.n[1], false);
        const bool truss = el.type == ElemType::Truss;
        ImU32 col = truss ? P.truss : P.bar;
        if (drawDeformed && !showUndeformed) continue;
        if (drawDeformed || (showRes && view != ResultView::Contour)) col = withAlpha(col, 120);
        dl->AddLine(a, b, col, truss ? 1.6f : 2.6f);
        if (el.kf > 0) {
            // ký hiệu nền đàn hồi: các vạch ngắn phía dưới
            const ImVec2 d = norm(b - a), nrm(d.y, -d.x);
            const ImVec2 down = nrm.y > 0 ? nrm : nrm * -1.0f;
            for (int k = 1; k < 4; ++k) {
                const ImVec2 p = a + (b - a) * (k / 4.0f);
                dl->AddLine(p + down * 3, p + down * 9, withAlpha(P.support, 160), 1.0f);
            }
        }
    }
    // hệ số sử dụng (cửa sổ kiểm tra thiết kế)
    if (showDesign && showUtil && util.size() == m.elements.size()) {
        for (size_t e = 0; e < m.elements.size(); ++e) {
            if (util[e] < 0 || !m.elements[e].isBar()) continue;
            const ImVec2 a = nodePos(m.elements[e].n[0], false), b = nodePos(m.elements[e].n[1], false);
            const double u = util[e];
            const ImU32 c = u <= 0.7 ? IM_COL32(60, 190, 100, 255) : u <= 1.0 ? IM_COL32(240, 180, 40, 255) : IM_COL32(240, 70, 60, 255);
            dl->AddLine(a, b, c, 6.0f);
            if (showValues) {
                char buf[16];
                std::snprintf(buf, sizeof buf, "%.2f", u);
                label(dl, (a + b) * 0.5f, buf, c);
            }
        }
    }
    // ---- hình biến dạng
    if (drawDeformed) {
        for (size_t e = 0; e < m.elements.size(); ++e) {
            const Element& el = m.elements[e];
            if (el.isBar()) {
                const BarResult& br = rs->bars[e];
                const Node& na = m.nodes[el.n[0]];
                const Node& nb = m.nodes[el.n[1]];
                const double L = std::hypot(nb.x - na.x, nb.y - na.y), c = (nb.x - na.x) / L, s = (nb.y - na.y) / L;
                std::vector<ImVec2> pts;
                for (size_t k = 0; k < br.x.size(); ++k) {
                    const double x = na.x + c * br.x[k] + (c * br.u[k] - s * br.v[k]) * S;
                    const double y = na.y + s * br.x[k] + (s * br.u[k] + c * br.v[k]) * S;
                    pts.push_back(cam.toScreen(x, y));
                }
                dl->AddPolyline(pts.data(), (int)pts.size(), P.deformed, 2.0f);
            } else if (!contourOn) {
                ImVec2 pts[8], ol[24];
                for (int k = 0; k < el.nodeCount(); ++k) pts[k] = nodePos(el.n[k], true);
                dl->AddPolyline(ol, planeOutline(el, pts, ol), P.deformed, 1.2f, ImDrawFlags_Closed);
            }
        }
        // chuyển vị lớn nhất
        int imax = -1;
        double best = 0;
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            const double u = std::hypot(rs->U[3 * i], rs->U[3 * i + 1]);
            if (u > best) {
                best = u;
                imax = (int)i;
            }
        }
        if (imax >= 0 && showValues) {
            char b[96];
            std::snprintf(b, sizeof b, "δmax = %.3f mm (nút %d)", best * 1e3, imax + 1);
            const ImVec2 p = nodePos(imax, true);
            dl->AddCircle(p, 5, P.deformed, 0, 2);
            label(dl, p + ImVec2(8, -10), b, P.deformed, ImVec2(0, 1));
        }
    }

    // ---- vẽ một đường biểu đồ dọc theo dãy điểm gốc (kiểu giáo trình: tô, gạch tung độ, viền, giá trị)
    // xs: tọa độ dọc thanh để nhận biết bước nhảy (hai điểm cùng x); sign: −1 với M (vẽ về phía thớ căng)
    auto drawCurve = [&](const std::vector<ImVec2>& base, ImVec2 nrm, const std::vector<double>& xs,
                         const std::vector<double>& vals, float sign, double kpx, double vmax, bool labels,
                         int fillAlpha) {
        const size_t n = base.size();
        if (n < 2) return;
        const double tiny = std::max(vmax * 1e-6, 1e-12);
        const ImU32 pf = withAlpha(P.posFill, fillAlpha), nf = withAlpha(P.negFill, fillAlpha);
        std::vector<ImVec2> off(n);
        for (size_t k = 0; k < n; ++k) off[k] = base[k] + nrm * (float)(sign * vals[k] * kpx);
        for (size_t k = 0; k + 1 < n; ++k) {
            const double v0 = vals[k], v1 = vals[k + 1];
            if ((v0 >= 0) == (v1 >= 0) || std::abs(v0) < tiny || std::abs(v1) < tiny) {
                const bool pos = (std::abs(v0) > std::abs(v1) ? v0 : v1) >= 0;
                dl->AddQuadFilled(base[k], off[k], off[k + 1], base[k + 1], pos ? pf : nf);
            } else {
                const float t = (float)(v0 / (v0 - v1));
                const ImVec2 z = base[k] + (base[k + 1] - base[k]) * t;
                dl->AddTriangleFilled(base[k], off[k], z, v0 >= 0 ? pf : nf);
                dl->AddTriangleFilled(z, off[k + 1], base[k + 1], v1 >= 0 ? pf : nf);
            }
        }
        const float spacing = 9.0f;
        float acc = spacing;
        for (size_t k = 0; k < n; ++k) {
            if (k > 0) acc += len(base[k] - base[k - 1]);
            if (acc >= spacing || k == 0 || k + 1 == n) {
                acc = 0;
                dl->AddLine(base[k], off[k], withAlpha(vals[k] >= 0 ? P.posLine : P.negLine, 110), 1.0f);
            }
        }
        for (size_t k = 0; k + 1 < n; ++k)
            dl->AddLine(off[k], off[k + 1], vals[k] + vals[k + 1] >= 0 ? P.posLine : P.negLine, 1.8f);
        if (!labels || !showValues) return;
        auto put = [&](size_t k) {
            if (std::abs(vals[k]) < vmax * 1e-4) return;
            const ImVec2 dir = nrm * (float)(sign * (vals[k] >= 0 ? 1 : -1));
            label(dl, off[k] + dir * 9, fmtVal(vals[k]), vals[k] >= 0 ? P.posLine : P.negLine);
        };
        put(0);
        put(n - 1);
        size_t kx = 0;
        for (size_t k = 1; k + 1 < n; ++k)
            if (std::abs(vals[k]) > std::abs(vals[kx])) kx = k;
        if (kx > 0 && kx + 1 < n && std::abs(vals[kx]) > 1.02 * std::max(std::abs(vals[0]), std::abs(vals[n - 1])) &&
            len(base[kx] - base[0]) > 30 && len(base[kx] - base[n - 1]) > 30)
            put(kx);
        for (size_t k = 1; k + 2 < n; ++k)
            if (xs[k] == xs[k + 1] && std::abs(vals[k] - vals[k + 1]) > vmax * 0.02) {
                put(k);
                put(k + 1);
            }
    };
    auto barFrame = [&](const Element& el, double& c, double& s) {
        const Node& na = m.nodes[el.n[0]];
        const Node& nb = m.nodes[el.n[1]];
        const double L = std::hypot(nb.x - na.x, nb.y - na.y);
        c = (nb.x - na.x) / L;
        s = (nb.y - na.y) / L;
    };
    auto basePoints = [&](const Element& el, const std::vector<double>& xs) {
        double c, s;
        barFrame(el, c, s);
        const Node& na = m.nodes[el.n[0]];
        std::vector<ImVec2> b(xs.size());
        for (size_t k = 0; k < xs.size(); ++k) b[k] = cam.toScreen(na.x + c * xs[k], na.y + s * xs[k]);
        return b;
    };

    // ---- biểu đồ nội lực
    const bool diagram = showRes && (view == ResultView::N || view == ResultView::Q || view == ResultView::M ||
                                     view == ResultView::Sigma);
    const bool envMode = diagram && showEnvelope && view != ResultView::Sigma;
    if (envMode) updateEnvelope();
    if (diagram) {
        const int qi = view == ResultView::N ? 0 : view == ResultView::Q ? 1 : 2;
        auto valuesOf = [&](size_t e, std::vector<double>& out) {
            const BarResult& br = rs->bars[e];
            switch (view) {
            case ResultView::N: out = br.N; break;
            case ResultView::Q: out = br.Q; break;
            case ResultView::M: out = br.M; break;
            default: {
                const Element& el = m.elements[e];
                const Section& sec = m.sections[el.sec];
                out.resize(br.x.size());
                for (size_t k = 0; k < br.x.size(); ++k) {
                    const double sN = br.N[k] / sec.A;
                    const bool bend = sec.I > 0 && el.type == ElemType::Frame;
                    const double top = sN - (bend ? br.M[k] / sec.I * sec.top() : 0.0);
                    const double bot = sN + (bend ? br.M[k] / sec.I * sec.bot() : 0.0);
                    out[k] = (std::abs(top) > std::abs(bot) ? top : bot) * 1e-3; // MPa
                }
            }
            }
        };
        double vmax = 0;
        std::vector<double> vals;
        for (size_t e = 0; e < m.elements.size(); ++e) {
            if (!m.elements[e].isBar()) continue;
            if (envMode) {
                if (env.size() != m.elements.size()) continue;
                for (double v : env[e].max[qi]) vmax = std::max(vmax, std::abs(v));
                for (double v : env[e].min[qi]) vmax = std::max(vmax, std::abs(v));
            } else {
                valuesOf(e, vals);
                for (double v : vals) vmax = std::max(vmax, std::abs(v));
            }
        }
        if (vmax > 1e-9) {
            const double kpx = 55.0 * diagramFactor / vmax;
            const float sign = view == ResultView::M ? -1.0f : 1.0f;
            for (size_t e = 0; e < m.elements.size(); ++e) {
                const Element& el = m.elements[e];
                if (!el.isBar()) continue;
                double c, s;
                barFrame(el, c, s);
                const ImVec2 nrm((float)-s, (float)-c); // trục y cục bộ trên màn hình
                if (envMode) {
                    if (env.size() != m.elements.size() || env[e].x.empty()) continue;
                    const auto base = basePoints(el, env[e].x);
                    if (el.type == ElemType::Truss) {
                        const ImVec2 a = nodePos(el.n[0], false), b = nodePos(el.n[1], false);
                        dl->AddLine(a, b, withAlpha(P.textDim, 160), 2.0f);
                        if (showValues)
                            label(dl, (a + b) * 0.5f, fmtVal(env[e].max[qi][0]) + " / " + fmtVal(env[e].min[qi][0]), P.text);
                        continue;
                    }
                    drawCurve(base, nrm, env[e].x, env[e].max[qi], sign, kpx, vmax, true, 55);
                    drawCurve(base, nrm, env[e].x, env[e].min[qi], sign, kpx, vmax, true, 55);
                    continue;
                }
                const BarResult& br = rs->bars[e];
                valuesOf(e, vals);
                if (el.type == ElemType::Truss && (view == ResultView::N || view == ResultView::Sigma)) {
                    // thanh dàn: lực dọc không đổi → tô màu thanh (xanh: kéo, đỏ: nén) và ghi giá trị
                    bool constant = true;
                    for (double x : vals) constant &= std::abs(x - vals[0]) < 1e-6 * vmax;
                    if (constant) {
                        const ImVec2 a = nodePos(el.n[0], false), b = nodePos(el.n[1], false);
                        const double v = vals[0];
                        const bool zero = std::abs(v) < vmax * 1e-6;
                        dl->AddLine(a, b, zero ? withAlpha(P.textDim, 160) : v > 0 ? P.posLine : P.negLine,
                                    zero ? 1.5f : 2.0f + 4.0f * (float)(std::abs(v) / vmax) * diagramFactor);
                        if (showValues) label(dl, (a + b) * 0.5f, zero ? "0" : fmtVal(v), zero ? P.textDim : v > 0 ? P.posLine : P.negLine);
                        continue;
                    }
                }
                drawCurve(basePoints(el, br.x), nrm, br.x, vals, sign, kpx, vmax, true, (int)(P.posFill >> 24));
            }
        }
        if (envMode && showValues) {
            const std::string t = std::string("Biểu đồ bao qua ") + std::to_string(envSets.size()) +
                                  (result.sets.size() > model.cases.size() ? " tổ hợp" : " trường hợp tải");
            label(dl, ImVec2(cam.c0.x + 10, cam.c0.y + 30), t, P.text, ImVec2(0, 0.5f));
        }
    }

    // ---- đường ảnh hưởng
    if (view == ResultView::Influence) {
        updateInfluence();
        if (il.ok) {
            const double vmax = std::max(std::abs(il.vmax), std::abs(il.vmin));
            const double kpx = vmax > 0 ? 70.0 * diagramFactor / vmax : 0;
            size_t i = 0;
            while (i < il.pts.size()) {
                size_t j = i;
                while (j < il.pts.size() && il.pts[j].elem == il.pts[i].elem) ++j;
                const Element& el = m.elements[il.pts[i].elem];
                double c, s;
                barFrame(el, c, s);
                // tung độ dương vẽ lên trên
                ImVec2 nrm((float)-s, (float)-c);
                if (c < -1e-9) nrm = nrm * -1.0f;
                std::vector<ImVec2> base;
                std::vector<double> xs, vals;
                for (size_t k = i; k < j; ++k) {
                    base.push_back(cam.toScreen(il.pts[k].X, il.pts[k].Y));
                    xs.push_back(il.pts[k].a);
                    vals.push_back(il.pts[k].value);
                }
                if (vmax > 0) drawCurve(base, nrm, xs, vals, 1.0f, kpx, vmax, true, (int)(P.posFill >> 24));
                i = j;
            }
            // đoàn tải ở vị trí bất lợi nhất
            if (trainShow) {
                const LoadTrain& tr = currentTrain();
                const TrainResult t = moveTrain(il, tr);
                auto pathY = [&](double X) {
                    double best = 1e300, y = 0;
                    for (const auto& p : il.pts)
                        if (std::abs(p.X - X) < best) {
                            best = std::abs(p.X - X);
                            y = p.Y;
                        }
                    return y;
                };
                for (int which = 0; which < 2 && t.ok; ++which) {
                    const double x0 = which == 0 ? t.xMax : t.xMin;
                    const bool rev = which == 0 ? t.revMax : t.revMin;
                    const ImU32 col = which == 0 ? P.posLine : P.negLine;
                    for (const auto& a : tr.axles) {
                        if (a.P == 0) continue;
                        const double X = rev ? x0 - a.offset : x0 + a.offset;
                        const ImVec2 p = cam.toScreen(X, pathY(X));
                        const float lift = which == 0 ? 0.0f : 0.0f;
                        arrow(dl, p - ImVec2(0, 46 + lift), p - ImVec2(0, 4), withAlpha(col, which == 0 ? 255 : 150), 2.0f, 8);
                        if (showValues) label(dl, p - ImVec2(0, 54), fmtVal(a.P), col, ImVec2(0.5f, 1), 0.8f);
                    }
                }
                if (t.ok) {
                    char b[160];
                    std::snprintf(b, sizeof b, "%s:  Smax = %s,  Smin = %s", tr.name.c_str(), fmtVal(t.Smax).c_str(), fmtVal(t.Smin).c_str());
                    label(dl, ImVec2(cam.c0.x + 10, cam.c0.y + 50), b, P.text, ImVec2(0, 0.5f));
                }
            }
            // đánh dấu đại lượng khảo sát
            if (ilTarget.kind == IlKind::Force && ilTarget.elem < (int)m.elements.size()) {
                const Element& el = m.elements[ilTarget.elem];
                double c, s;
                barFrame(el, c, s);
                const double x = std::clamp(ilTarget.x, 0.0, m.barLength(ilTarget.elem));
                const ImVec2 p = cam.toScreen(m.nodes[el.n[0]].x + c * x, m.nodes[el.n[0]].y + s * x);
                const ImVec2 nrm((float)-s, (float)-c);
                dl->AddLine(p - nrm * 12, p + nrm * 12, P.select, 3.0f);
                label(dl, p - nrm * 20, "k", P.select);
            } else if (ilTarget.node < (int)m.nodes.size()) {
                dl->AddCircle(nodePos(ilTarget.node, false), 9, P.select, 0, 2.5f);
            }
            char b[200];
            std::snprintf(b, sizeof b, "%s   |   max %s, min %s   |   Ω+ = %s, Ω− = %s", ilTargetName(m, ilTarget).c_str(),
                          fmtVal(il.vmax).c_str(), fmtVal(il.vmin).c_str(), fmtVal(il.areaPos).c_str(),
                          fmtVal(il.areaNeg).c_str());
            label(dl, ImVec2(cam.c0.x + 10, cam.c0.y + 30), b, P.text, ImVec2(0, 0.5f));
        } else {
            label(dl, ImVec2(cam.c0.x + 10, cam.c0.y + 30), "Đường ảnh hưởng: " + il.error, IM_COL32(200, 30, 30, 255),
                  ImVec2(0, 0.5f));
        }
    }

    // ---- dạng riêng (ổn định / dao động)
    if (view == ResultView::Mode) {
        const bool fresh = eig.ok && eigVersion == modelVersion && eigMode >= 0 && eigMode < (int)eig.modes.size();
        if (fresh) {
            const EigenMode& md = eig.modes[eigMode];
            const Model& sm = eig.sub;
            double amp = 0.1 * modelSize * deformFactor;
            if (animateMode) amp *= std::sin(ImGui::GetTime() * 2 * kPi * 0.6);
            auto pos = [&](int n) {
                return cam.toScreen(sm.nodes[n].x + md.U[3 * n] * amp, sm.nodes[n].y + md.U[3 * n + 1] * amp);
            };
            for (size_t e = 0; e < m.elements.size(); ++e) {
                const Element& el = m.elements[e];
                if (el.isBar()) {
                    if (e >= eig.barChain.size()) continue;
                    std::vector<ImVec2> pts;
                    for (int n : eig.barChain[e]) pts.push_back(pos(n));
                    dl->AddPolyline(pts.data(), (int)pts.size(), P.deformed, 2.2f);
                } else {
                    ImVec2 pts[8], ol[24];
                    for (int k = 0; k < el.nodeCount(); ++k) pts[k] = pos(el.n[k]);
                    dl->AddPolyline(ol, planeOutline(el, pts, ol), P.deformed, 1.2f, ImDrawFlags_Closed);
                }
            }
            char b[200];
            if (eig.kind == EigenKind::Mechanism)
                std::snprintf(b, sizeof b, "Dạng biến hình (cơ cấu): các thanh chuyển động như vật rắn, không biến dạng");
            else if (eig.kind == EigenKind::Buckling)
                std::snprintf(b, sizeof b, "Dạng mất ổn định %d:  λcr = %.4g  (tải trọng tới hạn = λcr × bộ tải \"%s\")",
                              eigMode + 1, md.lambda,
                              eigOpt.loadSet < (int)result.sets.size() ? result.sets[eigOpt.loadSet].name.c_str() : "");
            else
                std::snprintf(b, sizeof b, "Dạng dao động %d:  ω = %.4g rad/s,  f = %.4g Hz,  T = %.4g s", eigMode + 1,
                              md.omega, md.freq, md.period);
            label(dl, ImVec2(cam.c0.x + 10, cam.c0.y + 30), b, P.text, ImVec2(0, 0.5f));
        } else {
            label(dl, ImVec2(cam.c0.x + 10, cam.c0.y + 30),
                  eig.ok || eigVersion < 0 ? "Chưa có kết quả trị riêng — mở Tính toán → Ổn định & dao động riêng"
                                           : "Lỗi trị riêng: " + eig.error,
                  IM_COL32(200, 120, 0, 255), ImVec2(0, 0.5f));
        }
    }

    // ---- gối tựa
    auto connectedDir = [&](int node) {
        ImVec2 d(0, 0);
        for (const auto& el : m.elements) {
            if (!el.isBar()) continue;
            if (el.n[0] == node) d = d + norm(nodePos(el.n[1], false) - nodePos(node, false));
            if (el.n[1] == node) d = d + norm(nodePos(el.n[0], false) - nodePos(node, false));
        }
        return d;
    };
    if (showSupports) {
        const float s = 11;
        for (const auto& sp : m.supports) {
            if (sp.node < 0 || sp.node >= (int)m.nodes.size()) continue;
            const ImVec2 p = nodePos(sp.node, false);
            const ImU32 c = P.support;
            const bool fx = sp.fix[0], fy = sp.fix[1], fr = sp.fix[2];
            auto hatch = [&](ImVec2 a, ImVec2 b, ImVec2 out) {
                dl->AddLine(a, b, c, 2.0f);
                const int n = 5;
                for (int i = 0; i <= n; ++i) {
                    const ImVec2 q = a + (b - a) * (i / (float)n);
                    dl->AddLine(q, q + out * 6 + norm(b - a) * -5, c, 1.0f);
                }
            };
            if (fx && fy && fr) {
                // ngàm: mặt ngàm vuông góc với hướng thanh
                ImVec2 d = connectedDir(sp.node);
                if (len(d) < 0.1f) d = ImVec2(0, -1);
                d = norm(d);
                ImVec2 out = std::abs(d.x) > std::abs(d.y) ? ImVec2(d.x > 0 ? -1.f : 1.f, 0) : ImVec2(0, d.y > 0 ? -1.f : 1.f);
                const ImVec2 t(-out.y, out.x);
                hatch(p - t * (s * 1.3f), p + t * (s * 1.3f), out);
            } else if (fx && fy && !fr) {
                const ImVec2 l = p + ImVec2(-s * 0.8f, s * 1.4f), r = p + ImVec2(s * 0.8f, s * 1.4f);
                dl->AddTriangle(p, l, r, c, 1.6f);
                hatch(l - ImVec2(4, 0), r + ImVec2(4, 0), ImVec2(0, 1));
            } else if (!fx && fy && !fr) {
                const ImVec2 l = p + ImVec2(-s * 0.8f, s * 1.3f), r = p + ImVec2(s * 0.8f, s * 1.3f);
                dl->AddTriangle(p, l, r, c, 1.6f);
                dl->AddCircle(l + ImVec2(4, 3.5f), 3, c, 0, 1.2f);
                dl->AddCircle(r + ImVec2(-4, 3.5f), 3, c, 0, 1.2f);
                hatch(l + ImVec2(-4, 7), r + ImVec2(4, 7), ImVec2(0, 1));
            } else if (fx && !fy && !fr) {
                const ImVec2 u = p + ImVec2(-s * 1.3f, -s * 0.8f), d = p + ImVec2(-s * 1.3f, s * 0.8f);
                dl->AddTriangle(p, u, d, c, 1.6f);
                dl->AddCircle(u + ImVec2(-3.5f, 4), 3, c, 0, 1.2f);
                dl->AddCircle(d + ImVec2(-3.5f, -4), 3, c, 0, 1.2f);
                hatch(u + ImVec2(-7, -4), d + ImVec2(-7, 4), ImVec2(-1, 0));
            } else if (fx || fy || fr) {
                // liên kết khác (ngàm trượt...): ô vuông kèm ký hiệu
                std::string code;
                if (fx) code += "x";
                if (fy) code += "y";
                if (fr) code += "φ";
                dl->AddRect(p - ImVec2(7, 7), p + ImVec2(7, 7), c, 0, 1.8f);
                label(dl, p + ImVec2(10, 10), code, c, ImVec2(0, 0), 0.8f);
            }
            // gối đàn hồi
            auto spring = [&](ImVec2 a, ImVec2 dir) {
                const ImVec2 t(-dir.y, dir.x);
                ImVec2 prev = a;
                for (int i = 1; i <= 8; ++i) {
                    ImVec2 q = a + dir * (i * 3.0f) + t * (i == 8 ? 0.0f : (i % 2 ? 4.0f : -4.0f));
                    dl->AddLine(prev, q, c, 1.3f);
                    prev = q;
                }
                hatch(prev - t * 7, prev + t * 7, dir);
            };
            if (!fx && sp.k[0] > 0) spring(p + ImVec2(-4, 0), ImVec2(-1, 0));
            if (!fy && sp.k[1] > 0) spring(p + ImVec2(0, 4), ImVec2(0, 1));
            if (!fr && sp.k[2] > 0) {
                dl->PathArcTo(p, 9, 0, 1.6f * kPi);
                dl->PathStroke(c, 1.4f);
            }
        }
    }

    // ---- tải trọng của trường hợp đang soạn
    if (showLoads && activeCase < (int)m.cases.size() && !(showRes && view == ResultView::Contour) &&
        view != ResultView::Influence && view != ResultView::Mode && !envMode) {
        const LoadCase& lc = m.cases[activeCase];
        const bool faint = showRes && view != ResultView::Deformed;
        const ImU32 c = faint ? withAlpha(P.load, 80) : P.load;
        const bool showValues = this->showValues && !faint;
        double qmax = 0, fmax = 0;
        for (const auto& l : lc.elem)
            if (l.kind == ElemLoadKind::Distributed) qmax = std::max({qmax, std::abs(l.q1), std::abs(l.q2)});
            else if (l.kind == ElemLoadKind::Point) fmax = std::max(fmax, std::abs(l.q1));
        for (const auto& l : lc.nodal) fmax = std::max(fmax, std::hypot(l.F[0], l.F[1]));
        for (const auto& l : lc.nodal) {
            if (l.node < 0 || l.node >= (int)m.nodes.size()) continue;
            const ImVec2 p = nodePos(l.node, false);
            const double F = std::hypot(l.F[0], l.F[1]);
            if (F > 0) {
                const ImVec2 d = norm(ImVec2((float)l.F[0], (float)-l.F[1]));
                const float L = 22 + 28 * (float)(F / std::max(fmax, 1e-12));
                arrow(dl, p - d * L, p - d * 3, c, 2.0f, 9);
                if (showValues) label(dl, p - d * (L + 10), fmtVal(F) + " kN", c);
            }
            if (l.F[2] != 0) {
                momentArc(dl, p, 15, l.F[2] > 0, c, 1.8f);
                if (showValues) label(dl, p + ImVec2(17, -17), fmtVal(std::abs(l.F[2])) + " kNm", c, ImVec2(0, 1));
            }
        }
        for (const auto& l : lc.elem) {
            if (l.elem < 0 || l.elem >= (int)m.elements.size() || !m.elements[l.elem].isBar()) continue;
            const Element& el = m.elements[l.elem];
            const ImVec2 a = nodePos(el.n[0], false), b = nodePos(el.n[1], false);
            const ImVec2 t = norm(b - a), loc(t.y, -t.x); // loc: y cục bộ (màn hình)
            const ImVec2 locY(-t.y, t.x);
            (void)loc;
            ImVec2 dir;
            switch (l.dir) {
            case LoadDir::LocalX: dir = t; break;
            case LoadDir::LocalY: dir = locY; break;
            case LoadDir::GlobalX: dir = ImVec2(1, 0); break;
            case LoadDir::GlobalY: dir = ImVec2(0, -1); break;
            }
            const double Lw = m.barLength(l.elem);
            switch (l.kind) {
            case ElemLoadKind::Distributed: {
                const int n = std::max(3, (int)(len(b - a) / 22));
                std::vector<ImVec2> tails;
                const bool along = l.dir == LoadDir::LocalX || (std::abs(dir.x * t.y - dir.y * t.x) < 0.2f);
                // tải dọc trục: vẽ lệch sang bên để không trùng thanh
                const ImVec2 shift = along ? ImVec2(-t.y, t.x) * -10.0f : ImVec2(0, 0);
                for (int i = 0; i <= n; ++i) {
                    const float r = (float)i / n;
                    const double q = l.q1 + (l.q2 - l.q1) * r;
                    const ImVec2 p = a + (b - a) * r + shift;
                    const float Lq = (float)(6 + 30 * std::abs(q) / std::max(qmax, 1e-12));
                    const ImVec2 d = dir * (float)(q >= 0 ? 1 : -1);
                    if (along) {
                        if (i < n) arrow(dl, p, p + d * 14, c, 1.2f, 6);
                        continue;
                    }
                    const ImVec2 tail = p - d * Lq;
                    tails.push_back(tail);
                    if (std::abs(q) > 1e-12) arrow(dl, tail, p, c, 1.2f, 6);
                }
                if (!tails.empty()) dl->AddPolyline(tails.data(), (int)tails.size(), c, 1.4f);
                if (showValues) {
                    std::string s = l.q1 == l.q2 ? fmtVal(std::abs(l.q1)) : fmtVal(l.q1) + " → " + fmtVal(l.q2);
                    s += " kN/m";
                    if (l.projected) s += " (h.chiếu)";
                    const ImVec2 mid = tails.empty() ? (a + b) * 0.5f + shift * 2.0f : tails[tails.size() / 2];
                    label(dl, mid + ImVec2(0, -9), s, c, ImVec2(0.5f, 1));
                }
                break;
            }
            case ElemLoadKind::Point: {
                const ImVec2 p = a + (b - a) * (float)(l.a / Lw);
                const ImVec2 d = dir * (float)(l.q1 >= 0 ? 1 : -1);
                const float L = 22 + 28 * (float)(std::abs(l.q1) / std::max(fmax, 1e-12));
                arrow(dl, p - d * L, p - d * 2, c, 2.0f, 9);
                if (showValues) label(dl, p - d * (L + 10), fmtVal(std::abs(l.q1)) + " kN", c);
                break;
            }
            case ElemLoadKind::Moment: {
                const ImVec2 p = a + (b - a) * (float)(l.a / Lw);
                momentArc(dl, p, 14, l.q1 > 0, c, 1.8f);
                if (showValues) label(dl, p + ImVec2(16, -16), fmtVal(std::abs(l.q1)) + " kNm", c, ImVec2(0, 1));
                break;
            }
            case ElemLoadKind::LackOfFit: {
                char buf[64];
                std::snprintf(buf, sizeof buf, "ΔL=%+.3g mm", l.q1 * 1e3);
                label(dl, (a + b) * 0.5f + locY * 14, buf, c);
                break;
            }
            case ElemLoadKind::Temperature: {
                char buf[96];
                std::snprintf(buf, sizeof buf, "t=%+.0f°C Δt=%+.0f°C", l.q1, l.q2);
                label(dl, (a + b) * 0.5f + locY * 14, buf, c);
                break;
            }
            }
        }
        for (const auto& d : lc.disp) {
            if (d.node < 0 || d.node >= (int)m.nodes.size()) continue;
            char buf[64];
            if (d.dof == RZ) std::snprintf(buf, sizeof buf, "Δφ=%.4g rad", d.value);
            else std::snprintf(buf, sizeof buf, "Δ%s=%.4g mm", d.dof == UX ? "x" : "y", d.value * 1e3);
            label(dl, nodePos(d.node, false) + ImVec2(12, 22), buf, c, ImVec2(0, 0));
        }
    }

    // ---- phản lực
    // khoảng ký hiệu gối chiếm phía dưới nút (px), để mũi tên phản lực đứng không đè lên
    auto supportDepth = [&](int node) {
        const Support* sp = m.supportAt(node);
        if (!sp) return 4.0f;
        const float s = 11;
        const bool fx = sp->fix[0], fy = sp->fix[1], fr = sp->fix[2];
        float d = 8;
        if (fx && fy && fr) {
            const ImVec2 dir = connectedDir(node);
            d = len(dir) < 0.1f || std::abs(dir.y) >= std::abs(dir.x) ? (dir.y > 0 ? s * 1.3f + 2 : 8.0f) : s * 1.3f + 2;
        } else if (fx && fy) d = s * 1.4f + 7;
        else if (fy && !fx && !fr) d = s * 1.3f + 13;
        else if (fx && !fy && !fr) d = s * 0.8f + 10;
        if (!fy && sp->k[1] > 0) d = std::max(d, 34.0f);
        return d;
    };
    if (showRes && showReactions && view != ResultView::Contour && view != ResultView::Influence && view != ResultView::Mode && !envMode) {
        const ImU32 c = P.reaction;
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            const ImVec2 p = nodePos((int)i, false);
            const double rx = rs->R[3 * i], ry = rs->R[3 * i + 1], rm = rs->R[3 * i + 2];
            const double tol = 1e-9;
            if (std::abs(rx) > tol) {
                const float sx = rx > 0 ? 1.f : -1.f;
                arrow(dl, p + ImVec2(-sx * 46, 14), p + ImVec2(-sx * 16, 14), c, 2.0f, 8);
                if (showValues) label(dl, p + ImVec2(-sx * 50, 14), fmtVal(std::abs(rx)), c, ImVec2(sx > 0 ? 1.f : 0.f, 0.5f));
            }
            if (std::abs(ry) > tol) {
                // nằm dưới ký hiệu gối, trên đường thẳng đứng qua tâm gối
                const float y0 = supportDepth((int)i) + 5, y1 = y0 + 30;
                if (ry > 0) arrow(dl, p + ImVec2(0, y1), p + ImVec2(0, y0), c, 2.0f, 8);
                else arrow(dl, p + ImVec2(0, y0), p + ImVec2(0, y1), c, 2.0f, 8);
                if (showValues) label(dl, p + ImVec2(0, y1 + 4), fmtVal(std::abs(ry)), c, ImVec2(0.5f, 0.f));
            }
            if (std::abs(rm) > tol) {
                momentArc(dl, p, 24, rm > 0, c, 1.8f);
                if (showValues) label(dl, p + ImVec2(-26, -26), fmtVal(std::abs(rm)), c, ImVec2(1, 1));
            }
        }
    }

    // ---- nút
    std::vector<char> trussOnly(m.nodes.size(), 0), anyBar(m.nodes.size(), 0);
    for (const auto& el : m.elements)
        for (int k = 0; k < el.nodeCount(); ++k) {
            if (el.type == ElemType::Truss && !anyBar[el.n[k]]) trussOnly[el.n[k]] = 1;
            if (el.type == ElemType::Frame) trussOnly[el.n[k]] = 0, anyBar[el.n[k]] = 1;
        }
    std::vector<char> planeNode(m.nodes.size(), 0), hinged(m.nodes.size(), 0);
    for (const auto& el : m.elements) {
        if (el.isPlane())
            for (int k = 0; k < el.nodeCount(); ++k) planeNode[el.n[k]] = 1;
        if (el.type == ElemType::Frame) {
            if (el.hingeI) hinged[el.n[0]] = 1;
            if (el.hingeJ) hinged[el.n[1]] = 1;
        }
    }
    for (size_t i = 0; i < m.nodes.size(); ++i) {
        const ImVec2 p = nodePos((int)i, false);
        if (hinged[i]) { // khớp đầu thanh: vẽ tại tâm nút
            dl->AddCircleFilled(p, 4.2f, P.hinge);
            dl->AddCircle(p, 4.2f, P.bar, 0, 1.5f);
        } else if (trussOnly[i]) {
            dl->AddCircleFilled(p, 3.8f, P.hinge);
            dl->AddCircle(p, 3.8f, P.truss, 0, 1.4f);
        } else if (!planeNode[i] || anyBar[i] || showNodeIds || cmd) {
            dl->AddCircleFilled(p, planeNode[i] ? 1.8f : 2.6f, P.node);
        }
        if (selNodes.count((int)i)) dl->AddCircle(p, 6.5f, P.select, 0, 2.2f);
        if (showNodeIds) label(dl, p + ImVec2(5, -5), std::to_string(i + 1), P.textDim, ImVec2(0, 1), 0.8f);
    }
    // phần tử được chọn / số hiệu
    for (size_t e = 0; e < m.elements.size(); ++e) {
        const Element& el = m.elements[e];
        const int nc = el.nodeCount(), ncc = el.cornerCount();
        ImVec2 pts[8], cen(0, 0);
        for (int k = 0; k < nc; ++k) pts[k] = nodePos(el.n[k], false);
        for (int k = 0; k < ncc; ++k) cen = cen + pts[k] * (1.0f / ncc);
        if (selElems.count((int)e)) {
            if (el.isBar()) dl->AddLine(pts[0], pts[1], withAlpha(P.select, 200), 5.0f);
            else {
                ImVec2 ol[24];
                dl->AddPolyline(ol, planeOutline(el, pts, ol), P.select, 2.5f, ImDrawFlags_Closed);
            }
        }
        if (showElemIds) {
            const std::string s = std::to_string(e + 1);
            label(dl, cen, s, P.textDim, ImVec2(0.5f, 0.5f), 0.8f);
            const ImVec2 ts = ImGui::CalcTextSize(s.c_str()) * 0.8f;
            dl->AddRect(cen - ts * 0.5f - ImVec2(3, 1), cen + ts * 0.5f + ImVec2(3, 1), P.textDim, 2);
        }
    }

    // ---------------------------------------------------------------- tương tác kiểu AutoCAD
    const ImVec2 mp = io.MousePos;
    double wx, wy;
    cam.toWorld(mp, wx, wy);
    hoverNode = hoverElem = -1;
    auto pickNode = [&](float radius) {
        int best = -1;
        float bd = radius;
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            const float d = len(nodePos((int)i, false) - mp);
            if (d < bd) {
                bd = d;
                best = (int)i;
            }
        }
        return best;
    };
    auto pickElem = [&]() {
        int best = -1;
        float bd = 6;
        for (size_t e = 0; e < m.elements.size(); ++e) {
            const Element& el = m.elements[e];
            if (!el.isBar()) continue;
            const ImVec2 a = nodePos(el.n[0], false), b = nodePos(el.n[1], false), ab = b - a;
            const float L2 = ab.x * ab.x + ab.y * ab.y;
            if (L2 < 1e-6f) continue;
            const float t = std::clamp(((mp.x - a.x) * ab.x + (mp.y - a.y) * ab.y) / L2, 0.0f, 1.0f);
            const float d = len(a + ab * t - mp);
            if (d < bd) {
                bd = d;
                best = (int)e;
            }
        }
        if (best >= 0) return best;
        for (size_t e = 0; e < m.elements.size(); ++e) {
            const Element& el = m.elements[e];
            if (!el.isPlane()) continue;
            int olIdx[8];
            const int nc = el.outline(olIdx);
            bool inside = false;
            for (int i = 0, j = nc - 1; i < nc; j = i++) {
                const Node& a = m.nodes[el.n[olIdx[i]]];
                const Node& b = m.nodes[el.n[olIdx[j]]];
                if ((a.y > wy) != (b.y > wy) && wx < (b.x - a.x) * (wy - a.y) / (b.y - a.y) + a.x) inside = !inside;
            }
            if (inside) return (int)e;
        }
        return -1;
    };

    // ---- bắt điểm (OSNAP / lưới / ORTHO)
    const bool pointMode = cmd && cmd->wantsPoint();
    const bool selectMode = !cmd || cmd->wantsSelection();
    double bx = 0, by = 0;
    const bool hasBase = cmd && cmd->base(bx, by);
    SnapPoint sp;
    sp.x = wx;
    sp.y = wy;
    if (pointMode && hovered) {
        const float ap = 12; // ô bắt điểm (px)
        struct Cand {
            SnapPoint p;
            float d;
            int pri;
        };
        std::vector<Cand> cands;
        auto consider = [&](SnapPoint p, int pri) {
            const float d = len(cam.toScreen(p.x, p.y) - mp);
            if (d <= ap) cands.push_back({p, d, pri});
        };
        if (osnap) {
            for (size_t i = 0; i < m.nodes.size(); ++i) {
                SnapPoint p;
                p.x = m.nodes[i].x;
                p.y = m.nodes[i].y;
                p.kind = SnapKind::Node;
                p.node = (int)i;
                consider(p, 0);
            }
            std::vector<int> near;
            for (size_t e = 0; e < m.elements.size(); ++e) {
                const Element& el = m.elements[e];
                if (!el.isBar()) continue;
                const Node& A = m.nodes[el.n[0]];
                const Node& B = m.nodes[el.n[1]];
                SnapPoint mid;
                mid.x = 0.5 * (A.x + B.x);
                mid.y = 0.5 * (A.y + B.y);
                mid.kind = SnapKind::Mid;
                mid.elem = (int)e;
                mid.t = 0.5;
                consider(mid, 1);
                const double dx = B.x - A.x, dy = B.y - A.y, L2 = dx * dx + dy * dy;
                if (L2 <= 0) continue;
                const double t = std::clamp(((wx - A.x) * dx + (wy - A.y) * dy) / L2, 0.0, 1.0);
                if (len(cam.toScreen(A.x + t * dx, A.y + t * dy) - mp) <= ap) near.push_back((int)e);
                if (hasBase) {
                    const double tp = ((bx - A.x) * dx + (by - A.y) * dy) / L2;
                    if (tp > 1e-6 && tp < 1 - 1e-6) {
                        SnapPoint pp;
                        pp.x = A.x + tp * dx;
                        pp.y = A.y + tp * dy;
                        pp.kind = SnapKind::Perpendicular;
                        pp.elem = (int)e;
                        pp.t = tp;
                        consider(pp, 3);
                    }
                }
                if (t > 1e-6 && t < 1 - 1e-6) {
                    SnapPoint np;
                    np.x = A.x + t * dx;
                    np.y = A.y + t * dy;
                    np.kind = SnapKind::Nearest;
                    np.elem = (int)e;
                    np.t = t;
                    consider(np, 4);
                }
            }
            // giao điểm của các thanh gần con trỏ
            for (size_t i = 0; i < near.size(); ++i)
                for (size_t j = i + 1; j < near.size(); ++j) {
                    const Element& e1 = m.elements[near[i]];
                    const Element& e2 = m.elements[near[j]];
                    const Node &A = m.nodes[e1.n[0]], &B = m.nodes[e1.n[1]], &C = m.nodes[e2.n[0]], &D = m.nodes[e2.n[1]];
                    const double r1x = B.x - A.x, r1y = B.y - A.y, r2x = D.x - C.x, r2y = D.y - C.y;
                    const double den = r1x * r2y - r1y * r2x;
                    if (std::abs(den) < 1e-12) continue;
                    const double t = ((C.x - A.x) * r2y - (C.y - A.y) * r2x) / den;
                    const double u = ((C.x - A.x) * r1y - (C.y - A.y) * r1x) / den;
                    if (t <= 1e-9 || t >= 1 - 1e-9 || u < -1e-9 || u > 1 + 1e-9) continue;
                    SnapPoint ip;
                    ip.x = A.x + t * r1x;
                    ip.y = A.y + t * r1y;
                    ip.kind = SnapKind::Intersection;
                    ip.elem = near[i];
                    ip.t = t;
                    ip.elem2 = near[j];
                    ip.t2 = u;
                    consider(ip, 2);
                }
        }
        if (!cands.empty()) {
            std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) {
                return a.pri != b.pri ? a.pri < b.pri : a.d < b.d;
            });
            sp = cands[0].p;
        } else {
            if (ortho && hasBase) {
                if (std::abs(wx - bx) >= std::abs(wy - by)) sp.y = by;
                else sp.x = bx;
            }
            if (snapGrid && gridStep > 0) {
                const bool lockY = ortho && hasBase && sp.y == by, lockX = ortho && hasBase && sp.x == bx;
                if (!lockX) sp.x = std::round(sp.x / gridStep) * gridStep;
                if (!lockY) sp.y = std::round(sp.y / gridStep) * gridStep;
                sp.kind = SnapKind::Grid;
            }
        }
    }
    cursor = sp;
    cursorInCanvas = hovered;
    if (hovered) {
        cursorWX = pointMode ? sp.x : wx;
        cursorWY = pointMode ? sp.y : wy;
    }

    // ---- zoom, pan
    if (hovered) {
        if (io.MouseWheel != 0 && !io.KeyShift) {
            const double f = std::pow(1.15, io.MouseWheel);
            cam.scale = std::clamp(cam.scale * f, 1e-3, 1e7);
            const ImVec2 c = (cam.c0 + cam.c1) * 0.5f;
            cam.cx = wx - (mp.x - c.x) / cam.scale;
            cam.cy = wy + (mp.y - c.y) / cam.scale;
        }
        if (io.MouseWheelH != 0 || (io.MouseWheel != 0 && io.KeyShift))
            cam.cx -= (io.MouseWheelH + (io.KeyShift ? io.MouseWheel : 0)) * 40 / cam.scale;
        if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Middle)) fitRequested = true;
    }
    if (ImGui::IsItemActive() &&
        (ImGui::IsMouseDragging(ImGuiMouseButton_Right, 3) || ImGui::IsMouseDragging(ImGuiMouseButton_Middle, 1))) {
        cam.cx -= io.MouseDelta.x / cam.scale;
        cam.cy += io.MouseDelta.y / cam.scale;
    }

    // ---- đối tượng dưới con trỏ (chế độ chọn)
    if (hovered && selectMode && !g_in.boxing) {
        hoverNode = pickNode(8);
        if (hoverNode < 0) hoverElem = pickElem();
    }
    if (hoverNode >= 0) dl->AddCircle(nodePos(hoverNode, false), 7, P.hover, 0, 2);
    else if (hoverElem >= 0) {
        const Element& el = m.elements[hoverElem];
        ImVec2 pts[8], ol[24];
        for (int k = 0; k < el.nodeCount(); ++k) pts[k] = nodePos(el.n[k], false);
        if (el.isBar()) dl->AddLine(pts[0], pts[1], withAlpha(P.hover, 170), 4.5f);
        else dl->AddPolyline(ol, planeOutline(el, pts, ol), P.hover, 2.0f, ImDrawFlags_Closed);
    }

    // ---- grip của đối tượng đang chọn (khi không có lệnh)
    int hoverGrip = -1;
    if (!cmd) {
        for (int n : selectionNodes()) {
            const ImVec2 p = nodePos(n, false);
            const bool hv = hovered && std::abs(p.x - mp.x) <= 6 && std::abs(p.y - mp.y) <= 6;
            if (hv) hoverGrip = n;
            dl->AddRectFilled(p - ImVec2(4.5f, 4.5f), p + ImVec2(4.5f, 4.5f), hv ? P.gripHover : P.grip);
            dl->AddRect(p - ImVec2(4.5f, 4.5f), p + ImVec2(4.5f, 4.5f), withAlpha(P.bg, 255));
        }
    }

    // ---- xem trước của lệnh
    if (cmd && hovered) cmd->preview(*this, dl, sp);

    // ---- chuột trái
    const bool add = !io.KeyShift; // AutoCAD: nhấp để thêm vào tập chọn, Shift+nhấp để bỏ ra
    auto finishBox = [&](ImVec2 end) {
        const bool window = end.x >= g_in.start.x; // trái→phải: cửa sổ; phải→trái: cắt ngang
        const ImVec2 r0(std::min(g_in.start.x, end.x), std::min(g_in.start.y, end.y));
        const ImVec2 r1(std::max(g_in.start.x, end.x), std::max(g_in.start.y, end.y));
        auto inside = [&](ImVec2 p) { return p.x >= r0.x && p.x <= r1.x && p.y >= r0.y && p.y <= r1.y; };
        auto segHits = [&](ImVec2 a, ImVec2 b) {
            if (inside(a) || inside(b)) return true;
            // giao với 4 cạnh khung
            auto cross = [](ImVec2 p, ImVec2 q, ImVec2 r, ImVec2 s2) {
                auto o = [](ImVec2 a2, ImVec2 b2, ImVec2 c2) { return (b2.x - a2.x) * (c2.y - a2.y) - (b2.y - a2.y) * (c2.x - a2.x); };
                return (o(p, q, r) * o(p, q, s2) < 0) && (o(r, s2, p) * o(r, s2, q) < 0);
            };
            const ImVec2 c[4] = {r0, ImVec2(r1.x, r0.y), r1, ImVec2(r0.x, r1.y)};
            for (int k = 0; k < 4; ++k)
                if (cross(a, b, c[k], c[(k + 1) % 4])) return true;
            return false;
        };
        auto put = [&](std::set<int>& s2, int id) {
            if (add) s2.insert(id);
            else s2.erase(id);
        };
        for (size_t i = 0; i < m.nodes.size(); ++i)
            if (inside(nodePos((int)i, false))) put(selNodes, (int)i);
        for (size_t e = 0; e < m.elements.size(); ++e) {
            const Element& el = m.elements[e];
            int olIdx[8];
            const int nc = el.outline(olIdx);
            bool all = true, any = false;
            for (int k = 0; k < nc; ++k) {
                const ImVec2 a = nodePos(el.n[olIdx[k]], false), b = nodePos(el.n[olIdx[(k + 1) % nc]], false);
                all &= inside(a);
                if (el.isBar() && k == 1) break;
                any |= segHits(a, b);
            }
            if (window ? all : (all || any)) put(selElems, (int)e);
        }
        g_in.boxing = false;
    };
    if (hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        if (pointMode) {
            cmd->onPoint(*this, sp);
            if (cmd && cmd->done) cmd.reset();
        } else if (g_in.boxing) {
            finishBox(mp); // kiểu AutoCAD: nhấp – di chuột – nhấp
        } else if (!cmd && hoverGrip >= 0) {
            // nhấp grip: kéo giãn nút
            struct GripCmd : Command {
                int node;
                double x0, y0;
                GripCmd(int n, double x, double y) : node(n), x0(x), y0(y) { name = "** KÉO GIÃN **"; }
                bool wantsPoint() const override { return true; }
                bool base(double& x, double& y) const override {
                    x = x0;
                    y = y0;
                    return true;
                }
                std::string prompt(App&) override { return "Chỉ định điểm mới cho nút:"; }
                void onPoint(App& a, const SnapPoint& p) override {
                    a.snapshot();
                    a.model.nodes[node] = {p.x, p.y};
                    double X0, Y0, X1, Y1;
                    a.model.bounds(X0, Y0, X1, Y1);
                    mergeNodes(a.model, 1e-6 * std::max({1.0, X1 - X0, Y1 - Y0}));
                    a.selNodes.clear();
                    a.selElems.clear();
                    a.changed();
                    done = true;
                }
                void preview(App& a, ImDrawList* d, const SnapPoint& c) override {
                    for (const auto& el : a.model.elements) {
                        if (!el.isBar()) continue;
                        for (int k = 0; k < 2; ++k)
                            if (el.n[k] == node) {
                                const Node& o = a.model.nodes[el.n[1 - k]];
                                d->AddLine(a.cam.toScreen(o.x, o.y), a.cam.toScreen(c.x, c.y), a.pal.hover, 1.5f);
                            }
                    }
                    const ImVec2 p = a.cam.toScreen(x0, y0);
                    d->AddRectFilled(ImVec2(p.x - 4.5f, p.y - 4.5f), ImVec2(p.x + 4.5f, p.y + 4.5f), a.pal.gripHot);
                }
            };
            cmd = std::make_unique<GripCmd>(hoverGrip, m.nodes[hoverGrip].x, m.nodes[hoverGrip].y);
        } else if (selectMode) {
            g_in.leftDown = true;
            g_in.start = mp;
            if (hoverNode >= 0) {
                if (add) selNodes.insert(hoverNode);
                else selNodes.erase(hoverNode);
                g_in.leftDown = false;
            } else if (hoverElem >= 0) {
                if (add) selElems.insert(hoverElem);
                else selElems.erase(hoverElem);
                g_in.leftDown = false;
            } else {
                g_in.boxing = true; // bắt đầu khung chọn
            }
        }
    }
    if (g_in.leftDown && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
        g_in.leftDown = false;
        if (g_in.boxing && len(mp - g_in.start) > 5) finishBox(mp); // kéo – thả
    }
    if (g_in.boxing && !selectMode) g_in.boxing = false;
    if (g_in.boxing) {
        const bool window = mp.x >= g_in.start.x;
        const ImU32 c = window ? P.windowSel : P.crossingSel;
        dl->AddRectFilled(g_in.start, mp, withAlpha(c, 40));
        if (window) dl->AddRect(g_in.start, mp, c, 0, 1.2f);
        else {
            // viền nét đứt cho chọn cắt ngang
            const ImVec2 r0(std::min(g_in.start.x, mp.x), std::min(g_in.start.y, mp.y));
            const ImVec2 r1(std::max(g_in.start.x, mp.x), std::max(g_in.start.y, mp.y));
            auto dash = [&](ImVec2 a, ImVec2 b) {
                const float L = len(b - a);
                const ImVec2 d = norm(b - a);
                for (float t = 0; t < L; t += 8) dl->AddLine(a + d * t, a + d * std::min(t + 4, L), c, 1.2f);
            };
            dash(r0, ImVec2(r1.x, r0.y));
            dash(ImVec2(r1.x, r0.y), r1);
            dash(r1, ImVec2(r0.x, r1.y));
            dash(ImVec2(r0.x, r1.y), r0);
        }
    }

    // ---- chuột phải (không kéo): Enter khi đang có lệnh, menu ngữ cảnh khi rảnh
    if (hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
        g_in.rightDown = true;
        g_in.start = mp;
    }
    if (g_in.rightDown && ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
        g_in.rightDown = false;
        if (len(mp - g_in.start) < 4) {
            if (g_in.boxing) g_in.boxing = false;
            else if (cmd) {
                cmd->onEnter(*this);
                if (cmd && cmd->done) cmd.reset();
            } else {
                if (hoverNode >= 0 && !selNodes.count(hoverNode)) selNodes.insert(hoverNode);
                else if (hoverElem >= 0 && !selElems.count(hoverElem)) selElems.insert(hoverElem);
                ImGui::OpenPopup("ctx");
            }
        }
    }
    if (ImGui::BeginPopup("ctx")) {
        const bool sel = !selNodes.empty() || !selElems.empty();
        if (!lastCommand.empty() && ImGui::MenuItem(("Lặp lại " + lastCommand).c_str(), "Enter")) startCommand(lastCommand);
        ImGui::Separator();
        if (ImGui::MenuItem("Xóa", "Delete", false, sel)) startCommand("E");
        if (ImGui::MenuItem("Di chuyển", "M", false, sel)) startCommand("M");
        if (ImGui::MenuItem("Sao chép", "CO", false, sel)) startCommand("CO");
        if (ImGui::MenuItem("Xoay", "RO", false, sel)) startCommand("RO");
        if (ImGui::MenuItem("Đối xứng", "MI", false, sel)) startCommand("MI");
        if (ImGui::MenuItem("Chia thanh", "DIV", false, !selElems.empty())) startCommand("DIV");
        ImGui::Separator();
        if (ImGui::MenuItem("Ngàm", "NGAM", false, !selNodes.empty())) startCommand("NGAM");
        if (ImGui::MenuItem("Gối cố định", "KHOP", false, !selNodes.empty())) startCommand("KHOP");
        if (ImGui::MenuItem("Gối di động", "DIDONG", false, !selNodes.empty())) startCommand("DIDONG");
        if (ImGui::MenuItem("Tải nút...", "TAINUT", false, !selNodes.empty())) startCommand("TAINUT");
        if (ImGui::MenuItem("Tải phân bố...", "TAIPB", false, !selElems.empty())) startCommand("TAIPB");
        ImGui::Separator();
        if (ImGui::MenuItem("Bỏ chọn tất cả", "Esc", false, sel)) {
            selNodes.clear();
            selElems.clear();
        }
        if (ImGui::MenuItem("Thuộc tính", "Ctrl+1")) showPropsPanel = true;
        if (ImGui::MenuItem("Vừa khung nhìn", "Z E")) fitRequested = true;
        ImGui::EndPopup();
    }

    // ---- dấu bắt điểm
    if (pointMode && hovered && sp.kind != SnapKind::None && sp.kind != SnapKind::Grid) {
        const ImVec2 p = cam.toScreen(sp.x, sp.y);
        const float r = 7;
        const ImU32 c = P.snap;
        switch (sp.kind) {
        case SnapKind::Node: dl->AddRect(p - ImVec2(r, r), p + ImVec2(r, r), c, 0, 0, 2.0f); break;
        case SnapKind::Mid: dl->AddTriangle(p + ImVec2(0, -r), p + ImVec2(-r, r * 0.8f), p + ImVec2(r, r * 0.8f), c, 2.0f); break;
        case SnapKind::Intersection:
            dl->AddLine(p - ImVec2(r, r), p + ImVec2(r, r), c, 2.0f);
            dl->AddLine(p + ImVec2(-r, r), p + ImVec2(r, -r), c, 2.0f);
            break;
        case SnapKind::Perpendicular:
            dl->AddLine(p + ImVec2(-r, r), p + ImVec2(r, r), c, 2.0f);
            dl->AddLine(p + ImVec2(-r, r), p + ImVec2(-r, -r), c, 2.0f);
            dl->AddLine(p + ImVec2(-r, 0), p + ImVec2(0, 0), c, 2.0f);
            dl->AddLine(p + ImVec2(0, 0), p + ImVec2(0, r), c, 2.0f);
            break;
        case SnapKind::Nearest:
            dl->AddLine(p + ImVec2(-r, -r), p + ImVec2(r, -r), c, 2.0f);
            dl->AddLine(p + ImVec2(r, -r), p + ImVec2(-r, r), c, 2.0f);
            dl->AddLine(p + ImVec2(-r, r), p + ImVec2(r, r), c, 2.0f);
            dl->AddLine(p + ImVec2(r, r), p + ImVec2(-r, -r), c, 2.0f);
            break;
        default: break;
        }
        label(dl, p + ImVec2(12, 12), snapName(sp.kind), P.snap, ImVec2(0, 0), 0.85f);
    }

    // ---- con trỏ chữ thập
    if (hovered && !ImGui::IsPopupOpen("ctx")) {
        ImGui::SetMouseCursor(ImGuiMouseCursor_None);
        const ImVec2 cp = pointMode ? cam.toScreen(sp.x, sp.y) : mp;
        const ImU32 c = withAlpha(P.crosshair, 200);
        dl->AddLine(ImVec2(cam.c0.x, cp.y), ImVec2(cam.c1.x, cp.y), withAlpha(P.crosshair, 70), 1.0f);
        dl->AddLine(ImVec2(cp.x, cam.c0.y), ImVec2(cp.x, cam.c1.y), withAlpha(P.crosshair, 70), 1.0f);
        dl->AddLine(ImVec2(cp.x - 40, cp.y), ImVec2(cp.x + 40, cp.y), c, 1.0f);
        dl->AddLine(ImVec2(cp.x, cp.y - 40), ImVec2(cp.x, cp.y + 40), c, 1.0f);
        if (selectMode) dl->AddRect(cp - ImVec2(5, 5), cp + ImVec2(5, 5), c);
    }

    // ---- nhập liệu động cạnh con trỏ
    if (dynInput && hovered && cmd) {
        std::string t1 = cmd->prompt(*this);
        std::string t2 = cmdLine.empty() ? "" : "› " + cmdLine;
        if (hasBase && pointMode) {
            char b[96];
            std::snprintf(b, sizeof b, "L = %.3f   ∠ %.1f°", std::hypot(sp.x - bx, sp.y - by),
                          std::atan2(sp.y - by, sp.x - bx) * 180 / kPi);
            t2 = t2.empty() ? b : t2 + "   " + b;
        }
        const ImVec2 at = mp + ImVec2(22, 22);
        const ImVec2 s1 = ImGui::CalcTextSize(t1.c_str()), s2 = ImGui::CalcTextSize(t2.c_str());
        const float w = std::max(s1.x, s2.x) + 12, hh = s1.y + (t2.empty() ? 0 : s2.y) + 8;
        dl->AddRectFilled(at, at + ImVec2(w, hh), withAlpha(darkTheme ? IM_COL32(20, 24, 30, 255) : IM_COL32(255, 255, 240, 255), 225), 3);
        dl->AddRect(at, at + ImVec2(w, hh), withAlpha(P.textDim, 160), 3);
        dl->AddText(at + ImVec2(6, 4), P.textDim, t1.c_str());
        if (!t2.empty()) dl->AddText(at + ImVec2(6, 4 + s1.y), P.text, t2.c_str());
    }

    // ---- chú giải phổ màu
    if (contourOn) {
        const float bw = 18, bh = std::min(260.0f, sz.y - 60);
        const ImVec2 o(cam.c1.x - 120, cam.c0.y + 34);
        label(dl, o + ImVec2(bw, -8), contourFieldName(contour, m.planeMode == PlaneMode::Axisymmetric), P.text, ImVec2(0.5f, 1), 0.9f);
        const int nb = 64;
        for (int i = 0; i < nb; ++i) {
            const float y0 = o.y + bh * i / nb, y1 = o.y + bh * (i + 1) / nb;
            dl->AddRectFilled(ImVec2(o.x, y0), ImVec2(o.x + bw, y1 + 0.5f), colormapU32(1.0 - (i + 0.5) / nb));
        }
        dl->AddRect(o, o + ImVec2(bw, bh), P.textDim);
        const int nt = 8;
        for (int i = 0; i <= nt; ++i) {
            const double v = cmax - (cmax - cmin) * i / nt;
            const float y = o.y + bh * i / nt;
            dl->AddLine(ImVec2(o.x + bw, y), ImVec2(o.x + bw + 4, y), P.textDim);
            label(dl, ImVec2(o.x + bw + 7, y), fmtVal(v), P.text, ImVec2(0, 0.5f), 0.85f);
        }
    }

    // ---- biểu tượng UCS (góc trái dưới)
    {
        const ImVec2 o(cam.c0.x + 26, cam.c1.y - 26);
        dl->AddLine(o, o + ImVec2(42, 0), P.axisX, 2.0f);
        dl->AddTriangleFilled(o + ImVec2(48, 0), o + ImVec2(40, -4), o + ImVec2(40, 4), P.axisX);
        dl->AddLine(o, o + ImVec2(0, -42), P.axisY, 2.0f);
        dl->AddTriangleFilled(o + ImVec2(0, -48), o + ImVec2(-4, -40), o + ImVec2(4, -40), P.axisY);
        dl->AddRect(o - ImVec2(4, 4), o + ImVec2(4, 4), P.textDim);
        dl->AddText(o + ImVec2(52, -8), P.axisX, "X");
        dl->AddText(o + ImVec2(-4, -66), P.axisY, "Y");
    }
    if (!fmTitle.empty()) label(dl, ImVec2(cam.c0.x + 10, cam.c0.y + 30), fmTitle, P.select, ImVec2(0, 0.5f));
    if (!result.ok && !result.error.empty())
        label(dl, ImVec2(cam.c0.x + 10, cam.c0.y + 12), "⚠ " + result.error, IM_COL32(230, 60, 50, 255), ImVec2(0, 0.5f), 0.9f);
    else if (needAnalyze && !autoAnalyze)
        label(dl, ImVec2(cam.c0.x + 10, cam.c0.y + 12), "Mô hình đã thay đổi — nhấn F5 để tính lại", IM_COL32(220, 140, 0, 255),
              ImVec2(0, 0.5f), 0.9f);
    if (drawDeformed) {
        char b[64];
        std::snprintf(b, sizeof b, "phóng đại biến dạng ×%.0f", S);
        label(dl, ImVec2(cam.c1.x - 10, cam.c1.y - 12), b, P.textDim, ImVec2(1, 0.5f), 0.85f);
    }
    dl->PopClipRect();

    // ---- chú thích khi rê chuột
    if (hovered && !cmd && !g_in.boxing && !ImGui::IsMouseDragging(ImGuiMouseButton_Right)) {
        if (hoverNode >= 0) {
            ImGui::BeginTooltip();
            const Node& n = m.nodes[hoverNode];
            ImGui::Text("Nút %d  (%.3f; %.3f)", hoverNode + 1, n.x, n.y);
            if (rs) {
                ImGui::Text("ux = %.4f mm, uy = %.4f mm, φ = %.3e rad", rs->U[3 * hoverNode] * 1e3,
                            rs->U[3 * hoverNode + 1] * 1e3, rs->U[3 * hoverNode + 2]);
                const double* R = &rs->R[3 * hoverNode];
                if (R[0] != 0 || R[1] != 0 || R[2] != 0)
                    ImGui::Text("Phản lực: Rx = %.3f kN, Ry = %.3f kN, M = %.3f kNm", R[0], R[1], R[2]);
                if (!rs->nodeStressCount.empty() && rs->nodeStressCount[hoverNode] > 0) {
                    const auto& s = rs->nodeStress[hoverNode];
                    const StressInfo si = stressInfo(s[0], s[1], s[2], s[3]);
                    ImGui::Text("σx = %.3f, σy = %.3f, τxy = %.3f MPa", si.sx * 1e-3, si.sy * 1e-3, si.txy * 1e-3);
                    ImGui::Text("σ1 = %.3f, σ2 = %.3f MPa, góc %.1f°, σvM = %.3f MPa", si.s1 * 1e-3, si.s2 * 1e-3,
                                si.angle, si.vonMises * 1e-3);
                }
            }
            ImGui::EndTooltip();
        } else if (hoverElem >= 0) {
            const Element& el = m.elements[hoverElem];
            ImGui::BeginTooltip();
            ImGui::Text("Phần tử %d — %s", hoverElem + 1, elemTypeName(el.type));
            if (el.isBar()) {
                const Node& a = m.nodes[el.n[0]];
                const Node& b = m.nodes[el.n[1]];
                const double L = std::hypot(b.x - a.x, b.y - a.y);
                ImGui::Text("Nút %d → %d, L = %.3f m, %s", el.n[0] + 1, el.n[1] + 1, L,
                            el.sec < (int)m.sections.size() ? m.sections[el.sec].name.c_str() : "?");
                if (rs) {
                    const BarResult& br = rs->bars[hoverElem];
                    const ImVec2 pa = nodePos(el.n[0], false), pb = nodePos(el.n[1], false), ab = pb - pa;
                    const float t = std::clamp(((mp.x - pa.x) * ab.x + (mp.y - pa.y) * ab.y) / (ab.x * ab.x + ab.y * ab.y),
                                               0.0f, 1.0f);
                    const double x = t * L;
                    size_t k = 0;
                    for (size_t i = 0; i < br.x.size(); ++i)
                        if (std::abs(br.x[i] - x) < std::abs(br.x[k] - x)) k = i;
                    ImGui::Separator();
                    ImGui::Text("Tại x = %.3f m:", br.x[k]);
                    ImGui::Text("N = %.3f kN   Q = %.3f kN   M = %.3f kNm", br.N[k], br.Q[k], br.M[k]);
                    ImGui::Text("Chuyển vị cục bộ: u = %.4f mm, v = %.4f mm", br.u[k] * 1e3, br.v[k] * 1e3);
                    const auto [mn, mx] = std::minmax_element(br.M.begin(), br.M.end());
                    ImGui::Text("Mmax = %.3f, Mmin = %.3f kNm", *mx, *mn);
                }
            } else if (rs) {
                for (int k = 0; k < el.nodeCount(); ++k) {
                    const auto& s = rs->elemStress[hoverElem][k];
                    const StressInfo si = stressInfo(s[0], s[1], s[2], s[3]);
                    if (m.planeMode == PlaneMode::Axisymmetric)
                        ImGui::Text("Nút %d: σr %.3f  σz %.3f  τrz %.3f  σθ %.3f  σvM %.3f MPa", el.n[k] + 1, si.sx * 1e-3,
                                    si.sy * 1e-3, si.txy * 1e-3, si.sz * 1e-3, si.vonMises * 1e-3);
                    else
                        ImGui::Text("Nút %d: σx %.3f  σy %.3f  τxy %.3f  σvM %.3f MPa", el.n[k] + 1, si.sx * 1e-3,
                                    si.sy * 1e-3, si.txy * 1e-3, si.vonMises * 1e-3);
                }
            }
            ImGui::EndTooltip();
        }
    }
    ImGui::End();
}
