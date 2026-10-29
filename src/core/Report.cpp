#include "core/Report.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <sstream>

#include "core/PostProcess.h"
#include "core/I18n.h"

namespace fem {

namespace {

#ifdef __GNUC__
std::string fmt(const char* f, ...) __attribute__((format(printf, 1, 2)));
#endif
std::string fmt(const char* f, ...) {
    char b[1024];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(b, sizeof b, f, ap);
    va_end(ap);
    return b;
}

std::string esc(const std::string& s) {
    std::string r;
    for (char c : s) {
        if (c == '<') r += "&lt;";
        else if (c == '>') r += "&gt;";
        else if (c == '&') r += "&amp;";
        else r += c;
    }
    return r;
}

std::string val(double v) {
    const double a = std::abs(v);
    return fmt(a >= 100 ? "%.1f" : a >= 1 ? "%.2f" : a == 0 ? "0" : "%.3f", v);
}

// Khung vẽ SVG: tọa độ thực → điểm ảnh
struct Svg {
    double x0, y0, x1, y1, k, pad = 40;
    int W = 820, H = 0;
    std::ostringstream o;
    explicit Svg(const Model& m, double margin) {
        m.bounds(x0, y0, x1, y1);
        const double size = std::max({x1 - x0, y1 - y0, 1.0});
        x0 -= margin * size;
        x1 += margin * size;
        y0 -= margin * size;
        y1 += margin * size;
        k = (W - 2 * pad) / std::max(x1 - x0, 1e-9);
        H = (int)std::clamp((y1 - y0) * k + 2 * pad, 160.0, 900.0);
        k = std::min(k, (H - 2 * pad) / std::max(y1 - y0, 1e-9));
    }
    double X(double x) const { return pad + (x - x0) * k; }
    double Y(double y) const { return H - pad - (y - y0) * k; }
    void line(double xa, double ya, double xb, double yb, const char* col, double w, const char* extra = "") {
        o << fmt("<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' stroke='%s' stroke-width='%.1f' %s/>\n", X(xa), Y(ya), X(xb), Y(yb),
                 col, w, extra);
    }
    void text(double x, double y, const std::string& t, const char* col, int size = 11, const char* anchor = "middle") {
        o << fmt("<text x='%.1f' y='%.1f' fill='%s' font-size='%d' text-anchor='%s'>%s</text>\n", X(x), Y(y), col, size, anchor,
                 esc(t).c_str());
    }
    std::string str() const {
        return fmt("<svg xmlns='http://www.w3.org/2000/svg' width='%d' height='%d' viewBox='0 0 %d %d'>", W, H, W, H) +
               "<rect width='100%' height='100%' fill='white'/>\n" + o.str() + "</svg>\n";
    }
};

void drawStructure(Svg& s, const Model& m, bool faint) {
    const char* col = faint ? "#9aa3ad" : "#1d2530";
    for (const auto& el : m.elements) {
        int ol[8];
        const int nc = el.outline(ol);
        if (el.isBar()) {
            const Node& a = m.nodes[el.n[0]];
            const Node& b = m.nodes[el.n[1]];
            s.line(a.x, a.y, b.x, b.y, el.type == ElemType::Truss ? "#3a6ea5" : col, faint ? 1.2 : 2.2);
        } else {
            for (int k = 0; k < nc; ++k) {
                const Node& a = m.nodes[el.n[ol[k]]];
                const Node& b = m.nodes[el.n[ol[(k + 1) % nc]]];
                s.line(a.x, a.y, b.x, b.y, "#b8c2cf", 0.8);
            }
        }
    }
    for (const auto& sp : m.supports) {
        if (sp.empty() || sp.node < 0 || sp.node >= (int)m.nodes.size()) continue;
        const double X = s.X(m.nodes[sp.node].x), Y = s.Y(m.nodes[sp.node].y);
        if (sp.fix[0] && sp.fix[1] && sp.fix[2])
            s.o << fmt("<path d='M%.1f %.1f h-12 h24' stroke='#2a8a45' stroke-width='2.5' fill='none'/>"
                       "<path d='M%.1f %.1f l-5 7 M%.1f %.1f l-5 7 M%.1f %.1f l-5 7 M%.1f %.1f l-5 7' stroke='#2a8a45'/>\n",
                       X, Y, X - 10, Y, X - 3, Y, X + 4, Y, X + 11, Y);
        else
            s.o << fmt("<path d='M%.1f %.1f l-8 13 h16 z' stroke='#2a8a45' stroke-width='1.6' fill='none'/>\n", X, Y);
    }
}

std::string diagramSvg(const Model& m, const LoadSetResult& rs, BarQty q) {
    Svg s(m, 0.12);
    drawStructure(s, m, true);
    double vmax = 0;
    for (size_t e = 0; e < m.elements.size(); ++e)
        if (m.elements[e].isBar())
            for (double v : barValues(rs.bars[e], q)) vmax = std::max(vmax, std::abs(v));
    if (vmax <= 1e-12) return "";
    const double sc = 50.0 / s.k / vmax; // biểu đồ lớn nhất ≈ 50 px
    for (size_t e = 0; e < m.elements.size(); ++e) {
        if (!m.elements[e].isBar()) continue;
        const DiagramGeom g = diagramGeometry(m, rs, (int)e, q, sc);
        const size_t n = g.base.size();
        for (size_t k = 0; k + 1 < n; ++k) {
            const bool pos = g.values[k] + g.values[k + 1] >= 0;
            s.o << fmt("<polygon points='%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f' fill='%s' fill-opacity='0.25' stroke='none'/>\n",
                       s.X(g.base[k][0]), s.Y(g.base[k][1]), s.X(g.off[k][0]), s.Y(g.off[k][1]), s.X(g.off[k + 1][0]),
                       s.Y(g.off[k + 1][1]), s.X(g.base[k + 1][0]), s.Y(g.base[k + 1][1]), pos ? "#2f6fd6" : "#d9443a");
            s.line(g.off[k][0], g.off[k][1], g.off[k + 1][0], g.off[k + 1][1], pos ? "#1f55b5" : "#b8302a", 1.6);
        }
        auto put = [&](size_t k) {
            if (std::abs(g.values[k]) < 1e-4 * vmax) return;
            s.text(g.off[k][0], g.off[k][1], val(g.values[k]), g.values[k] >= 0 ? "#1f55b5" : "#b8302a", 10);
        };
        if (n >= 2) {
            put(0);
            put(n - 1);
            size_t kx = 0;
            for (size_t k = 1; k + 1 < n; ++k)
                if (std::abs(g.values[k]) > std::abs(g.values[kx])) kx = k;
            if (kx > 0 && kx + 1 < n && std::abs(g.values[kx]) > 1.02 * std::max(std::abs(g.values[0]), std::abs(g.values[n - 1])))
                put(kx);
        }
    }
    return s.str();
}

std::string deformedSvg(const Model& m, const LoadSetResult& rs) {
    Svg s(m, 0.12);
    drawStructure(s, m, true);
    double umax = 0;
    int imax = 0;
    for (size_t i = 0; i < m.nodes.size(); ++i) {
        const double u = std::hypot(rs.U[3 * i], rs.U[3 * i + 1]);
        if (u > umax) {
            umax = u;
            imax = (int)i;
        }
    }
    if (umax <= 0) return "";
    double x0, y0, x1, y1;
    m.bounds(x0, y0, x1, y1);
    const double S = 0.08 * std::max({x1 - x0, y1 - y0, 1.0}) / umax;
    for (size_t e = 0; e < m.elements.size(); ++e) {
        const Element& el = m.elements[e];
        if (el.isBar()) {
            const BarResult& b = rs.bars[e];
            const Node& a = m.nodes[el.n[0]];
            const Node& q = m.nodes[el.n[1]];
            const double L = std::hypot(q.x - a.x, q.y - a.y), c = (q.x - a.x) / L, sn = (q.y - a.y) / L;
            std::string pts;
            for (size_t k = 0; k < b.x.size(); ++k)
                pts += fmt("%.1f,%.1f ", s.X(a.x + c * b.x[k] + (c * b.u[k] - sn * b.v[k]) * S),
                           s.Y(a.y + sn * b.x[k] + (sn * b.u[k] + c * b.v[k]) * S));
            s.o << "<polyline points='" << pts << "' fill='none' stroke='#c2268f' stroke-width='1.8'/>\n";
        } else {
            std::string pts;
            int ol[8];
            const int nc = el.outline(ol);
            for (int k = 0; k < nc; ++k) {
                const int i = el.n[ol[k]];
                pts += fmt("%.1f,%.1f ", s.X(m.nodes[i].x + rs.U[3 * i] * S), s.Y(m.nodes[i].y + rs.U[3 * i + 1] * S));
            }
            s.o << "<polygon points='" << pts << "' fill='none' stroke='#c2268f' stroke-width='0.8'/>\n";
        }
    }
    s.text(m.nodes[imax].x + rs.U[3 * imax] * S, m.nodes[imax].y + rs.U[3 * imax + 1] * S,
           trf("δmax = %.3f mm (nút %d)", umax * 1e3, imax + 1), "#c2268f", 11, "start");
    return s.str();
}

} // namespace

std::string htmlReport(const Model& m, const AnalysisResult& r, const HtmlReportOptions& opt) {
    std::ostringstream o;
    o << "<!doctype html><html lang='" << languageCode() << "'><head><meta charset='utf-8'><title>" << esc(m.title.empty() ? tr("Báo cáo FEM2D") : m.title)
      << "</title><style>"
         "body{font-family:-apple-system,Segoe UI,Arial,sans-serif;margin:24px auto;max-width:900px;color:#1d2530}"
         "h1{font-size:22px;border-bottom:2px solid #0a6ed1;padding-bottom:6px}h2{font-size:18px;margin-top:28px;color:#0a5aa8}"
         "h3{font-size:15px;margin-top:18px}table{border-collapse:collapse;font-size:12px;margin:6px 0}"
         "td,th{border:1px solid #c8d0da;padding:3px 8px;text-align:right}th{background:#eef2f7}"
         "td:first-child,th:first-child{text-align:center}.note{color:#667;font-size:12px}"
         "svg{border:1px solid #dde3ea;margin:6px 0;max-width:100%;height:auto}"
         "@media print{h2{page-break-before:always}h2.first{page-break-before:avoid}}"
         "</style></head><body>\n";
    o << "<h1>" << esc(m.title.empty() ? tr("Báo cáo tính toán kết cấu") : m.title) << "</h1>\n";
    o << tr("<p class='note'>Lập bằng FEM2D — phương pháp phần tử hữu hạn (phương pháp chuyển vị). Đơn vị: kN, m; chuyển vị mm; ứng suất MPa. Quy ước: N &gt; 0 kéo; M &gt; 0 căng thớ dưới, biểu đồ M vẽ về phía thớ căng.</p>\n");
    // ---- dữ liệu
    o << tr("<h2 class='first'>1. Sơ đồ tính và dữ liệu</h2>\n");
    {
        Svg s(m, 0.1);
        drawStructure(s, m, false);
        for (size_t i = 0; i < m.nodes.size() && m.nodes.size() <= 120; ++i)
            s.text(m.nodes[i].x, m.nodes[i].y, std::to_string(i + 1), "#556", 9, "start");
        o << s.str();
    }
    int nb = 0, np = 0;
    for (const auto& e : m.elements) (e.isBar() ? nb : np)++;
    o << trf("<p>Số nút: %zu; số thanh: %d; số phần tử phẳng: %d; số gối: %zu; số phương trình: %d.</p>\n", m.nodes.size(), nb,
             np, m.supports.size(), r.numEquations);
    o << tr("<h3>Vật liệu</h3><table><tr><th>#</th><th>Tên</th><th>E (MPa)</th><th>ν</th><th>γ (kN/m³)</th></tr>\n");
    for (size_t i = 0; i < m.materials.size(); ++i) {
        const auto& x = m.materials[i];
        o << fmt("<tr><td>%zu</td><td style='text-align:left'>%s</td><td>%.0f</td><td>%.2f</td><td>%.1f</td></tr>\n", i + 1,
                 esc(x.name).c_str(), x.E / 1e3, x.nu, x.gamma);
    }
    o << "</table>\n";
    if (nb) {
        o << tr("<h3>Tiết diện</h3><table><tr><th>#</th><th>Tên</th><th>A (cm²)</th><th>I (cm⁴)</th><th>h (cm)</th></tr>\n");
        for (size_t i = 0; i < m.sections.size(); ++i) {
            const auto& x = m.sections[i];
            o << fmt("<tr><td>%zu</td><td style='text-align:left'>%s</td><td>%.2f</td><td>%.1f</td><td>%.1f</td></tr>\n", i + 1,
                     esc(x.name).c_str(), x.A * 1e4, x.I * 1e8, x.h * 100);
        }
        o << "</table>\n";
    }
    o << tr("<h3>Trường hợp tải và tổ hợp</h3><ul>\n");
    for (const auto& c : m.cases)
        o << "<li>" << esc(c.name) << trf(": %zu tải nút, %zu tải trên thanh%s</li>\n", c.nodal.size(), c.elem.size(),
                                           c.selfWeight != 0 ? tr(", có trọng lượng bản thân") : "");
    for (const auto& c : m.combos) {
        o << "<li>" << tr("Tổ hợp") << " " << esc(c.name) << ": ";
        for (size_t k = 0; k < c.terms.size(); ++k)
            o << (k ? " + " : "") << fmt("%.3g×", c.terms[k].second)
              << esc(c.terms[k].first < (int)m.cases.size() ? m.cases[c.terms[k].first].name : "?");
        o << "</li>\n";
    }
    o << "</ul>\n";
    if (!r.ok) {
        o << "<p><b>" << tr("Lỗi tính toán:") << "</b> " << esc(r.error) << "</p></body></html>";
        return o.str();
    }
    for (const auto& w : r.warnings) o << "<p class='note'>⚠ " << esc(w) << "</p>\n";
    // ---- kết quả từng bộ tải
    int sec = 2;
    for (const auto& s : r.sets) {
        o << fmt("<h2>%d. %s: ", sec++, s.combo ? tr("Tổ hợp") : tr("Trường hợp tải")) << esc(s.name) << "</h2>\n";
        o << trf("<p>Tổng tải ΣFx = %.3f, ΣFy = %.3f kN; tổng phản lực ΣRx = %.3f, ΣRy = %.3f kN.</p>\n", s.sumFx, s.sumFy, s.sumRx,
                 s.sumRy);
        if (nb) {
            if (opt.diagramM) {
                const std::string d = diagramSvg(m, s, BarQty::M);
                if (!d.empty()) o << tr("<h3>Biểu đồ mômen M (kNm)</h3>") << d;
            }
            if (opt.diagramQ) {
                const std::string d = diagramSvg(m, s, BarQty::Q);
                if (!d.empty()) o << tr("<h3>Biểu đồ lực cắt Q (kN)</h3>") << d;
            }
            if (opt.diagramN) {
                const std::string d = diagramSvg(m, s, BarQty::N);
                if (!d.empty()) o << tr("<h3>Biểu đồ lực dọc N (kN)</h3>") << d;
            }
        }
        if (opt.deformed) {
            const std::string d = deformedSvg(m, s);
            if (!d.empty()) o << tr("<h3>Hình biến dạng</h3>") << d;
        }
        if (!opt.tables) continue;
        o << tr("<h3>Phản lực gối</h3><table><tr><th>Nút</th><th>Rx (kN)</th><th>Ry (kN)</th><th>M (kNm)</th></tr>\n");
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            bool any = false;
            for (int c = 0; c < 3; ++c) any |= r.dofState[3 * i + c] == DofSupport || s.R[3 * i + c] != 0;
            if (any) o << fmt("<tr><td>%zu</td><td>%.3f</td><td>%.3f</td><td>%.3f</td></tr>\n", i + 1, s.R[3 * i], s.R[3 * i + 1], s.R[3 * i + 2]);
        }
        o << tr("</table>\n<h3>Chuyển vị nút</h3><table><tr><th>Nút</th><th>ux (mm)</th><th>uy (mm)</th><th>φz (rad)</th></tr>\n");
        for (size_t i = 0; i < m.nodes.size() && (int)i < opt.maxRows; ++i)
            o << fmt("<tr><td>%zu</td><td>%.4f</td><td>%.4f</td><td>%.3e</td></tr>\n", i + 1, s.U[3 * i] * 1e3, s.U[3 * i + 1] * 1e3,
                     s.U[3 * i + 2]);
        o << "</table>\n";
        if (nb) {
            o << tr("<h3>Nội lực đầu thanh</h3><table><tr><th>Thanh</th><th>Ni</th><th>Qi</th><th>Mi</th><th>Nj</th><th>Qj</th><th>Mj</th><th>M nhịp max</th></tr>\n");
            int rows = 0;
            for (size_t e = 0; e < m.elements.size() && rows < opt.maxRows; ++e) {
                if (!m.elements[e].isBar()) continue;
                const BarResult& b = s.bars[e];
                if (b.x.empty()) continue;
                const size_t n = b.x.size() - 1;
                size_t kx = 0;
                for (size_t k = 0; k <= n; ++k)
                    if (std::abs(b.M[k]) > std::abs(b.M[kx])) kx = k;
                o << fmt("<tr><td>%zu</td><td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%.2f</td></tr>\n",
                         e + 1, b.N[0], b.Q[0], b.M[0], b.N[n], b.Q[n], b.M[n], b.M[kx]);
                ++rows;
            }
            o << "</table>\n";
        }
        if (np) {
            const bool axi = m.planeMode == PlaneMode::Axisymmetric;
            o << (axi ? tr("<h3>Ứng suất tại nút (trung bình, MPa)</h3><table><tr><th>Nút</th><th>σr</th><th>σz</th><th>τrz</th><th>σθ</th>")
                      : tr("<h3>Ứng suất tại nút (trung bình, MPa)</h3><table><tr><th>Nút</th><th>σx</th><th>σy</th><th>τxy</th><th>σz</th>"))
              << "<th>σ1</th><th>σ2</th><th>σvM</th></tr>\n";
            int rows = 0;
            for (size_t i = 0; i < m.nodes.size() && rows < opt.maxRows; ++i) {
                if (s.nodeStressCount[i] == 0) continue;
                const auto& v = s.nodeStress[i];
                const StressInfo si = stressInfo(v[0], v[1], v[2], v[3]);
                o << fmt("<tr><td>%zu</td><td>%.3f</td><td>%.3f</td><td>%.3f</td><td>%.3f</td><td>%.3f</td><td>%.3f</td><td>%.3f</td></tr>\n",
                         i + 1, si.sx * 1e-3, si.sy * 1e-3, si.txy * 1e-3, si.sz * 1e-3, si.s1 * 1e-3, si.s2 * 1e-3,
                         si.vonMises * 1e-3);
                ++rows;
            }
            o << "</table>\n";
        }
    }
    o << "</body></html>\n";
    return o.str();
}

} // namespace fem
