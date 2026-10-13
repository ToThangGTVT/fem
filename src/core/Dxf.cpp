#include "core/Dxf.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>

#include "core/Mesher.h"

namespace fem {

namespace {

struct Writer {
    std::ostringstream o;
    double unit = 1;
    void pair(int code, const std::string& v) { o << code << "\n" << v << "\n"; }
    void pair(int code, double v) {
        char b[48];
        std::snprintf(b, sizeof b, "%.6f", v);
        pair(code, std::string(b));
    }
    void line(const char* layer, double x1, double y1, double x2, double y2) {
        pair(0, "LINE");
        pair(8, layer);
        pair(10, x1 * unit);
        pair(20, y1 * unit);
        pair(30, 0.0);
        pair(11, x2 * unit);
        pair(21, y2 * unit);
        pair(31, 0.0);
    }
    void text(const char* layer, double x, double y, double h, const std::string& s, int halign = 1, int valign = 2) {
        pair(0, "TEXT");
        pair(8, layer);
        pair(10, x * unit);
        pair(20, y * unit);
        pair(30, 0.0);
        pair(40, h * unit);
        pair(1, s);
        pair(72, std::to_string(halign));
        pair(11, x * unit);
        pair(21, y * unit);
        pair(31, 0.0);
        pair(73, std::to_string(valign));
    }
};

std::string fmtNum(double v) {
    char b[32];
    const double a = std::abs(v);
    std::snprintf(b, sizeof b, a >= 100 ? "%.1f" : a >= 1 ? "%.2f" : "%.3f", v);
    return b;
}

} // namespace

std::string exportDxf(const Model& m, const LoadSetResult* rs, const DxfExportOptions& opt) {
    Writer w;
    w.unit = opt.unit;
    double x0, y0, x1, y1;
    m.bounds(x0, y0, x1, y1);
    const double size = std::max({x1 - x0, y1 - y0, 1.0});
    const double th = opt.textHeight > 0 ? opt.textHeight : 0.015 * size;
    struct Layer {
        const char* name;
        int color;
    };
    const Layer layers[] = {{"FEM_THANH", 7},   {"FEM_DAN", 5},      {"FEM_PHANTU", 8},      {"FEM_GOI", 3},
                            {"FEM_BIEUDO_DUONG", 5}, {"FEM_BIEUDO_AM", 1}, {"FEM_GIATRI", 2}, {"FEM_BIENDANG", 6},
                            {"FEM_SOHIEU", 9}};
    w.pair(0, "SECTION");
    w.pair(2, "HEADER");
    w.pair(9, "$ACADVER");
    w.pair(1, "AC1009");
    w.pair(9, "$INSUNITS");
    w.pair(70, opt.unit == 1000 ? "4" : "6"); // 4: mm, 6: m
    w.pair(0, "ENDSEC");
    w.pair(0, "SECTION");
    w.pair(2, "TABLES");
    w.pair(0, "TABLE");
    w.pair(2, "LAYER");
    w.pair(70, std::to_string(std::size(layers)));
    for (const auto& l : layers) {
        w.pair(0, "LAYER");
        w.pair(2, l.name);
        w.pair(70, "0");
        w.pair(62, std::to_string(l.color));
        w.pair(6, "CONTINUOUS");
    }
    w.pair(0, "ENDTAB");
    w.pair(0, "ENDSEC");
    w.pair(0, "SECTION");
    w.pair(2, "ENTITIES");

    // ---- sơ đồ
    for (size_t e = 0; e < m.elements.size(); ++e) {
        const Element& el = m.elements[e];
        const int nc = el.nodeCount();
        if (el.isBar()) {
            const Node& a = m.nodes[el.n[0]];
            const Node& b = m.nodes[el.n[1]];
            w.line(el.type == ElemType::Truss ? "FEM_DAN" : "FEM_THANH", a.x, a.y, b.x, b.y);
        } else {
            for (int k = 0; k < nc; ++k) {
                const Node& a = m.nodes[el.n[k]];
                const Node& b = m.nodes[el.n[(k + 1) % nc]];
                w.line("FEM_PHANTU", a.x, a.y, b.x, b.y);
            }
        }
    }
    // gối tựa: tam giác (khớp) hoặc gạch ngàm
    const double gs = 0.025 * size;
    for (const auto& s : m.supports) {
        if (s.node < 0 || s.node >= (int)m.nodes.size() || s.empty()) continue;
        const Node& p = m.nodes[s.node];
        if (s.fix[0] && s.fix[1] && s.fix[2]) {
            w.line("FEM_GOI", p.x - gs, p.y, p.x + gs, p.y);
            for (int i = 0; i <= 4; ++i) {
                const double x = p.x - gs + 0.5 * gs * i;
                w.line("FEM_GOI", x, p.y, x - 0.4 * gs, p.y - 0.5 * gs);
            }
        } else {
            w.line("FEM_GOI", p.x, p.y, p.x - 0.6 * gs, p.y - gs);
            w.line("FEM_GOI", p.x, p.y, p.x + 0.6 * gs, p.y - gs);
            w.line("FEM_GOI", p.x - 0.6 * gs, p.y - gs, p.x + 0.6 * gs, p.y - gs);
            if (!s.fix[0]) w.line("FEM_GOI", p.x - 0.8 * gs, p.y - 1.3 * gs, p.x + 0.8 * gs, p.y - 1.3 * gs);
        }
    }
    if (!rs) {
        for (size_t i = 0; i < m.nodes.size(); ++i)
            w.text("FEM_SOHIEU", m.nodes[i].x + 0.3 * th, m.nodes[i].y + 0.3 * th, 0.7 * th, std::to_string(i + 1), 0, 1);
    }

    // ---- biểu đồ nội lực
    if (rs && opt.diagram && rs->bars.size() == m.elements.size()) {
        const double sc = opt.scale > 0 ? opt.scale : autoDiagramScale(m, *rs, opt.qty);
        double vmax = 0;
        for (size_t e = 0; e < m.elements.size(); ++e)
            if (m.elements[e].isBar())
                for (double v : barValues(rs->bars[e], opt.qty)) vmax = std::max(vmax, std::abs(v));
        for (size_t e = 0; e < m.elements.size(); ++e) {
            if (!m.elements[e].isBar()) continue;
            const DiagramGeom g = diagramGeometry(m, *rs, (int)e, opt.qty, sc);
            const size_t n = g.base.size();
            if (n < 2) continue;
            for (size_t k = 0; k + 1 < n; ++k) {
                const bool pos = g.values[k] + g.values[k + 1] >= 0;
                w.line(pos ? "FEM_BIEUDO_DUONG" : "FEM_BIEUDO_AM", g.off[k][0], g.off[k][1], g.off[k + 1][0], g.off[k + 1][1]);
            }
            // gạch tung độ
            const size_t step = std::max<size_t>(1, n / 12);
            for (size_t k = 0; k < n; k += step)
                w.line(g.values[k] >= 0 ? "FEM_BIEUDO_DUONG" : "FEM_BIEUDO_AM", g.base[k][0], g.base[k][1], g.off[k][0], g.off[k][1]);
            w.line(g.values[n - 1] >= 0 ? "FEM_BIEUDO_DUONG" : "FEM_BIEUDO_AM", g.base[n - 1][0], g.base[n - 1][1], g.off[n - 1][0],
                   g.off[n - 1][1]);
            if (opt.values) {
                auto put = [&](size_t k) {
                    if (std::abs(g.values[k]) < 1e-4 * vmax) return;
                    w.text("FEM_GIATRI", g.off[k][0], g.off[k][1], th, fmtNum(g.values[k]));
                };
                put(0);
                put(n - 1);
                size_t kx = 0;
                for (size_t k = 1; k + 1 < n; ++k)
                    if (std::abs(g.values[k]) > std::abs(g.values[kx])) kx = k;
                if (kx > 0 && kx + 1 < n && std::abs(g.values[kx]) > 1.02 * std::max(std::abs(g.values[0]), std::abs(g.values[n - 1])))
                    put(kx);
            }
        }
    }
    // ---- hình biến dạng
    if (rs && opt.deformed && rs->U.size() == 3 * m.nodes.size()) {
        double umax = 0;
        for (size_t i = 0; i < m.nodes.size(); ++i) umax = std::max(umax, std::hypot(rs->U[3 * i], rs->U[3 * i + 1]));
        const double S = opt.deformScale > 0 ? opt.deformScale : (umax > 0 ? 0.08 * size / umax : 1);
        for (size_t e = 0; e < m.elements.size(); ++e) {
            const Element& el = m.elements[e];
            if (el.isBar() && e < rs->bars.size()) {
                const BarResult& b = rs->bars[e];
                const Node& a = m.nodes[el.n[0]];
                const Node& q = m.nodes[el.n[1]];
                const double L = std::hypot(q.x - a.x, q.y - a.y), c = (q.x - a.x) / L, s = (q.y - a.y) / L;
                for (size_t k = 0; k + 1 < b.x.size(); ++k) {
                    auto P = [&](size_t i, double& X, double& Y) {
                        X = a.x + c * b.x[i] + (c * b.u[i] - s * b.v[i]) * S;
                        Y = a.y + s * b.x[i] + (s * b.u[i] + c * b.v[i]) * S;
                    };
                    double X0, Y0, X1, Y1;
                    P(k, X0, Y0);
                    P(k + 1, X1, Y1);
                    w.line("FEM_BIENDANG", X0, Y0, X1, Y1);
                }
            } else if (el.isPlane()) {
                const int nc = el.nodeCount();
                for (int k = 0; k < nc; ++k) {
                    const int i = el.n[k], j = el.n[(k + 1) % nc];
                    w.line("FEM_BIENDANG", m.nodes[i].x + rs->U[3 * i] * S, m.nodes[i].y + rs->U[3 * i + 1] * S,
                           m.nodes[j].x + rs->U[3 * j] * S, m.nodes[j].y + rs->U[3 * j + 1] * S);
                }
            }
        }
    }
    w.pair(0, "ENDSEC");
    w.pair(0, "EOF");
    return w.o.str();
}

// ============================================================================ nhập
namespace {

struct Pair {
    int code;
    std::string value;
};

std::vector<Pair> readPairs(const std::string& text) {
    std::vector<Pair> out;
    std::istringstream in(text);
    std::string a, b;
    auto trim = [](std::string& s) {
        while (!s.empty() && (s.back() == '\r' || s.back() == ' ' || s.back() == '\t')) s.pop_back();
        size_t i = 0;
        while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
        s = s.substr(i);
    };
    while (std::getline(in, a) && std::getline(in, b)) {
        trim(a);
        trim(b);
        out.push_back({std::atoi(a.c_str()), b});
    }
    return out;
}

struct Seg {
    std::string layer;
    double x1, y1, x2, y2;
};

std::vector<Seg> readSegments(const std::string& text) {
    const std::vector<Pair> p = readPairs(text);
    std::vector<Seg> segs;
    size_t i = 0;
    auto num = [](const std::string& s) { return std::strtod(s.c_str(), nullptr); };
    while (i < p.size()) {
        if (p[i].code != 0) {
            ++i;
            continue;
        }
        const std::string type = p[i].value;
        ++i;
        if (type == "LINE") {
            Seg s{"", 0, 0, 0, 0};
            for (; i < p.size() && p[i].code != 0; ++i) {
                switch (p[i].code) {
                case 8: s.layer = p[i].value; break;
                case 10: s.x1 = num(p[i].value); break;
                case 20: s.y1 = num(p[i].value); break;
                case 11: s.x2 = num(p[i].value); break;
                case 21: s.y2 = num(p[i].value); break;
                default: break;
                }
            }
            segs.push_back(s);
        } else if (type == "LWPOLYLINE") {
            std::string layer;
            bool closed = false;
            std::vector<std::array<double, 2>> v;
            for (; i < p.size() && p[i].code != 0; ++i) {
                if (p[i].code == 8) layer = p[i].value;
                else if (p[i].code == 70) closed = (std::atoi(p[i].value.c_str()) & 1) != 0;
                else if (p[i].code == 10) v.push_back({num(p[i].value), 0});
                else if (p[i].code == 20 && !v.empty()) v.back()[1] = num(p[i].value);
            }
            for (size_t k = 0; k + 1 < v.size(); ++k) segs.push_back({layer, v[k][0], v[k][1], v[k + 1][0], v[k + 1][1]});
            if (closed && v.size() > 2) segs.push_back({layer, v.back()[0], v.back()[1], v[0][0], v[0][1]});
        } else if (type == "POLYLINE") {
            std::string layer;
            bool closed = false;
            for (; i < p.size() && p[i].code != 0; ++i) {
                if (p[i].code == 8) layer = p[i].value;
                else if (p[i].code == 70) closed = (std::atoi(p[i].value.c_str()) & 1) != 0;
            }
            std::vector<std::array<double, 2>> v;
            while (i < p.size() && p[i].code == 0 && p[i].value == "VERTEX") {
                ++i;
                std::array<double, 2> q{0, 0};
                for (; i < p.size() && p[i].code != 0; ++i) {
                    if (p[i].code == 10) q[0] = num(p[i].value);
                    else if (p[i].code == 20) q[1] = num(p[i].value);
                }
                v.push_back(q);
            }
            for (size_t k = 0; k + 1 < v.size(); ++k) segs.push_back({layer, v[k][0], v[k][1], v[k + 1][0], v[k + 1][1]});
            if (closed && v.size() > 2) segs.push_back({layer, v.back()[0], v.back()[1], v[0][0], v[0][1]});
        }
    }
    return segs;
}

} // namespace

std::vector<std::string> dxfLayers(const std::string& text) {
    std::set<std::string> s;
    for (const auto& g : readSegments(text)) s.insert(g.layer);
    return std::vector<std::string>(s.begin(), s.end());
}

int splitBarsAtNodes(Model& m, double tol) {
    int splits = 0;
    for (size_t node = 0; node < m.nodes.size(); ++node) {
        const Node P = m.nodes[node];
        for (size_t e = 0; e < m.elements.size(); ++e) {
            const Element& el = m.elements[e];
            if (!el.isBar() || el.n[0] == (int)node || el.n[1] == (int)node) continue;
            const Node& A = m.nodes[el.n[0]];
            const Node& B = m.nodes[el.n[1]];
            const double dx = B.x - A.x, dy = B.y - A.y, L2 = dx * dx + dy * dy;
            if (L2 <= 0) continue;
            const double t = ((P.x - A.x) * dx + (P.y - A.y) * dy) / L2;
            if (t <= 1e-9 || t >= 1 - 1e-9) continue;
            if (std::hypot(A.x + t * dx - P.x, A.y + t * dy - P.y) > tol) continue;
            const int nn = splitBarAt(m, (int)e, t);
            // nút mới trùng với nút đang xét: nối thẳng vào nút đó
            for (auto& x : m.elements)
                for (int k = 0; k < x.nodeCount(); ++k)
                    if (x.n[k] == nn) x.n[k] = (int)node;
            m.nodes.pop_back();
            ++splits;
        }
    }
    return splits;
}

DxfImportResult importDxf(const std::string& text, Model& m, const DxfImportOptions& o) {
    DxfImportResult r;
    const std::vector<Seg> segs = readSegments(text);
    std::set<std::string> layers;
    for (const auto& s : segs) layers.insert(s.layer);
    r.layers.assign(layers.begin(), layers.end());
    if (segs.empty()) {
        r.error = "Bản vẽ không có đoạn thẳng (LINE / LWPOLYLINE / POLYLINE).";
        return r;
    }
    if (m.materials.empty()) m.materials.push_back(Material{});
    if (m.sections.empty()) m.sections.push_back(Section{});
    double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
    for (const auto& s : segs) {
        x0 = std::min({x0, s.x1, s.x2});
        x1 = std::max({x1, s.x1, s.x2});
        y0 = std::min({y0, s.y1, s.y2});
        y1 = std::max({y1, s.y1, s.y2});
    }
    const double tol = 1e-6 * std::max({1.0, (x1 - x0) * o.unit, (y1 - y0) * o.unit});
    const size_t n0 = m.nodes.size(), e0 = m.elements.size();
    for (const auto& s : segs) {
        if (!o.layer.empty() && s.layer != o.layer) continue;
        const double ax = s.x1 * o.unit, ay = s.y1 * o.unit, bx = s.x2 * o.unit, by = s.y2 * o.unit;
        if (std::hypot(bx - ax, by - ay) <= tol) continue;
        ++r.lines;
        Element el;
        el.type = o.type;
        el.n = {m.addNode(ax, ay, tol), m.addNode(bx, by, tol), -1, -1};
        el.mat = std::clamp(o.mat, 0, (int)m.materials.size() - 1);
        el.sec = std::clamp(o.sec, 0, (int)m.sections.size() - 1);
        if (m.findElement(el.type, el.n) < 0) m.elements.push_back(el);
    }
    r.splits = splitBarsAtNodes(m, tol);
    mergeNodes(m, tol);
    r.nodes = (int)(m.nodes.size() - n0);
    r.elements = (int)(m.elements.size() - e0);
    r.ok = r.lines > 0;
    if (!r.ok) r.error = "Không có đoạn thẳng nào trên lớp \"" + o.layer + "\".";
    return r;
}

} // namespace fem
