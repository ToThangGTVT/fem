#include "core/Mesher.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <set>

namespace fem {

std::vector<int> meshPatch(Model& m, const double cx[4], const double cy[4], int nx, int ny, ElemType type, int mat,
                           double t, bool incompatible) {
    nx = std::max(nx, 1);
    ny = std::max(ny, 1);
    double scale = 0;
    for (int i = 0; i < 4; ++i) scale = std::max({scale, std::abs(cx[i]), std::abs(cy[i])});
    const double tol = 1e-9 * std::max(scale, 1.0);
    std::vector<int> grid((nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            const double u = (double)i / nx, v = (double)j / ny;
            const double N[4] = {(1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v};
            double x = 0, y = 0;
            for (int k = 0; k < 4; ++k) {
                x += N[k] * cx[k];
                y += N[k] * cy[k];
            }
            grid[j * (nx + 1) + i] = m.addNode(x, y, tol);
        }
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int a = grid[j * (nx + 1) + i], b = grid[j * (nx + 1) + i + 1];
            const int c = grid[(j + 1) * (nx + 1) + i + 1], d = grid[(j + 1) * (nx + 1) + i];
            Element e;
            e.mat = mat;
            e.t = t;
            e.incompatible = incompatible;
            if (type == ElemType::Quad4) {
                e.type = ElemType::Quad4;
                e.n = {a, b, c, d};
                m.elements.push_back(e);
            } else {
                e.type = ElemType::Tri3;
                // đường chéo xen kẽ để lưới đối xứng hơn
                if ((i + j) % 2 == 0) {
                    e.n = {a, b, c, -1};
                    m.elements.push_back(e);
                    e.n = {a, c, d, -1};
                    m.elements.push_back(e);
                } else {
                    e.n = {a, b, d, -1};
                    m.elements.push_back(e);
                    e.n = {b, c, d, -1};
                    m.elements.push_back(e);
                }
            }
        }
    return grid;
}

std::vector<int> meshRectangle(Model& m, double x0, double y0, double w, double h, int nx, int ny, ElemType type,
                               int mat, double t, bool incompatible) {
    const double cx[4] = {x0, x0 + w, x0 + w, x0};
    const double cy[4] = {y0, y0, y0 + h, y0 + h};
    return meshPatch(m, cx, cy, nx, ny, type, mat, t, incompatible);
}

void subdivideBar(Model& m, int e, int n) {
    if (n < 2 || e < 0 || e >= (int)m.elements.size() || !m.elements[e].isBar()) return;
    const Element orig = m.elements[e];
    const Node a = m.nodes[orig.n[0]], b = m.nodes[orig.n[1]];
    const double L = std::hypot(b.x - a.x, b.y - a.y);
    std::vector<int> nodes = {orig.n[0]};
    for (int k = 1; k < n; ++k) {
        m.nodes.push_back({a.x + (b.x - a.x) * k / n, a.y + (b.y - a.y) * k / n});
        nodes.push_back((int)m.nodes.size() - 1);
    }
    nodes.push_back(orig.n[1]);
    std::vector<int> ids = {e};
    for (int k = 0; k < n; ++k) {
        Element piece = orig;
        piece.n[0] = nodes[k];
        piece.n[1] = nodes[k + 1];
        piece.hingeI = (k == 0) && orig.hingeI;
        piece.hingeJ = (k == n - 1) && orig.hingeJ;
        if (k == 0) m.elements[e] = piece;
        else {
            m.elements.push_back(piece);
            ids.push_back((int)m.elements.size() - 1);
        }
    }
    const double seg = L / n;
    for (auto& lc : m.cases) {
        std::vector<ElemLoad> added;
        for (auto& l : lc.elem) {
            if (l.elem != e) continue;
            switch (l.kind) {
            case ElemLoadKind::Distributed: {
                const double q1 = l.q1, q2 = l.q2;
                for (int k = 0; k < n; ++k) {
                    ElemLoad p = l;
                    p.elem = ids[k];
                    p.q1 = q1 + (q2 - q1) * k / n;
                    p.q2 = q1 + (q2 - q1) * (k + 1) / n;
                    if (k == 0) l = p;
                    else added.push_back(p);
                }
                break;
            }
            case ElemLoadKind::Point:
            case ElemLoadKind::Moment: {
                const int k = std::clamp((int)std::floor(l.a / seg), 0, n - 1);
                l.elem = ids[k];
                l.a -= k * seg;
                break;
            }
            case ElemLoadKind::Temperature:
                for (int k = 1; k < n; ++k) {
                    ElemLoad p = l;
                    p.elem = ids[k];
                    added.push_back(p);
                }
                break;
            case ElemLoadKind::LackOfFit: // ΔL chia đều cho các đoạn
                l.q1 /= n;
                for (int k = 1; k < n; ++k) {
                    ElemLoad p = l;
                    p.elem = ids[k];
                    added.push_back(p);
                }
                break;
            }
        }
        lc.elem.insert(lc.elem.end(), added.begin(), added.end());
    }
}

int applyEdgeLoad(Model& m, int caseIndex, const std::vector<int>& nodes, const EdgeLoadSpec& spec) {
    if (caseIndex < 0 || caseIndex >= (int)m.cases.size()) return 0;
    const std::set<int> sel(nodes.begin(), nodes.end());
    struct Edge {
        int a, b;      // theo chiều ngược kim đồng hồ của phần tử
        double t;
        int count = 0;
    };
    std::map<std::pair<int, int>, Edge> edges;
    for (const auto& el : m.elements) {
        if (!el.isPlane()) continue;
        const int nc = el.nodeCount();
        double area2 = 0;
        for (int k = 0; k < nc; ++k) {
            const Node& p = m.nodes[el.n[k]];
            const Node& q = m.nodes[el.n[(k + 1) % nc]];
            area2 += p.x * q.y - q.x * p.y;
        }
        for (int k = 0; k < nc; ++k) {
            int a = el.n[k], b = el.n[(k + 1) % nc];
            if (area2 < 0) std::swap(a, b);
            auto key = std::make_pair(std::min(a, b), std::max(a, b));
            Edge& e = edges[key];
            e.a = a;
            e.b = b;
            e.t = el.t;
            e.count++;
        }
    }
    std::map<int, std::array<double, 2>> F;
    int used = 0;
    for (const auto& [key, e] : edges) {
        if (e.count != 1 || !sel.count(e.a) || !sel.count(e.b)) continue; // chỉ cạnh biên
        const Node& A = m.nodes[e.a];
        const Node& B = m.nodes[e.b];
        const double dx = B.x - A.x, dy = B.y - A.y, L = std::hypot(dx, dy);
        if (L <= 0) continue;
        ++used;
        if (!spec.pressure) {
            for (int n : {e.a, e.b}) {
                F[n][0] += spec.qx * L / 2;
                F[n][1] += spec.qy * L / 2;
            }
        } else {
            // pháp tuyến ngoài của cạnh (phần tử ngược chiều kim đồng hồ): (dy, −dx)/L
            const double nx = dy / L, ny = -dx / L;
            double pa = spec.p0 + spec.ky * A.y, pb = spec.p0 + spec.ky * B.y;
            if (spec.clipNegative) {
                pa = std::max(pa, 0.0);
                pb = std::max(pb, 0.0);
            }
            const double fa = L * (2 * pa + pb) / 6 * e.t, fb = L * (pa + 2 * pb) / 6 * e.t;
            F[e.a][0] -= fa * nx;
            F[e.a][1] -= fa * ny;
            F[e.b][0] -= fb * nx;
            F[e.b][1] -= fb * ny;
        }
    }
    for (const auto& [n, f] : F)
        if (f[0] != 0 || f[1] != 0) m.cases[caseIndex].nodal.push_back({n, {f[0], f[1], 0}});
    return used;
}

int splitBarAt(Model& m, int e, double t) {
    if (e < 0 || e >= (int)m.elements.size() || !m.elements[e].isBar() || t <= 0 || t >= 1) return -1;
    const Element orig = m.elements[e];
    const Node a = m.nodes[orig.n[0]], b = m.nodes[orig.n[1]];
    const double L = std::hypot(b.x - a.x, b.y - a.y), cut = t * L;
    m.nodes.push_back({a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t});
    const int nn = (int)m.nodes.size() - 1;
    Element p1 = orig, p2 = orig;
    p1.n[1] = nn;
    p1.hingeJ = false;
    p2.n[0] = nn;
    p2.hingeI = false;
    m.elements[e] = p1;
    m.elements.push_back(p2);
    const int e2 = (int)m.elements.size() - 1;
    for (auto& lc : m.cases) {
        std::vector<ElemLoad> added;
        for (auto& l : lc.elem) {
            if (l.elem != e) continue;
            switch (l.kind) {
            case ElemLoadKind::Distributed: {
                const double qt = l.q1 + (l.q2 - l.q1) * t;
                ElemLoad r = l;
                r.elem = e2;
                r.q1 = qt;
                l.q2 = qt;
                added.push_back(r);
                break;
            }
            case ElemLoadKind::Point:
            case ElemLoadKind::Moment:
                if (l.a > cut) {
                    l.elem = e2;
                    l.a -= cut;
                }
                break;
            case ElemLoadKind::Temperature: {
                ElemLoad r = l;
                r.elem = e2;
                added.push_back(r);
                break;
            }
            case ElemLoadKind::LackOfFit: {
                ElemLoad r = l;
                r.elem = e2;
                r.q1 = l.q1 * (1 - t);
                l.q1 *= t;
                added.push_back(r);
                break;
            }
            }
        }
        lc.elem.insert(lc.elem.end(), added.begin(), added.end());
    }
    return nn;
}

int mergeNodes(Model& m, double tol) {
    const int nn = (int)m.nodes.size();
    std::vector<int> order(nn), rep(nn);
    for (int i = 0; i < nn; ++i) order[i] = rep[i] = i;
    std::function<int(int)> find = [&](int i) { return rep[i] == i ? i : (rep[i] = find(rep[i])); };
    std::sort(order.begin(), order.end(), [&](int a, int b) { return m.nodes[a].x < m.nodes[b].x; });
    for (int ia = 0; ia < nn; ++ia) {
        const int a = order[ia];
        for (int ib = ia + 1; ib < nn; ++ib) {
            const int b = order[ib];
            if (m.nodes[b].x - m.nodes[a].x > tol) break;
            if (std::abs(m.nodes[b].y - m.nodes[a].y) > tol) continue;
            const int ra = find(a), rb = find(b);
            if (ra != rb) rep[std::max(ra, rb)] = std::min(ra, rb); // giữ nút có chỉ số nhỏ hơn
        }
    }
    int merged = 0;
    for (int i = 0; i < nn; ++i) {
        rep[i] = find(i);
        if (rep[i] != i) ++merged;
    }
    if (!merged) return 0;
    for (auto& el : m.elements)
        for (int k = 0; k < el.nodeCount(); ++k) el.n[k] = rep[el.n[k]];
    // gối: gộp điều kiện biên vào nút giữ lại
    for (auto& s : m.supports) {
        const int r = rep[s.node];
        if (r == s.node) continue;
        Support& t = m.supportFor(r);
        for (int c = 0; c < 3; ++c) {
            t.fix[c] = t.fix[c] || s.fix[c];
            t.k[c] = std::max(t.k[c], s.k[c]);
        }
    }
    for (auto& lc : m.cases) {
        for (auto& l : lc.nodal) l.node = rep[l.node];
        for (auto& d : lc.disp) d.node = rep[d.node];
    }
    // phần tử suy biến (nút lặp) hoặc trùng lặp
    std::vector<int> dead;
    std::set<std::vector<int>> seen;
    for (size_t e = 0; e < m.elements.size(); ++e) {
        const Element& el = m.elements[e];
        std::vector<int> key(el.n.begin(), el.n.begin() + el.nodeCount());
        std::sort(key.begin(), key.end());
        const bool degenerate = std::adjacent_find(key.begin(), key.end()) != key.end();
        key.push_back((int)el.type);
        if (degenerate || !seen.insert(key).second) dead.push_back((int)e);
    }
    m.removeElements(dead);
    std::vector<int> drop;
    for (int i = 0; i < nn; ++i)
        if (rep[i] != i) drop.push_back(i);
    m.supports.erase(std::remove_if(m.supports.begin(), m.supports.end(),
                                    [&](const Support& s) { return rep[s.node] != s.node; }),
                     m.supports.end());
    m.removeNodes(drop);
    return merged;
}

} // namespace fem

