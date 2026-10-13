#include "core/Mesher.h"

#include <cstdio>

#include "core/SectionGeom.h"

#include <algorithm>
#include <array>
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
    // nút trên lưới nửa bước (I, J) ∈ [0, 2nx]×[0, 2ny]; nút góc phần tử ở chỉ số chẵn, nút giữa cạnh có chỉ số lẻ
    std::map<std::pair<int, int>, int> made;
    auto node = [&](int I, int J) {
        auto it = made.find({I, J});
        if (it != made.end()) return it->second;
        const double u = 0.5 * I / nx, v = 0.5 * J / ny;
        const double N[4] = {(1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v};
        double x = 0, y = 0;
        for (int k = 0; k < 4; ++k) {
            x += N[k] * cx[k];
            y += N[k] * cy[k];
        }
        const int id = m.addNode(x, y, tol);
        made[{I, J}] = id;
        return id;
    };
    std::vector<int> grid((nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) grid[j * (nx + 1) + i] = node(2 * i, 2 * j);
    const bool quad = type == ElemType::Quad4 || type == ElemType::Quad8;
    const bool second = type == ElemType::Tri6 || type == ElemType::Quad8;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            // góc theo tọa độ nửa bước: a(0,0) b(2,0) c(2,2) d(0,2)
            const int I = 2 * i, J = 2 * j;
            const std::array<std::pair<int, int>, 4> C = {{{I, J}, {I + 2, J}, {I + 2, J + 2}, {I, J + 2}}};
            Element e;
            e.mat = mat;
            e.t = t;
            e.incompatible = incompatible;
            e.n.fill(-1);
            auto add = [&](std::initializer_list<int> corners) {
                std::vector<std::pair<int, int>> cs;
                for (int k : corners) cs.push_back(C[k]);
                const int nc = (int)cs.size();
                for (int k = 0; k < nc; ++k) e.n[k] = node(cs[k].first, cs[k].second);
                if (second)
                    for (int k = 0; k < nc; ++k) {
                        const auto& p = cs[k];
                        const auto& q = cs[(k + 1) % nc];
                        e.n[nc + k] = node((p.first + q.first) / 2, (p.second + q.second) / 2);
                    }
                m.elements.push_back(e);
            };
            if (quad) {
                e.type = type;
                add({0, 1, 2, 3});
            } else {
                e.type = second ? ElemType::Tri6 : ElemType::Tri3;
                // đường chéo xen kẽ để lưới đối xứng hơn
                if ((i + j) % 2 == 0) {
                    add({0, 1, 2});
                    add({0, 2, 3});
                } else {
                    add({0, 1, 3});
                    add({1, 2, 3});
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
    const bool axi = m.planeMode == PlaneMode::Axisymmetric;
    struct Edge {
        int a, b, mid = -1; // theo chiều ngược kim đồng hồ của phần tử; mid: nút giữa cạnh (bậc hai)
        double t;
        int count = 0;
    };
    std::map<std::pair<int, int>, Edge> edges;
    for (const auto& el : m.elements) {
        if (!el.isPlane()) continue;
        const int nc = el.cornerCount();
        const bool second = el.nodeCount() > nc;
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
            e.mid = second ? el.n[nc + k] : -1;
            e.t = el.t;
            e.count++;
        }
    }
    std::map<int, std::array<double, 2>> F;
    int used = 0;
    // Gauss 3 điểm trên [0, 1]
    const double gs[3] = {0.5 - 0.5 * std::sqrt(0.6), 0.5, 0.5 + 0.5 * std::sqrt(0.6)};
    const double gw[3] = {5.0 / 18, 8.0 / 18, 5.0 / 18};
    for (const auto& [key, e] : edges) {
        if (e.count != 1 || !sel.count(e.a) || !sel.count(e.b)) continue; // chỉ cạnh biên
        if (e.mid >= 0 && !sel.count(e.mid)) continue;
        const Node& A = m.nodes[e.a];
        const Node& B = m.nodes[e.b];
        if (std::hypot(B.x - A.x, B.y - A.y) <= 0) continue;
        ++used;
        const Node M = e.mid >= 0 ? m.nodes[e.mid] : Node{(A.x + B.x) / 2, (A.y + B.y) / 2};
        const int ids[3] = {e.a, e.b, e.mid};
        const int nn = e.mid >= 0 ? 3 : 2;
        for (int g = 0; g < 3; ++g) {
            const double s = gs[g];
            // hàm dạng trên cạnh: tuyến tính (A, B) hoặc bậc hai (A, B, M)
            double N[3], dN[3];
            if (nn == 2) {
                N[0] = 1 - s, N[1] = s, dN[0] = -1, dN[1] = 1;
            } else {
                N[0] = (1 - s) * (1 - 2 * s), N[1] = s * (2 * s - 1), N[2] = 4 * s * (1 - s);
                dN[0] = 4 * s - 3, dN[1] = 4 * s - 1, dN[2] = 4 - 8 * s;
            }
            const double px[3] = {A.x, B.x, M.x}, py[3] = {A.y, B.y, M.y};
            double x = 0, y = 0, tx = 0, ty = 0;
            for (int k = 0; k < nn; ++k) {
                x += N[k] * px[k];
                y += N[k] * py[k];
                tx += dN[k] * px[k];
                ty += dN[k] * py[k];
            }
            const double J = std::hypot(tx, ty);
            if (J <= 0) continue;
            const double w = gw[g] * J * (axi ? 2 * 3.14159265358979323846 * std::max(x, 0.0) : 1.0);
            double fx, fy;
            if (!spec.pressure) {
                fx = spec.qx, fy = spec.qy;
            } else if (!axi) {
                double p = spec.p0 + spec.ky * y;
                if (spec.clipNegative) p = std::max(p, 0.0);
                fx = -p * e.t * ty / J, fy = p * e.t * tx / J;
            } else {
                // pháp tuyến ngoài (phần tử ngược chiều kim đồng hồ): (ty, −tx)/J; áp lực dương nén vào vật thể
                double p = spec.p0 + spec.ky * y;
                if (spec.clipNegative) p = std::max(p, 0.0);
                fx = -p * ty / J, fy = p * tx / J;
            }
            for (int k = 0; k < nn; ++k) {
                F[ids[k]][0] += N[k] * fx * w;
                F[ids[k]][1] += N[k] * fy * w;
            }
        }
    }
    // bỏ phần dư do làm tròn (ví dụ nút trên trục đối xứng: lực đúng bằng 0)
    double fmax = 0;
    for (const auto& [n, f] : F) fmax = std::max({fmax, std::abs(f[0]), std::abs(f[1])});
    for (const auto& [n, f] : F)
        if (std::abs(f[0]) > 1e-12 * fmax || std::abs(f[1]) > 1e-12 * fmax)
            m.cases[caseIndex].nodal.push_back({n, {f[0], f[1], 0}});
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

int taperBar(Model& m, int e, int secJ, int n) {
    if (e < 0 || e >= (int)m.elements.size() || !m.elements[e].isBar() || secJ < 0 || secJ >= (int)m.sections.size() || n < 1)
        return 0;
    const Section si = m.sections[m.elements[e].sec], sj = m.sections[secJ];
    const bool sameShape = si.shape.kind == sj.shape.kind && si.shape.kind != ShapeKind::None &&
                           si.shape.kind != ShapeKind::Polygon;
    const int firstElem = (int)m.elements.size();
    if (n > 1) subdivideBar(m, e, n);
    std::vector<int> pieces = {e};
    for (int k = firstElem; k < (int)m.elements.size(); ++k) pieces.push_back(k);
    for (int k = 0; k < n; ++k) {
        const double t = (k + 0.5) / n; // nội suy tại giữa đoạn
        Section s;
        if (sameShape) {
            SectionShape sh = si.shape;
            for (int i = 0; i < 6; ++i) sh.p[i] = si.shape.p[i] + t * (sj.shape.p[i] - si.shape.p[i]);
            s = sectionFromShape("", sh);
        } else {
            s.A = si.A + t * (sj.A - si.A);
            s.I = si.I + t * (sj.I - si.I);
            s.h = si.h + t * (sj.h - si.h);
        }
        char b[160];
        std::snprintf(b, sizeof b, "%s→%s %d/%d", si.name.c_str(), sj.name.c_str(), k + 1, n);
        s.name = b;
        m.sections.push_back(s);
        m.elements[pieces[k]].sec = (int)m.sections.size() - 1;
    }
    return n;
}

std::vector<int> makeArch(Model& m, double x0, double y0, double L, double f, int n, int shape, int supports, int mat, int sec) {
    n = std::max(2, n + (n % 2)); // số đoạn chẵn để có nút đỉnh vòm
    std::vector<int> nodes;
    const double R = (L * L / 4 + f * f) / (2 * f); // bán kính cung tròn
    for (int i = 0; i <= n; ++i) {
        const double x = L * i / n;
        double y;
        if (shape == 0) y = 4 * f * x * (L - x) / (L * L);
        else y = std::sqrt(std::max(0.0, R * R - (x - L / 2) * (x - L / 2))) - (R - f);
        nodes.push_back(m.addNode(x0 + x, y0 + y));
    }
    for (int i = 0; i < n; ++i) {
        Element el;
        el.type = ElemType::Frame;
        el.n = {nodes[i], nodes[i + 1], -1, -1};
        el.mat = mat;
        el.sec = sec;
        if (supports == 1 && i == n / 2 - 1) el.hingeJ = true; // khớp đỉnh vòm
        m.elements.push_back(el);
    }
    const bool fixRot = supports == 2;
    m.supportFor(nodes.front()).fix = {true, true, fixRot};
    m.supportFor(nodes.back()).fix = {true, true, fixRot};
    return nodes;
}

} // namespace fem

