#include "core/PostProcess.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace fem {

const std::vector<double>& barValues(const BarResult& r, BarQty q) {
    switch (q) {
    case BarQty::N: return r.N;
    case BarQty::Q: return r.Q;
    case BarQty::M: return r.M;
    }
    return r.M;
}

double barValueAt(const BarResult& r, BarQty q, double x, bool right) {
    const std::vector<double>& v = barValues(r, q);
    const size_t n = r.x.size();
    if (n == 0) return 0;
    const double eps = 1e-9 * std::max(1.0, r.x.back());
    size_t k = std::lower_bound(r.x.begin(), r.x.end(), x - eps) - r.x.begin();
    if (k >= n) return v[n - 1];
    if (std::abs(r.x[k] - x) <= eps) {
        if (right && k + 1 < n && std::abs(r.x[k + 1] - r.x[k]) <= eps) return v[k + 1];
        return v[k];
    }
    if (k == 0) return v[0];
    const double t = (x - r.x[k - 1]) / (r.x[k] - r.x[k - 1]);
    return v[k - 1] + t * (v[k] - v[k - 1]);
}

std::vector<BarEnvelope> envelope(const Model& m, const AnalysisResult& r, const std::vector<int>& sets) {
    std::vector<BarEnvelope> out(m.elements.size());
    if (!r.ok || sets.empty()) return out;
    for (size_t e = 0; e < m.elements.size(); ++e) {
        if (!m.elements[e].isBar()) continue;
        // các điểm: hợp các điểm mẫu của mọi bộ kết quả, giữ cả hai phía bước nhảy
        std::vector<std::pair<double, bool>> pts;
        for (int s : sets) {
            if (s < 0 || s >= (int)r.sets.size() || r.sets[s].bars.size() <= e) continue;
            const auto& x = r.sets[s].bars[e].x;
            for (size_t k = 0; k < x.size(); ++k) {
                const bool left = k + 1 < x.size() && x[k + 1] == x[k];
                pts.push_back({x[k], !left});
            }
        }
        std::sort(pts.begin(), pts.end());
        // gộp các điểm quá gần nhau
        std::vector<std::pair<double, bool>> u;
        for (const auto& p : pts)
            if (u.empty() || p.first - u.back().first > 1e-9 || p.second != u.back().second) u.push_back(p);
        BarEnvelope& env = out[e];
        env.x.resize(u.size());
        for (int q = 0; q < 3; ++q) {
            env.max[q].assign(u.size(), -1e300);
            env.min[q].assign(u.size(), 1e300);
        }
        for (size_t k = 0; k < u.size(); ++k) {
            env.x[k] = u[k].first;
            for (int s : sets) {
                if (s < 0 || s >= (int)r.sets.size() || r.sets[s].bars.size() <= e) continue;
                for (int q = 0; q < 3; ++q) {
                    const double v = barValueAt(r.sets[s].bars[e], (BarQty)q, u[k].first, u[k].second);
                    env.max[q][k] = std::max(env.max[q][k], v);
                    env.min[q][k] = std::min(env.min[q][k], v);
                }
            }
        }
    }
    return out;
}

InfluenceLine influenceLine(const Model& m, const std::vector<int>& path, const IlTarget& t, int pointsPerBar) {
    InfluenceLine il;
    const int nn = (int)m.nodes.size(), ne = (int)m.elements.size();
    if (path.empty()) {
        il.error = "Chưa chọn đường chạy tải (các thanh mà lực P = 1 di chuyển trên đó).";
        return il;
    }
    for (int e : path)
        if (e < 0 || e >= ne || !m.elements[e].isBar()) {
            il.error = "Đường chạy tải chỉ gồm các thanh.";
            return il;
        }
    if (t.kind == IlKind::Force) {
        if (t.elem < 0 || t.elem >= ne || !m.elements[t.elem].isBar()) {
            il.error = "Tiết diện khảo sát phải nằm trên một thanh.";
            return il;
        }
    } else if (t.node < 0 || t.node >= nn || t.dof < 0 || t.dof > 2) {
        il.error = "Nút khảo sát không hợp lệ.";
        return il;
    }

    const double tx = t.kind == IlKind::Force ? std::clamp(t.x, 0.0, m.barLength(t.elem)) : 0.0;
    Model mc = m;
    mc.cases.clear();
    mc.combos.clear();
    pointsPerBar = std::max(pointsPerBar, 2);
    // Mỗi vị trí lực là một trường hợp tải. Tại tiết diện khảo sát nằm trên đường chạy, nội lực có bước nhảy
    // nên tạo hai điểm: lực ngay bên trái (tiết diện ở bên phải lực) và ngay bên phải tiết diện.
    struct Pos {
        int elem;
        double a;
        int side; // 0: thường, 1: lực ở bên trái tiết diện, 2: lực ở bên phải tiết diện
    };
    std::vector<Pos> pos;
    for (int e : path) {
        const double L = m.barLength(e);
        std::vector<double> as;
        for (int k = 0; k <= pointsPerBar; ++k) as.push_back(L * k / pointsPerBar);
        const bool cut = t.kind == IlKind::Force && e == t.elem;
        if (cut) as.push_back(tx);
        std::sort(as.begin(), as.end());
        as.erase(std::unique(as.begin(), as.end(), [&](double u, double v) { return std::abs(u - v) < 1e-9 * L; }), as.end());
        for (double a : as) {
            if (cut && std::abs(a - tx) < 1e-9 * L) {
                pos.push_back({e, tx, 1});
                pos.push_back({e, tx, 2});
            } else {
                pos.push_back({e, a, 0});
            }
        }
    }
    std::vector<int> caseOf;
    for (size_t i = 0; i < pos.size(); ++i) {
        if (pos[i].side == 2) {
            caseOf.push_back(caseOf.back()); // cùng trường hợp tải với điểm trước
            continue;
        }
        LoadCase lc("P=1");
        ElemLoad l;
        l.elem = pos[i].elem;
        l.kind = ElemLoadKind::Point;
        l.dir = LoadDir::GlobalY;
        l.q1 = -1.0;
        l.a = pos[i].a;
        lc.elem.push_back(l);
        mc.cases.push_back(lc);
        caseOf.push_back((int)mc.cases.size() - 1);
    }
    const AnalysisResult res = analyze(mc, 16);
    if (!res.ok) {
        il.error = res.error;
        return il;
    }
    if (t.kind == IlKind::Reaction) {
        const int d = 3 * t.node + t.dof;
        const Support* s = m.supportAt(t.node);
        if (res.dofState[d] != DofSupport && !(s && s->k[t.dof] > 0)) {
            il.error = "Nút đã chọn không có liên kết theo phương này.";
            return il;
        }
    }
    for (size_t i = 0; i < pos.size(); ++i) {
        const LoadSetResult& s = res.sets[caseOf[i]];
        const Element& el = m.elements[pos[i].elem];
        const Node& a = m.nodes[el.n[0]];
        const Node& b = m.nodes[el.n[1]];
        const double L = m.barLength(pos[i].elem), r = pos[i].a / L;
        double v = 0;
        switch (t.kind) {
        case IlKind::Reaction: v = s.R[3 * t.node + t.dof]; break;
        case IlKind::Displacement: v = s.U[3 * t.node + t.dof]; break;
        case IlKind::Force: {
            const BarResult& br = s.bars[t.elem];
            if (pos[i].side == 0) {
                v = barValueAt(br, t.qty, tx, true);
                break;
            }
            // giá trị hai phía của lực P = 1 đặt đúng tại tiết diện
            const double Lt = m.barLength(t.elem), eps = 1e-9 * Lt;
            const double c = (m.nodes[m.elements[t.elem].n[1]].x - m.nodes[m.elements[t.elem].n[0]].x) / Lt;
            const double sn = (m.nodes[m.elements[t.elem].n[1]].y - m.nodes[m.elements[t.elem].n[0]].y) / Lt;
            // bước nhảy (phải − trái) khi đi qua lực P = 1 hướng xuống: ΔQ = −cos, ΔN = +sin, ΔM = 0
            const double jump = t.qty == BarQty::Q ? -c : t.qty == BarQty::N ? sn : 0.0;
            double vRight, vLeft;
            if (tx <= eps) {
                vRight = barValueAt(br, t.qty, tx, true);
                vLeft = vRight - jump;
            } else if (tx >= Lt - eps) {
                vLeft = barValueAt(br, t.qty, tx, false);
                vRight = vLeft + jump;
            } else {
                vRight = barValueAt(br, t.qty, tx, true);
                vLeft = barValueAt(br, t.qty, tx, false);
            }
            // lực ở bên trái tiết diện → tiết diện nằm bên phải lực
            v = pos[i].side == 1 ? vRight : vLeft;
            break;
        }
        }
        il.pts.push_back({pos[i].elem, pos[i].a, a.x + (b.x - a.x) * r, a.y + (b.y - a.y) * r, v});
    }
    il.vmax = il.vmin = il.pts[0].value;
    for (const auto& p : il.pts) {
        il.vmax = std::max(il.vmax, p.value);
        il.vmin = std::min(il.vmin, p.value);
    }
    // diện tích theo hình chiếu ngang, tách phần dương/âm tại điểm đổi dấu
    for (size_t i = 0; i + 1 < il.pts.size(); ++i) {
        const IlPoint &p0 = il.pts[i], &p1 = il.pts[i + 1];
        if (p0.elem != p1.elem) continue;
        const double dx = std::abs(p1.X - p0.X), v0 = p0.value, v1 = p1.value;
        if ((v0 >= 0) == (v1 >= 0)) {
            const double a = 0.5 * (v0 + v1) * dx;
            (a >= 0 ? il.areaPos : il.areaNeg) += a;
        } else {
            const double z = v0 / (v0 - v1) * dx;
            const double a0 = 0.5 * v0 * z, a1 = 0.5 * v1 * (dx - z);
            (a0 >= 0 ? il.areaPos : il.areaNeg) += a0;
            (a1 >= 0 ? il.areaPos : il.areaNeg) += a1;
        }
    }
    il.ok = true;
    return il;
}

std::string ilTargetName(const Model& m, const IlTarget& t) {
    char b[160];
    static const char* dofs[] = {"x", "y", "φ"};
    const int d = std::clamp(t.dof, 0, 2);
    switch (t.kind) {
    case IlKind::Reaction:
        std::snprintf(b, sizeof b, "đah R%s tại nút %d", d == RZ ? "M" : dofs[d], t.node + 1);
        break;
    case IlKind::Displacement:
        std::snprintf(b, sizeof b, "đah %s%s tại nút %d", d == RZ ? "" : "u", dofs[d], t.node + 1);
        break;
    case IlKind::Force: {
        static const char* q[] = {"N", "Q", "M"};
        const double L = (t.elem >= 0 && t.elem < (int)m.elements.size() && m.elements[t.elem].isBar()) ? m.barLength(t.elem) : 0;
        std::snprintf(b, sizeof b, "đah %s tại tiết diện x = %.3g m của thanh %d", q[(int)t.qty], std::clamp(t.x, 0.0, L),
                      t.elem + 1);
        break;
    }
    }
    return b;
}

const std::vector<LoadTrain>& trainPresets() {
    static const std::vector<LoadTrain> v = {
        {"Lực đơn vị P = 1", {{1, 0}}, 0},
        {"Xe H30 (22TCN 18-79): 60 + 120 + 120 kN", {{60, 0}, {120, 6.0}, {120, 7.6}}, 0},
        {"Xe tải HL-93 (TCVN 11823): 35 + 145 + 145 kN, làn 9,3 kN/m", {{35, 0}, {145, 4.3}, {145, 8.6}}, 9.3},
        {"Xe hai trục HL-93: 2 × 110 kN, làn 9,3 kN/m", {{110, 0}, {110, 1.2}}, 9.3},
        {"НК-80: 4 × 200 kN, cách 1,2 m", {{200, 0}, {200, 1.2}, {200, 2.4}, {200, 3.6}}, 0},
    };
    return v;
}

double ilValueAtX(const InfluenceLine& il, double X) {
    // các đoạn theo thứ tự X tăng dần (đường chạy đã sắp từ trái sang phải)
    std::vector<std::pair<double, double>> pts;
    for (const auto& p : il.pts) pts.push_back({p.X, p.value});
    std::stable_sort(pts.begin(), pts.end(), [](auto& a, auto& b) { return a.first < b.first; });
    if (pts.empty() || X < pts.front().first - 1e-12 || X > pts.back().first + 1e-12) return 0;
    for (size_t k = 0; k + 1 < pts.size(); ++k) {
        const double x0 = pts[k].first, x1 = pts[k + 1].first;
        if (X >= x0 && X <= x1) {
            if (x1 - x0 <= 1e-12) continue; // bước nhảy: lấy đoạn kế tiếp
            return pts[k].second + (X - x0) / (x1 - x0) * (pts[k + 1].second - pts[k].second);
        }
    }
    return pts.back().second;
}

TrainResult moveTrain(const InfluenceLine& il, const LoadTrain& t) {
    TrainResult r;
    if (!il.ok || il.pts.empty() || t.axles.empty()) return r;
    double X0 = 1e300, X1 = -1e300;
    for (const auto& p : il.pts) {
        X0 = std::min(X0, p.X);
        X1 = std::max(X1, p.X);
    }
    double len = 0;
    for (const auto& a : t.axles) len = std::max(len, a.offset);
    // tải làn: phủ phần cùng dấu
    r.laneMax = t.lane * std::max(0.0, il.areaPos);
    r.laneMin = t.lane * std::min(0.0, il.areaNeg);
    // vị trí thử: mỗi trục đặt tại mỗi đỉnh đường ảnh hưởng (lời giải tối ưu cho đường gãy khúc) + lưới đều
    std::vector<double> cand;
    for (const auto& p : il.pts)
        for (const auto& a : t.axles) {
            cand.push_back(p.X - a.offset);
            cand.push_back(p.X + a.offset); // khi quay đầu
        }
    const int nGrid = 600;
    for (int i = 0; i <= nGrid; ++i) cand.push_back(X0 - len + (X1 - X0 + 2 * len) * i / nGrid);
    r.Smax = -1e300;
    r.Smin = 1e300;
    for (double x : cand)
        for (int rev = 0; rev < 2; ++rev) {
            double S = 0;
            for (const auto& a : t.axles) {
                const double xa = rev ? x - a.offset : x + a.offset;
                S += a.P * ilValueAtX(il, xa);
                // tại bước nhảy: thử thêm cả hai phía
            }
            if (S > r.Smax) {
                r.Smax = S;
                r.xMax = x;
                r.revMax = rev;
            }
            if (S < r.Smin) {
                r.Smin = S;
                r.xMin = x;
                r.revMin = rev;
            }
        }
    r.Smax = std::max(r.Smax, 0.0) + r.laneMax;
    r.Smin = std::min(r.Smin, 0.0) + r.laneMin;
    r.ok = true;
    return r;
}

DiagramGeom diagramGeometry(const Model& m, const LoadSetResult& s, int e, BarQty qty, double scale) {
    DiagramGeom g;
    if (e < 0 || e >= (int)m.elements.size() || !m.elements[e].isBar() || e >= (int)s.bars.size()) return g;
    const Element& el = m.elements[e];
    const Node& a = m.nodes[el.n[0]];
    const Node& b = m.nodes[el.n[1]];
    const double L = std::hypot(b.x - a.x, b.y - a.y), c = (b.x - a.x) / L, sn = (b.y - a.y) / L;
    const double sign = qty == BarQty::M ? -1.0 : 1.0;
    const BarResult& r = s.bars[e];
    const std::vector<double>& v = barValues(r, qty);
    for (size_t k = 0; k < r.x.size(); ++k) {
        const double bx = a.x + c * r.x[k], by = a.y + sn * r.x[k], d = sign * v[k] * scale;
        g.base.push_back({bx, by});
        g.off.push_back({bx - sn * d, by + c * d}); // trục y cục bộ (−s, c)
        g.values.push_back(v[k]);
    }
    return g;
}

double autoDiagramScale(const Model& m, const LoadSetResult& s, BarQty qty, double fraction) {
    double vmax = 0;
    for (size_t e = 0; e < m.elements.size() && e < s.bars.size(); ++e)
        if (m.elements[e].isBar())
            for (double v : barValues(s.bars[e], qty)) vmax = std::max(vmax, std::abs(v));
    double x0, y0, x1, y1;
    m.bounds(x0, y0, x1, y1);
    const double size = std::max({x1 - x0, y1 - y0, 1.0});
    return vmax > 0 ? fraction * size / vmax : 0;
}

} // namespace fem

