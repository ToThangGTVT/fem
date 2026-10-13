#include "core/ForceMethod.h"

#include "core/PostProcess.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>

namespace fem {

namespace {

std::string fmt(const char* f, ...) __attribute__((format(printf, 1, 2)));
std::string fmt(const char* f, ...) {
    char b[512];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(b, sizeof b, f, ap);
    va_end(ap);
    return b;
}

// Phân tích LU với chọn trụ theo hàng; giải nhiều vế phải.
struct LU {
    int n = 0;
    std::vector<double> a;
    std::vector<int> piv;
    bool factor(const std::vector<double>& A, int N) {
        n = N;
        a = A;
        piv.resize(n);
        for (int i = 0; i < n; ++i) piv[i] = i;
        double scale = 0;
        for (double v : a) scale = std::max(scale, std::abs(v));
        for (int k = 0; k < n; ++k) {
            int p = k;
            for (int i = k + 1; i < n; ++i)
                if (std::abs(a[i * n + k]) > std::abs(a[p * n + k])) p = i;
            if (std::abs(a[p * n + k]) <= 1e-13 * scale) return false;
            if (p != k) {
                for (int j = 0; j < n; ++j) std::swap(a[p * n + j], a[k * n + j]);
                std::swap(piv[p], piv[k]);
            }
            for (int i = k + 1; i < n; ++i) {
                const double f = a[i * n + k] /= a[k * n + k];
                for (int j = k + 1; j < n; ++j) a[i * n + j] -= f * a[k * n + j];
            }
        }
        return true;
    }
    std::vector<double> solve(const std::vector<double>& b) const {
        std::vector<double> x(n);
        for (int i = 0; i < n; ++i) x[i] = b[piv[i]];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j) x[i] -= a[i * n + j] * x[j];
        for (int i = n - 1; i >= 0; --i) {
            for (int j = i + 1; j < n; ++j) x[i] -= a[i * n + j] * x[j];
            x[i] /= a[i * n + i];
        }
        return x;
    }
};

const char* dofTxt(int d) { return d == UX ? "x" : d == UY ? "y" : "M"; }

// Các tải trên thanh và tải nút của một bộ tải
void collectLoads(const Model& m, int set, std::vector<BarLoads>& bl, std::vector<double>& P, std::vector<std::string>& warn) {
    const int nn = (int)m.nodes.size(), ne = (int)m.elements.size(), nc = (int)m.cases.size();
    bl.assign(ne, BarLoads{});
    P.assign(3 * nn, 0.0);
    std::vector<std::pair<int, double>> terms;
    if (set < nc) terms = {{set, 1.0}};
    else if (set - nc < (int)m.combos.size()) terms = m.combos[set - nc].terms;
    double sw = 0;
    bool warned = false;
    for (const auto& [ci, f] : terms) {
        if (ci < 0 || ci >= nc) continue;
        const LoadCase& lc = m.cases[ci];
        sw += f * lc.selfWeight;
        for (const auto& l : lc.nodal)
            if (l.node >= 0 && l.node < nn)
                for (int c = 0; c < 3; ++c) P[3 * l.node + c] += f * l.F[c];
        for (const auto& l : lc.elem) {
            if (l.elem < 0 || l.elem >= ne || !m.elements[l.elem].isBar()) continue;
            if (l.kind == ElemLoadKind::Temperature || l.kind == ElemLoadKind::LackOfFit) {
                if (!warned) warn.push_back("Tác dụng nhiệt độ / sai số chế tạo không xét trong lời giải phương pháp lực.");
                warned = true;
                continue;
            }
            addBarLoad(barProps(m, l.elem), l, f, bl[l.elem]);
        }
        if (!lc.disp.empty()) warn.push_back("Chuyển vị cưỡng bức không xét trong lời giải phương pháp lực.");
    }
    if (sw != 0)
        for (int e = 0; e < ne; ++e)
            if (m.elements[e].isBar())
                addBarSelfWeight(barProps(m, e), sw * m.materials[m.elements[e].mat].gamma * m.sections[m.elements[e].sec].A, bl[e]);
}

// Ma trận cân bằng: ẩn là N, Mi, Mj của các thanh và phản lực; phương trình là cân bằng nút
struct EqSys {
    int R = 0, C = 0;
    std::vector<int> row;   // theo bậc tự do
    std::vector<std::vector<double>> col; // cột đầy đủ (R phần tử)
    std::vector<int> bar, kind, dof;      // mô tả cột
};

bool endRigid(const Element& el, int k) { return el.type == ElemType::Frame && !(k == 0 ? el.hingeI : el.hingeJ); }

EqSys buildEq(const Model& m) {
    EqSys q;
    const int nn = (int)m.nodes.size();
    std::vector<int> rigid(nn, 0);
    std::vector<char> used(nn, 0);
    for (const auto& el : m.elements) {
        if (!el.isBar()) continue;
        for (int k = 0; k < 2; ++k) {
            used[el.n[k]] = 1;
            if (endRigid(el, k)) rigid[el.n[k]]++;
        }
    }
    q.row.assign(3 * nn, -1);
    for (int i = 0; i < nn; ++i) {
        if (!used[i]) continue;
        q.row[3 * i] = q.R++;
        q.row[3 * i + 1] = q.R++;
        const Support* s = m.supportAt(i);
        if (rigid[i] > 0 || (s && s->fix[2])) q.row[3 * i + 2] = q.R++;
    }
    auto newCol = [&](int bar, int kind, int dof) {
        q.col.emplace_back(q.R, 0.0);
        q.bar.push_back(bar);
        q.kind.push_back(kind);
        q.dof.push_back(dof);
        return (int)q.col.size() - 1;
    };
    auto put = [&](int c, int node, int d, double v) {
        const int r = q.row[3 * node + d];
        if (r >= 0) q.col[c][r] += v;
    };
    for (size_t e = 0; e < m.elements.size(); ++e) {
        const Element& el = m.elements[e];
        if (!el.isBar()) continue;
        const Node& a = m.nodes[el.n[0]];
        const Node& b = m.nodes[el.n[1]];
        const double L = std::hypot(b.x - a.x, b.y - a.y), c = (b.x - a.x) / L, s = (b.y - a.y) / L;
        int k = newCol((int)e, 0, -1); // N (kéo): thanh kéo nút i về phía j
        put(k, el.n[0], 0, c);
        put(k, el.n[0], 1, s);
        put(k, el.n[1], 0, -c);
        put(k, el.n[1], 1, -s);
        for (int end = 0; end < 2; ++end) {
            if (!endRigid(el, end)) continue;
            // mômen đầu thanh = 1 (ngược kim đồng hồ, nút tác dụng lên thanh), lực cắt cân bằng V = 1/L
            k = newCol((int)e, 1 + end, -1);
            put(k, el.n[0], 0, s / L);
            put(k, el.n[0], 1, -c / L);
            put(k, el.n[1], 0, -s / L);
            put(k, el.n[1], 1, c / L);
            put(k, el.n[end], 2, -1);
        }
    }
    for (const auto& sp : m.supports) {
        if (sp.node < 0 || sp.node >= nn || !used[sp.node]) continue;
        for (int d = 0; d < 3; ++d)
            if (sp.fix[d] && q.row[3 * sp.node + d] >= 0) {
                const int k = newCol(-1, 3, 3 * sp.node + d);
                q.col[k][q.row[3 * sp.node + d]] = 1;
            }
    }
    q.C = (int)q.col.size();
    return q;
}

// Chọn cột trụ (ẩn cơ bản) bằng khử Gauss trên các cột đã chuẩn hóa; các cột còn lại là ẩn thừa
void choosePivots(const EqSys& q, std::vector<int>& piv, std::vector<int>& red) {
    const int R = q.R, C = q.C;
    std::vector<double> A((size_t)R * C);
    for (int j = 0; j < C; ++j) {
        double mx = 0;
        for (int i = 0; i < R; ++i) mx = std::max(mx, std::abs(q.col[j][i]));
        for (int i = 0; i < R; ++i) A[(size_t)i * C + j] = mx > 0 ? q.col[j][i] / mx : 0;
    }
    piv.clear();
    red.clear();
    int r = 0;
    for (int c = 0; c < C; ++c) {
        int best = -1;
        double bv = 1e-9;
        for (int i = r; i < R; ++i)
            if (std::abs(A[(size_t)i * C + c]) > bv) {
                bv = std::abs(A[(size_t)i * C + c]);
                best = i;
            }
        if (best < 0) {
            red.push_back(c);
            continue;
        }
        for (int k = 0; k < C; ++k) std::swap(A[(size_t)best * C + k], A[(size_t)r * C + k]);
        for (int i = 0; i < R; ++i) {
            if (i == r) continue;
            const double f = A[(size_t)i * C + c] / A[(size_t)r * C + c];
            if (f == 0) continue;
            for (int k = c; k < C; ++k) A[(size_t)i * C + k] -= f * A[(size_t)r * C + k];
        }
        piv.push_back(c);
        ++r;
    }
}

// Lưới điểm dọc mỗi thanh và biểu đồ riêng (tải trong nhịp trên dầm đơn giản)
struct BarGrid {
    std::vector<BarResult> part; // biểu đồ phần riêng (N, Q, M do tải trong nhịp, hai đầu khớp)
    std::vector<BarProps> props;
    std::vector<BarLoads> loads;
    std::vector<double> EI, EA, L;
};

// Trạng thái nội lực của thanh: N, Mi, Mj (ẩn) và hệ số của phần riêng
struct BarState {
    double N = 0, Mi = 0, Mj = 0, pf = 0;
};

BarGrid particular(const Model& m, const std::vector<BarLoads>& bl, std::vector<double>& cpart) {
    BarGrid g;
    const int ne = (int)m.elements.size();
    g.part.resize(ne);
    g.props.resize(ne);
    g.loads = bl;
    g.EI.assign(ne, 0);
    g.EA.assign(ne, 0);
    g.L.assign(ne, 0);
    cpart.assign(3 * m.nodes.size(), 0.0);
    for (int e = 0; e < ne; ++e) {
        const Element& el = m.elements[e];
        if (!el.isBar()) continue;
        BarProps p = barProps(m, e);
        g.EA[e] = p.EA;
        g.EI[e] = el.type == ElemType::Frame ? p.EI : 0;
        g.L[e] = p.L;
        p.relI = p.relJ = true;
        p.kf = p.kG = 0;
        Mat6 K;
        Vec6 fg;
        barGlobal(p, bl[e], K, fg);
        for (int k = 0; k < 2; ++k)
            for (int c = 0; c < 2; ++c) cpart[3 * el.n[k] + c] += fg[3 * k + c];
        barRecover(p, bl[e], Vec6{}, 40, g.part[e]);
        g.props[e] = p;
    }
    return g;
}

// Trạng thái từng thanh từ véc tơ ẩn q
std::vector<BarState> barStates(const Model& m, const EqSys& q, const std::vector<double>& qv, bool withPart) {
    std::vector<BarState> st(m.elements.size());
    for (int j = 0; j < q.C; ++j) {
        if (q.kind[j] == 0) st[q.bar[j]].N = qv[j];
        else if (q.kind[j] == 1) st[q.bar[j]].Mi = qv[j];
        else if (q.kind[j] == 2) st[q.bar[j]].Mj = qv[j];
    }
    for (auto& b : st) b.pf = withPart ? 1.0 : 0.0;
    return st;
}

void evalState(const BarGrid& g, int e, const BarState& s, double x, bool after, double& N, double& Q, double& M) {
    const double L = g.L[e];
    N = s.N;
    Q = (s.Mi + s.Mj) / L;
    M = -s.Mi + (s.Mi + s.Mj) * x / L;
    if (s.pf != 0) {
        double n, qq, mm;
        barForcesAt(g.props[e], g.loads[e], g.part[e], x, after, n, qq, mm);
        N += s.pf * n;
        Q += s.pf * qq;
        M += s.pf * mm;
    }
}

// Biểu đồ để vẽ
void stateDiagrams(const Model& m, const EqSys& q, const std::vector<double>& qv, const std::vector<BarState>& st,
                   const BarGrid& g, LoadSetResult& out) {
    const int ne = (int)m.elements.size(), nn = (int)m.nodes.size();
    out.bars.assign(ne, BarResult{});
    out.U.assign(3 * nn, 0.0);
    out.R.assign(3 * nn, 0.0);
    out.elemStress.assign(ne, {});
    out.nodeStress.assign(nn, {0, 0, 0, 0});
    out.nodeStressCount.assign(nn, 0);
    for (int j = 0; j < q.C; ++j)
        if (q.kind[j] == 3) out.R[q.dof[j]] = qv[j];
    for (int e = 0; e < ne; ++e) {
        if (!m.elements[e].isBar()) continue;
        const BarResult& p = g.part[e];
        BarResult& b = out.bars[e];
        b.x = p.x;
        const size_t n = p.x.size();
        b.N.resize(n);
        b.Q.resize(n);
        b.M.resize(n);
        b.u.assign(n, 0);
        b.v.assign(n, 0);
        for (size_t k = 0; k < n; ++k) {
            const bool after = !(k + 1 < n && p.x[k + 1] == p.x[k]);
            evalState(g, e, st[e], p.x[k], after, b.N[k], b.Q[k], b.M[k]);
        }
        b.endForces = {-b.N.front(), b.Q.front(), -b.M.front(), b.N.back(), -b.Q.back(), b.M.back()};
    }
}

// ∫ (Mᵃ·Mᵇ/EI + Nᵃ·Nᵇ/EA) dx trên thanh e: Gauss 3 điểm trong mỗi khoảng giữa các điểm mẫu (hàm trơn trong khoảng)
double integrate(const BarGrid& g, int e, const BarState& a, const BarState& b, double* mPart = nullptr) {
    static const double gx[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
    static const double gw[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
    const auto& xs = g.part[e].x;
    double sm = 0, sn = 0;
    for (size_t k = 0; k + 1 < xs.size(); ++k) {
        const double x0 = xs[k], h = xs[k + 1] - x0;
        if (h <= 0) continue;
        for (int i = 0; i < 3; ++i) {
            const double x = x0 + 0.5 * h * (gx[i] + 1), w = 0.5 * h * gw[i];
            double Na, Qa, Ma, Nb, Qb, Mb;
            evalState(g, e, a, x, true, Na, Qa, Ma);
            evalState(g, e, b, x, true, Nb, Qb, Mb);
            if (g.EI[e] > 0) sm += w * Ma * Mb / g.EI[e];
            if (g.EA[e] > 0) sn += w * Na * Nb / g.EA[e];
        }
    }
    if (mPart) *mPart = sm;
    return sm + sn;
}

std::string validateBars(const Model& m) {
    if (!m.hasBars()) return "Chỉ áp dụng cho hệ thanh.";
    if (m.hasPlane()) return "Mô hình có phần tử phẳng — lời giải theo phương pháp lực chỉ áp dụng cho hệ thanh.";
    for (const auto& el : m.elements)
        if (el.kf > 0 || el.kG > 0) return "Không áp dụng cho thanh trên nền đàn hồi.";
    for (const auto& s : m.supports)
        for (int c = 0; c < 3; ++c)
            if (!s.fix[c] && s.k[c] > 0) return "Không áp dụng cho gối đàn hồi.";
    return {};
}

} // namespace

ForceMethodResult forceMethod(const Model& m, int set) {
    ForceMethodResult r;
    r.error = validateBars(m);
    if (!r.error.empty()) return r;
    const int nsets = (int)(m.cases.size() + m.combos.size());
    if (set < 0 || set >= nsets) {
        r.error = "Chưa chọn bộ tải.";
        return r;
    }
    const int nn = (int)m.nodes.size();
    std::vector<BarLoads> bl;
    std::vector<double> P, cpart;
    collectLoads(m, set, bl, P, r.warnings);
    const EqSys q = buildEq(m);
    std::vector<int> piv, red;
    choosePivots(q, piv, red);
    if ((int)piv.size() < q.R) {
        r.error = "Hệ biến hình — không giải được bằng phương pháp lực.";
        return r;
    }
    r.n = (int)red.size();
    const BarGrid g = particular(m, bl, cpart);

    // ma trận các cột cơ bản
    const int R = q.R;
    std::vector<double> Bp((size_t)R * R);
    for (int i = 0; i < R; ++i)
        for (int k = 0; k < R; ++k) Bp[(size_t)i * R + k] = q.col[piv[k]][i];
    LU lu;
    if (!lu.factor(Bp, R)) {
        r.error = "Không lập được hệ cơ bản tĩnh định.";
        return r;
    }
    auto solveState = [&](const std::vector<double>& rhs, int redCol, double redVal) {
        std::vector<double> b = rhs;
        if (redCol >= 0)
            for (int i = 0; i < R; ++i) b[i] -= q.col[redCol][i] * redVal;
        const std::vector<double> x = lu.solve(b);
        std::vector<double> qv(q.C, 0.0);
        for (int k = 0; k < R; ++k) qv[piv[k]] = x[k];
        if (redCol >= 0) qv[redCol] = redVal;
        return qv;
    };
    // trạng thái tải: B q = −(P + c_part)
    std::vector<double> rhsP(R, 0.0);
    for (int d = 0; d < 3 * nn; ++d) {
        if (q.row[d] >= 0) rhsP[q.row[d]] -= P[d] + cpart[d];
        else if (P[d] != 0 && d % 3 == RZ)
            r.warnings.push_back(fmt("Mômen tập trung tại nút %d (không có thanh nối cứng) bị bỏ qua.", d / 3 + 1));
    }
    const std::vector<double> qP = solveState(rhsP, -1, 0);
    const std::vector<BarState> sP = barStates(m, q, qP, true);
    stateDiagrams(m, q, qP, sP, g, r.load);
    r.load.name = "M⁰_P (hệ cơ bản)";
    std::vector<std::vector<double>> qU(r.n);
    std::vector<std::vector<BarState>> sU(r.n);
    r.unit.resize(r.n);
    for (int i = 0; i < r.n; ++i) {
        qU[i] = solveState(std::vector<double>(R, 0.0), red[i], 1.0);
        sU[i] = barStates(m, q, qU[i], false);
        stateDiagrams(m, q, qU[i], sU[i], g, r.unit[i]);
        r.unit[i].name = fmt("M̄%d", i + 1);
        const int c = red[i];
        if (q.kind[c] == 3) r.redundants.push_back(fmt("X%d: phản lực R%s tại gối nút %d (bỏ liên kết)", i + 1, dofTxt(q.dof[c] % 3), q.dof[c] / 3 + 1));
        else if (q.kind[c] == 0) r.redundants.push_back(fmt("X%d: lực dọc trong thanh %d (cắt thanh)", i + 1, q.bar[c] + 1));
        else {
            const Element& el = m.elements[q.bar[c]];
            r.redundants.push_back(fmt("X%d: mômen đầu %s thanh %d (đặt khớp tại nút %d)", i + 1, q.kind[c] == 1 ? "i" : "j",
                                       q.bar[c] + 1, el.n[q.kind[c] - 1] + 1));
        }
    }
    // hệ số δᵢₖ, Δᵢₚ
    const int ne = (int)m.elements.size();
    r.delta.assign(r.n, std::vector<double>(r.n, 0.0));
    r.deltaP.assign(r.n, 0.0);
    for (int i = 0; i < r.n; ++i) {
        for (int k = i; k < r.n; ++k) {
            double s = 0;
            for (int e = 0; e < ne; ++e)
                if (m.elements[e].isBar()) s += integrate(g, e, sU[i][e], sU[k][e]);
            r.delta[i][k] = r.delta[k][i] = s;
        }
        double s = 0;
        for (int e = 0; e < ne; ++e)
            if (m.elements[e].isBar()) s += integrate(g, e, sU[i][e], sP[e]);
        r.deltaP[i] = s;
    }
    // giải hệ phương trình chính tắc
    r.X.assign(r.n, 0.0);
    if (r.n > 0) {
        std::vector<double> D((size_t)r.n * r.n);
        for (int i = 0; i < r.n; ++i)
            for (int k = 0; k < r.n; ++k) D[(size_t)i * r.n + k] = r.delta[i][k];
        LU ld;
        if (!ld.factor(D, r.n)) {
            r.error = "Hệ phương trình chính tắc suy biến.";
            return r;
        }
        std::vector<double> rhs(r.n);
        for (int i = 0; i < r.n; ++i) rhs[i] = -r.deltaP[i];
        r.X = ld.solve(rhs);
    }
    std::vector<double> qF = qP;
    for (int i = 0; i < r.n; ++i)
        for (int j = 0; j < q.C; ++j) qF[j] += r.X[i] * qU[i][j];
    stateDiagrams(m, q, qF, barStates(m, q, qF, true), g, r.final);
    r.final.name = "M (cuối cùng)";

    // đối chiếu với phần tử hữu hạn
    const AnalysisResult fe = analyze(m);
    if (fe.ok && set < (int)fe.sets.size()) {
        double mmax = 1e-12, dmax = 0;
        for (int e = 0; e < ne; ++e) {
            if (!m.elements[e].isBar()) continue;
            const BarResult& b = fe.sets[set].bars[e];
            const BarState st = barStates(m, q, qF, true)[e];
            for (size_t k = 0; k < b.x.size(); ++k) {
                const bool after = !(k + 1 < b.x.size() && b.x[k + 1] == b.x[k]);
                double N, Q, M;
                evalState(g, e, st, b.x[k], after, N, Q, M);
                mmax = std::max(mmax, std::abs(b.M[k]));
                dmax = std::max(dmax, std::abs(b.M[k] - M));
            }
        }
        r.maxDiffM = dmax;
        (void)mmax;
        r.final.U = fe.sets[set].U; // chuyển vị lấy từ lời giải PTHH để vẽ biến dạng nếu cần
    }
    // dữ liệu cho Maxwell–Mohr
    r.rowOf = q.row;
    r.pivots = piv;
    r.redundantCols = red;
    r.Bfull.assign(q.C, {});
    for (int j = 0; j < q.C; ++j) r.Bfull[j] = q.col[j];
    r.barOf = q.bar;
    r.kindOf = q.kind;
    r.loadSet = set;
    r.qFinal = qF;
    r.ok = true;
    return r;
}

MohrResult mohrDisplacement(const Model& m, const ForceMethodResult& fm, int node, int dof) {
    MohrResult res;
    if (!fm.ok) {
        res.error = "Chưa có lời giải phương pháp lực.";
        return res;
    }
    const int nn = (int)m.nodes.size();
    if (node < 0 || node >= nn || dof < 0 || dof > 2) {
        res.error = "Nút không hợp lệ.";
        return res;
    }
    const int row = fm.rowOf[3 * node + dof];
    if (row < 0) {
        res.error = "Nút không có phương trình cân bằng theo phương này (ví dụ góc xoay tại nút khớp).";
        return res;
    }
    // dựng lại hệ với lưới và trạng thái đơn vị
    const EqSys q = buildEq(m);
    const int R = q.R;
    std::vector<double> Bp((size_t)R * R);
    for (int i = 0; i < R; ++i)
        for (int k = 0; k < R; ++k) Bp[(size_t)i * R + k] = q.col[fm.pivots[k]][i];
    LU lu;
    if (!lu.factor(Bp, R)) {
        res.error = "Không lập được hệ cơ bản.";
        return res;
    }
    std::vector<double> rhs(R, 0.0);
    rhs[row] = -1.0; // lực đơn vị P = 1 tác dụng tại nút theo phương dof
    const std::vector<double> x = lu.solve(rhs);
    std::vector<double> qv(q.C, 0.0);
    for (int k = 0; k < R; ++k) qv[fm.pivots[k]] = x[k];
    std::vector<BarLoads> bl;
    std::vector<double> P, cpart;
    std::vector<std::string> warn;
    collectLoads(m, fm.loadSet, bl, P, warn);
    const BarGrid g = particular(m, bl, cpart);
    const std::vector<BarState> su = barStates(m, q, qv, false), sf = barStates(m, q, fm.qFinal, true);
    stateDiagrams(m, q, qv, su, g, res.unitState);
    res.unitState.name = "M̄_k (lực đơn vị)";
    for (size_t e = 0; e < m.elements.size(); ++e) {
        if (!m.elements[e].isBar()) continue;
        double mPart = 0;
        const double v = integrate(g, (int)e, su[e], sf[e], &mPart);
        if (std::abs(v) > 1e-15)
            res.terms.push_back(fmt("Thanh %zu: ∫M̄M/EI = %.4e, ∫N̄N/EA = %.4e", e + 1, mPart, v - mPart));
        res.value += v;
    }
    if (!fm.final.U.empty()) res.fem = fm.final.U[3 * node + dof];
    res.ok = true;
    return res;
}

MatrixMethodResult matrixMethod(const Model& m, int set, int maxEq) {
    MatrixMethodResult r;
    if (!m.hasBars() || m.hasPlane()) {
        r.error = "Chỉ áp dụng cho hệ thanh.";
        return r;
    }
    const int nn = (int)m.nodes.size(), ne = (int)m.elements.size(), ndof = 3 * nn;
    const AnalysisResult fe = analyze(m);
    if (!fe.ok) {
        r.error = fe.error;
        return r;
    }
    if (set < 0 || set >= (int)fe.sets.size()) {
        r.error = "Chưa chọn bộ tải.";
        return r;
    }
    std::vector<int> eq(ndof, -1);
    int neq = 0;
    for (int d = 0; d < ndof; ++d)
        if (fe.dofState[d] == DofFree) {
            eq[d] = neq++;
            r.dofLabels.push_back(fmt("%s%d", d % 3 == UX ? "u" : d % 3 == UY ? "v" : "φ", d / 3 + 1));
        }
    if (neq > maxEq) {
        r.error = fmt("Mô hình có %d ẩn chuyển vị (> %d) — quá lớn để trình bày dạng ma trận.", neq, maxEq);
        return r;
    }
    r.K.assign(neq, std::vector<double>(neq, 0.0));
    r.P.assign(neq, 0.0);
    std::vector<BarLoads> bl;
    std::vector<double> P;
    std::vector<std::string> warn;
    collectLoads(m, set, bl, P, warn);
    for (int e = 0; e < ne; ++e) {
        const Element& el = m.elements[e];
        if (!el.isBar()) continue;
        const BarProps p = barProps(m, e);
        Mat6 Kl, Kg;
        Vec6 fg;
        barLocalStiffness(p, Kl);
        barGlobal(p, bl[e], Kg, fg);
        MatrixMethodResult::Elem info;
        info.e = e;
        info.L = p.L;
        info.c = p.c;
        info.s = p.s;
        info.kLocal.assign(6, std::vector<double>(6));
        info.kGlobal.assign(6, std::vector<double>(6));
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                info.kLocal[i][j] = Kl[i][j];
                info.kGlobal[i][j] = Kg[i][j];
            }
        for (int k = 0; k < 2; ++k)
            for (int c = 0; c < 3; ++c) {
                const int d = 3 * el.n[k] + c;
                info.dofs.push_back(eq[d] >= 0 ? r.dofLabels[eq[d]] : "—");
            }
        for (int i = 0; i < 6; ++i) {
            const int di = 3 * el.n[i / 3] + i % 3, ei = eq[di];
            if (ei < 0) continue;
            r.P[ei] += fg[i];
            for (int j = 0; j < 6; ++j) {
                const int dj = 3 * el.n[j / 3] + j % 3, ej = eq[dj];
                if (ej >= 0) r.K[ei][ej] += Kg[i][j];
                else r.P[ei] -= Kg[i][j] * fe.sets[set].U[dj]; // chuyển vị gối đã biết
            }
        }
        r.elems.push_back(info);
    }
    for (const auto& s : m.supports)
        for (int c = 0; c < 3; ++c)
            if (!s.fix[c] && s.k[c] > 0 && eq[3 * s.node + c] >= 0) r.K[eq[3 * s.node + c]][eq[3 * s.node + c]] += s.k[c];
    for (int d = 0; d < ndof; ++d)
        if (eq[d] >= 0) r.P[eq[d]] += P[d];
    std::vector<double> Kd((size_t)neq * neq);
    for (int i = 0; i < neq; ++i)
        for (int j = 0; j < neq; ++j) Kd[(size_t)i * neq + j] = r.K[i][j];
    LU lu;
    if (neq > 0 && !lu.factor(Kd, neq)) {
        r.error = "Ma trận độ cứng suy biến.";
        return r;
    }
    r.Z = neq > 0 ? lu.solve(r.P) : std::vector<double>{};
    r.ok = true;
    return r;
}

} // namespace fem
