#include "core/Dynamics.h"

#include <algorithm>
#include <cmath>

namespace fem {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kGrav = 9.81;

struct DenseLU {
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
                if (std::abs(a[(size_t)i * n + k]) > std::abs(a[(size_t)p * n + k])) p = i;
            if (std::abs(a[(size_t)p * n + k]) <= 1e-13 * scale) return false;
            if (p != k) {
                for (int j = 0; j < n; ++j) std::swap(a[(size_t)p * n + j], a[(size_t)k * n + j]);
                std::swap(piv[p], piv[k]);
            }
            const double d = a[(size_t)k * n + k];
            for (int i = k + 1; i < n; ++i) {
                double& f = a[(size_t)i * n + k];
                f /= d;
                if (f == 0) continue;
                const double* rk = &a[(size_t)k * n];
                double* ri = &a[(size_t)i * n];
                for (int j = k + 1; j < n; ++j) ri[j] -= f * rk[j];
            }
        }
        return true;
    }
    std::vector<double> solve(const std::vector<double>& b) const {
        std::vector<double> x(n);
        for (int i = 0; i < n; ++i) x[i] = b[piv[i]];
        for (int i = 0; i < n; ++i) {
            const double* ri = &a[(size_t)i * n];
            double s = x[i];
            for (int j = 0; j < i; ++j) s -= ri[j] * x[j];
            x[i] = s;
        }
        for (int i = n - 1; i >= 0; --i) {
            const double* ri = &a[(size_t)i * n];
            double s = x[i];
            for (int j = i + 1; j < n; ++j) s -= ri[j] * x[j];
            x[i] = s / ri[i];
        }
        return x;
    }
};

std::vector<double> matvec(const std::vector<double>& A, const std::vector<double>& x) {
    const size_t n = x.size();
    std::vector<double> y(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        const double* r = &A[i * n];
        double s = 0;
        for (size_t j = 0; j < n; ++j) s += r[j] * x[j];
        y[i] = s;
    }
    return y;
}

// Mô hình đã chia nhỏ dùng để phân tích tĩnh với lực nút cho trước
Model linearSub(const DenseSystem& sys) {
    Model sm = sys.sub;
    sm.secondOrder = false;
    for (auto& e : sm.elements) e.kfCompOnly = false;
    return sm;
}

// Phân tích tĩnh mô hình chia nhỏ với các véc tơ lực nút (bậc tự do tự do) → kết quả từng véc tơ, quy về mô hình gốc
std::vector<LoadSetResult> staticOnSub(const Model& orig, const DenseSystem& sys, const std::vector<std::vector<double>>& F,
                                       const std::vector<std::string>& names, std::string& err) {
    Model sm = linearSub(sys);
    sm.cases.clear();
    sm.combos.clear();
    const int ndof = (int)sys.eq.size();
    for (size_t k = 0; k < F.size(); ++k) {
        LoadCase lc(k < names.size() ? names[k] : "F");
        for (int node = 0; node < ndof / 3; ++node) {
            NodalLoad l{node, {0, 0, 0}};
            bool any = false;
            for (int c = 0; c < 3; ++c) {
                const int e = sys.eq[3 * node + c];
                if (e >= 0 && F[k][e] != 0) {
                    l.F[c] = F[k][e];
                    any = true;
                }
            }
            if (any) lc.nodal.push_back(l);
        }
        sm.cases.push_back(lc);
    }
    std::vector<LoadSetResult> out;
    if (sm.cases.empty()) return out;
    const AnalysisResult r = analyze(sm, 8);
    if (!r.ok) {
        err = r.error;
        return out;
    }
    for (const auto& s : r.sets) out.push_back(mergeToOriginal(orig, sys, s));
    return out;
}

// u (bậc tự do tự do) của lời giải tĩnh một bộ tải trên mô hình chia nhỏ
bool staticFree(const DenseSystem& sys, int set, std::vector<double>& u, LoadSetResult* full, const Model& orig, std::string& err) {
    const Model sm = linearSub(sys);
    const AnalysisResult r = analyze(sm, 8);
    if (!r.ok) {
        err = r.error;
        return false;
    }
    if (set < 0 || set >= (int)r.sets.size()) {
        err = "Chưa chọn bộ tải.";
        return false;
    }
    u.assign(sys.neq, 0.0);
    for (size_t d = 0; d < sys.eq.size(); ++d)
        if (sys.eq[d] >= 0) u[sys.eq[d]] = r.sets[set].U[d];
    if (full) *full = mergeToOriginal(orig, sys, r.sets[set]);
    return true;
}

// Tổ hợp phản ứng các dạng: √(Σᵢ Σⱼ ρᵢⱼ vᵢ vⱼ)
struct Combiner {
    std::vector<std::vector<double>> rho;
    double operator()(const std::vector<double>& v) const {
        double s = 0;
        const size_t n = v.size();
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j) s += rho[i][j] * v[i] * v[j];
        return std::sqrt(std::max(0.0, s));
    }
};

} // namespace

// ============================================================================ ghép về mô hình gốc
LoadSetResult mergeToOriginal(const Model& orig, const DenseSystem& sys, const LoadSetResult& s) {
    LoadSetResult o;
    o.name = s.name;
    o.combo = s.combo;
    const int nn = (int)orig.nodes.size(), ne = (int)orig.elements.size();
    o.U.assign(s.U.begin(), s.U.begin() + std::min<size_t>(s.U.size(), 3 * nn));
    o.R.assign(s.R.begin(), s.R.begin() + std::min<size_t>(s.R.size(), 3 * nn));
    o.U.resize(3 * nn, 0.0);
    o.R.resize(3 * nn, 0.0);
    o.sumFx = s.sumFx;
    o.sumFy = s.sumFy;
    o.sumRx = s.sumRx;
    o.sumRy = s.sumRy;
    o.foundRx = s.foundRx;
    o.foundRy = s.foundRy;
    o.bars.assign(ne, BarResult{});
    o.elemStress.assign(ne, {});
    for (int e = 0; e < ne; ++e) {
        if (e < (int)s.elemStress.size()) o.elemStress[e] = s.elemStress[e];
        if (!orig.elements[e].isBar() || e >= (int)sys.pieces.size()) continue;
        BarResult& b = o.bars[e];
        double off = 0;
        const auto& pcs = sys.pieces[e];
        for (size_t k = 0; k < pcs.size(); ++k) {
            const BarResult& p = s.bars[pcs[k]];
            for (size_t i = 0; i < p.x.size(); ++i) {
                if (k > 0 && i == 0) continue; // điểm nối trùng với điểm cuối đoạn trước
                b.x.push_back(off + p.x[i]);
                b.N.push_back(p.N[i]);
                b.Q.push_back(p.Q[i]);
                b.M.push_back(p.M[i]);
                b.u.push_back(p.u[i]);
                b.v.push_back(p.v[i]);
            }
            off += sys.sub.barLength(pcs[k]);
        }
        if (!pcs.empty()) {
            const BarResult& f = s.bars[pcs.front()];
            const BarResult& l = s.bars[pcs.back()];
            b.endForces = {f.endForces[0], f.endForces[1], f.endForces[2], l.endForces[3], l.endForces[4], l.endForces[5]};
            b.uLocal = {f.uLocal[0], f.uLocal[1], f.uLocal[2], l.uLocal[3], l.uLocal[4], l.uLocal[5]};
        }
    }
    o.nodeStress.assign(s.nodeStress.begin(), s.nodeStress.begin() + std::min<size_t>(s.nodeStress.size(), nn));
    o.nodeStress.resize(nn, {0, 0, 0, 0});
    o.nodeStressCount.assign(s.nodeStressCount.begin(), s.nodeStressCount.begin() + std::min<size_t>(s.nodeStressCount.size(), nn));
    o.nodeStressCount.resize(nn, 0);
    return o;
}

// ============================================================================ phổ phản ứng
void soilParams(int soil, double& S, double& TB, double& TC, double& TD) {
    // Phổ đàn hồi loại 1 (TCVN 9386:2012, bảng 3.2)
    static const double tbl[5][4] = {{1.00, 0.15, 0.4, 2.0},
                                     {1.20, 0.15, 0.5, 2.0},
                                     {1.15, 0.20, 0.6, 2.0},
                                     {1.35, 0.20, 0.8, 2.0},
                                     {1.40, 0.15, 0.5, 2.0}};
    soil = std::clamp(soil, 0, 4);
    S = tbl[soil][0];
    TB = tbl[soil][1];
    TC = tbl[soil][2];
    TD = tbl[soil][3];
}

double spectrumSd(const SpectrumParams& p, double T) {
    double S, TB, TC, TD;
    soilParams(p.soil, S, TB, TC, TD);
    const double ag = p.agR * p.gammaI * kGrav, q = p.q;
    if (T <= TB) return ag * S * (2.0 / 3.0 + T / TB * (2.5 / q - 2.0 / 3.0));
    if (T <= TC) return ag * S * 2.5 / q;
    if (T <= TD) return std::max(ag * S * 2.5 / q * TC / T, p.beta * ag);
    return std::max(ag * S * 2.5 / q * TC * TD / (T * T), p.beta * ag);
}

SpectrumResult spectrumAnalysis(const Model& m, const EigenOptions& eo, const SpectrumParams& sp) {
    SpectrumResult res;
    EigenOptions o = eo;
    o.kind = EigenKind::Vibration;
    o.modes = std::max(1, eo.modes);
    const EigenResult er = eigenAnalysis(m, o);
    if (!er.ok) {
        res.error = er.error;
        return res;
    }
    DenseSystem sys;
    std::string err;
    if (!buildDenseSystem(m, o, sys, err)) {
        res.error = err;
        return res;
    }
    const int neq = sys.neq, ndof = (int)sys.eq.size();
    const int dirDof = sp.direction == 0 ? UX : UY;
    std::vector<double> r(neq, 0.0);
    for (int d = 0; d < ndof; ++d)
        if (sys.eq[d] >= 0 && d % 3 == dirDof) r[sys.eq[d]] = 1.0;
    const std::vector<double> Mr = matvec(sys.B, r);
    for (int i = 0; i < neq; ++i) res.totalMass += r[i] * Mr[i];
    if (res.totalMass <= 0) {
        res.error = "Không có khối lượng theo phương đã chọn.";
        return res;
    }
    std::vector<std::vector<double>> F;
    std::vector<std::string> names;
    std::vector<double> omegas;
    for (size_t i = 0; i < er.modes.size(); ++i) {
        const EigenMode& md = er.modes[i];
        std::vector<double> phi(neq, 0.0);
        for (int d = 0; d < ndof; ++d)
            if (sys.eq[d] >= 0) phi[sys.eq[d]] = md.U[d];
        const std::vector<double> Mphi = matvec(sys.B, phi);
        double Mn = 0, Ln = 0;
        for (int k = 0; k < neq; ++k) {
            Mn += phi[k] * Mphi[k];
            Ln += r[k] * Mphi[k];
        }
        ModalInfo mi;
        mi.omega = md.omega;
        mi.f = md.freq;
        mi.T = md.period;
        mi.gamma = Ln / Mn;
        mi.mEff = Ln * Ln / Mn;
        mi.ratio = mi.mEff / res.totalMass;
        mi.Sd = spectrumSd(sp, mi.T);
        res.sumRatio += mi.ratio;
        // lực quán tính tương đương: F = Γ·Sd·M·φ
        std::vector<double> f(neq);
        for (int k = 0; k < neq; ++k) f[k] = mi.gamma * mi.Sd * Mphi[k];
        F.push_back(f);
        char b[32];
        std::snprintf(b, sizeof b, "Dạng %zu", i + 1);
        names.push_back(b);
        omegas.push_back(md.omega);
        res.modes.push_back(mi);
    }
    if (res.sumRatio < 0.9)
        res.warnings.push_back("Tổng khối lượng hữu hiệu < 90% (TCVN 9386 mục 4.3.3.3) — nên tăng số dạng dao động.");
    res.modal = staticOnSub(m, sys, F, names, err);
    if (res.modal.size() != F.size()) {
        res.error = err.empty() ? "Lỗi phân tích tĩnh các dạng." : err;
        return res;
    }
    // hệ số tương quan
    const size_t nm = res.modal.size();
    Combiner comb;
    comb.rho.assign(nm, std::vector<double>(nm, 0.0));
    for (size_t i = 0; i < nm; ++i)
        for (size_t j = 0; j < nm; ++j) {
            if (!sp.cqc) {
                comb.rho[i][j] = i == j ? 1.0 : 0.0;
                continue;
            }
            const double z = sp.zeta, rr = omegas[j] / omegas[i];
            comb.rho[i][j] = 8 * z * z * (1 + rr) * std::pow(rr, 1.5) /
                             ((1 - rr * rr) * (1 - rr * rr) + 4 * z * z * rr * (1 + rr) * (1 + rr));
        }
    // lực cắt đáy từng dạng (tổng phản lực theo phương đang xét)
    std::vector<double> vb(nm);
    for (size_t i = 0; i < nm; ++i) {
        double s = 0;
        for (size_t d = dirDof; d < res.modal[i].R.size(); d += 3) s += res.modal[i].R[d];
        // kể phản lực nền nếu có
        s += sp.direction == 0 ? res.modal[i].foundRx : res.modal[i].foundRy;
        vb[i] = s;
        res.modes[i].baseShear = std::abs(s);
    }
    res.baseShear = comb(vb);
    // tổ hợp trường kết quả
    LoadSetResult& c = res.combined;
    c = res.modal[0];
    c.name = std::string("Động đất ") + (sp.direction == 0 ? "X" : "Y") + (sp.cqc ? " (CQC)" : " (SRSS)");
    std::vector<double> v(nm);
    auto combineArr = [&](auto get, auto set, size_t n) {
        for (size_t k = 0; k < n; ++k) {
            for (size_t i = 0; i < nm; ++i) v[i] = get(res.modal[i], k);
            set(c, k, comb(v));
        }
    };
    combineArr([](const LoadSetResult& s, size_t k) { return s.U[k]; }, [](LoadSetResult& s, size_t k, double x) { s.U[k] = x; }, c.U.size());
    combineArr([](const LoadSetResult& s, size_t k) { return s.R[k]; }, [](LoadSetResult& s, size_t k, double x) { s.R[k] = x; }, c.R.size());
    for (size_t e = 0; e < c.bars.size(); ++e) {
        BarResult& b = c.bars[e];
        for (size_t k = 0; k < b.x.size(); ++k) {
            for (size_t i = 0; i < nm; ++i) v[i] = res.modal[i].bars[e].N[k];
            b.N[k] = comb(v);
            for (size_t i = 0; i < nm; ++i) v[i] = res.modal[i].bars[e].Q[k];
            b.Q[k] = comb(v);
            for (size_t i = 0; i < nm; ++i) v[i] = res.modal[i].bars[e].M[k];
            b.M[k] = comb(v);
            for (size_t i = 0; i < nm; ++i) v[i] = res.modal[i].bars[e].u[k];
            b.u[k] = comb(v);
            for (size_t i = 0; i < nm; ++i) v[i] = res.modal[i].bars[e].v[k];
            b.v[k] = comb(v);
        }
    }
    for (size_t n = 0; n < c.nodeStress.size(); ++n)
        for (int q = 0; q < 4; ++q) {
            for (size_t i = 0; i < nm; ++i) v[i] = res.modal[i].nodeStress[n][q];
            c.nodeStress[n][q] = comb(v);
        }
    for (size_t e = 0; e < c.elemStress.size(); ++e)
        for (int k = 0; k < 4; ++k)
            for (int q = 0; q < 4; ++q) {
                for (size_t i = 0; i < nm; ++i) v[i] = res.modal[i].elemStress[e][k][q];
                c.elemStress[e][k][q] = comb(v);
            }
    c.sumRx = sp.direction == 0 ? res.baseShear : 0;
    c.sumRy = sp.direction == 1 ? res.baseShear : 0;
    res.ok = true;
    return res;
}

// ============================================================================ dao động cưỡng bức điều hòa
HarmonicResult harmonicAnalysis(const Model& m, const EigenOptions& eo, int loadSet, double theta) {
    HarmonicResult res;
    res.theta = theta;
    EigenOptions o = eo;
    o.kind = EigenKind::Vibration;
    DenseSystem sys;
    std::string err;
    if (!buildDenseSystem(m, o, sys, err)) {
        res.error = err;
        return res;
    }
    std::vector<double> ust;
    if (!staticFree(sys, loadSet, ust, &res.stat, m, err)) {
        res.error = err;
        return res;
    }
    const int n = sys.neq;
    const std::vector<double> F = matvec(sys.K, ust);
    std::vector<double> A = sys.K;
    for (size_t i = 0; i < A.size(); ++i) A[i] -= theta * theta * sys.B[i];
    DenseLU lu;
    if (!lu.factor(A, n)) {
        res.error = "Cộng hưởng: tần số lực kích thích trùng tần số riêng (K − θ²M suy biến).";
        return res;
    }
    const std::vector<double> u = lu.solve(F);
    double umax = 0, smax = 0;
    for (size_t d = 0; d < sys.eq.size(); ++d) {
        if (sys.eq[d] < 0 || d % 3 == RZ) continue;
        umax = std::max(umax, std::abs(u[sys.eq[d]]));
        smax = std::max(smax, std::abs(ust[sys.eq[d]]));
    }
    res.dynFactor = smax > 0 ? umax / smax : 0;
    const auto amp = staticOnSub(m, sys, {matvec(sys.K, u)}, {"Biên độ"}, err);
    if (amp.empty()) {
        res.error = err;
        return res;
    }
    res.amp = amp[0];
    res.amp.name = "Biên độ dao động cưỡng bức";
    o.modes = 6;
    const EigenResult er = eigenAnalysis(m, o);
    if (er.ok)
        for (const auto& md : er.modes) res.omegas.push_back(md.omega);
    res.ok = true;
    return res;
}

// ============================================================================ theo thời gian (Newmark)
double timeFunction(const TimeHistoryParams& p, double t) {
    switch (p.fn) {
    case TimeFn::Step: return t > 0 ? p.amp : 0.0;
    case TimeFn::Rect: return (t > 0 && t <= p.td) ? p.amp : 0.0;
    case TimeFn::HalfSine: return (t >= 0 && t <= p.td) ? p.amp * std::sin(kPi * t / p.td) : 0.0;
    case TimeFn::Harmonic:
    case TimeFn::GroundSine: return p.amp * std::sin(2 * kPi * p.freq * t);
    case TimeFn::GroundRecord: {
        const auto& r = p.record;
        if (r.empty() || t < r.front().first || t > r.back().first) return 0.0;
        for (size_t k = 0; k + 1 < r.size(); ++k)
            if (t >= r[k].first && t <= r[k + 1].first) {
                const double h = r[k + 1].first - r[k].first;
                return h > 0 ? p.amp * (r[k].second + (t - r[k].first) / h * (r[k + 1].second - r[k].second)) : p.amp * r[k].second;
            }
        return 0.0;
    }
    }
    return 0.0;
}

TimeHistoryResult timeHistory(const Model& m, const EigenOptions& eo, const TimeHistoryParams& p) {
    TimeHistoryResult res;
    EigenOptions o = eo;
    o.kind = EigenKind::Vibration;
    o.modes = 2;
    auto sys = std::make_shared<DenseSystem>();
    std::string err;
    if (!buildDenseSystem(m, o, *sys, err)) {
        res.error = err;
        return res;
    }
    res.sys = sys;
    const int n = sys->neq;
    const int steps = (int)std::ceil(p.tEnd / p.dt);
    if (p.dt <= 0 || steps < 1) {
        res.error = "Bước thời gian không hợp lệ.";
        return res;
    }
    if ((double)n * n * steps > 6e9) {
        res.error = "Quá nhiều phép tính (số ẩn² × số bước) — hãy giảm số đoạn chia, thời gian hoặc tăng bước thời gian.";
        return res;
    }
    // cản Rayleigh theo hai tần số riêng đầu
    const EigenResult er = eigenAnalysis(m, o);
    if (!er.ok) {
        res.error = er.error;
        return res;
    }
    const double w1 = er.modes[0].omega, w2 = er.modes.size() > 1 ? er.modes[1].omega : 0.0;
    res.omega1 = w1;
    if (w2 > 0) {
        res.alpha = 2 * p.zeta * w1 * w2 / (w1 + w2);
        res.betaR = 2 * p.zeta / (w1 + w2);
    } else {
        res.alpha = 2 * p.zeta * w1;
    }
    if (p.dt > 0.2 * 2 * kPi / std::max(w1, 1e-9) * 0.25)
        res.warnings.push_back("Bước thời gian lớn so với chu kỳ riêng thứ nhất — nên giảm Δt để kết quả chính xác.");
    // dạng tải
    std::vector<double> Fpat(n, 0.0);
    const bool ground = p.fn == TimeFn::GroundSine || p.fn == TimeFn::GroundRecord;
    if (ground) {
        std::vector<double> rr(n, 0.0);
        const int dirDof = p.groundDir == 0 ? UX : UY;
        for (size_t d = 0; d < sys->eq.size(); ++d)
            if (sys->eq[d] >= 0 && (int)(d % 3) == dirDof) rr[sys->eq[d]] = 1.0;
        Fpat = matvec(sys->B, rr);
        for (double& v : Fpat) v = -v; // F = −M·r·a_g
    } else {
        std::vector<double> ust;
        if (!staticFree(*sys, p.loadSet, ust, nullptr, m, err)) {
            res.error = err;
            return res;
        }
        Fpat = matvec(sys->K, ust);
    }
    // Newmark gia tốc trung bình
    const double gam = 0.5, bet = 0.25, dt = p.dt;
    const double a1 = 1 / (bet * dt * dt), a2 = 1 / (bet * dt), a3 = 1 / (2 * bet) - 1;
    const double b1 = gam / (bet * dt), b2 = gam / bet - 1, b3 = dt * (gam / (2 * bet) - 1);
    const auto& K = sys->K;
    const auto& M = sys->B;
    std::vector<double> C(K.size()), Kh(K.size());
    for (size_t i = 0; i < K.size(); ++i) {
        C[i] = res.alpha * M[i] + res.betaR * K[i];
        Kh[i] = K[i] + a1 * M[i] + b1 * C[i];
    }
    DenseLU lu;
    if (!lu.factor(Kh, n)) {
        res.error = "Ma trận hiệu dụng suy biến.";
        return res;
    }
    std::vector<double> u(n, 0.0), v(n, 0.0), a(n, 0.0), w(n), z(n);
    res.t.reserve(steps + 1);
    res.u.reserve(steps + 1);
    res.maxAbs.assign(sys->eq.size(), 0.0);
    auto store = [&](double t) {
        res.t.push_back(t);
        res.u.emplace_back(u.begin(), u.end());
        for (size_t d = 0; d < sys->eq.size(); ++d)
            if (sys->eq[d] >= 0) res.maxAbs[d] = std::max(res.maxAbs[d], std::abs(u[sys->eq[d]]));
    };
    store(0);
    for (int s = 1; s <= steps; ++s) {
        const double t = s * dt, f = timeFunction(p, t);
        for (int i = 0; i < n; ++i) {
            w[i] = a1 * u[i] + a2 * v[i] + a3 * a[i];
            z[i] = b1 * u[i] + b2 * v[i] + b3 * a[i];
        }
        const std::vector<double> Mw = matvec(M, w), Cz = matvec(C, z);
        std::vector<double> rhs(n);
        for (int i = 0; i < n; ++i) rhs[i] = f * Fpat[i] + Mw[i] + Cz[i];
        const std::vector<double> un = lu.solve(rhs);
        for (int i = 0; i < n; ++i) {
            const double an = a1 * (un[i] - u[i]) - a2 * v[i] - a3 * a[i];
            v[i] += dt * ((1 - gam) * a[i] + gam * an);
            a[i] = an;
            u[i] = un[i];
        }
        store(t);
    }
    res.ok = true;
    return res;
}

std::vector<double> historyOf(const TimeHistoryResult& r, int node, int dof) {
    std::vector<double> h(r.t.size(), 0.0);
    if (!r.sys) return h;
    const int d = 3 * node + dof;
    if (d < 0 || d >= (int)r.sys->eq.size() || r.sys->eq[d] < 0) return h;
    const int e = r.sys->eq[d];
    for (size_t k = 0; k < r.u.size(); ++k) h[k] = r.u[k][e];
    return h;
}

LoadSetResult timeSnapshot(const Model& m, const TimeHistoryResult& r, int k) {
    if (!r.ok || !r.sys || k < 0 || k >= (int)r.u.size()) return LoadSetResult{};
    std::vector<double> u(r.u[k].begin(), r.u[k].end());
    std::string err;
    char b[64];
    std::snprintf(b, sizeof b, "t = %.3f s", r.t[k]);
    auto out = staticOnSub(m, *r.sys, {matvec(r.sys->K, u)}, {b}, err);
    return out.empty() ? LoadSetResult{} : out[0];
}

} // namespace fem
