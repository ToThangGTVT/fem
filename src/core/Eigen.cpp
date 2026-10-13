#include "core/Eigen.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>

#include "core/Analysis.h"
#include "core/Elements.h"
#include "core/Mesher.h"

namespace fem {

namespace {

constexpr double kG = 9.81; // m/s²
constexpr double kPi = 3.14159265358979323846;

// Householder đưa về ba đường chéo (thuật toán tred2, theo EISPACK/JAMA)
void tred2(int n, std::vector<double>& Vm, std::vector<double>& d, std::vector<double>& e) {
    auto V = [&](int i, int j) -> double& { return Vm[(size_t)i * n + j]; };
    for (int j = 0; j < n; ++j) d[j] = V(n - 1, j);
    for (int i = n - 1; i > 0; --i) {
        double scale = 0, h = 0;
        for (int k = 0; k < i; ++k) scale += std::abs(d[k]);
        if (scale == 0) {
            e[i] = d[i - 1];
            for (int j = 0; j < i; ++j) {
                d[j] = V(i - 1, j);
                V(i, j) = 0;
                V(j, i) = 0;
            }
        } else {
            for (int k = 0; k < i; ++k) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            double f = d[i - 1], g = std::sqrt(h);
            if (f > 0) g = -g;
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for (int j = 0; j < i; ++j) e[j] = 0;
            for (int j = 0; j < i; ++j) {
                f = d[j];
                V(j, i) = f;
                g = e[j] + V(j, j) * f;
                for (int k = j + 1; k <= i - 1; ++k) {
                    g += V(k, j) * d[k];
                    e[k] += V(k, j) * f;
                }
                e[j] = g;
            }
            f = 0;
            for (int j = 0; j < i; ++j) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            const double hh = f / (h + h);
            for (int j = 0; j < i; ++j) e[j] -= hh * d[j];
            for (int j = 0; j < i; ++j) {
                f = d[j];
                g = e[j];
                for (int k = j; k <= i - 1; ++k) V(k, j) -= (f * e[k] + g * d[k]);
                d[j] = V(i - 1, j);
                V(i, j) = 0;
            }
        }
        d[i] = h;
    }
    for (int i = 0; i < n - 1; ++i) {
        V(n - 1, i) = V(i, i);
        V(i, i) = 1;
        const double h = d[i + 1];
        if (h != 0) {
            for (int k = 0; k <= i; ++k) d[k] = V(k, i + 1) / h;
            for (int j = 0; j <= i; ++j) {
                double g = 0;
                for (int k = 0; k <= i; ++k) g += V(k, i + 1) * V(k, j);
                for (int k = 0; k <= i; ++k) V(k, j) -= g * d[k];
            }
        }
        for (int k = 0; k <= i; ++k) V(k, i + 1) = 0;
    }
    for (int j = 0; j < n; ++j) {
        d[j] = V(n - 1, j);
        V(n - 1, j) = 0;
    }
    V(n - 1, n - 1) = 1;
    e[0] = 0;
}

// QL ẩn cho ma trận ba đường chéo (tql2)
void tql2(int n, std::vector<double>& Vm, std::vector<double>& d, std::vector<double>& e) {
    auto V = [&](int i, int j) -> double& { return Vm[(size_t)i * n + j]; };
    for (int i = 1; i < n; ++i) e[i - 1] = e[i];
    e[n - 1] = 0;
    double f = 0, tst1 = 0;
    const double eps = std::numeric_limits<double>::epsilon();
    for (int l = 0; l < n; ++l) {
        tst1 = std::max(tst1, std::abs(d[l]) + std::abs(e[l]));
        int m = l;
        while (m < n) {
            if (std::abs(e[m]) <= eps * tst1) break;
            ++m;
        }
        if (m == n) m = n - 1;
        if (m > l) {
            int iter = 0;
            do {
                if (++iter > 60) break;
                double g = d[l];
                double p = (d[l + 1] - g) / (2 * e[l]);
                double r = std::hypot(p, 1.0);
                if (p < 0) r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                const double dl1 = d[l + 1];
                double h = g - d[l];
                for (int i = l + 2; i < n; ++i) d[i] -= h;
                f += h;
                p = d[m];
                double c = 1, c2 = c, c3 = c;
                const double el1 = e[l + 1];
                double s = 0, s2 = 0;
                for (int i = m - 1; i >= l; --i) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = std::hypot(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for (int k = 0; k < n; ++k) {
                        h = V(k, i + 1);
                        V(k, i + 1) = s * V(k, i) + c * h;
                        V(k, i) = c * V(k, i) - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
            } while (std::abs(e[l]) > eps * tst1);
        }
        d[l] += f;
        e[l] = 0;
    }
    // sắp xếp tăng dần
    for (int i = 0; i < n - 1; ++i) {
        int k = i;
        double p = d[i];
        for (int j = i + 1; j < n; ++j)
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < n; ++j) std::swap(V(j, i), V(j, k));
        }
    }
}

// Tải thẳng đứng tương đương tại nút của một bộ tải (không kể trọng lượng bản thân)
std::vector<double> verticalNodalLoads(const Model& m, int set) {
    std::vector<double> Fy(m.nodes.size(), 0.0);
    std::vector<std::pair<int, double>> terms;
    const int nc = (int)m.cases.size();
    if (set < nc) terms = {{set, 1.0}};
    else if (set - nc < (int)m.combos.size()) terms = m.combos[set - nc].terms;
    std::vector<BarLoads> bl(m.elements.size());
    for (const auto& [ci, f] : terms) {
        if (ci < 0 || ci >= nc) continue;
        for (const auto& l : m.cases[ci].nodal)
            if (l.node >= 0 && l.node < (int)m.nodes.size()) Fy[l.node] += f * l.F[1];
        for (const auto& l : m.cases[ci].elem)
            if (l.elem >= 0 && l.elem < (int)m.elements.size() && m.elements[l.elem].isBar() &&
                l.kind != ElemLoadKind::Temperature)
                addBarLoad(barProps(m, l.elem), l, f, bl[l.elem]);
    }
    for (size_t e = 0; e < m.elements.size(); ++e) {
        if (bl[e].empty()) continue;
        Mat6 K;
        Vec6 fg;
        barGlobal(barProps(m, (int)e), bl[e], K, fg);
        Fy[m.elements[e].n[0]] += fg[1];
        Fy[m.elements[e].n[1]] += fg[4];
    }
    return Fy;
}

} // namespace

void symmetricEigen(int n, std::vector<double>& A, std::vector<double>& d, std::vector<double>& V) {
    V = A;
    d.assign(n, 0.0);
    std::vector<double> e(n, 0.0);
    if (n == 0) return;
    if (n == 1) {
        d[0] = A[0];
        V[0] = 1;
        return;
    }
    tred2(n, V, d, e);
    tql2(n, V, d, e);
}

bool buildDenseSystem(const Model& m, const EigenOptions& opt, DenseSystem& s, std::string& err) {
    s = DenseSystem{};
    const int nsets = (int)(m.cases.size() + m.combos.size());
    if (m.elements.empty()) {
        err = "Mô hình chưa có phần tử.";
        return false;
    }
    // ---- chia nhỏ các thanh
    Model& sub = s.sub;
    sub = m;
    const int nsub = std::max(1, opt.subdiv);
    auto& pieces = s.pieces;
    pieces.assign(m.elements.size(), {});
    s.barChain.assign(m.elements.size(), {});
    for (size_t e = 0; e < m.elements.size(); ++e) {
        if (!m.elements[e].isBar()) continue;
        const int n0 = m.elements[e].n[0], n1 = m.elements[e].n[1];
        const int firstNode = (int)sub.nodes.size(), firstElem = (int)sub.elements.size();
        if (nsub > 1) subdivideBar(sub, (int)e, nsub);
        auto& chain = s.barChain[e];
        chain.push_back(n0);
        for (int k = firstNode; k < (int)sub.nodes.size(); ++k) chain.push_back(k);
        chain.push_back(n1);
        pieces[e].push_back((int)e);
        for (int k = firstElem; k < (int)sub.elements.size(); ++k) pieces[e].push_back(k);
    }

    // ---- lực dọc cho bài toán ổn định
    std::vector<double> Nsub(sub.elements.size(), 0.0);
    if (opt.kind == EigenKind::Buckling) {
        if (opt.loadSet < 0 || opt.loadSet >= nsets) {
            err = "Chưa chọn bộ tải gây lực dọc.";
            return false;
        }
        const AnalysisResult st = analyze(sub, 4);
        if (!st.ok) {
            err = st.error;
            return false;
        }
        const LoadSetResult& ls = st.sets[opt.loadSet];
        for (size_t e = 0; e < sub.elements.size(); ++e) {
            const auto& b = ls.bars[e];
            if (b.N.empty()) continue;
            double sum = 0;
            for (double v : b.N) sum += v;
            Nsub[e] = sum / b.N.size();
        }
        s.memberN.assign(m.elements.size(), 0.0);
        for (size_t e = 0; e < m.elements.size(); ++e) {
            if (pieces[e].empty()) continue;
            double sum = 0;
            for (int k : pieces[e]) sum += Nsub[k];
            s.memberN[e] = sum / pieces[e].size();
        }
    }

    // ---- ma trận phần tử
    const int nn = (int)sub.nodes.size(), ndof = 3 * nn;
    struct Blk {
        std::vector<int> dofs;
        std::vector<double> K, B;
    };
    std::vector<Blk> blk;
    std::vector<char> hasStiff(ndof, 0);
    bool planeWarn = false;
    for (size_t e = 0; e < sub.elements.size(); ++e) {
        const Element& el = sub.elements[e];
        if (el.mat < 0 || el.mat >= (int)sub.materials.size()) {
            err = "Phần tử " + std::to_string(e + 1) + " chưa gán vật liệu.";
            return false;
        }
        const Material& mat = sub.materials[el.mat];
        const double rho = mat.gamma / kG * opt.selfMassFactor;
        Blk b;
        if (el.isBar()) {
            if (el.sec < 0 || el.sec >= (int)sub.sections.size()) {
                err = "Thanh " + std::to_string(e + 1) + " chưa gán tiết diện.";
                return false;
            }
            const BarProps p = barProps(sub, (int)e);
            Mat6 K, KG, M;
            barEigenMatrices(p, Nsub[e], rho * sub.sections[el.sec].A, K, KG, M);
            for (int k = 0; k < 2; ++k)
                for (int c = 0; c < 3; ++c) b.dofs.push_back(3 * el.n[k] + c);
            b.K.resize(36);
            b.B.resize(36);
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) {
                    b.K[i * 6 + j] = K[i][j];
                    b.B[i * 6 + j] = opt.kind == EigenKind::Buckling ? -KG[i][j] : M[i][j];
                }
        } else {
            double D[3][3], x[4], y[4];
            planeD(mat, sub.planeMode, D);
            const int nc = el.nodeCount();
            for (int k = 0; k < nc; ++k) {
                x[k] = sub.nodes[el.n[k]].x;
                y[k] = sub.nodes[el.n[k]].y;
                for (int c = 0; c < 2; ++c) b.dofs.push_back(3 * el.n[k] + c);
            }
            const int n = 2 * nc;
            b.K.assign(n * n, 0.0);
            b.B.assign(n * n, 0.0);
            bool ok;
            if (el.type == ElemType::Tri3) {
                double K[6][6];
                ok = tri3Stiffness(x, y, D, el.t, K);
                for (int i = 0; i < 6; ++i)
                    for (int j = 0; j < 6; ++j) b.K[i * 6 + j] = K[i][j];
            } else {
                double K[8][8];
                ok = quad4Stiffness(x, y, D, el.t, el.incompatible, K);
                for (int i = 0; i < 8; ++i)
                    for (int j = 0; j < 8; ++j) b.K[i * 8 + j] = K[i][j];
            }
            if (!ok) {
                err = "Phần tử " + std::to_string(e + 1) + " bị suy biến.";
                return false;
            }
            if (opt.kind == EigenKind::Vibration) {
                double area2 = 0;
                for (int k = 0; k < nc; ++k) area2 += x[k] * y[(k + 1) % nc] - x[(k + 1) % nc] * y[k];
                const double mn = rho * el.t * std::abs(area2) / 2 / nc; // khối lượng tập trung
                for (int i = 0; i < n; ++i) b.B[i * n + i] = mn;
            } else {
                planeWarn = true;
            }
        }
        const int n = (int)b.dofs.size();
        double kmax = 0;
        for (int i = 0; i < n; ++i) kmax = std::max(kmax, std::abs(b.K[i * n + i]));
        for (int i = 0; i < n; ++i)
            if (b.K[i * n + i] > 1e-12 * kmax) hasStiff[b.dofs[i]] = 1;
        blk.push_back(std::move(b));
    }
    if (planeWarn)
        s.warnings.push_back("Phần tử phẳng không được kể độ cứng hình học — chỉ các thanh tham gia mất ổn định.");

    // ---- điều kiện biên & đánh số
    std::vector<char> fixed(ndof, 0);
    std::vector<double> spring(ndof, 0.0);
    for (const auto& sp : sub.supports)
        for (int c = 0; c < 3; ++c) {
            if (sp.node < 0 || sp.node >= nn) continue;
            if (sp.fix[c]) fixed[3 * sp.node + c] = 1;
            else if (sp.k[c] > 0) {
                spring[3 * sp.node + c] += sp.k[c];
                hasStiff[3 * sp.node + c] = 1;
            }
        }
    std::vector<int> eq(ndof, -1);
    int neq = 0;
    for (int d = 0; d < ndof; ++d)
        if (!fixed[d] && hasStiff[d]) eq[d] = neq++;
    s.neq = neq;
    if (neq == 0) {
        err = "Không có bậc tự do tự do.";
        return false;
    }
    if (neq > 2500) {
        err = "Mô hình quá lớn cho bộ giải trị riêng dạng ma trận đầy (" + std::to_string(neq) +
                    " phương trình > 2500). Hãy giảm số đoạn chia hoặc lưới.";
        return false;
    }
    const size_t N = (size_t)neq;
    s.K.assign(N * N, 0.0);
    s.B.assign(N * N, 0.0);
    auto& K = s.K;
    auto& B = s.B;
    for (const auto& b : blk) {
        const int n = (int)b.dofs.size();
        for (int i = 0; i < n; ++i) {
            const int ei = eq[b.dofs[i]];
            if (ei < 0) continue;
            for (int j = 0; j < n; ++j) {
                const int ej = eq[b.dofs[j]];
                if (ej < 0) continue;
                K[ei * N + ej] += b.K[i * n + j];
                B[ei * N + ej] += b.B[i * n + j];
            }
        }
    }
    for (int d = 0; d < ndof; ++d)
        if (eq[d] >= 0 && spring[d] > 0) K[eq[d] * N + eq[d]] += spring[d];
    if (opt.kind == EigenKind::Vibration && opt.massSet >= 0 && opt.massSet < nsets) {
        const std::vector<double> Fy = verticalNodalLoads(sub, opt.massSet);
        for (int i = 0; i < nn; ++i) {
            const double mass = std::abs(Fy[i]) / kG * opt.massSetFactor;
            for (int c = 0; c < 2; ++c)
                if (eq[3 * i + c] >= 0) B[eq[3 * i + c] * N + eq[3 * i + c]] += mass;
        }
    }

    s.eq = eq;
    s.spring = spring;
    return true;
}

EigenResult eigenAnalysis(const Model& m, const EigenOptions& opt) {
    const auto t0 = std::chrono::steady_clock::now();
    EigenResult res;
    res.kind = opt.kind;
    DenseSystem sys;
    {
        std::string err;
        if (!buildDenseSystem(m, opt, sys, err)) {
            res.error = err;
            return res;
        }
    }
    res.sub = sys.sub;
    res.barChain = sys.barChain;
    res.memberN = sys.memberN;
    res.warnings = sys.warnings;
    res.numEquations = sys.neq;
    const int nn = (int)res.sub.nodes.size(), ndof = 3 * nn, neq = sys.neq;
    const size_t N = (size_t)neq;
    const std::vector<int>& eq = sys.eq;
    std::vector<double>& K = sys.K;
    std::vector<double>& B = sys.B;
    // ---- K = L·Lᵀ (Cholesky)
    std::vector<double>& L = K; // phân tích tại chỗ, phần tam giác dưới
    double kmax = 0;
    for (size_t i = 0; i < N; ++i) kmax = std::max(kmax, L[i * N + i]);
    for (size_t j = 0; j < N; ++j) {
        double s = L[j * N + j];
        for (size_t k = 0; k < j; ++k) s -= L[j * N + k] * L[j * N + k];
        if (!(s > 1e-13 * kmax)) {
            int dof = 0;
            for (int d = 0; d < ndof; ++d)
                if (eq[d] == (int)j) dof = d;
            res.error = "Hệ biến hình: ma trận độ cứng suy biến tại nút " + std::to_string(dof / 3 + 1) + " (" +
                        dofName(dof % 3) + ").";
            return res;
        }
        const double ljj = std::sqrt(s);
        L[j * N + j] = ljj;
        for (size_t i = j + 1; i < N; ++i) {
            double t = L[i * N + j];
            for (size_t k = 0; k < j; ++k) t -= L[i * N + k] * L[j * N + k];
            L[i * N + j] = t / ljj;
        }
    }
    // ---- C = L⁻¹ B L⁻ᵀ
    // X = L⁻¹ B (giải tiến theo từng cột), rồi C = L⁻¹ Xᵀ
    auto forward = [&](std::vector<double>& M) {
        for (size_t c = 0; c < N; ++c)
            for (size_t i = 0; i < N; ++i) {
                double v = M[i * N + c];
                for (size_t k = 0; k < i; ++k) v -= L[i * N + k] * M[k * N + c];
                M[i * N + c] = v / L[i * N + i];
            }
    };
    forward(B);
    std::vector<double> C(N * N);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j) C[i * N + j] = B[j * N + i];
    forward(C);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = i + 1; j < N; ++j) {
            const double v = 0.5 * (C[i * N + j] + C[j * N + i]);
            C[i * N + j] = C[j * N + i] = v;
        }
    std::vector<double> mu, V;
    symmetricEigen(neq, C, mu, V);

    // μ lớn nhất ↔ λ = 1/μ nhỏ nhất
    double mumax = 0;
    for (double v : mu) mumax = std::max(mumax, std::abs(v));
    std::vector<int> order;
    for (int i = neq - 1; i >= 0; --i)
        if (mu[i] > 1e-10 * mumax) order.push_back(i);
    if (order.empty()) {
        res.error = opt.kind == EigenKind::Buckling
                        ? "Không có thanh chịu nén dưới bộ tải đã chọn — kết cấu không mất ổn định."
                        : "Không có khối lượng (kiểm tra γ của vật liệu hoặc bộ tải quy khối lượng).";
        return res;
    }
    const int nm = std::min((int)order.size(), std::max(1, opt.modes));
    for (int k = 0; k < nm; ++k) {
        const int i = order[k];
        EigenMode md;
        md.lambda = 1.0 / mu[i];
        if (opt.kind == EigenKind::Vibration) {
            md.omega = std::sqrt(md.lambda);
            md.freq = md.omega / (2 * kPi);
            md.period = md.freq > 0 ? 1 / md.freq : 0;
        }
        // φ = L⁻ᵀ y
        std::vector<double> phi(N);
        for (size_t r = 0; r < N; ++r) phi[r] = V[r * N + i];
        for (size_t r = N; r-- > 0;) {
            double v = phi[r];
            for (size_t k2 = r + 1; k2 < N; ++k2) v -= L[k2 * N + r] * phi[k2];
            phi[r] = v / L[r * N + r];
        }
        md.U.assign(ndof, 0.0);
        for (int d = 0; d < ndof; ++d)
            if (eq[d] >= 0) md.U[d] = phi[eq[d]];
        // chuẩn hóa: chuyển vị thẳng lớn nhất = 1, dấu dương
        double big = 0;
        int ib = 0;
        for (int nd = 0; nd < nn; ++nd)
            for (int c = 0; c < 2; ++c)
                if (std::abs(md.U[3 * nd + c]) > big) {
                    big = std::abs(md.U[3 * nd + c]);
                    ib = 3 * nd + c;
                }
        if (big == 0)
            for (int d = 0; d < ndof; ++d)
                if (std::abs(md.U[d]) > big) {
                    big = std::abs(md.U[d]);
                    ib = d;
                }
        const double sc = big > 0 ? (md.U[ib] > 0 ? 1 : -1) / big : 1;
        for (double& v : md.U) v *= sc;
        res.modes.push_back(std::move(md));
    }
    res.ok = true;
    res.timeMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return res;
}

double effectiveLengthFactor(const Model& m, const EigenResult& r, int e) {
    if (!r.ok || r.kind != EigenKind::Buckling || r.modes.empty() || e < 0 || e >= (int)m.elements.size()) return 0;
    const Element& el = m.elements[e];
    if (!el.isBar() || e >= (int)r.memberN.size()) return 0;
    const double Ncr = -r.memberN[e] * r.modes[0].lambda; // lực nén tới hạn (dương)
    if (Ncr <= 0) return 0;
    const double EI = m.materials[el.mat].E * m.sections[el.sec].I;
    if (EI <= 0) return 0;
    return kPi / m.barLength(e) * std::sqrt(EI / Ncr);
}

} // namespace fem
