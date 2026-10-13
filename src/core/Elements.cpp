#include "core/Elements.h"

#include <algorithm>
#include <cmath>

namespace fem {

namespace {

// Nghịch đảo ma trận vuông n×n (lưu theo hàng) bằng Gauss–Jordan. Trả về false nếu suy biến.
bool invert(std::vector<double>& A, int n) {
    std::vector<double> inv(n * n, 0.0);
    for (int i = 0; i < n; ++i) inv[i * n + i] = 1.0;
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int r = col + 1; r < n; ++r)
            if (std::abs(A[r * n + col]) > std::abs(A[piv * n + col])) piv = r;
        if (std::abs(A[piv * n + col]) < 1e-300) return false;
        if (piv != col)
            for (int k = 0; k < n; ++k) {
                std::swap(A[piv * n + k], A[col * n + k]);
                std::swap(inv[piv * n + k], inv[col * n + k]);
            }
        const double d = A[col * n + col];
        for (int k = 0; k < n; ++k) {
            A[col * n + k] /= d;
            inv[col * n + k] /= d;
        }
        for (int r = 0; r < n; ++r) {
            if (r == col) continue;
            const double f = A[r * n + col];
            if (f == 0.0) continue;
            for (int k = 0; k < n; ++k) {
                A[r * n + k] -= f * A[col * n + k];
                inv[r * n + k] -= f * inv[col * n + k];
            }
        }
    }
    A.swap(inv);
    return true;
}

std::vector<int> released(const BarProps& p) {
    std::vector<int> r;
    if (p.relI) r.push_back(2);
    if (p.relJ) r.push_back(5);
    return r;
}

bool isReleased(const std::vector<int>& R, int i) { return std::find(R.begin(), R.end(), i) != R.end(); }

// Krr⁻¹ cho các bậc tự do được giải phóng.
std::vector<double> releasedInverse(const Mat6& K, const std::vector<int>& R) {
    const int n = (int)R.size();
    std::vector<double> Krr(n * n);
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) Krr[a * n + b] = K[R[a]][R[b]];
    invert(Krr, n);
    return Krr;
}

void rotate(const BarProps& p, const Vec6& in, Vec6& out, bool toLocal) {
    const double c = p.c, s = toLocal ? p.s : -p.s;
    for (int k = 0; k < 2; ++k) {
        out[3 * k + 0] = c * in[3 * k + 0] + s * in[3 * k + 1];
        out[3 * k + 1] = -s * in[3 * k + 0] + c * in[3 * k + 1];
        out[3 * k + 2] = in[3 * k + 2];
    }
}

// Hàm dạng Hermite cho chuyển vị ngang của thanh
void hermite(double L, double x, double N[4]) {
    const double xi = x / L, xi2 = xi * xi, xi3 = xi2 * xi;
    N[0] = 1 - 3 * xi2 + 2 * xi3;
    N[1] = L * (xi - 2 * xi2 + xi3);
    N[2] = 3 * xi2 - 2 * xi3;
    N[3] = L * (xi3 - xi2);
}

double hermiteV(const BarProps& p, const Vec6& uL, double x) {
    double N[4];
    hermite(p.L, x, N);
    return N[0] * uL[1] + N[1] * uL[2] + N[2] * uL[4] + N[3] * uL[5];
}

const double kGauss4x[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
const double kGauss4w[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
const double kGauss3x[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
const double kGauss3w[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};

// Nội lực tại tiết diện x, xét cân bằng phần thanh bên trái [0, x].
// after = true: kể cả lực tập trung đặt đúng tại x (giá trị bên phải của bước nhảy).
void forcesAt(const BarProps& p, const BarLoads& l, const Vec6& fe, const Vec6& uL, double x, bool after, double& N,
              double& Q, double& M) {
    const double L = p.L;
    const double dqx = (l.qx2 - l.qx1) / L, dqy = (l.qy2 - l.qy1) / L;
    N = -fe[0] - (l.qx1 * x + dqx * x * x / 2);
    Q = fe[1] + l.qy1 * x + dqy * x * x / 2;
    M = -fe[2] + fe[1] * x + l.qy1 * x * x / 2 + dqy * x * x * x / 6;
    const double eps = 1e-9 * L;
    auto included = [&](double a) { return a < x - eps || (after && std::abs(a - x) <= eps); };
    for (const auto& pt : l.points)
        if (included(pt.a)) {
            N -= pt.Px;
            Q += pt.Py;
            M += pt.Py * (x - pt.a);
        }
    for (const auto& mo : l.moments)
        if (included(mo.a)) M -= mo.M;
    if (p.kf > 0 && x > 0) {
        // phản lực nền p(s) = −kf·v(s)
        double I0 = 0, I1 = 0;
        for (int g = 0; g < 4; ++g) {
            const double s = 0.5 * x * (kGauss4x[g] + 1), w = 0.5 * x * kGauss4w[g];
            const double ps = -p.kf * hermiteV(p, uL, s);
            I0 += w * ps;
            I1 += w * ps * (x - s);
        }
        Q += I0;
        M += I1;
    }
}

} // namespace

BarProps barProps(const Model& m, int e) {
    const Element& el = m.elements[e];
    const Node& a = m.nodes[el.n[0]];
    const Node& b = m.nodes[el.n[1]];
    BarProps p;
    const double dx = b.x - a.x, dy = b.y - a.y;
    p.L = std::hypot(dx, dy);
    if (p.L > 0) {
        p.c = dx / p.L;
        p.s = dy / p.L;
    }
    const Material& mat = m.materials[el.mat];
    const Section& sec = m.sections[el.sec];
    p.EA = mat.E * sec.A;
    p.EI = mat.E * sec.I;
    p.alpha = mat.alpha;
    p.h = sec.h;
    p.truss = el.type == ElemType::Truss;
    if (p.truss) {
        p.relI = p.relJ = true;
        p.kf = 0;
        p.hasBending = p.EI > 0;
        // Thanh dàn không có I: dùng EI giả định chỉ để quy tải trong nhịp về nút như dầm đơn giản.
        if (!p.hasBending) p.EI = std::max(p.EA * p.L * p.L * 1e-6, 1e-12);
    } else {
        p.relI = el.hingeI;
        p.relJ = el.hingeJ;
        p.kf = el.kf;
    }
    return p;
}

void addBarLoad(const BarProps& p, const ElemLoad& l, double factor, BarLoads& out) {
    double dx = 0, dy = 0; // véc tơ đơn vị phương tải trong hệ cục bộ
    switch (l.dir) {
    case LoadDir::LocalX: dx = 1; break;
    case LoadDir::LocalY: dy = 1; break;
    case LoadDir::GlobalX: dx = p.c; dy = -p.s; break;
    case LoadDir::GlobalY: dx = p.s; dy = p.c; break;
    }
    switch (l.kind) {
    case ElemLoadKind::Distributed: {
        double f = factor;
        if (l.projected) {
            if (l.dir == LoadDir::GlobalY) f *= std::abs(p.c);
            else if (l.dir == LoadDir::GlobalX) f *= std::abs(p.s);
        }
        out.qx1 += l.q1 * f * dx;
        out.qx2 += l.q2 * f * dx;
        out.qy1 += l.q1 * f * dy;
        out.qy2 += l.q2 * f * dy;
        break;
    }
    case ElemLoadKind::Point: {
        const double a = std::clamp(l.a, 0.0, p.L);
        out.points.push_back({a, l.q1 * factor * dx, l.q1 * factor * dy});
        break;
    }
    case ElemLoadKind::Moment: out.moments.push_back({std::clamp(l.a, 0.0, p.L), l.q1 * factor}); break;
    case ElemLoadKind::Temperature:
        out.T += l.q1 * factor;
        out.dT += l.q2 * factor;
        break;
    }
}

void addBarSelfWeight(const BarProps& p, double w, BarLoads& out) {
    // tải −w theo phương Y tổng thể, trên một đơn vị chiều dài thanh
    out.qx1 += -w * p.s;
    out.qx2 += -w * p.s;
    out.qy1 += -w * p.c;
    out.qy2 += -w * p.c;
}

void barLocalStiffness(const BarProps& p, Mat6& K) {
    for (auto& r : K) r.fill(0.0);
    const double L = p.L;
    const double k = p.EA / L;
    K[0][0] = k;
    K[0][3] = -k;
    K[3][0] = -k;
    K[3][3] = k;
    const double a = 12 * p.EI / (L * L * L), b = 6 * p.EI / (L * L), c = 4 * p.EI / L, d = 2 * p.EI / L;
    const int idx[4] = {1, 2, 4, 5};
    const double kb[4][4] = {{a, b, -a, b}, {b, c, -b, d}, {-a, -b, a, -b}, {b, d, -b, c}};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) K[idx[i]][idx[j]] = kb[i][j];
    if (p.kf > 0) {
        // ma trận nền Winkler tương thích: kf ∫ NᵀN dx
        const double f = p.kf * L / 420.0;
        const double kw[4][4] = {{156, 22 * L, 54, -13 * L},
                                 {22 * L, 4 * L * L, 13 * L, -3 * L * L},
                                 {54, 13 * L, 156, -22 * L},
                                 {-13 * L, -3 * L * L, -22 * L, 4 * L * L}};
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) K[idx[i]][idx[j]] += f * kw[i][j];
    }
}

void barFixedEndLoads(const BarProps& p, const BarLoads& l, Vec6& f) {
    f.fill(0.0);
    const double L = p.L;
    f[0] += L * (2 * l.qx1 + l.qx2) / 6;
    f[3] += L * (l.qx1 + 2 * l.qx2) / 6;
    f[1] += L * (7 * l.qy1 + 3 * l.qy2) / 20;
    f[2] += L * L * (3 * l.qy1 + 2 * l.qy2) / 60;
    f[4] += L * (3 * l.qy1 + 7 * l.qy2) / 20;
    f[5] += -L * L * (2 * l.qy1 + 3 * l.qy2) / 60;
    for (const auto& pt : l.points) {
        const double xi = pt.a / L;
        double N[4];
        hermite(L, pt.a, N);
        f[0] += pt.Px * (1 - xi);
        f[3] += pt.Px * xi;
        f[1] += pt.Py * N[0];
        f[2] += pt.Py * N[1];
        f[4] += pt.Py * N[2];
        f[5] += pt.Py * N[3];
    }
    for (const auto& mo : l.moments) {
        const double xi = mo.a / L;
        f[1] += mo.M * (-6 * xi + 6 * xi * xi) / L;
        f[2] += mo.M * (1 - 4 * xi + 3 * xi * xi);
        f[4] += mo.M * (6 * xi - 6 * xi * xi) / L;
        f[5] += mo.M * (3 * xi * xi - 2 * xi);
    }
    if (l.T != 0) {
        const double N0 = p.EA * p.alpha * l.T;
        f[0] -= N0;
        f[3] += N0;
    }
    if (l.dT != 0 && p.h > 0) {
        const double M0 = p.EI * p.alpha * l.dT / p.h;
        f[2] -= M0;
        f[5] += M0;
    }
}

void barGlobal(const BarProps& p, const BarLoads& l, Mat6& Kg, Vec6& fg) {
    Mat6 K;
    Vec6 f;
    barLocalStiffness(p, K);
    barFixedEndLoads(p, l, f);
    const std::vector<int> R = released(p);
    if (!R.empty()) {
        const int n = (int)R.size();
        const std::vector<double> Kinv = releasedInverse(K, R);
        Mat6 Kc = K;
        Vec6 fc = f;
        for (int i = 0; i < 6; ++i) {
            if (isReleased(R, i)) continue;
            for (int a = 0; a < n; ++a)
                for (int b = 0; b < n; ++b) {
                    const double kia = K[i][R[a]] * Kinv[a * n + b];
                    fc[i] -= kia * f[R[b]];
                    for (int j = 0; j < 6; ++j)
                        if (!isReleased(R, j)) Kc[i][j] -= kia * K[R[b]][j];
                }
        }
        for (int r : R) {
            fc[r] = 0;
            for (int j = 0; j < 6; ++j) Kc[r][j] = Kc[j][r] = 0;
        }
        K = Kc;
        f = fc;
    }
    // Kg = Tᵀ K T, fg = Tᵀ f; T: ma trận quay tổng thể → cục bộ. Trước hết KT = K·T.
    Mat6 KT;
    for (int i = 0; i < 6; ++i)
        for (int blk = 0; blk < 2; ++blk) {
            const double a = K[i][3 * blk], b = K[i][3 * blk + 1];
            KT[i][3 * blk] = a * p.c - b * p.s;
            KT[i][3 * blk + 1] = a * p.s + b * p.c;
            KT[i][3 * blk + 2] = K[i][3 * blk + 2];
        }
    for (int j = 0; j < 6; ++j) {
        Vec6 col, out;
        for (int i = 0; i < 6; ++i) col[i] = KT[i][j];
        rotate(p, col, out, false);
        for (int i = 0; i < 6; ++i) Kg[i][j] = out[i];
    }
    rotate(p, f, fg, false);
}

void barRecover(const BarProps& p, const BarLoads& l, const Vec6& uGlobal, int samples, BarResult& out) {
    Mat6 K;
    Vec6 f, uL;
    barLocalStiffness(p, K);
    barFixedEndLoads(p, l, f);
    rotate(p, uGlobal, uL, true);
    const std::vector<int> R = released(p);
    if (!R.empty()) {
        const int n = (int)R.size();
        const std::vector<double> Kinv = releasedInverse(K, R);
        std::vector<double> rhs(n);
        for (int a = 0; a < n; ++a) {
            rhs[a] = f[R[a]];
            for (int j = 0; j < 6; ++j)
                if (!isReleased(R, j)) rhs[a] -= K[R[a]][j] * uL[j];
        }
        for (int a = 0; a < n; ++a) {
            double v = 0;
            for (int b = 0; b < n; ++b) v += Kinv[a * n + b] * rhs[b];
            uL[R[a]] = v;
        }
    }
    Vec6 fe;
    for (int i = 0; i < 6; ++i) {
        double s = -f[i];
        for (int j = 0; j < 6; ++j) s += K[i][j] * uL[j];
        fe[i] = s;
    }
    for (int r : R) fe[r] = 0.0; // khớp: mômen đầu thanh bằng 0 (loại sai số làm tròn)
    out.endForces = fe;
    out.uLocal = uL;

    // Các điểm lấy mẫu: chia đều + hai phía của mỗi lực tập trung
    const double L = p.L, eps = 1e-9 * L;
    std::vector<std::pair<double, bool>> xs;
    std::vector<double> special;
    for (const auto& pt : l.points)
        if (pt.a > eps && pt.a < L - eps) special.push_back(pt.a);
    for (const auto& mo : l.moments)
        if (mo.a > eps && mo.a < L - eps) special.push_back(mo.a);
    for (int i = 0; i <= samples; ++i) {
        const double x = L * i / samples;
        bool clash = false;
        for (double a : special) clash |= std::abs(a - x) <= eps;
        if (!clash) xs.push_back({x, i != samples});
    }
    for (double a : special) {
        xs.push_back({a, false});
        xs.push_back({a, true});
    }
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());

    const size_t n = xs.size();
    out.x.resize(n);
    out.N.resize(n);
    out.Q.resize(n);
    out.M.resize(n);
    out.u.resize(n);
    out.v.resize(n);
    for (size_t k = 0; k < n; ++k) {
        out.x[k] = xs[k].first;
        forcesAt(p, l, fe, uL, xs[k].first, xs[k].second, out.N[k], out.Q[k], out.M[k]);
    }

    // Chuyển vị dọc thanh: tích phân độ cong κ = M/EI + κ_nhiệt và biến dạng dọc ε = N/EA + αt
    const double kappaT = (p.h > 0) ? p.alpha * l.dT / p.h : 0.0;
    const double epsT = p.alpha * l.T;
    double theta = uL[2];
    out.u[0] = uL[0];
    out.v[0] = uL[1];
    for (size_t k = 0; k + 1 < n; ++k) {
        const double x0 = out.x[k], x1 = out.x[k + 1], h = x1 - x0;
        if (h <= 0) {
            out.u[k + 1] = out.u[k];
            out.v[k + 1] = out.v[k];
            continue;
        }
        double ik = 0, iks = 0, ie = 0;
        for (int g = 0; g < 3; ++g) {
            const double s = x0 + 0.5 * h * (kGauss3x[g] + 1), w = 0.5 * h * kGauss3w[g];
            double Ns, Qs, Ms;
            forcesAt(p, l, fe, uL, s, true, Ns, Qs, Ms);
            const double kap = Ms / p.EI + kappaT;
            ik += w * kap;
            iks += w * kap * (x1 - s);
            ie += w * (Ns / p.EA + epsT);
        }
        out.v[k + 1] = out.v[k] + theta * h + iks;
        theta += ik;
        out.u[k + 1] = out.u[k] + ie;
    }
    const double ev = out.v[n - 1] - uL[4], eu = out.u[n - 1] - uL[3];
    for (size_t k = 0; k < n; ++k) {
        const double r = out.x[k] / L;
        if (p.hasBending) out.v[k] -= ev * r;
        else out.v[k] = uL[1] + (uL[4] - uL[1]) * r;
        out.u[k] -= eu * r;
    }
}

void barFoundationForce(const BarProps& p, const BarResult& r, double& Fx, double& Fy) {
    Fx = Fy = 0;
    if (p.kf <= 0) return;
    double fy = 0;
    for (int g = 0; g < 4; ++g) {
        const double s = 0.5 * p.L * (kGauss4x[g] + 1);
        fy += 0.5 * p.L * kGauss4w[g] * (-p.kf * hermiteV(p, r.uLocal, s));
    }
    Fx = -p.s * fy;
    Fy = p.c * fy;
}

// ---------------------------------------------------------------- phần tử phẳng

void planeD(const Material& mat, PlaneMode mode, double D[3][3]) {
    const double E = mat.E, nu = mat.nu;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) D[i][j] = 0;
    if (mode == PlaneMode::Stress) {
        const double c = E / (1 - nu * nu);
        D[0][0] = D[1][1] = c;
        D[0][1] = D[1][0] = c * nu;
        D[2][2] = c * (1 - nu) / 2;
    } else {
        const double c = E / ((1 + nu) * (1 - 2 * nu));
        D[0][0] = D[1][1] = c * (1 - nu);
        D[0][1] = D[1][0] = c * nu;
        D[2][2] = c * (1 - 2 * nu) / 2;
    }
}

namespace {

void tri3B(const double x[3], const double y[3], double B[3][6], double& A2) {
    const double b[3] = {y[1] - y[2], y[2] - y[0], y[0] - y[1]};
    const double c[3] = {x[2] - x[1], x[0] - x[2], x[1] - x[0]};
    A2 = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 6; ++j) B[i][j] = 0;
    if (A2 == 0) return;
    for (int i = 0; i < 3; ++i) {
        B[0][2 * i] = b[i] / A2;
        B[1][2 * i + 1] = c[i] / A2;
        B[2][2 * i] = c[i] / A2;
        B[2][2 * i + 1] = b[i] / A2;
    }
}

const double kXiN[4] = {-1, 1, 1, -1};
const double kEtaN[4] = {-1, -1, 1, 1};

struct QuadPoint {
    double B[3][8];
    double Ba[3][4];
    double detJ;
};

// Tính ma trận B (và B của mode không tương thích) tại (ξ, η). Trả về false nếu Jacobian suy biến.
bool quadPoint(const double x[4], const double y[4], double xi, double eta, const double J0inv[2][2], double detJ0,
               bool inc, QuadPoint& q) {
    double dNxi[4], dNeta[4];
    for (int i = 0; i < 4; ++i) {
        dNxi[i] = 0.25 * kXiN[i] * (1 + eta * kEtaN[i]);
        dNeta[i] = 0.25 * kEtaN[i] * (1 + xi * kXiN[i]);
    }
    double J[2][2] = {{0, 0}, {0, 0}};
    for (int i = 0; i < 4; ++i) {
        J[0][0] += dNxi[i] * x[i];
        J[0][1] += dNxi[i] * y[i];
        J[1][0] += dNeta[i] * x[i];
        J[1][1] += dNeta[i] * y[i];
    }
    q.detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    if (q.detJ == 0 || q.detJ * detJ0 <= 0) return false;
    const double Ji[2][2] = {{J[1][1] / q.detJ, -J[0][1] / q.detJ}, {-J[1][0] / q.detJ, J[0][0] / q.detJ}};
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 8; ++c) q.B[r][c] = 0;
        for (int c = 0; c < 4; ++c) q.Ba[r][c] = 0;
    }
    for (int i = 0; i < 4; ++i) {
        const double dx = Ji[0][0] * dNxi[i] + Ji[0][1] * dNeta[i];
        const double dy = Ji[1][0] * dNxi[i] + Ji[1][1] * dNeta[i];
        q.B[0][2 * i] = dx;
        q.B[1][2 * i + 1] = dy;
        q.B[2][2 * i] = dy;
        q.B[2][2 * i + 1] = dx;
    }
    if (inc) {
        // P1 = 1−ξ², P2 = 1−η²; đạo hàm theo J0 và nhân detJ0/detJ (Taylor) để thỏa mãn patch test
        const double r = detJ0 / q.detJ;
        const double dP[2][2] = {{-2 * xi, 0}, {0, -2 * eta}};
        for (int k = 0; k < 2; ++k) {
            const double dx = r * (J0inv[0][0] * dP[k][0] + J0inv[0][1] * dP[k][1]);
            const double dy = r * (J0inv[1][0] * dP[k][0] + J0inv[1][1] * dP[k][1]);
            q.Ba[0][2 * k] = dx;
            q.Ba[1][2 * k + 1] = dy;
            q.Ba[2][2 * k] = dy;
            q.Ba[2][2 * k + 1] = dx;
        }
    }
    return true;
}

bool quadCenter(const double x[4], const double y[4], double J0inv[2][2], double& detJ0) {
    double J[2][2] = {{0, 0}, {0, 0}};
    for (int i = 0; i < 4; ++i) {
        J[0][0] += 0.25 * kXiN[i] * x[i];
        J[0][1] += 0.25 * kXiN[i] * y[i];
        J[1][0] += 0.25 * kEtaN[i] * x[i];
        J[1][1] += 0.25 * kEtaN[i] * y[i];
    }
    detJ0 = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    if (detJ0 == 0) return false;
    J0inv[0][0] = J[1][1] / detJ0;
    J0inv[0][1] = -J[0][1] / detJ0;
    J0inv[1][0] = -J[1][0] / detJ0;
    J0inv[1][1] = J[0][0] / detJ0;
    return true;
}

const double kG2 = 0.5773502691896258; // 1/√3

// Lắp các ma trận Kuu (8×8), Kua (8×4), Kaa (4×4)
bool quadMatrices(const double x[4], const double y[4], const double D[3][3], double t, bool inc,
                  std::vector<double>& Kuu, std::vector<double>& Kua, std::vector<double>& Kaa) {
    double J0inv[2][2], detJ0;
    if (!quadCenter(x, y, J0inv, detJ0)) return false;
    Kuu.assign(64, 0.0);
    Kua.assign(32, 0.0);
    Kaa.assign(16, 0.0);
    for (int g = 0; g < 4; ++g) {
        QuadPoint q;
        if (!quadPoint(x, y, kXiN[g] * kG2, kEtaN[g] * kG2, J0inv, detJ0, inc, q)) return false;
        const double dv = t * std::abs(q.detJ);
        double DB[3][8], DBa[3][4];
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 8; ++c) DB[r][c] = D[r][0] * q.B[0][c] + D[r][1] * q.B[1][c] + D[r][2] * q.B[2][c];
            for (int c = 0; c < 4; ++c)
                DBa[r][c] = D[r][0] * q.Ba[0][c] + D[r][1] * q.Ba[1][c] + D[r][2] * q.Ba[2][c];
        }
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j)
                Kuu[i * 8 + j] += dv * (q.B[0][i] * DB[0][j] + q.B[1][i] * DB[1][j] + q.B[2][i] * DB[2][j]);
            for (int j = 0; j < 4; ++j)
                Kua[i * 4 + j] += dv * (q.B[0][i] * DBa[0][j] + q.B[1][i] * DBa[1][j] + q.B[2][i] * DBa[2][j]);
        }
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                Kaa[i * 4 + j] += dv * (q.Ba[0][i] * DBa[0][j] + q.Ba[1][i] * DBa[1][j] + q.Ba[2][i] * DBa[2][j]);
    }
    return true;
}

} // namespace

bool tri3Stiffness(const double x[3], const double y[3], const double D[3][3], double t, double K[6][6]) {
    double B[3][6], A2;
    tri3B(x, y, B, A2);
    if (std::abs(A2) < 1e-14) return false;
    const double vol = t * std::abs(A2) / 2;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            double s = 0;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) s += B[a][i] * D[a][b] * B[b][j];
            K[i][j] = vol * s;
        }
    return true;
}

void tri3Stress(const double x[3], const double y[3], const double D[3][3], const double u[6], double sig[3]) {
    double B[3][6], A2;
    tri3B(x, y, B, A2);
    double eps[3] = {0, 0, 0};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 6; ++c) eps[r] += B[r][c] * u[c];
    for (int r = 0; r < 3; ++r) sig[r] = D[r][0] * eps[0] + D[r][1] * eps[1] + D[r][2] * eps[2];
}

bool quad4Stiffness(const double x[4], const double y[4], const double D[3][3], double t, bool inc, double K[8][8]) {
    std::vector<double> Kuu, Kua, Kaa;
    if (!quadMatrices(x, y, D, t, inc, Kuu, Kua, Kaa)) return false;
    if (inc) {
        if (!invert(Kaa, 4)) return false;
        // K = Kuu − Kua Kaa⁻¹ Kau
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j) {
                double s = 0;
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b) s += Kua[i * 4 + a] * Kaa[a * 4 + b] * Kua[j * 4 + b];
                Kuu[i * 8 + j] -= s;
            }
    }
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) K[i][j] = Kuu[i * 8 + j];
    return true;
}

void quad4Stress(const double x[4], const double y[4], const double D[3][3], double t, bool inc, const double u[8],
                 double sig[4][3]) {
    std::vector<double> Kuu, Kua, Kaa;
    double alpha[4] = {0, 0, 0, 0};
    if (!quadMatrices(x, y, D, t, inc, Kuu, Kua, Kaa)) {
        for (int i = 0; i < 4; ++i) sig[i][0] = sig[i][1] = sig[i][2] = 0;
        return;
    }
    if (inc && invert(Kaa, 4)) {
        double rhs[4] = {0, 0, 0, 0};
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 8; ++i) rhs[a] += Kua[i * 4 + a] * u[i];
        for (int a = 0; a < 4; ++a) {
            double s = 0;
            for (int b = 0; b < 4; ++b) s += Kaa[a * 4 + b] * rhs[b];
            alpha[a] = -s;
        }
    }
    double J0inv[2][2], detJ0;
    quadCenter(x, y, J0inv, detJ0);
    double sg[4][3];
    for (int g = 0; g < 4; ++g) {
        QuadPoint q;
        quadPoint(x, y, kXiN[g] * kG2, kEtaN[g] * kG2, J0inv, detJ0, inc, q);
        double eps[3] = {0, 0, 0};
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 8; ++c) eps[r] += q.B[r][c] * u[c];
            if (inc)
                for (int c = 0; c < 4; ++c) eps[r] += q.Ba[r][c] * alpha[c];
        }
        for (int r = 0; r < 3; ++r) sg[g][r] = D[r][0] * eps[0] + D[r][1] * eps[1] + D[r][2] * eps[2];
    }
    // ngoại suy từ điểm Gauss ra nút
    const double s3 = std::sqrt(3.0);
    for (int i = 0; i < 4; ++i) {
        const double r = kXiN[i] * s3, s = kEtaN[i] * s3;
        for (int c = 0; c < 3; ++c) {
            double v = 0;
            for (int g = 0; g < 4; ++g) v += 0.25 * (1 + r * kXiN[g]) * (1 + s * kEtaN[g]) * sg[g][c];
            sig[i][c] = v;
        }
    }
}

void tri3BodyForce(const double x[3], const double y[3], double t, double bx, double by, double f[6]) {
    const double A = std::abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])) / 2;
    for (int i = 0; i < 3; ++i) {
        f[2 * i] = bx * t * A / 3;
        f[2 * i + 1] = by * t * A / 3;
    }
}

void quad4BodyForce(const double x[4], const double y[4], double t, double bx, double by, double f[8]) {
    for (int i = 0; i < 8; ++i) f[i] = 0;
    for (int g = 0; g < 4; ++g) {
        const double xi = kXiN[g] * kG2, eta = kEtaN[g] * kG2;
        double J[2][2] = {{0, 0}, {0, 0}}, N[4];
        for (int i = 0; i < 4; ++i) {
            N[i] = 0.25 * (1 + xi * kXiN[i]) * (1 + eta * kEtaN[i]);
            const double dxi = 0.25 * kXiN[i] * (1 + eta * kEtaN[i]);
            const double deta = 0.25 * kEtaN[i] * (1 + xi * kXiN[i]);
            J[0][0] += dxi * x[i];
            J[0][1] += dxi * y[i];
            J[1][0] += deta * x[i];
            J[1][1] += deta * y[i];
        }
        const double dv = t * std::abs(J[0][0] * J[1][1] - J[0][1] * J[1][0]);
        for (int i = 0; i < 4; ++i) {
            f[2 * i] += N[i] * bx * dv;
            f[2 * i + 1] += N[i] * by * dv;
        }
    }
}

StressInfo stressInfo(double sx, double sy, double txy, double sz) {
    StressInfo s;
    s.sx = sx;
    s.sy = sy;
    s.txy = txy;
    s.sz = sz;
    const double c = 0.5 * (sx + sy), R = std::hypot(0.5 * (sx - sy), txy);
    s.s1 = c + R;
    s.s2 = c - R;
    s.angle = 0.5 * std::atan2(2 * txy, sx - sy) * 180.0 / 3.14159265358979323846;
    s.vonMises = std::sqrt(std::max(0.0, sx * sx + sy * sy + sz * sz - sx * sy - sy * sz - sz * sx + 3 * txy * txy));
    return s;
}

} // namespace fem
