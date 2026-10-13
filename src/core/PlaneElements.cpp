#include "core/PlaneElements.h"

#include <algorithm>
#include <cmath>

#include "core/Elements.h"

namespace fem {

namespace {

constexpr double kPi = 3.14159265358979323846;

struct GP {
    double xi, eta, w;
};

std::vector<GP> gaussRule(ElemType t) {
    switch (t) {
    case ElemType::Tri3: return {{1.0 / 3, 1.0 / 3, 0.5}};
    case ElemType::Tri6: {
        const double a = 0.445948490915965, b = 0.091576213509771;
        const double w1 = 0.223381589678011 / 2, w2 = 0.109951743655322 / 2;
        return {{a, a, w1}, {1 - 2 * a, a, w1}, {a, 1 - 2 * a, w1}, {b, b, w2}, {1 - 2 * b, b, w2}, {b, 1 - 2 * b, w2}};
    }
    case ElemType::Quad4: {
        const double g = 0.5773502691896258;
        return {{-g, -g, 1}, {g, -g, 1}, {g, g, 1}, {-g, g, 1}};
    }
    case ElemType::Quad8: {
        const double g = std::sqrt(0.6), p[3] = {-g, 0, g}, w[3] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
        std::vector<GP> r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.push_back({p[i], p[j], w[i] * w[j]});
        return r;
    }
    default: return {};
    }
}

void axisymD(const Material& m, double D[4][4]) {
    const double E = m.E, nu = m.nu, c = E / ((1 + nu) * (1 - 2 * nu));
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) D[i][j] = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) D[i][j] = c * (i == j ? 1 - nu : nu);
    D[3][3] = c * (1 - 2 * nu) / 2;
}

// Đánh giá tại (ξ, η): N, dN/dx, dN/dy, |detJ|, r. Trả về false nếu Jacobian suy biến.
bool evalPoint(ElemType t, int n, const double* x, const double* y, double xi, double eta, double N[8], double dx[8],
               double dy[8], double& detJ, double& r) {
    double dxi[8], deta[8];
    planeShape(t, xi, eta, N, dxi, deta);
    double J00 = 0, J01 = 0, J10 = 0, J11 = 0;
    r = 0;
    for (int i = 0; i < n; ++i) {
        J00 += dxi[i] * x[i];
        J01 += dxi[i] * y[i];
        J10 += deta[i] * x[i];
        J11 += deta[i] * y[i];
        r += N[i] * x[i];
    }
    detJ = J00 * J11 - J01 * J10;
    if (std::abs(detJ) < 1e-300) return false;
    for (int i = 0; i < n; ++i) {
        dx[i] = (J11 * dxi[i] - J01 * deta[i]) / detJ;
        dy[i] = (-J10 * dxi[i] + J00 * deta[i]) / detJ;
    }
    return true;
}

bool useLegacy(const Element& el, PlaneMode mode) {
    return mode != PlaneMode::Axisymmetric && (el.type == ElemType::Tri3 || el.type == ElemType::Quad4);
}

} // namespace

bool isTriangle(ElemType t) { return t == ElemType::Tri3 || t == ElemType::Tri6; }

void planeShape(ElemType t, double xi, double eta, double N[8], double dxi[8], double deta[8]) {
    for (int i = 0; i < 8; ++i) N[i] = dxi[i] = deta[i] = 0;
    switch (t) {
    case ElemType::Tri3:
        N[0] = 1 - xi - eta;
        N[1] = xi;
        N[2] = eta;
        dxi[0] = -1, dxi[1] = 1;
        deta[0] = -1, deta[2] = 1;
        break;
    case ElemType::Tri6: {
        const double L1 = 1 - xi - eta, L2 = xi, L3 = eta;
        const double a1 = -1, a2 = 1, a3 = 0, b1 = -1, b2 = 0, b3 = 1; // dL/dξ, dL/dη
        N[0] = L1 * (2 * L1 - 1);
        N[1] = L2 * (2 * L2 - 1);
        N[2] = L3 * (2 * L3 - 1);
        N[3] = 4 * L1 * L2;
        N[4] = 4 * L2 * L3;
        N[5] = 4 * L3 * L1;
        dxi[0] = (4 * L1 - 1) * a1, deta[0] = (4 * L1 - 1) * b1;
        dxi[1] = (4 * L2 - 1) * a2, deta[1] = (4 * L2 - 1) * b2;
        dxi[2] = (4 * L3 - 1) * a3, deta[2] = (4 * L3 - 1) * b3;
        dxi[3] = 4 * (a1 * L2 + L1 * a2), deta[3] = 4 * (b1 * L2 + L1 * b2);
        dxi[4] = 4 * (a2 * L3 + L2 * a3), deta[4] = 4 * (b2 * L3 + L2 * b3);
        dxi[5] = 4 * (a3 * L1 + L3 * a1), deta[5] = 4 * (b3 * L1 + L3 * b1);
        break;
    }
    case ElemType::Quad4: {
        static const double xs[4] = {-1, 1, 1, -1}, ys[4] = {-1, -1, 1, 1};
        for (int i = 0; i < 4; ++i) {
            N[i] = 0.25 * (1 + xi * xs[i]) * (1 + eta * ys[i]);
            dxi[i] = 0.25 * xs[i] * (1 + eta * ys[i]);
            deta[i] = 0.25 * ys[i] * (1 + xi * xs[i]);
        }
        break;
    }
    case ElemType::Quad8: {
        static const double xs[4] = {-1, 1, 1, -1}, ys[4] = {-1, -1, 1, 1};
        for (int i = 0; i < 4; ++i) {
            const double a = xs[i], b = ys[i];
            N[i] = 0.25 * (1 + a * xi) * (1 + b * eta) * (a * xi + b * eta - 1);
            dxi[i] = 0.25 * a * (1 + b * eta) * (2 * a * xi + b * eta);
            deta[i] = 0.25 * b * (1 + a * xi) * (a * xi + 2 * b * eta);
        }
        N[4] = 0.5 * (1 - xi * xi) * (1 - eta), dxi[4] = -xi * (1 - eta), deta[4] = -0.5 * (1 - xi * xi);
        N[5] = 0.5 * (1 + xi) * (1 - eta * eta), dxi[5] = 0.5 * (1 - eta * eta), deta[5] = -(1 + xi) * eta;
        N[6] = 0.5 * (1 - xi * xi) * (1 + eta), dxi[6] = -xi * (1 + eta), deta[6] = 0.5 * (1 - xi * xi);
        N[7] = 0.5 * (1 - xi) * (1 - eta * eta), dxi[7] = -0.5 * (1 - eta * eta), deta[7] = -(1 - xi) * eta;
        break;
    }
    default: break;
    }
}

void planeNodeNatural(ElemType t, int k, double& xi, double& eta) {
    static const double tx[6] = {0, 1, 0, 0.5, 0.5, 0}, ty[6] = {0, 0, 1, 0, 0.5, 0.5};
    static const double qx[8] = {-1, 1, 1, -1, 0, 1, 0, -1}, qy[8] = {-1, -1, 1, 1, -1, 0, 1, 0};
    if (isTriangle(t)) {
        xi = tx[k];
        eta = ty[k];
    } else {
        xi = qx[k];
        eta = qy[k];
    }
}

bool planeStiffness(const Element& el, const double* x, const double* y, const Material& mat, PlaneMode mode,
                    std::vector<double>& K) {
    const int n = el.nodeCount(), nd = 2 * n;
    K.assign((size_t)nd * nd, 0.0);
    if (useLegacy(el, mode)) {
        double D[3][3];
        planeD(mat, mode, D);
        if (el.type == ElemType::Tri3) {
            double k[6][6];
            if (!tri3Stiffness(x, y, D, el.t, k)) return false;
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) K[i * 6 + j] = k[i][j];
        } else {
            double k[8][8];
            if (!quad4Stiffness(x, y, D, el.t, el.incompatible, k)) return false;
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < 8; ++j) K[i * 8 + j] = k[i][j];
        }
        return true;
    }
    const bool ax = mode == PlaneMode::Axisymmetric;
    double D[4][4];
    int ns = 3;
    if (ax) {
        axisymD(mat, D);
        ns = 4;
    } else {
        double d3[3][3];
        planeD(mat, mode, d3);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) D[i][j] = i < 3 && j < 3 ? d3[i][j] : 0;
    }
    double sgn = 0;
    for (const GP& g : gaussRule(el.type)) {
        double N[8], dx[8], dy[8], detJ, r;
        if (!evalPoint(el.type, n, x, y, g.xi, g.eta, N, dx, dy, detJ, r)) return false;
        if (sgn == 0) sgn = detJ > 0 ? 1 : -1;
        else if (detJ * sgn <= 0) return false; // phần tử méo (Jacobian đổi dấu)
        if (ax && r <= 0) return false;
        std::vector<double> B((size_t)ns * nd, 0.0);
        for (int i = 0; i < n; ++i) {
            if (ax) {
                B[0 * nd + 2 * i] = dx[i];
                B[1 * nd + 2 * i + 1] = dy[i];
                B[2 * nd + 2 * i] = N[i] / r;
                B[3 * nd + 2 * i] = dy[i];
                B[3 * nd + 2 * i + 1] = dx[i];
            } else {
                B[0 * nd + 2 * i] = dx[i];
                B[1 * nd + 2 * i + 1] = dy[i];
                B[2 * nd + 2 * i] = dy[i];
                B[2 * nd + 2 * i + 1] = dx[i];
            }
        }
        const double w = g.w * std::abs(detJ) * (ax ? 2 * kPi * r : el.t);
        std::vector<double> DB((size_t)ns * nd, 0.0);
        for (int a = 0; a < ns; ++a)
            for (int c = 0; c < nd; ++c) {
                double s = 0;
                for (int b = 0; b < ns; ++b) s += D[a][b] * B[b * nd + c];
                DB[a * nd + c] = s;
            }
        for (int i = 0; i < nd; ++i)
            for (int j = 0; j < nd; ++j) {
                double s = 0;
                for (int a = 0; a < ns; ++a) s += B[a * nd + i] * DB[a * nd + j];
                K[(size_t)i * nd + j] += w * s;
            }
    }
    return true;
}

void planeStress(const Element& el, const double* x, const double* y, const Material& mat, PlaneMode mode, const double* u,
                 std::array<std::array<double, 4>, 8>& sig) {
    const int n = el.nodeCount();
    for (auto& s : sig) s = {0, 0, 0, 0};
    auto out = [&](double sx, double sy, double txy) {
        return std::array<double, 4>{sx, sy, txy, mode == PlaneMode::Strain ? mat.nu * (sx + sy) : 0.0};
    };
    if (useLegacy(el, mode)) {
        double D[3][3];
        planeD(mat, mode, D);
        if (el.type == ElemType::Tri3) {
            double s[3];
            tri3Stress(x, y, D, u, s);
            for (int k = 0; k < 3; ++k) sig[k] = out(s[0], s[1], s[2]);
        } else {
            double s[4][3];
            quad4Stress(x, y, D, el.t, el.incompatible, u, s);
            for (int k = 0; k < 4; ++k) sig[k] = out(s[k][0], s[k][1], s[k][2]);
        }
        return;
    }
    const bool ax = mode == PlaneMode::Axisymmetric;
    double D3[3][3], D4[4][4];
    if (ax) axisymD(mat, D4);
    else planeD(mat, mode, D3);
    double xmax = 1e-12;
    for (int i = 0; i < n; ++i) xmax = std::max(xmax, std::abs(x[i]) + std::abs(y[i]));
    for (int k = 0; k < n; ++k) {
        double xi, eta, N[8], dx[8], dy[8], detJ, r;
        planeNodeNatural(el.type, k, xi, eta);
        if (!evalPoint(el.type, n, x, y, xi, eta, N, dx, dy, detJ, r)) continue;
        double ex = 0, ey = 0, gxy = 0, eth = 0;
        for (int i = 0; i < n; ++i) {
            ex += dx[i] * u[2 * i];
            ey += dy[i] * u[2 * i + 1];
            gxy += dy[i] * u[2 * i] + dx[i] * u[2 * i + 1];
            eth += N[i] * u[2 * i];
        }
        if (ax) {
            eth = r > 1e-9 * xmax ? eth / r : ex; // trên trục đối xứng: εθ = εr
            const double e[4] = {ex, ey, eth, gxy};
            double s[4];
            for (int a = 0; a < 4; ++a) s[a] = D4[a][0] * e[0] + D4[a][1] * e[1] + D4[a][2] * e[2] + D4[a][3] * e[3];
            sig[k] = {s[0], s[1], s[3], s[2]}; // σr, σz, τrz, σθ
        } else {
            const double e[3] = {ex, ey, gxy};
            double s[3];
            for (int a = 0; a < 3; ++a) s[a] = D3[a][0] * e[0] + D3[a][1] * e[1] + D3[a][2] * e[2];
            sig[k] = out(s[0], s[1], s[2]);
        }
    }
}

void planeBodyForce(const Element& el, const double* x, const double* y, PlaneMode mode, double bx, double by,
                    std::vector<double>& f) {
    const int n = el.nodeCount();
    f.assign(2 * n, 0.0);
    if (useLegacy(el, mode)) {
        if (el.type == ElemType::Tri3) tri3BodyForce(x, y, el.t, bx, by, f.data());
        else quad4BodyForce(x, y, el.t, bx, by, f.data());
        return;
    }
    const bool ax = mode == PlaneMode::Axisymmetric;
    for (const GP& g : gaussRule(el.type == ElemType::Tri3 ? ElemType::Tri6 : el.type)) {
        double N[8], dx[8], dy[8], detJ, r;
        if (!evalPoint(el.type, n, x, y, g.xi, g.eta, N, dx, dy, detJ, r)) return;
        const double w = g.w * std::abs(detJ) * (ax ? 2 * kPi * r : el.t);
        for (int i = 0; i < n; ++i) {
            f[2 * i] += N[i] * bx * w;
            f[2 * i + 1] += N[i] * by * w;
        }
    }
}

void planeLumpedMass(const Element& el, const double* x, const double* y, PlaneMode mode, double rho, std::vector<double>& m) {
    const int n = el.nodeCount();
    m.assign(n, 0.0);
    const bool ax = mode == PlaneMode::Axisymmetric;
    double total = 0, diag = 0;
    const ElemType rule = el.type == ElemType::Tri3 ? ElemType::Tri6 : el.type == ElemType::Quad4 ? ElemType::Quad8 : el.type;
    for (const GP& g : gaussRule(rule)) {
        double N[8], dx[8], dy[8], detJ, r;
        if (!evalPoint(el.type, n, x, y, g.xi, g.eta, N, dx, dy, detJ, r)) return;
        const double w = g.w * std::abs(detJ) * (ax ? 2 * kPi * r : el.t) * rho;
        total += w;
        for (int i = 0; i < n; ++i) m[i] += w * N[i] * N[i];
    }
    for (double v : m) diag += v;
    if (diag > 0)
        for (double& v : m) v *= total / diag; // HRZ: giữ tổng khối lượng
}

} // namespace fem
