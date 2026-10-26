#include "core/Fluid.h"
#include "core/TriMesher.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <unordered_map>

namespace fem {

namespace {

constexpr double kPi = 3.14159265358979323846;

std::string fmt(const char* f, int a) {
    char b[256];
    std::snprintf(b, sizeof b, f, a);
    return b;
}

bool isDirichlet(FluidBcKind k) { return k == FluidBcKind::Wall || k == FluidBcKind::MovingWall || k == FluidBcKind::Inlet; }

} // namespace

const char* fluidSideName(FluidSide s) {
    switch (s) {
    case FluidSide::Left: return "Trái";
    case FluidSide::Right: return "Phải";
    case FluidSide::Bottom: return "Dưới";
    case FluidSide::Top: return "Trên";
    case FluidSide::Inner: return "Trong (vật cản)";
    default: return "?";
    }
}

const char* fluidBcName(FluidBcKind k) {
    switch (k) {
    case FluidBcKind::Wall: return "Thành cứng";
    case FluidBcKind::MovingWall: return "Thành di động";
    case FluidBcKind::Inlet: return "Cửa vào";
    case FluidBcKind::Outlet: return "Cửa ra (p = 0)";
    case FluidBcKind::Symmetry: return "Đối xứng";
    case FluidBcKind::Opening: return "Thông khí (p = 0)";
    }
    return "?";
}

const char* fluidFieldName(FluidField f) {
    switch (f) {
    case FluidField::Speed: return "|u| (m/s)";
    case FluidField::U: return "u (m/s)";
    case FluidField::V: return "v (m/s)";
    case FluidField::P: return "p (Pa)";
    case FluidField::Cp: return "Cp";
    case FluidField::Vorticity: return "ω (1/s)";
    case FluidField::PTotal: return "p toàn phần (Pa)";
    case FluidField::Alpha: return "α tỉ lệ nước";
    case FluidField::TurbK: return "k rối (m²/s²)";
    }
    return "?";
}

// ============================================================================ lưới
std::string buildFluidMesh(const Model& m, FluidMesh& M, double flowAngle) {
    M = FluidMesh{};
    M.ca = std::cos(flowAngle);
    M.sa = std::sin(flowAngle);
    for (const Node& n : m.nodes) {
        M.px.push_back(n.x);
        M.py.push_back(n.y);
    }
    // ---- ô
    M.cellStart.push_back(0);
    for (size_t e = 0; e < m.elements.size(); ++e) {
        const Element& el = m.elements[e];
        if (!el.isPlane()) continue;
        const int nc = el.cornerCount();
        int ids[4];
        for (int k = 0; k < nc; ++k) ids[k] = el.n[k];
        double A = 0, cx = 0, cy = 0;
        for (int k = 0; k < nc; ++k) {
            const Node& a = m.nodes[ids[k]];
            const Node& b = m.nodes[ids[(k + 1) % nc]];
            const double cr = a.x * b.y - b.x * a.y;
            A += cr;
            cx += (a.x + b.x) * cr;
            cy += (a.y + b.y) * cr;
        }
        A *= 0.5;
        if (std::abs(A) < 1e-14) return fmt("Phần tử %d có diện tích bằng 0.", (int)e + 1);
        cx /= 6 * A;
        cy /= 6 * A;
        if (A < 0) std::reverse(ids, ids + nc);
        for (int k = 0; k < nc; ++k) M.cellNodes.push_back(ids[k]);
        M.cellStart.push_back((int)M.cellNodes.size());
        M.cellElem.push_back((int)e);
        M.cx.push_back(cx);
        M.cy.push_back(cy);
        M.vol.push_back(std::abs(A));
    }
    M.nCells = (int)M.cellElem.size();
    if (M.nCells == 0) return "Mô hình chưa có phần tử phẳng: tạo lưới miền chất lỏng trước.";

    // ---- mặt: mỗi cạnh nối hai nút góc
    struct TmpFace {
        int a, b, owner, neigh;
    };
    std::vector<TmpFace> tf;
    std::unordered_map<long long, int> edge;
    const long long NN = (long long)m.nodes.size();
    for (int c = 0; c < M.nCells; ++c) {
        const int s = M.cellStart[c], n = M.cellStart[c + 1] - s;
        for (int k = 0; k < n; ++k) {
            const int a = M.cellNodes[s + k], b = M.cellNodes[s + (k + 1) % n];
            const long long key = (long long)std::min(a, b) * NN + std::max(a, b);
            auto it = edge.find(key);
            if (it == edge.end()) {
                edge[key] = (int)tf.size();
                tf.push_back({a, b, c, -1});
            } else {
                TmpFace& f = tf[it->second];
                if (f.neigh >= 0 || f.a == a)
                    return fmt("Phần tử %d chồng lấn lên phần tử khác (một cạnh thuộc hơn hai phần tử hoặc hai phần tử "
                               "trùng nhau).",
                               M.cellElem[c] + 1);
                f.neigh = c;
            }
        }
    }
    // mặt trong: owner < neigh, sắp theo (owner, neigh) để dùng tiền điều kiện DIC; mặt biên đứng sau
    for (TmpFace& f : tf)
        if (f.neigh >= 0 && f.owner > f.neigh) {
            std::swap(f.owner, f.neigh);
            std::swap(f.a, f.b);
        }
    std::stable_sort(tf.begin(), tf.end(), [](const TmpFace& p, const TmpFace& q) {
        const bool pi = p.neigh >= 0, qi = q.neigh >= 0;
        if (pi != qi) return pi;
        if (!pi) return p.owner < q.owner;
        return p.owner != q.owner ? p.owner < q.owner : p.neigh < q.neigh;
    });
    M.nFaces = (int)tf.size();
    M.nInternal = (int)std::count_if(tf.begin(), tf.end(), [](const TmpFace& f) { return f.neigh >= 0; });

    // ---- hình chữ nhật bao
    M.x0 = M.y0 = std::numeric_limits<double>::max();
    M.x1 = M.y1 = -std::numeric_limits<double>::max();
    for (int id : M.cellNodes) {
        M.x0 = std::min(M.x0, m.nodes[id].x);
        M.x1 = std::max(M.x1, m.nodes[id].x);
        M.y0 = std::min(M.y0, m.nodes[id].y);
        M.y1 = std::max(M.y1, m.nodes[id].y);
    }
    // hình chữ nhật bao trong hệ tọa độ dòng chảy (để nhận diện các đoạn biên)
    auto lx = [&](int n) { return M.ca * m.nodes[n].x + M.sa * m.nodes[n].y; };
    auto ly = [&](int n) { return -M.sa * m.nodes[n].x + M.ca * m.nodes[n].y; };
    double X0 = std::numeric_limits<double>::max(), X1 = -X0, Y0 = X0, Y1 = -X0;
    for (int id : M.cellNodes) {
        X0 = std::min(X0, lx(id)), X1 = std::max(X1, lx(id));
        Y0 = std::min(Y0, ly(id)), Y1 = std::max(Y1, ly(id));
    }
    const double tol = 1e-7 * std::max({X1 - X0, Y1 - Y0, 1e-9});

    for (const TmpFace& t : tf) {
        const Node& A = m.nodes[t.a];
        const Node& B = m.nodes[t.b];
        M.owner.push_back(t.owner);
        M.neigh.push_back(t.neigh);
        M.fa.push_back(t.a);
        M.fb.push_back(t.b);
        // cạnh a → b đi ngược chiều kim đồng hồ quanh owner ⇒ pháp tuyến ngoài (dy, −dx)
        const double sx = B.y - A.y, sy = -(B.x - A.x);
        M.sx.push_back(sx);
        M.sy.push_back(sy);
        const double fx = 0.5 * (A.x + B.x), fy = 0.5 * (A.y + B.y);
        M.fx.push_back(fx);
        M.fy.push_back(fy);
        if (t.neigh >= 0) {
            const double L = std::hypot(sx, sy);
            const double nx = sx / L, ny = sy / L;
            const double dP = std::abs((fx - M.cx[t.owner]) * nx + (fy - M.cy[t.owner]) * ny);
            const double dN = std::abs((M.cx[t.neigh] - fx) * nx + (M.cy[t.neigh] - fy) * ny);
            M.w.push_back(dP + dN > 0 ? dN / (dP + dN) : 0.5);
            M.side.push_back(-1);
        } else {
            M.w.push_back(1.0);
            FluidSide s = FluidSide::Inner;
            const double ax = lx(t.a), ay = ly(t.a), bx = lx(t.b), by = ly(t.b);
            if (std::abs(ax - X0) <= tol && std::abs(bx - X0) <= tol) s = FluidSide::Left;
            else if (std::abs(ax - X1) <= tol && std::abs(bx - X1) <= tol) s = FluidSide::Right;
            else if (std::abs(ay - Y0) <= tol && std::abs(by - Y0) <= tol) s = FluidSide::Bottom;
            else if (std::abs(ay - Y1) <= tol && std::abs(by - Y1) <= tol) s = FluidSide::Top;
            M.side.push_back((int)s);
            ++M.sideFaces[(int)s];
        }
    }
    // ---- mặt của từng ô
    std::vector<int> cnt(M.nCells + 1, 0);
    for (int f = 0; f < M.nFaces; ++f) {
        ++cnt[M.owner[f] + 1];
        if (M.neigh[f] >= 0) ++cnt[M.neigh[f] + 1];
    }
    for (int c = 0; c < M.nCells; ++c) cnt[c + 1] += cnt[c];
    M.cfStart = cnt;
    M.cfFace.assign(cnt.back(), 0);
    for (int f = 0; f < M.nFaces; ++f) {
        M.cfFace[cnt[M.owner[f]]++] = f;
        if (M.neigh[f] >= 0) M.cfFace[cnt[M.neigh[f]]++] = f;
    }
    return {};
}

// ============================================================================ bộ giải
namespace {

// hằng số k-ε chuẩn và hàm tường
constexpr double kCmu = 0.09, kC1 = 1.44, kC2 = 1.92, kSigK = 1.0, kSigE = 1.3, kKappa = 0.41, kE = 9.8;
constexpr double kYLam = 11.06; // y* chuyển tiếp lớp nhớt – lớp logarit (hàm tường "co giãn": y* ≥ 11,06)

bool isWall(FluidBcKind k) { return k == FluidBcKind::Wall || k == FluidBcKind::MovingWall; }

using Poly = std::vector<std::array<double, 2>>;

// Phần đa giác lồi thỏa a·x + b·y ≤ c (Sutherland–Hodgman)
Poly clipHalf(const Poly& P, double a, double b, double c) {
    Poly r;
    for (size_t i = 0; i < P.size(); ++i) {
        const auto& p = P[i];
        const auto& q = P[(i + 1) % P.size()];
        const double fp = a * p[0] + b * p[1] - c, fq = a * q[0] + b * q[1] - c;
        if (fp <= 0) r.push_back(p);
        if ((fp < 0) != (fq < 0) && fp != fq) {
            const double t = fp / (fp - fq);
            r.push_back({p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])});
        }
    }
    return r;
}

double polyArea(const Poly& P) {
    double A = 0;
    for (size_t i = 0; i < P.size(); ++i) A += P[i][0] * P[(i + 1) % P.size()][1] - P[(i + 1) % P.size()][0] * P[i][1];
    return 0.5 * std::abs(A);
}

} // namespace

bool FluidSolver::init(const Model& m, const FluidParams& p) {
    P_ = p;
    err_ = buildFluidMesh(m, M_, p.flowAngle);
    if (!err_.empty()) return false;
    if (P_.rho <= 0 || P_.mu <= 0 || (P_.vof && (P_.rho2 <= 0 || P_.mu2 <= 0))) {
        err_ = "Khối lượng riêng và độ nhớt phải > 0.";
        return false;
    }
    if (P_.vof) {
        // mặt thoáng luôn là bài toán theo thời gian; bước thời gian thay đổi nên dùng Euler ẩn
        P_.transient = true;
        P_.timeScheme = TimeScheme::Euler;
    }
    if (P_.transient && P_.dt <= 0) {
        err_ = "Bước thời gian Δt phải > 0.";
        return false;
    }
    skx_.assign(M_.nInternal, 0);
    sky_.assign(M_.nInternal, 0);
    for (int f = 0; f < M_.nInternal; ++f) {
        const int P = M_.owner[f], N = M_.neigh[f];
        const double w = M_.w[f];
        skx_[f] = M_.fx[f] - (w * M_.cx[P] + (1 - w) * M_.cx[N]);
        sky_[f] = M_.fy[f] - (w * M_.cy[P] + (1 - w) * M_.cy[N]);
    }
    // ma trận bình phương tối thiểu Σ w·d·dᵀ của từng ô, nghịch đảo sẵn
    lsq_.assign(4 * (size_t)M_.nCells, 0);
    for (int f = 0; f < M_.nFaces; ++f) {
        const int P = M_.owner[f], N = M_.neigh[f];
        const double dx = (N >= 0 ? M_.cx[N] : M_.fx[f]) - M_.cx[P], dy = (N >= 0 ? M_.cy[N] : M_.fy[f]) - M_.cy[P];
        const double w = 1.0 / (dx * dx + dy * dy);
        for (int c : {P, N}) {
            if (c < 0) continue;
            lsq_[4 * c] += w * dx * dx;
            lsq_[4 * c + 1] += w * dx * dy;
            lsq_[4 * c + 3] += w * dy * dy;
        }
    }
    for (int c = 0; c < M_.nCells; ++c) {
        double* G = &lsq_[4 * c];
        const double a = G[0], b = G[1], d = G[3], det = a * d - b * b;
        if (std::abs(det) < 1e-12 * (a * d + 1e-300)) {
            err_ = fmt("Ô %d không đủ lân cận để tính gradient.", c + 1);
            return false;
        }
        G[0] = d / det, G[1] = -b / det, G[2] = -b / det, G[3] = a / det;
    }
    // ma trận dựng lại véc tơ tại tâm ô từ các thành phần pháp tuyến trên mặt (cân bằng lực của VOF)
    recon_.assign(4 * (size_t)M_.nCells, 0);
    for (int f = 0; f < M_.nFaces; ++f) {
        const double L = std::hypot(M_.sx[f], M_.sy[f]);
        for (int c : {M_.owner[f], M_.neigh[f]}) {
            if (c < 0) continue;
            recon_[4 * c] += M_.sx[f] * M_.sx[f] / L;
            recon_[4 * c + 1] += M_.sx[f] * M_.sy[f] / L;
            recon_[4 * c + 3] += M_.sy[f] * M_.sy[f] / L;
        }
    }
    for (int c = 0; c < M_.nCells; ++c) {
        double* G = &recon_[4 * c];
        const double a = G[0], b = G[1], d = G[3], det = a * d - b * b;
        G[0] = d / det, G[1] = -b / det, G[2] = -b / det, G[3] = a / det;
    }
    P_.alphaU = std::clamp(P_.alphaU, 0.05, 1.0);
    P_.alphaP = std::clamp(P_.alphaP, 0.01, 1.0);
    P_.outerIters = std::max(P_.outerIters, 1);

    // phạm vi từng đoạn biên (để dựng phân bố parabol); có biên áp suất cho trước hay không; có đường ra hay không
    hasOutlet_ = false;
    bool anyOutflow = false;
    for (int s = 0; s < (int)FluidSide::Count; ++s) {
        sideMin_[s] = std::numeric_limits<double>::max();
        sideMax_[s] = -std::numeric_limits<double>::max();
    }
    double innerY0 = std::numeric_limits<double>::max(), innerY1 = -innerY0;
    for (int f = M_.nInternal; f < M_.nFaces; ++f) {
        const int s = M_.side[f];
        const bool vert = s == (int)FluidSide::Left || s == (int)FluidSide::Right;
        const int na = M_.fa[f], nb = M_.fb[f];
        const double lxa = M_.ca * M_.px[na] + M_.sa * M_.py[na], lya = -M_.sa * M_.px[na] + M_.ca * M_.py[na];
        const double lxb = M_.ca * M_.px[nb] + M_.sa * M_.py[nb], lyb = -M_.sa * M_.px[nb] + M_.ca * M_.py[nb];
        const double a = vert ? lya : lxa, b = vert ? lyb : lxb;
        sideMin_[s] = std::min({sideMin_[s], a, b});
        sideMax_[s] = std::max({sideMax_[s], a, b});
        if (s == (int)FluidSide::Inner) {
            innerY0 = std::min({innerY0, lya, lyb});
            innerY1 = std::max({innerY1, lya, lyb});
        }
        if (fixedPressure(f)) hasOutlet_ = true;
        if (P_.bc[s].kind == FluidBcKind::Outlet || P_.bc[s].kind == FluidBcKind::Opening) anyOutflow = true;
    }

    // vận tốc, chiều dài tham chiếu
    double Uin = 0, Uwall = 0, Lin = 0, Lwall = 0;
    for (int s = 0; s < (int)FluidSide::Count; ++s) {
        if (!M_.sideFaces[s]) continue;
        const FluidBc& b = P_.bc[s];
        const double U = std::hypot(b.u, b.v), L = sideMax_[s] - sideMin_[s];
        if (b.kind == FluidBcKind::Inlet && U * (b.parabolic ? 2.0 / 3 : 1) > Uin) {
            Uin = U * (b.parabolic ? 2.0 / 3 : 1);
            Lin = L;
        }
        if (b.kind == FluidBcKind::MovingWall && U > Uwall) {
            Uwall = U;
            Lwall = L;
        }
    }
    Uref_ = P_.Uref > 0 ? P_.Uref : Uin > 0 ? Uin : Uwall > 0 ? Uwall : 1.0;
    if (P_.Lref > 0) Lref_ = P_.Lref;
    else if (M_.sideFaces[(int)FluidSide::Inner] && innerY1 > innerY0) Lref_ = innerY1 - innerY0;
    else if (Uin > 0) Lref_ = Lin;
    else if (Uwall > 0) Lref_ = Lwall;
    else Lref_ = M_.y1 - M_.y0;

    const int nc = M_.nCells, nf = M_.nFaces;
    u_.assign(nc, 0);
    v_.assign(nc, 0);
    p_.assign(nc, 0);
    F_.assign(nf, 0);
    gpx_.assign(nc, 0);
    gpy_.assign(nc, 0);
    pgx_.assign(nc, 0);
    pgy_.assign(nc, 0);
    aP_.assign(nc, 1);
    at_.assign(nc, 0);
    Fn_.assign(nf, 0);
    Fnn_.assign(nf, 0);
    nOld_ = 0;
    dt_ = P_.dt;

    // ---- phân số thể tích nước ban đầu: diện tích phần ô nằm dưới mực nước (và trước cửa cống khi vỡ đập)
    alpha_.assign(nc, 1.0);
    if (P_.vof) {
        double X0 = std::numeric_limits<double>::max();
        for (int id : M_.cellNodes) X0 = std::min(X0, M_.ca * M_.px[id] + M_.sa * M_.py[id]);
        for (int c = 0; c < nc; ++c) {
            Poly P;
            for (int k = M_.cellStart[c]; k < M_.cellStart[c + 1]; ++k) P.push_back({M_.px[M_.cellNodes[k]], M_.py[M_.cellNodes[k]]});
            P = clipHalf(P, 0, 1, P_.initLevel);
            if (P_.initGate > 0) P = clipHalf(P, M_.ca, M_.sa, X0 + P_.initGate);
            alpha_[c] = P.size() >= 3 ? std::clamp(polyArea(P) / M_.vol[c], 0.0, 1.0) : 0.0;
        }
    }
    alphaF_.assign(nf, 1.0);
    for (int f = 0; f < M_.nInternal; ++f) alphaF_[f] = M_.w[f] * alpha_[M_.owner[f]] + (1 - M_.w[f]) * alpha_[M_.neigh[f]];
    k_.assign(nc, 0);
    eps_.assign(nc, 0);
    mut_.assign(nc, 0);
    muWall_.assign(nf, P_.mu);
    F_.assign(nf, 0); // updateProperties dùng F_ cho α trên biên (chiều dòng)
    updateProperties();

    std::vector<double> ub, vb;
    boundaryVelocity(ub, vb);
    double netIn = 0, inA = 0, inU = 0, inV = 0;
    for (int f = M_.nInternal; f < nf; ++f) {
        const FluidBcKind k = P_.bc[M_.side[f]].kind;
        if (isDirichlet(k)) F_[f] = ub[f] * M_.sx[f] + vb[f] * M_.sy[f];
        netIn -= F_[f];
        if (k == FluidBcKind::Inlet) {
            const double L = std::hypot(M_.sx[f], M_.sy[f]);
            if (ub[f] != 0 || vb[f] != 0) inA += L;
            inU += ub[f] * L;
            inV += vb[f] * L;
        }
    }
    // giá trị ban đầu: vận tốc trung bình cửa vào trên toàn miền (VOF: chỉ trong nước, trừ bài toán vỡ đập). Bắt đầu từ
    // trường bằng 0 thì ở số Reynolds lớn đường chéo phương trình động lượng gần bằng 0 và hiệu chỉnh áp suất vòng lặp
    // đầu tiên bùng nổ.
    if (inA > 0) {
        const double u0 = inU / inA, v0 = inV / inA;
        for (int c = 0; c < nc; ++c) {
            const bool fill = !P_.vof || (alpha_[c] > 0.5 && P_.initGate <= 0);
            u_[c] = fill ? u0 : 0;
            v_[c] = fill ? v0 : 0;
        }
        for (int f = 0; f < M_.nInternal; ++f) {
            const int P = M_.owner[f], N = M_.neigh[f];
            const double w = M_.w[f];
            F_[f] = (w * u_[P] + (1 - w) * u_[N]) * M_.sx[f] + (w * v_[P] + (1 - w) * v_[N]) * M_.sy[f];
        }
        for (int f = M_.nInternal; f < nf; ++f) {
            const FluidBcKind k = P_.bc[M_.side[f]].kind;
            if (k == FluidBcKind::Outlet || k == FluidBcKind::Opening)
                F_[f] = u_[M_.owner[f]] * M_.sx[f] + v_[M_.owner[f]] * M_.sy[f];
        }
    }
    un_ = unn_ = u_;
    vn_ = vnn_ = v_;
    if (!anyOutflow && std::abs(netIn) > 1e-9 * Uref_ * Lref_) {
        err_ = "Có dòng vào nhưng không có cửa ra: chất lỏng không nén cần ít nhất một biên \"Cửa ra\".";
        return false;
    }
    // ---- áp suất ban đầu (VOF): thủy tĩnh theo cột nước – không khí tại vị trí ô
    if (P_.vof) {
        const double g = P_.gravity, yTop = M_.y1;
        double X0 = std::numeric_limits<double>::max();
        for (int id : M_.cellNodes) X0 = std::min(X0, M_.ca * M_.px[id] + M_.sa * M_.py[id]);
        for (int c = 0; c < nc; ++c) {
            const double y = M_.cy[c];
            // mực nước tại cột này: nước ban đầu chỉ trong đoạn trước cửa cống (vỡ đập)
            const bool wet = P_.initGate <= 0 || M_.ca * M_.cx[c] + M_.sa * M_.cy[c] <= X0 + P_.initGate;
            const double level = wet ? P_.initLevel : -std::numeric_limits<double>::max();
            double pt = P_.rho2 * g * (yTop - std::max(y, std::min(level, yTop)));
            if (y < level) pt += P_.rho * g * (std::min(level, yTop) - y);
            p_[c] = pt + rhoC_[c] * g * y; // p_rgh = p − ρ·g·x = p + ρ·g·y
        }
    }
    initTurbulence();
    iter_ = 0;
    time_ = 0;
    converged_ = diverged_ = false;
    hist_.clear();
    return true;
}

bool FluidSolver::fixedPressure(int f) const {
    const FluidBcKind k = P_.bc[M_.side[f]].kind;
    return k == FluidBcKind::Opening || (k == FluidBcKind::Outlet && !P_.vof);
}

// áp suất (VOF: p_rgh) cho trước trên biên: áp suất toàn phần bằng 0
double FluidSolver::fixedPressureValue(int f) const {
    return P_.vof ? rhoC_[M_.owner[f]] * P_.gravity * M_.fy[f] : 0.0;
}

double FluidSolver::inletProfile(int f) const {
    const int s = M_.side[f];
    const FluidBc& b = P_.bc[s];
    if (!b.parabolic) return 1.0;
    const double L = sideMax_[s] - sideMin_[s];
    if (L <= 0) return 1.0;
    // giá trị trung bình của 4s(1 − s) trên mặt: lưu lượng vào đúng bằng 2/3 · Umax · L
    const bool vert = s == (int)FluidSide::Left || s == (int)FluidSide::Right;
    const double half = 0.5 * std::hypot(M_.sx[f], M_.sy[f]);
    const double c = vert ? -M_.sa * M_.fx[f] + M_.ca * M_.fy[f] : M_.ca * M_.fx[f] + M_.sa * M_.fy[f];
    const double s0 = (c - half - sideMin_[s]) / L, s1 = (c + half - sideMin_[s]) / L;
    auto F = [](double t) { return 2 * t * t - 4.0 / 3 * t * t * t; };
    return s1 > s0 ? (F(s1) - F(s0)) / (s1 - s0) : 4 * s0 * (1 - s0);
}

void FluidSolver::boundaryVelocity(std::vector<double>& ub, std::vector<double>& vb) const {
    ub.assign(M_.nFaces, 0);
    vb.assign(M_.nFaces, 0);
    for (int f = M_.nInternal; f < M_.nFaces; ++f) {
        const FluidBc& b = P_.bc[M_.side[f]];
        const int c = M_.owner[f];
        switch (b.kind) {
        case FluidBcKind::Wall: break;
        case FluidBcKind::MovingWall:
            ub[f] = b.u;
            vb[f] = b.v;
            break;
        case FluidBcKind::Inlet: {
            // VOF: nước vào dưới mực nước cửa vào, phần cao hơn là thành
            double k = P_.vof ? 1.0 : inletProfile(f);
            if (P_.vof) {
                const double ya = M_.py[M_.fa[f]], yb = M_.py[M_.fb[f]];
                const double lo = std::min(ya, yb), hi = std::max(ya, yb);
                k = hi - lo > 1e-12 ? std::clamp((P_.inletLevel - lo) / (hi - lo), 0.0, 1.0) : (M_.fy[f] <= P_.inletLevel ? 1.0 : 0.0);
            }
            ub[f] = b.u * k;
            vb[f] = b.v * k;
            break;
        }
        case FluidBcKind::Outlet:
        case FluidBcKind::Opening:
            ub[f] = u_[c];
            vb[f] = v_[c];
            break;
        case FluidBcKind::Symmetry: {
            const double L = std::hypot(M_.sx[f], M_.sy[f]), nx = M_.sx[f] / L, ny = M_.sy[f] / L;
            const double un = u_[c] * nx + v_[c] * ny;
            ub[f] = u_[c] - un * nx;
            vb[f] = v_[c] - un * ny;
            break;
        }
        }
    }
}

void FluidSolver::boundaryAlpha(std::vector<double>& ab) const {
    ab.assign(M_.nFaces, 1.0);
    if (!P_.vof) return;
    for (int f = M_.nInternal; f < M_.nFaces; ++f) {
        const int c = M_.owner[f];
        const FluidBcKind k = P_.bc[M_.side[f]].kind;
        ab[f] = alpha_[c];
        if (k == FluidBcKind::Inlet) {
            const double ya = M_.py[M_.fa[f]], yb = M_.py[M_.fb[f]];
            const double lo = std::min(ya, yb), hi = std::max(ya, yb);
            ab[f] = hi - lo > 1e-12 ? std::clamp((P_.inletLevel - lo) / (hi - lo), 0.0, 1.0) : (M_.fy[f] <= P_.inletLevel ? 1.0 : 0.0);
        } else if (k == FluidBcKind::Opening && F_[f] < 0)
            ab[f] = 0; // không khí đi vào qua biên thông khí
    }
}

void FluidSolver::updateProperties() {
    const int nc = M_.nCells, nf = M_.nFaces;
    rhoC_.resize(nc);
    muC_.resize(nc);
    rhoF_.resize(nf);
    if (!P_.vof) {
        std::fill(rhoC_.begin(), rhoC_.end(), P_.rho);
        std::fill(muC_.begin(), muC_.end(), P_.mu);
        std::fill(rhoF_.begin(), rhoF_.end(), P_.rho);
    } else {
        for (int c = 0; c < nc; ++c) {
            const double a = std::clamp(alpha_[c], 0.0, 1.0);
            rhoC_[c] = a * P_.rho + (1 - a) * P_.rho2;
            muC_[c] = a * P_.mu + (1 - a) * P_.mu2;
        }
        std::vector<double> ab;
        boundaryAlpha(ab);
        for (int f = 0; f < nf; ++f) {
            const double a = std::clamp(f < M_.nInternal ? alphaF_[f] : ab[f], 0.0, 1.0);
            rhoF_[f] = a * P_.rho + (1 - a) * P_.rho2;
        }
    }
    if (P_.turbulence == TurbModel::Laminar)
        for (int f = M_.nInternal; f < nf; ++f) muWall_[f] = muC_[M_.owner[f]];
}

void FluidSolver::boundaryPressure(const std::vector<double>& p, const std::vector<double>& gx, const std::vector<double>& gy,
                                   std::vector<double>& pb) const {
    pb.assign(M_.nFaces, 0);
    for (int f = M_.nInternal; f < M_.nFaces; ++f) {
        if (fixedPressure(f)) {
            pb[f] = fixedPressureValue(f);
            continue;
        }
        const int c = M_.owner[f];
        pb[f] = p[c] + gx[c] * (M_.fx[f] - M_.cx[c]) + gy[c] * (M_.fy[f] - M_.cy[c]);
    }
}

// Gradient bình phương tối thiểu có trọng số 1/|d|²: chính xác với trường tuyến tính trên mọi lưới (kể cả tam giác lệch).
// Mặt biên tham gia như một ô lân cận đặt tại trung điểm mặt với giá trị faceVal.
void FluidSolver::gradient(const std::vector<double>& phi, const std::vector<double>& faceVal, std::vector<double>& gx,
                           std::vector<double>& gy) const {
    const int nc = M_.nCells;
    gx.assign(nc, 0);
    gy.assign(nc, 0);
    for (int f = 0; f < M_.nFaces; ++f) {
        const int P = M_.owner[f], N = M_.neigh[f];
        const double dx = (N >= 0 ? M_.cx[N] : M_.fx[f]) - M_.cx[P], dy = (N >= 0 ? M_.cy[N] : M_.fy[f]) - M_.cy[P];
        const double w = 1.0 / (dx * dx + dy * dy);
        const double dphi = (N >= 0 ? phi[N] : faceVal[f]) - phi[P];
        gx[P] += w * dx * dphi;
        gy[P] += w * dy * dphi;
        if (N >= 0) {
            gx[N] += w * dx * dphi; // (−d)·(−Δφ)
            gy[N] += w * dy * dphi;
        }
    }
    for (int c = 0; c < nc; ++c) {
        const double* G = &lsq_[4 * c];
        const double a = gx[c], b = gy[c];
        gx[c] = G[0] * a + G[1] * b;
        gy[c] = G[2] * a + G[3] * b;
    }
}

void FluidSolver::gradientGG(const std::vector<double>& phi, const std::vector<double>& faceVal, std::vector<double>& gx,
                             std::vector<double>& gy) const {
    std::vector<double> lx, ly;
    gradient(phi, faceVal, lx, ly);
    gx.assign(M_.nCells, 0);
    gy.assign(M_.nCells, 0);
    for (int f = 0; f < M_.nInternal; ++f) {
        const int P = M_.owner[f], N = M_.neigh[f];
        const double w = M_.w[f];
        const double pf = w * phi[P] + (1 - w) * phi[N] + (w * lx[P] + (1 - w) * lx[N]) * skx_[f] +
                          (w * ly[P] + (1 - w) * ly[N]) * sky_[f];
        gx[P] += pf * M_.sx[f];
        gy[P] += pf * M_.sy[f];
        gx[N] -= pf * M_.sx[f];
        gy[N] -= pf * M_.sy[f];
    }
    for (int f = M_.nInternal; f < M_.nFaces; ++f) {
        gx[M_.owner[f]] += faceVal[f] * M_.sx[f];
        gy[M_.owner[f]] += faceVal[f] * M_.sy[f];
    }
    for (int c = 0; c < M_.nCells; ++c) {
        gx[c] /= M_.vol[c];
        gy[c] /= M_.vol[c];
    }
}

// u_f·S trên các mặt trong (lưu lượng thể tích), u_f nội suy có hiệu chỉnh độ lệch về trung điểm mặt
void FluidSolver::faceFlux(std::vector<double>& out) const {
    std::vector<double> ub, vb, gux, guy, gvx, gvy;
    boundaryVelocity(ub, vb);
    gradient(u_, ub, gux, guy);
    gradient(v_, vb, gvx, gvy);
    out.assign(M_.nInternal, 0);
    for (int f = 0; f < M_.nInternal; ++f) {
        const int P = M_.owner[f], N = M_.neigh[f];
        const double w = M_.w[f], rx = skx_[f], ry = sky_[f];
        const double uf = w * u_[P] + (1 - w) * u_[N] + (w * gux[P] + (1 - w) * gux[N]) * rx + (w * guy[P] + (1 - w) * guy[N]) * ry;
        const double vf = w * v_[P] + (1 - w) * v_[N] + (w * gvx[P] + (1 - w) * gvx[N]) * rx + (w * gvy[P] + (1 - w) * gvy[N]) * ry;
        out[f] = uf * M_.sx[f] + vf * M_.sy[f];
    }
}

void FluidSolver::gaussSeidel(std::vector<double>& x, const std::vector<double>& b, int sweeps) const {
    auto relax = [&](int c) {
        double s = b[c];
        for (int k = M_.cfStart[c]; k < M_.cfStart[c + 1]; ++k) {
            const int f = M_.cfFace[k];
            if (f >= M_.nInternal) continue;
            if (M_.owner[f] == c) s -= upper_[f] * x[M_.neigh[f]];
            else s -= lower_[f] * x[M_.owner[f]];
        }
        x[c] = s / diag_[c];
    };
    for (int s = 0; s < sweeps; ++s) {
        for (int c = 0; c < M_.nCells; ++c) relax(c);
        for (int c = M_.nCells - 1; c >= 0; --c) relax(c);
    }
}

// Gradient liên hợp có tiền điều kiện Cholesky không đầy đủ dạng đường chéo (DIC) cho ma trận đối xứng diag_/upper_.
int FluidSolver::pcg(std::vector<double>& x, const std::vector<double>& b, double relTol, int maxIt) const {
    const int n = M_.nCells, ni = M_.nInternal;
    std::vector<double> rD(diag_);
    for (int f = 0; f < ni; ++f) rD[M_.neigh[f]] -= upper_[f] * upper_[f] / rD[M_.owner[f]];
    for (double& d : rD) d = 1.0 / d;
    auto mult = [&](const std::vector<double>& a, std::vector<double>& y) {
        for (int c = 0; c < n; ++c) y[c] = diag_[c] * a[c];
        for (int f = 0; f < ni; ++f) {
            y[M_.owner[f]] += upper_[f] * a[M_.neigh[f]];
            y[M_.neigh[f]] += upper_[f] * a[M_.owner[f]];
        }
    };
    auto precond = [&](const std::vector<double>& r, std::vector<double>& z) {
        for (int c = 0; c < n; ++c) z[c] = rD[c] * r[c];
        for (int f = 0; f < ni; ++f) z[M_.neigh[f]] -= rD[M_.neigh[f]] * upper_[f] * z[M_.owner[f]];
        for (int f = ni - 1; f >= 0; --f) z[M_.owner[f]] -= rD[M_.owner[f]] * upper_[f] * z[M_.neigh[f]];
    };
    std::vector<double> r(n), z(n), pv(n), q(n);
    mult(x, q);
    double bn = 0, rn = 0;
    for (int c = 0; c < n; ++c) {
        r[c] = b[c] - q[c];
        bn += std::abs(b[c]);
        rn += std::abs(r[c]);
    }
    if (bn < 1e-300 || rn <= relTol * bn) return 0;
    precond(r, z);
    pv = z;
    double rz = std::inner_product(r.begin(), r.end(), z.begin(), 0.0);
    int it = 0;
    while (it < maxIt) {
        ++it;
        mult(pv, q);
        const double pq = std::inner_product(pv.begin(), pv.end(), q.begin(), 0.0);
        if (std::abs(pq) < 1e-300) break;
        const double alpha = rz / pq;
        rn = 0;
        for (int c = 0; c < n; ++c) {
            x[c] += alpha * pv[c];
            r[c] -= alpha * q[c];
            rn += std::abs(r[c]);
        }
        if (rn <= relTol * bn) break;
        precond(r, z);
        const double rz2 = std::inner_product(r.begin(), r.end(), z.begin(), 0.0);
        const double beta = rz2 / rz;
        rz = rz2;
        for (int c = 0; c < n; ++c) pv[c] = z[c] + beta * pv[c];
    }
    return it;
}

void FluidSolver::pressureSource() {
    std::vector<double> pb;
    // áp suất trên biên ngoại suy bằng gradient của vòng lặp trước
    boundaryPressure(p_, pgx_, pgy_, pb);
    gradientGG(p_, pb, pgx_, pgy_);
    gpx_ = pgx_;
    gpy_ = pgy_;
    if (P_.vof) {
        // ∇p_rgh + (g·x)·∇ρ dựng lại từ thành phần pháp tuyến trên mặt — cùng dạng với số hạng áp suất của Rhie–Chow:
        // ở trạng thái thủy tĩnh mọi mặt cân bằng nên véc tơ tại tâm ô bằng 0 đúng (không sinh dòng giả trong không khí)
        std::vector<double> sx(M_.nCells, 0), sy(M_.nCells, 0);
        for (int f = 0; f < M_.nFaces; ++f) {
            const int P = M_.owner[f], N = M_.neigh[f];
            double gn = 0; // thành phần pháp tuyến × |S|
            const double gh = -P_.gravity * M_.fy[f];
            if (N >= 0) {
                const double dx = M_.cx[N] - M_.cx[P], dy = M_.cy[N] - M_.cy[P];
                const double SS = M_.sx[f] * M_.sx[f] + M_.sy[f] * M_.sy[f];
                gn = SS / (dx * M_.sx[f] + dy * M_.sy[f]) * ((p_[N] - p_[P]) + gh * (rhoC_[N] - rhoC_[P]));
            } else if (fixedPressure(f)) {
                const double dx = M_.fx[f] - M_.cx[P], dy = M_.fy[f] - M_.cy[P];
                const double SS = M_.sx[f] * M_.sx[f] + M_.sy[f] * M_.sy[f];
                gn = SS / (dx * M_.sx[f] + dy * M_.sy[f]) * (fixedPressureValue(f) - p_[P]);
            } // thành, cửa vào: građien pháp tuyến bằng 0
            const double L = std::hypot(M_.sx[f], M_.sy[f]);
            sx[P] += M_.sx[f] / L * gn;
            sy[P] += M_.sy[f] / L * gn;
            if (N >= 0) {
                sx[N] += M_.sx[f] / L * gn;
                sy[N] += M_.sy[f] / L * gn;
            }
        }
        for (int c = 0; c < M_.nCells; ++c) {
            const double* G = &recon_[4 * c];
            gpx_[c] = G[0] * sx[c] + G[1] * sy[c];
            gpy_[c] = G[2] * sx[c] + G[3] * sy[c];
        }
    }
}

void FluidSolver::momentum(double& resU, double& resV) {
    const int nc = M_.nCells, ni = M_.nInternal, nf = M_.nFaces;
    const bool turb = P_.turbulence != TurbModel::Laminar;
    std::vector<double> ub, vb, gux, guy, gvx, gvy;
    boundaryVelocity(ub, vb);
    gradient(u_, ub, gux, guy);
    gradient(v_, vb, gvx, gvy);
    pressureSource();
    auto muEff = [&](int c) { return muC_[c] + mut_[c]; };

    diag_.assign(nc, 0);
    upper_.assign(ni, 0);
    lower_.assign(ni, 0);
    std::vector<double> bu(nc, 0), bv(nc, 0);
    const bool lud = P_.scheme == ConvScheme::LinearUpwind;
    // Barth–Jespersen: giá trị tái dựng trên mặt không vượt ra ngoài [min, max] của ô và các lân cận
    // (ngược dòng tuyến tính không sinh dao động ở số Reynolds lớn)
    auto limit = [&](const std::vector<double>& phi, const std::vector<double>& bval, std::vector<double>& gx,
                     std::vector<double>& gy) {
        std::vector<double> lo(phi), hi(phi), psi(nc, 1.0);
        for (int f = 0; f < nf; ++f) {
            const int P = M_.owner[f], N = M_.neigh[f];
            const double q = N >= 0 ? phi[N] : bval[f];
            lo[P] = std::min(lo[P], q), hi[P] = std::max(hi[P], q);
            if (N >= 0) lo[N] = std::min(lo[N], phi[P]), hi[N] = std::max(hi[N], phi[P]);
        }
        auto clip = [&](int c, double d) {
            if (d > 1e-12 * (std::abs(phi[c]) + 1e-30)) psi[c] = std::min(psi[c], (hi[c] - phi[c]) / d);
            else if (d < -1e-12 * (std::abs(phi[c]) + 1e-30)) psi[c] = std::min(psi[c], (lo[c] - phi[c]) / d);
        };
        for (int f = 0; f < nf; ++f) {
            const int P = M_.owner[f], N = M_.neigh[f];
            clip(P, gx[P] * (M_.fx[f] - M_.cx[P]) + gy[P] * (M_.fy[f] - M_.cy[P]));
            if (N >= 0) clip(N, gx[N] * (M_.fx[f] - M_.cx[N]) + gy[N] * (M_.fy[f] - M_.cy[N]));
        }
        for (int c = 0; c < nc; ++c) {
            gx[c] *= std::max(psi[c], 0.0);
            gy[c] *= std::max(psi[c], 0.0);
        }
    };
    std::vector<double> lux(gux), luy(guy), lvx(gvx), lvy(gvy); // gradient đã giới hạn, chỉ dùng cho đối lưu
    if (lud) {
        limit(u_, ub, lux, luy);
        limit(v_, vb, lvx, lvy);
    }
    for (int f = 0; f < ni; ++f) {
        const int P = M_.owner[f], N = M_.neigh[f];
        const double F = rhoF_[f] * F_[f]; // lưu lượng khối
        // đối lưu ngược dòng (ẩn), dạng bị chặn Σ F·(φf − φP): trừ sai số liên tục khỏi đường chéo để đường chéo
        // luôn bằng tổng hệ số lân cận (ô chưa có dòng ra ở các vòng lặp đầu không làm hệ suy biến); khi đã thỏa
        // phương trình liên tục, hai dạng trùng nhau
        diag_[P] += std::max(-F, 0.0);
        upper_[f] += std::min(F, 0.0);
        diag_[N] += std::max(F, 0.0);
        lower_[f] += std::min(-F, 0.0);
        // khuếch tán: phần trực giao ẩn
        const double w = M_.w[f];
        const double mu = w * muEff(P) + (1 - w) * muEff(N);
        const double dx = M_.cx[N] - M_.cx[P], dy = M_.cy[N] - M_.cy[P];
        const double SS = M_.sx[f] * M_.sx[f] + M_.sy[f] * M_.sy[f], dS = dx * M_.sx[f] + dy * M_.sy[f];
        const double Sd = SS / dS, Df = mu * Sd;
        diag_[P] += Df;
        diag_[N] += Df;
        upper_[f] -= Df;
        lower_[f] -= Df;
        // phần không trực giao (hiệu chỉnh trễ)
        const double kx = M_.sx[f] - Sd * dx, ky = M_.sy[f] - Sd * dy;
        const double fux = w * gux[P] + (1 - w) * gux[N], fuy = w * guy[P] + (1 - w) * guy[N];
        const double fvx = w * gvx[P] + (1 - w) * gvx[N], fvy = w * gvy[P] + (1 - w) * gvy[N];
        double cu = mu * (fux * kx + fuy * ky), cv = mu * (fvx * kx + fvy * ky);
        // độ nhớt thay đổi (chảy rối): thêm ∇·(μ_eff·∇uᵀ)
        if (turb) {
            cu += mu * (fux * M_.sx[f] + fvx * M_.sy[f]);
            cv += mu * (fuy * M_.sx[f] + fvy * M_.sy[f]);
        }
        bu[P] += cu;
        bu[N] -= cu;
        bv[P] += cv;
        bv[N] -= cv;
        // ngược dòng tuyến tính bậc hai (hiệu chỉnh trễ)
        if (lud) {
            const int U = F >= 0 ? P : N;
            const double rx = M_.fx[f] - M_.cx[U], ry = M_.fy[f] - M_.cy[U];
            const double du = F * (lux[U] * rx + luy[U] * ry), dv = F * (lvx[U] * rx + lvy[U] * ry);
            bu[P] -= du;
            bu[N] += du;
            bv[P] -= dv;
            bv[N] += dv;
        }
    }
    for (int f = ni; f < nf; ++f) {
        const int P = M_.owner[f];
        const FluidBcKind k = P_.bc[M_.side[f]].kind;
        const double F = rhoF_[f] * F_[f];
        const double dx = M_.fx[f] - M_.cx[P], dy = M_.fy[f] - M_.cy[P];
        const double SS = M_.sx[f] * M_.sx[f] + M_.sy[f] * M_.sy[f], dS = dx * M_.sx[f] + dy * M_.sy[f];
        if (k == FluidBcKind::Outlet || k == FluidBcKind::Opening) {
            // dòng ra: građien bằng 0, φb = φP (không đóng góp ở dạng bị chặn); dòng chảy ngược vào: φb = 0
            if (F < 0) diag_[P] -= F;
            continue;
        }
        // thành: độ nhớt hiệu dụng của hàm tường; cửa vào, đối xứng: độ nhớt của ô
        const double Db = (isWall(k) ? muWall_[f] : muEff(P)) * SS / dS;
        // vận tốc cho trước trên mặt (thành, cửa vào) hoặc thành phần tiếp tuyến của vận tốc ô (đối xứng): F·(φb − φP)
        bu[P] -= F * ub[f];
        bv[P] -= F * vb[f];
        diag_[P] -= F;
        diag_[P] += Db;
        bu[P] += Db * ub[f];
        bv[P] += Db * vb[f];
    }
    // áp suất, thời gian, giảm dư
    const bool bdf2 = P_.transient && P_.timeScheme == TimeScheme::BDF2 && nOld_ >= 2;
    for (int c = 0; c < nc; ++c) {
        bu[c] -= gpx_[c] * M_.vol[c];
        bv[c] -= gpy_[c] * M_.vol[c];
        at_[c] = 0;
        if (P_.transient) {
            const double r = rhoC_[c] * M_.vol[c] / dt_;
            at_[c] = bdf2 ? 1.5 * r : r;
            diag_[c] += at_[c];
            bu[c] += bdf2 ? r * (2 * un_[c] - 0.5 * unn_[c]) : r * un_[c];
            bv[c] += bdf2 ? r * (2 * vn_[c] - 0.5 * vnn_[c]) : r * vn_[c];
        }
        const double a = diag_[c];
        diag_[c] = a / P_.alphaU;
        bu[c] += (1 - P_.alphaU) * diag_[c] * u_[c];
        bv[c] += (1 - P_.alphaU) * diag_[c] * v_[c];
        aP_[c] = diag_[c];
    }
    // sai số dư chuẩn hóa: Σ|b − A·φ| / Σ|aP·φ|
    auto residual = [&](const std::vector<double>& x, const std::vector<double>& b) {
        std::vector<double> r(b);
        double sa = 0;
        for (int c = 0; c < nc; ++c) {
            r[c] -= diag_[c] * x[c];
            sa += std::abs(diag_[c] * x[c]);
        }
        for (int f = 0; f < ni; ++f) {
            r[M_.owner[f]] -= upper_[f] * x[M_.neigh[f]];
            r[M_.neigh[f]] -= lower_[f] * x[M_.owner[f]];
        }
        double sr = 0, sb = 0;
        for (double v : r) sr += std::abs(v);
        for (double v : b) sb += std::abs(v);
        return sr / std::max(sa + sb, 1e-300);
    };
    resU = residual(u_, bu);
    resV = residual(v_, bv);
    gaussSeidel(u_, bu, 2);
    gaussSeidel(v_, bv, 2);
}

double FluidSolver::pressureCorrection(const std::vector<double>& majumdar) {
    const int nc = M_.nCells, ni = M_.nInternal, nf = M_.nFaces;
    std::vector<double> D(nc), imb(nc, 0);
    for (int c = 0; c < nc; ++c) D[c] = M_.vol[c] / aP_[c];
    std::vector<double> ub, vb;
    boundaryVelocity(ub, vb);

    // ---- lưu lượng thể tích qua mặt theo Rhie–Chow; VOF: số hạng áp suất mặt gồm cả g·x·Δρ (cân bằng lực tại mặt)
    const bool bdf2 = P_.transient && P_.timeScheme == TimeScheme::BDF2 && nOld_ >= 2;
    std::vector<double> Df(ni), Sdf(ni), Fu;
    faceFlux(Fu);
    for (int f = 0; f < ni; ++f) {
        const int P = M_.owner[f], N = M_.neigh[f];
        const double w = M_.w[f];
        const double dx = M_.cx[N] - M_.cx[P], dy = M_.cy[N] - M_.cy[P];
        const double SS = M_.sx[f] * M_.sx[f] + M_.sy[f] * M_.sy[f], Sd = SS / (dx * M_.sx[f] + dy * M_.sy[f]);
        const double gx = w * gpx_[P] + (1 - w) * gpx_[N], gy = w * gpy_[P] + (1 - w) * gpy_[N];
        Df[f] = w * D[P] + (1 - w) * D[N];
        Sdf[f] = Sd;
        double dp = p_[N] - p_[P];
        if (P_.vof) dp += -P_.gravity * M_.fy[f] * (rhoC_[N] - rhoC_[P]);
        double F = Fu[f] - Df[f] * Sd * (dp - (gx * dx + gy * dy));
        // hiệu chỉnh Majumdar: nghiệm hội tụ không phụ thuộc hệ số giảm dư
        F += majumdar[f];
        // hiệu chỉnh Choi cho số hạng thời gian: Fn_ lưu (F − ū·S) của các bước trước
        if (P_.transient) {
            const double tau = w * at_[P] / aP_[P] + (1 - w) * at_[N] / aP_[N];
            F += bdf2 ? tau / 1.5 * (2 * Fn_[f] - 0.5 * Fnn_[f]) : tau * Fn_[f];
        }
        F_[f] = F;
    }
    for (int f = 0; f < ni; ++f) {
        imb[M_.owner[f]] += F_[f];
        imb[M_.neigh[f]] -= F_[f];
    }
    std::vector<double> SdOut(nf, 0);
    for (int f = ni; f < nf; ++f) {
        const int P = M_.owner[f];
        const FluidBcKind k = P_.bc[M_.side[f]].kind;
        if (fixedPressure(f)) {
            const double dx = M_.fx[f] - M_.cx[P], dy = M_.fy[f] - M_.cy[P];
            const double SS = M_.sx[f] * M_.sx[f] + M_.sy[f] * M_.sy[f];
            SdOut[f] = SS / (dx * M_.sx[f] + dy * M_.sy[f]);
            F_[f] = u_[P] * M_.sx[f] + v_[P] * M_.sy[f] -
                    D[P] * SdOut[f] * ((fixedPressureValue(f) - p_[P]) - (gpx_[P] * dx + gpy_[P] * dy));
        } else if (k == FluidBcKind::Outlet) {
            // VOF: cửa ra građien áp suất bằng 0 — lưu lượng theo vận tốc ô, không cho chảy ngược vào
            F_[f] = std::max(u_[P] * M_.sx[f] + v_[P] * M_.sy[f], 0.0);
        } else F_[f] = ub[f] * M_.sx[f] + vb[f] * M_.sy[f];
        imb[P] += F_[f];
    }
    double sumF = 0, sumImb = 0;
    for (int f = 0; f < nf; ++f) sumF += std::abs(F_[f]);
    for (double v : imb) sumImb += std::abs(v);
    contNorm_ = std::max(sumF, Uref_ * Lref_);
    const double resC = sumImb / contNorm_;

    // ---- phương trình hiệu chỉnh áp suất
    diag_.assign(nc, 0);
    upper_.assign(ni, 0);
    for (int f = 0; f < ni; ++f) {
        const double a = Df[f] * Sdf[f];
        diag_[M_.owner[f]] += a;
        diag_[M_.neigh[f]] += a;
        upper_[f] = -a;
    }
    for (int f = ni; f < nf; ++f)
        if (SdOut[f] > 0) diag_[M_.owner[f]] += D[M_.owner[f]] * SdOut[f];
    if (!hasOutlet_) diag_[0] *= 2; // miền kín: cố định mức áp suất (hệ vẫn tương thích, p'₀ = 0)
    std::vector<double> rhs(nc), pc(nc, 0);
    for (int c = 0; c < nc; ++c) rhs[c] = -imb[c];
    pcg(pc, rhs, P_.transient ? 1e-3 : 1e-2, 2000);

    // ---- hiệu chỉnh áp suất, vận tốc, lưu lượng
    std::vector<double> pcb(nf, 0), gx, gy;
    for (int f = ni; f < nf; ++f)
        if (SdOut[f] == 0) pcb[f] = pc[M_.owner[f]];
    gradient(pc, pcb, gx, gy);
    for (int c = 0; c < nc; ++c) {
        p_[c] += P_.alphaP * pc[c];
        u_[c] -= D[c] * gx[c];
        v_[c] -= D[c] * gy[c];
    }
    for (int f = 0; f < ni; ++f) F_[f] -= Df[f] * Sdf[f] * (pc[M_.neigh[f]] - pc[M_.owner[f]]);
    for (int f = ni; f < nf; ++f)
        if (SdOut[f] > 0) F_[f] += D[M_.owner[f]] * SdOut[f] * pc[M_.owner[f]];
    if (!hasOutlet_) {
        double mean = 0, V = 0;
        for (int c = 0; c < nc; ++c) {
            mean += p_[c] * M_.vol[c];
            V += M_.vol[c];
        }
        mean /= V;
        for (double& v : p_) v -= mean;
    }
    return resC;
}

// ---------------------------------------------------------------- chảy rối k-ε
void FluidSolver::initTurbulence() {
    const int nc = M_.nCells;
    if (P_.turbulence == TurbModel::Laminar) {
        std::fill(mut_.begin(), mut_.end(), 0.0);
        return;
    }
    // k, ε ở cửa vào theo cường độ rối và thang chiều dài (mặc định 7% bề rộng cửa vào)
    double Win = 0;
    for (int s = 0; s < (int)FluidSide::Count; ++s)
        if (M_.sideFaces[s] && P_.bc[s].kind == FluidBcKind::Inlet) Win = std::max(Win, sideMax_[s] - sideMin_[s]);
    const double ell = P_.turbLength > 0 ? P_.turbLength : 0.07 * (Win > 0 ? Win : Lref_);
    kIn_ = 1.5 * std::pow(P_.turbIntensity * Uref_, 2);
    kIn_ = std::max(kIn_, 1e-10 * Uref_ * Uref_);
    epsIn_ = std::pow(kCmu, 0.75) * std::pow(kIn_, 1.5) / ell;
    k_.assign(nc, kIn_);
    eps_.assign(nc, epsIn_);
    kn_ = k_;
    en_ = eps_;
    for (int c = 0; c < nc; ++c) mut_[c] = rhoC_[c] * kCmu * k_[c] * k_[c] / eps_[c];
    for (int f = M_.nInternal; f < M_.nFaces; ++f) muWall_[f] = muC_[M_.owner[f]];
}

void FluidSolver::solveTurbulence() {
    const int nc = M_.nCells, ni = M_.nInternal, nf = M_.nFaces;
    const double kMin = 1e-10 * Uref_ * Uref_, eMin = 1e-12 * Uref_ * Uref_ * Uref_ / Lref_;
    const double cmu25 = std::pow(kCmu, 0.25), cmu75 = std::pow(kCmu, 0.75);
    std::vector<double> ub, vb, gux, guy, gvx, gvy;
    boundaryVelocity(ub, vb);
    gradient(u_, ub, gux, guy);
    gradient(v_, vb, gvx, gvy);
    // sinh rối Kato–Launder: Pk = μt·S·Ω (không làm k tăng vọt tại điểm dừng như μt·S²)
    std::vector<double> Pk(nc);
    for (int c = 0; c < nc; ++c) {
        const double S = std::sqrt(2 * (gux[c] * gux[c] + gvy[c] * gvy[c]) + std::pow(guy[c] + gvx[c], 2));
        const double W = std::abs(gvx[c] - guy[c]);
        Pk[c] = mut_[c] * S * W;
    }
    // ---- hàm tường "co giãn" (y* ≥ 11,06): ứng suất tiếp trên thành, sinh rối và ε trong ô sát thành
    std::vector<double> wPk(nc, 0), wEps(nc, 0), wLen(nc, 0);
    for (int f = ni; f < nf; ++f) {
        if (!isWall(P_.bc[M_.side[f]].kind)) continue;
        const int c = M_.owner[f];
        const double L = std::hypot(M_.sx[f], M_.sy[f]), nx = M_.sx[f] / L, ny = M_.sy[f] / L;
        const double yP = std::abs((M_.fx[f] - M_.cx[c]) * nx + (M_.fy[f] - M_.cy[c]) * ny);
        const double k = std::max(k_[c], kMin), us = cmu25 * std::sqrt(k);
        const double mu = muC_[c], rho = rhoC_[c];
        const double ys = std::max(rho * us * yP / mu, kYLam);
        const double yEff = ys * mu / (rho * us); // khoảng cách hiệu dụng (thành "dịch ra" khi lưới quá mịn)
        double du = u_[c] - ub[f], dv = v_[c] - vb[f];
        const double dn = du * nx + dv * ny;
        du -= dn * nx;
        dv -= dn * ny;
        const double Ut = std::hypot(du, dv);
        const double muW = mu * ys * kKappa / std::log(kE * ys);
        muWall_[f] = muW;
        const double tauW = muW * Ut / yP;
        wPk[c] += tauW * us / (kKappa * yEff) * L;
        wEps[c] += cmu75 * std::pow(k, 1.5) / (kKappa * yEff) * L;
        wLen[c] += L;
    }
    std::vector<char> fixedE(nc, 0);
    for (int c = 0; c < nc; ++c)
        if (wLen[c] > 0) {
            Pk[c] = wPk[c] / wLen[c];
            eps_[c] = std::max(wEps[c] / wLen[c], eMin);
            fixedE[c] = 1;
        }
    // ---- lắp ráp và giải một phương trình vận chuyển ρ(∂φ/∂t + u·∇φ) = ∇·(Γ∇φ) + nguồn − sink·φ
    auto solve = [&](std::vector<double>& phi, const std::vector<double>& phiOld, double sigma, double phiIn,
                     const std::vector<double>& src, const std::vector<double>& sink, const std::vector<char>* fixed,
                     double phiMin) {
        diag_.assign(nc, 0);
        upper_.assign(ni, 0);
        lower_.assign(ni, 0);
        std::vector<double> b(nc, 0);
        for (int f = 0; f < ni; ++f) {
            const int P = M_.owner[f], N = M_.neigh[f];
            const double F = rhoF_[f] * F_[f];
            diag_[P] += std::max(-F, 0.0);
            upper_[f] += std::min(F, 0.0);
            diag_[N] += std::max(F, 0.0);
            lower_[f] += std::min(-F, 0.0);
            const double w = M_.w[f];
            const double G = w * (muC_[P] + mut_[P] / sigma) + (1 - w) * (muC_[N] + mut_[N] / sigma);
            const double dx = M_.cx[N] - M_.cx[P], dy = M_.cy[N] - M_.cy[P];
            const double SS = M_.sx[f] * M_.sx[f] + M_.sy[f] * M_.sy[f];
            const double Df = G * SS / (dx * M_.sx[f] + dy * M_.sy[f]);
            diag_[P] += Df;
            diag_[N] += Df;
            upper_[f] -= Df;
            lower_[f] -= Df;
        }
        for (int f = ni; f < nf; ++f) {
            const int P = M_.owner[f];
            const FluidBcKind k = P_.bc[M_.side[f]].kind;
            const double F = rhoF_[f] * F_[f];
            if (k == FluidBcKind::Inlet) {
                // giá trị cho trước ở cửa vào
                const double dx = M_.fx[f] - M_.cx[P], dy = M_.fy[f] - M_.cy[P];
                const double SS = M_.sx[f] * M_.sx[f] + M_.sy[f] * M_.sy[f];
                const double Db = (muC_[P] + mut_[P] / sigma) * SS / (dx * M_.sx[f] + dy * M_.sy[f]);
                diag_[P] += Db - F;
                b[P] += Db * phiIn - F * phiIn;
            } else if ((k == FluidBcKind::Outlet || k == FluidBcKind::Opening) && F < 0) {
                diag_[P] -= F; // chảy ngược vào: mang giá trị cửa vào
                b[P] -= F * phiIn;
            }
            // thành, đối xứng, dòng ra: thông lượng bằng 0 / građien bằng 0
        }
        for (int c = 0; c < nc; ++c) {
            const double V = M_.vol[c];
            diag_[c] += sink[c] * V;
            b[c] += src[c] * V;
            if (P_.transient) {
                const double r = rhoC_[c] * V / dt_;
                diag_[c] += r;
                b[c] += r * phiOld[c];
            }
            const double a = diag_[c] / 0.7; // giảm dư 0,7
            b[c] += (a - diag_[c]) * phi[c];
            diag_[c] = a;
        }
        for (int sweep = 0; sweep < 3; ++sweep)
            for (int dir = 0; dir < 2; ++dir)
                for (int i = 0; i < nc; ++i) {
                    const int c = dir == 0 ? i : nc - 1 - i;
                    if (fixed && (*fixed)[c]) continue;
                    double s = b[c];
                    for (int q = M_.cfStart[c]; q < M_.cfStart[c + 1]; ++q) {
                        const int f = M_.cfFace[q];
                        if (f >= ni) continue;
                        if (M_.owner[f] == c) s -= upper_[f] * phi[M_.neigh[f]];
                        else s -= lower_[f] * phi[M_.owner[f]];
                    }
                    phi[c] = std::max(s / diag_[c], phiMin);
                }
    };
    std::vector<double> src(nc), sink(nc);
    // phương trình k: nguồn Pk, tiêu tán ρ·ε = (ρ·ε/k)·k (ẩn)
    for (int c = 0; c < nc; ++c) {
        src[c] = Pk[c];
        sink[c] = rhoC_[c] * eps_[c] / std::max(k_[c], kMin);
    }
    solve(k_, kn_, kSigK, kIn_, src, sink, nullptr, kMin);
    // phương trình ε: nguồn C1·(ε/k)·Pk, tiêu tán C2·ρ·ε²/k; ô sát thành giữ giá trị hàm tường
    for (int c = 0; c < nc; ++c) {
        const double ek = eps_[c] / std::max(k_[c], kMin);
        src[c] = kC1 * ek * Pk[c];
        sink[c] = kC2 * rhoC_[c] * ek;
    }
    solve(eps_, en_, kSigE, epsIn_, src, sink, &fixedE, eMin);
    for (int c = 0; c < nc; ++c) mut_[c] = std::min(rhoC_[c] * kCmu * k_[c] * k_[c] / eps_[c], 1e5 * muC_[c]);
}

// ---------------------------------------------------------------- vận chuyển α (VOF)
double FluidSolver::courantRate() const {
    std::vector<double> out(M_.nCells, 0);
    for (int f = 0; f < M_.nFaces; ++f) {
        const double F = F_[f];
        if (F > 0) out[M_.owner[f]] += F;
        else if (M_.neigh[f] >= 0) out[M_.neigh[f]] -= F;
    }
    double r = 0;
    for (int c = 0; c < M_.nCells; ++c) r = std::max(r, out[c] / M_.vol[c]);
    return r;
}

// HRIC (Muzaferija & Perić): giá trị α trên mặt pha trộn giữa ngược dòng và xuôi dòng (nén mặt phân cách), hiệu chỉnh
// theo số Courant và góc giữa mặt phân cách với mặt ô; tích phân hiện, chia bước con để Courant ≤ 0,25.
void FluidSolver::advectAlpha(double dt) {
    const int ni = M_.nInternal, nf = M_.nFaces;
    const int nSub = std::max(1, (int)std::ceil(courantRate() * dt / 0.25));
    const double h = dt / nSub;
    std::vector<double> sumAf(nf, 0), ab, gx, gy, flux(nf);
    for (int s = 0; s < nSub; ++s) {
        boundaryAlpha(ab);
        gradient(alpha_, ab, gx, gy);
        for (int f = 0; f < ni; ++f) {
            const double F = F_[f];
            const int D = F >= 0 ? M_.owner[f] : M_.neigh[f], A = F >= 0 ? M_.neigh[f] : M_.owner[f];
            const double aD = alpha_[D], aA = alpha_[A];
            const double dx = M_.cx[A] - M_.cx[D], dy = M_.cy[A] - M_.cy[D];
            const double aUU = std::clamp(aA - 2 * (gx[D] * dx + gy[D] * dy), 0.0, 1.0);
            const double den = aA - aUU;
            double af = aD;
            if (std::abs(den) > 1e-9) {
                const double tD = (aD - aUU) / den;
                double tf = tD;
                if (tD >= 0 && tD <= 1) tf = tD < 0.5 ? 2 * tD : 1.0;
                const double Co = std::abs(F) * h / M_.vol[D];
                if (Co > 0.7) tf = tD;
                else if (Co > 0.3) tf = tD + (tf - tD) * (0.7 - Co) / 0.4;
                const double gn = std::hypot(gx[D], gy[D]), dn = std::hypot(dx, dy);
                const double cosT = gn > 1e-12 ? std::abs(gx[D] * dx + gy[D] * dy) / (gn * dn) : 1.0;
                const double sq = std::sqrt(cosT);
                tf = tf * sq + tD * (1 - sq);
                af = std::clamp(tf * den + aUU, 0.0, 1.0);
            }
            flux[f] = F * af;
            sumAf[f] += af;
        }
        for (int f = ni; f < nf; ++f) {
            const double af = F_[f] >= 0 ? alpha_[M_.owner[f]] : ab[f];
            flux[f] = F_[f] * af;
            sumAf[f] += af;
        }
        for (int f = 0; f < nf; ++f) {
            alpha_[M_.owner[f]] -= h * flux[f] / M_.vol[M_.owner[f]];
            if (M_.neigh[f] >= 0) alpha_[M_.neigh[f]] += h * flux[f] / M_.vol[M_.neigh[f]];
        }
        for (double& a : alpha_) a = std::clamp(a, 0.0, 1.0);
    }
    for (int f = 0; f < nf; ++f) alphaF_[f] = sumAf[f] / nSub;
}

double FluidSolver::waterVolume() const {
    double V = 0;
    for (int c = 0; c < M_.nCells; ++c) V += alpha_[c] * M_.vol[c];
    return V;
}

double FluidSolver::pressureAt(int c) const {
    return P_.vof ? p_[c] - rhoC_[c] * P_.gravity * M_.cy[c] : p_[c] + hydrostatic(M_.cy[c]);
}

bool FluidSolver::simpleIteration(double& resU, double& resV, double& resC) {
    // hiệu chỉnh Majumdar (1 − α)·(F − ū·S), lấy theo nghiệm trước khi giải động lượng
    std::vector<double> mj;
    faceFlux(mj);
    for (int f = 0; f < M_.nInternal; ++f) mj[f] = (1 - P_.alphaU) * (F_[f] - mj[f]);
    momentum(resU, resV);
    resC = pressureCorrection(mj);
    if (P_.turbulence != TurbModel::Laminar) solveTurbulence();
    if (!std::isfinite(resU) || !std::isfinite(resV) || !std::isfinite(resC) || resC > 1e6) return false;
    for (int c = 0; c < M_.nCells; ++c)
        if (!std::isfinite(u_[c]) || !std::isfinite(v_[c]) || !std::isfinite(p_[c])) return false;
    return true;
}

void FluidSolver::record(double resU, double resV, double resC) {
    FluidHistory h;
    h.iter = iter_;
    h.time = time_;
    h.resU = resU;
    h.resV = resV;
    h.resC = resC;
    const FluidForces fo = forces();
    h.Cd = fo.Cd;
    h.Cl = fo.Cl;
    hist_.push_back(h);
}

bool FluidSolver::step() {
    if (diverged_ || M_.nCells == 0) return false;
    double ru = 0, rv = 0, rc = 0;
    if (!P_.transient) {
        if (converged_) return true;
        if (!simpleIteration(ru, rv, rc)) {
            diverged_ = true;
            return false;
        }
        ++iter_;
        record(ru, rv, rc);
        converged_ = iter_ > 5 && std::max({ru, rv, rc}) < P_.tol;
        return true;
    }
    // ---- bước thời gian: VOF tự chọn Δt theo số Courant của dòng chảy và của sóng trọng lực √(g·H)
    dt_ = P_.dt;
    if (P_.vof) {
        double hMin = std::numeric_limits<double>::max();
        for (double V : M_.vol) hMin = std::min(hMin, std::sqrt(V));
        const double wave = std::sqrt(P_.gravity * std::max(M_.y1 - M_.y0, 1e-9)) / hMin;
        dt_ = std::min(P_.dt, P_.maxCourant / std::max(courantRate(), wave));
    }
    // ---- lưu nghiệm cũ, (F − ū·S) cũ cho hiệu chỉnh Choi
    const int ni = M_.nInternal;
    unn_.swap(un_);
    vnn_.swap(vn_);
    Fnn_.swap(Fn_);
    un_ = u_;
    vn_ = v_;
    kn_ = k_;
    en_ = eps_;
    faceFlux(Fn_);
    for (int f = 0; f < ni; ++f) Fn_[f] = F_[f] - Fn_[f];
    nOld_ = std::min(nOld_ + 1, 2);
    time_ += dt_;
    for (int k = 0; k < P_.outerIters; ++k) {
        if (!simpleIteration(ru, rv, rc)) {
            diverged_ = true;
            return false;
        }
        if (k > 0 && std::max({ru, rv, rc}) < P_.outerTol) break;
    }
    // VOF: vận chuyển α theo lưu lượng đã bảo toàn khối lượng của bước này, cập nhật ρ, μ cho bước sau
    if (P_.vof) {
        advectAlpha(dt_);
        updateProperties();
    }
    ++iter_;
    record(ru, rv, rc);
    return true;
}

// ============================================================================ hậu xử lý
std::vector<double> FluidSolver::cellField(FluidField f) const {
    const int nc = M_.nCells;
    std::vector<double> r(nc, 0);
    const double q = 0.5 * P_.rho * Uref_ * Uref_;
    if (f == FluidField::Vorticity) {
        std::vector<double> ub, vb, gux, guy, gvx, gvy;
        boundaryVelocity(ub, vb);
        gradient(u_, ub, gux, guy);
        gradient(v_, vb, gvx, gvy);
        for (int c = 0; c < nc; ++c) r[c] = gvx[c] - guy[c];
        return r;
    }
    for (int c = 0; c < nc; ++c) switch (f) {
        case FluidField::Speed: r[c] = std::hypot(u_[c], v_[c]); break;
        case FluidField::U: r[c] = u_[c]; break;
        case FluidField::V: r[c] = v_[c]; break;
        case FluidField::P: r[c] = P_.vof ? pressureAt(c) : p_[c]; break;
        case FluidField::Cp: r[c] = (P_.vof ? pressureAt(c) : p_[c]) / q; break;
        case FluidField::PTotal: r[c] = pressureAt(c); break;
        case FluidField::Alpha: r[c] = alpha_[c]; break;
        case FluidField::TurbK: r[c] = k_[c]; break;
        default: break;
        }
    return r;
}

std::vector<double> FluidSolver::nodeField(FluidField f, int nodeCount) const {
    if (f == FluidField::PTotal && !P_.vof) {
        // phần thủy tĩnh tuyến tính: cộng đúng tại tọa độ nút thay vì lấy trung bình từ tâm ô
        std::vector<double> r = nodeField(FluidField::P, nodeCount);
        for (size_t n = 0; n < r.size(); ++n) r[n] += hydrostatic(M_.py[n]);
        return r;
    }
    const std::vector<double> cv = cellField(f);
    nodeCount = std::min(nodeCount, (int)M_.px.size());
    const std::vector<double>& nx = M_.px;
    const std::vector<double>& ny = M_.py;
    std::vector<double> sum(nodeCount, 0), wsum(nodeCount, 0);
    for (int c = 0; c < M_.nCells; ++c)
        for (int k = M_.cellStart[c]; k < M_.cellStart[c + 1]; ++k) {
            const int n = M_.cellNodes[k];
            if (n >= nodeCount) continue;
            const double w = 1.0 / std::max(std::hypot(nx[n] - M_.cx[c], ny[n] - M_.cy[c]), 1e-12);
            sum[n] += w * cv[c];
            wsum[n] += w;
        }
    std::vector<double> r(nodeCount, std::numeric_limits<double>::quiet_NaN());
    for (int n = 0; n < nodeCount; ++n)
        if (wsum[n] > 0) r[n] = sum[n] / wsum[n];
    if (f == FluidField::Speed || f == FluidField::U || f == FluidField::V) {
        // nút trên biên có vận tốc cho trước: ưu tiên thành cứng > thành di động > cửa vào
        std::vector<int> prio(nodeCount, 0);
        for (int fc = M_.nInternal; fc < M_.nFaces; ++fc) {
            const FluidBc& b = P_.bc[M_.side[fc]];
            const int pr = b.kind == FluidBcKind::Wall ? 3 : b.kind == FluidBcKind::MovingWall ? 2 : b.kind == FluidBcKind::Inlet ? 1 : 0;
            if (!pr) continue;
            for (int n : {M_.fa[fc], M_.fb[fc]}) {
                if (n >= nodeCount || prio[n] >= pr) continue;
                prio[n] = pr;
                double k = 1;
                if (b.kind == FluidBcKind::Inlet && b.parabolic) {
                    const int s = M_.side[fc];
                    const bool vert = s == (int)FluidSide::Left || s == (int)FluidSide::Right;
                    const double L = sideMax_[s] - sideMin_[s];
                    const double c = vert ? -M_.sa * nx[n] + M_.ca * ny[n] : M_.ca * nx[n] + M_.sa * ny[n];
                    const double t = L > 0 ? (c - sideMin_[s]) / L : 0.5;
                    k = 4 * t * (1 - t);
                }
                const double uu = b.kind == FluidBcKind::Wall ? 0 : b.u * k, vv = b.kind == FluidBcKind::Wall ? 0 : b.v * k;
                r[n] = f == FluidField::U ? uu : f == FluidField::V ? vv : std::hypot(uu, vv);
            }
        }
    }
    return r;
}

int FluidSolver::locate(double x, double y, int hint) const {
    auto inside = [&](int c, int& exitFace) {
        double worst = 0;
        exitFace = -1;
        for (int k = M_.cfStart[c]; k < M_.cfStart[c + 1]; ++k) {
            const int f = M_.cfFace[k];
            const double sg = M_.owner[f] == c ? 1 : -1;
            const double L = std::hypot(M_.sx[f], M_.sy[f]);
            const double d = sg * ((x - M_.fx[f]) * M_.sx[f] + (y - M_.fy[f]) * M_.sy[f]) / L;
            if (d > worst + 1e-12 * L) {
                worst = d;
                exitFace = f;
            }
        }
        return exitFace < 0;
    };
    if (hint >= 0 && hint < M_.nCells) {
        int c = hint;
        for (int s = 0; s < 64; ++s) {
            int ef;
            if (inside(c, ef)) return c;
            const int nb = M_.owner[ef] == c ? M_.neigh[ef] : M_.owner[ef];
            if (nb < 0) break;
            c = nb;
        }
    }
    if (x < M_.x0 || x > M_.x1 || y < M_.y0 || y > M_.y1) return -1;
    for (int c = 0; c < M_.nCells; ++c) {
        int ef;
        if (inside(c, ef)) return c;
    }
    return -1;
}

bool FluidSolver::probe(double x, double y, double& u, double& v, double& p) const {
    const int c = locate(x, y);
    if (c < 0) return false;
    // tái dựng tuyến tính trong ô: φ(x) = φc + ∇φc·(x − xc)
    std::vector<double> ub, vb, pb, gx, gy;
    boundaryVelocity(ub, vb);
    const double rx = x - M_.cx[c], ry = y - M_.cy[c];
    gradient(u_, ub, gx, gy);
    u = u_[c] + gx[c] * rx + gy[c] * ry;
    gradient(v_, vb, gx, gy);
    v = v_[c] + gx[c] * rx + gy[c] * ry;
    p = p_[c] + pgx_[c] * rx + pgy_[c] * ry;
    if (P_.vof) p -= rhoC_[c] * P_.gravity * y; // p = p_rgh + ρ·g·x
    return true;
}

FluidForces FluidSolver::forces() const {
    FluidForces r;
    if (!M_.sideFaces[(int)FluidSide::Inner]) return r;
    r.valid = true;
    std::vector<double> ub, vb;
    boundaryVelocity(ub, vb);
    for (int f = M_.nInternal; f < M_.nFaces; ++f) {
        if (M_.side[f] != (int)FluidSide::Inner) continue;
        const int c = M_.owner[f];
        const double rx = M_.fx[f] - M_.cx[c], ry = M_.fy[f] - M_.cy[c];
        // phần thủy tĩnh: một pha ρ·g·(yFree − y); VOF: tách theo áp suất thủy tĩnh của nước ρ·g·(y_mặt − y) không có
        // nghĩa khi mặt nước thay đổi, nên toàn bộ áp suất tính vào phần động
        const double ph = P_.vof ? 0.0 : hydrostatic(M_.fy[f]);
        const double pb = p_[c] + pgx_[c] * rx + pgy_[c] * ry + (P_.vof ? -rhoC_[c] * P_.gravity * M_.fy[f] : ph);
        // S hướng ra khỏi chất lỏng (vào vật cản): áp suất đẩy vật theo chiều S
        r.Fpx += pb * M_.sx[f];
        r.Fpy += pb * M_.sy[f];
        r.Fhx += ph * M_.sx[f];
        r.Fhy += ph * M_.sy[f];
        const double L = std::hypot(M_.sx[f], M_.sy[f]), nx = M_.sx[f] / L, ny = M_.sy[f] / L;
        const double dn = std::abs(rx * nx + ry * ny);
        double du = u_[c] - ub[f], dv = v_[c] - vb[f];
        const double dnrm = du * nx + dv * ny;
        du -= dnrm * nx;
        dv -= dnrm * ny;
        if (P_.bc[(int)FluidSide::Inner].kind == FluidBcKind::Symmetry) du = dv = 0;
        r.Fvx += muWall_[f] * L * du / dn;
        r.Fvy += muWall_[f] * L * dv / dn;
    }
    r.Fx = r.Fpx + r.Fvx;
    r.Fy = r.Fpy + r.Fvy;
    r.Fd = M_.ca * r.Fx + M_.sa * r.Fy;
    r.Fl = -M_.sa * r.Fx + M_.ca * r.Fy;
    const double q = 0.5 * P_.rho * Uref_ * Uref_ * Lref_;
    const double dx = r.Fx - r.Fhx, dy = r.Fy - r.Fhy;
    r.Cd = (M_.ca * dx + M_.sa * dy) / q;
    r.Cl = (-M_.sa * dx + M_.ca * dy) / q;
    return r;
}

std::vector<std::array<double, 2>> FluidSolver::innerCp() const {
    std::vector<std::array<double, 2>> r;
    double xc = 0, yc = 0, Lsum = 0;
    for (int f = M_.nInternal; f < M_.nFaces; ++f) {
        if (M_.side[f] != (int)FluidSide::Inner) continue;
        const double L = std::hypot(M_.sx[f], M_.sy[f]);
        xc += M_.fx[f] * L;
        yc += M_.fy[f] * L;
        Lsum += L;
    }
    if (Lsum <= 0) return r;
    xc /= Lsum;
    yc /= Lsum;
    const double q = 0.5 * P_.rho * Uref_ * Uref_;
    for (int f = M_.nInternal; f < M_.nFaces; ++f) {
        if (M_.side[f] != (int)FluidSide::Inner) continue;
        const int c = M_.owner[f];
        const double pb = p_[c] + pgx_[c] * (M_.fx[f] - M_.cx[c]) + pgy_[c] * (M_.fy[f] - M_.cy[c]) -
                          (P_.vof ? rhoC_[c] * P_.gravity * M_.fy[f] : 0.0);
        // góc tính từ điểm đón dòng (phía −x), dương ở phía trên
        const double dx = M_.fx[f] - xc, dy = M_.fy[f] - yc;
        const double th = std::atan2(-M_.sa * dx + M_.ca * dy, -(M_.ca * dx + M_.sa * dy)) * 180 / kPi;
        r.push_back({th, pb / q});
    }
    std::sort(r.begin(), r.end());
    return r;
}

double FluidSolver::strouhal() const {
    if (!P_.transient || hist_.size() < 20) return 0;
    const double tEnd = hist_.back().time, tStart = 0.5 * tEnd;
    double mean = 0;
    int n = 0;
    for (const auto& h : hist_)
        if (h.time >= tStart) {
            mean += h.Cl;
            ++n;
        }
    if (n < 10) return 0;
    mean /= n;
    double amp = 0;
    for (const auto& h : hist_)
        if (h.time >= tStart) amp = std::max(amp, std::abs(h.Cl - mean));
    if (amp < 1e-6) return 0;
    std::vector<double> up; // thời điểm Cl cắt giá trị trung bình theo chiều đi lên
    for (size_t i = 1; i < hist_.size(); ++i) {
        if (hist_[i - 1].time < tStart) continue;
        const double a = hist_[i - 1].Cl - mean, b = hist_[i].Cl - mean;
        if (a < 0 && b >= 0) up.push_back(hist_[i - 1].time + (hist_[i].time - hist_[i - 1].time) * (-a) / (b - a));
    }
    if (up.size() < 3) return 0;
    const double T = (up.back() - up.front()) / (up.size() - 1);
    return T > 0 ? Lref_ / (T * Uref_) : 0;
}

double FluidSolver::massImbalance() const {
    std::vector<double> imb(M_.nCells, 0);
    double s = 0;
    for (int f = 0; f < M_.nFaces; ++f) {
        imb[M_.owner[f]] += F_[f];
        if (M_.neigh[f] >= 0) imb[M_.neigh[f]] -= F_[f];
        s += std::abs(F_[f]);
    }
    double r = 0;
    for (double v : imb) r += std::abs(v);
    return s > 0 ? r / s : 0;
}

std::vector<std::array<double, 2>> FluidSolver::defaultSeeds(int n) const {
    std::vector<std::array<double, 2>> s;
    n = std::max(n, 2);
    // hình chữ nhật bao trong hệ tọa độ dòng chảy
    double X0 = std::numeric_limits<double>::max(), X1 = -X0, Y0 = X0, Y1 = -X0;
    for (int id : M_.cellNodes) {
        const double x = M_.ca * M_.px[id] + M_.sa * M_.py[id], y = -M_.sa * M_.px[id] + M_.ca * M_.py[id];
        X0 = std::min(X0, x), X1 = std::max(X1, x), Y0 = std::min(Y0, y), Y1 = std::max(Y1, y);
    }
    auto add = [&](double x, double y) { s.push_back({M_.ca * x - M_.sa * y, M_.sa * x + M_.ca * y}); };
    const double eps = 1e-3 * std::max(X1 - X0, Y1 - Y0);
    bool inlet = false;
    for (int sd = 0; sd < (int)FluidSide::Inner; ++sd) {
        if (!M_.sideFaces[sd] || P_.bc[sd].kind != FluidBcKind::Inlet) continue;
        inlet = true;
        for (int i = 0; i < n; ++i) {
            const double t = sideMin_[sd] + (sideMax_[sd] - sideMin_[sd]) * (i + 0.5) / n;
            switch ((FluidSide)sd) {
            case FluidSide::Left: add(X0 + eps, t); break;
            case FluidSide::Right: add(X1 - eps, t); break;
            case FluidSide::Bottom: add(t, Y0 + eps); break;
            default: add(t, Y1 - eps); break;
            }
        }
    }
    if (!inlet) {
        // không có cửa vào (hộp kín): các điểm trên đường thẳng qua giữa miền
        for (int i = 0; i < n; ++i) add(0.5 * (X0 + X1), Y0 + (Y1 - Y0) * (i + 0.5) / n);
    }
    return s;
}

std::vector<std::vector<std::array<double, 2>>> FluidSolver::streamlines(const std::vector<std::array<double, 2>>& seeds,
                                                                       int maxSteps) const {
    std::vector<std::vector<std::array<double, 2>>> out;
    if (M_.nCells == 0) return out;
    std::vector<double> ub, vb, gux, guy, gvx, gvy;
    boundaryVelocity(ub, vb);
    gradient(u_, ub, gux, guy);
    gradient(v_, vb, gvx, gvy);
    double umax = 0;
    for (int c = 0; c < M_.nCells; ++c) umax = std::max(umax, std::hypot(u_[c], v_[c]));
    if (umax <= 0) return out;
    auto vel = [&](int c, double x, double y, double& a, double& b) {
        const double rx = x - M_.cx[c], ry = y - M_.cy[c];
        a = u_[c] + gux[c] * rx + guy[c] * ry;
        b = v_[c] + gvx[c] * rx + gvy[c] * ry;
    };
    auto trace = [&](double x, double y, int c, double dir) {
        std::vector<std::array<double, 2>> pts{{x, y}};
        const double x0 = x, y0 = y;
        double len = 0;
        for (int s = 0; s < maxSteps; ++s) {
            const double h = 0.35 * std::sqrt(M_.vol[c]);
            double a, b;
            vel(c, x, y, a, b);
            double sp = std::hypot(a, b);
            if (sp < 1e-4 * umax) break;
            const double xm = x + 0.5 * h * dir * a / sp, ym = y + 0.5 * h * dir * b / sp;
            const int cm = locate(xm, ym, c);
            if (cm < 0) break;
            vel(cm, xm, ym, a, b);
            sp = std::hypot(a, b);
            if (sp < 1e-4 * umax) break;
            const double xn = x + h * dir * a / sp, yn = y + h * dir * b / sp;
            const int cn = locate(xn, yn, cm);
            if (cn < 0) break;
            x = xn;
            y = yn;
            c = cn;
            len += h;
            pts.push_back({x, y});
            // đường dòng khép kín (xoáy)
            if (len > 20 * h && std::hypot(x - x0, y - y0) < 0.6 * h) {
                pts.push_back({x0, y0});
                break;
            }
        }
        return pts;
    };
    for (const auto& sd : seeds) {
        const int c = locate(sd[0], sd[1]);
        if (c < 0) continue;
        auto fw = trace(sd[0], sd[1], c, 1);
        auto bw = trace(sd[0], sd[1], c, -1);
        std::reverse(bw.begin(), bw.end());
        if (!bw.empty()) bw.pop_back();
        bw.insert(bw.end(), fw.begin(), fw.end());
        if (bw.size() > 2) out.push_back(std::move(bw));
    }
    return out;
}

// ============================================================================ bài toán từ miền chảy khai báo
std::string buildFlowProblem(const Model& m, Model& mesh, FluidParams& p, std::string* note) {
    const FlowDomain& f = m.flow;
    if (!f.active) return "Chưa khai báo dòng chảy (Trang chủ › Dòng chảy).";
    const double W = f.width();
    if (W <= 0 || f.length <= 0) return "Miền chảy có kích thước bằng 0.";
    if (f.rho <= 0 || f.mu <= 0) return "Khối lượng riêng và độ nhớt phải > 0.";
    double nx, ny;
    f.direction(nx, ny);
    double c[4][2];
    f.corners(c);
    DomainSpec ds;
    for (auto& q : c) ds.outer.push_back({q[0], q[1]});
    if (polygonArea(ds.outer) < 0) std::reverse(ds.outer.begin(), ds.outer.end());

    // vật cản: các vòng thanh khép kín; vòng nằm ngoài miền chảy thì bỏ qua
    std::vector<int> bars;
    for (size_t e = 0; e < m.elements.size(); ++e)
        if (m.elements[e].isBar()) bars.push_back((int)e);
    int ignored = 0, outside = 0;
    const auto loops = closedLoops(m, bars, &ignored);
    double minSize = 1e300;
    for (const auto& L : loops) {
        Polygon2 P;
        for (int id : L) P.push_back({m.nodes[id].x, m.nodes[id].y});
        int in = 0;
        for (const auto& v : P) in += pointInPolygon(ds.outer, v[0], v[1]);
        // mọi đỉnh nằm ngoài nhưng cạnh vẫn có thể cắt ngang miền (tường chắn ngang)
        auto crosses = [&](const std::array<double, 2>& a, const std::array<double, 2>& b) {
            for (size_t j = 0; j < ds.outer.size(); ++j) {
                const auto& c0 = ds.outer[j];
                const auto& d0 = ds.outer[(j + 1) % ds.outer.size()];
                auto o = [](const std::array<double, 2>& p, const std::array<double, 2>& q, const std::array<double, 2>& r) {
                    return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
                };
                if ((o(c0, d0, a) > 0) != (o(c0, d0, b) > 0) && (o(a, b, c0) > 0) != (o(a, b, d0) > 0)) return true;
            }
            return false;
        };
        for (size_t i = 0; in == 0 && i < P.size(); ++i)
            if (crosses(P[i], P[(i + 1) % P.size()])) in = 1;
        if (in == 0) {
            ++outside;
            continue;
        }
        minSize = std::min(minSize, std::sqrt(std::abs(polygonArea(P))));
        ds.holes.push_back(P);
    }
    std::string msg;
    if (ignored) msg += std::to_string(ignored) + " nhóm thanh không khép kín bị bỏ qua (khép kín đường vẽ bằng C). ";
    if (outside) msg += std::to_string(outside) + " đường khép kín nằm ngoài miền chảy bị bỏ qua. ";
    // kích thước ô: tự chọn theo bề rộng dòng chảy và kích thước vật cản
    double h = f.h > 0 ? f.h : W / 40;
    if (f.h <= 0 && !ds.holes.empty()) h = std::min(h, minSize / 12);
    ds.h = h;
    ds.hFar = f.hFar > 0 ? std::max(f.hFar, h) : std::max(h, W / 12);
    ds.grading = 0.2;
    mesh = Model{};
    Material mat;
    mat.name = "Chất lỏng";
    mat.E = 1;
    mesh.materials.push_back(mat);
    mesh.cases.push_back(LoadCase("TH1"));
    std::string err;
    const int n = meshDomain(mesh, ds, 0, &err);
    if (n <= 0) return "Không tạo được lưới: " + err;
    {
        // vật cản chắn ngang toàn bộ dòng chảy chia miền thành nhiều phần rời nhau
        FluidMesh fm;
        buildFluidMesh(mesh, fm);
        std::vector<char> seen(fm.nCells, 0);
        std::vector<int> st{0};
        seen[0] = 1;
        int cnt = 1;
        while (!st.empty()) {
            const int c = st.back();
            st.pop_back();
            for (int k = fm.cfStart[c]; k < fm.cfStart[c + 1]; ++k) {
                const int fc = fm.cfFace[k], o = fm.owner[fc] == c ? fm.neigh[fc] : fm.owner[fc];
                if (o >= 0 && !seen[o]) {
                    seen[o] = 1;
                    ++cnt;
                    st.push_back(o);
                }
            }
        }
        if (cnt < fm.nCells) return "Vật cản chắn kín dòng chảy: miền chảy bị chia thành nhiều phần rời nhau.";
    }
    char b[160];
    std::snprintf(b, sizeof b, "Lưới tam giác %d ô (ô sát vật cản %.3g m, ô xa %.3g m), %zu vật cản.", n, ds.h, ds.hFar,
                  ds.holes.size());
    msg = b + std::string(msg.empty() ? "" : " ") + msg;
    if (note) *note = msg;

    // điều kiện biên trong hệ tọa độ dòng chảy: trái = cửa vào, phải = cửa ra, dưới/trên = hai biên bên
    p.rho = f.rho;
    p.mu = f.mu;
    p.flowAngle = std::atan2(ny, nx);
    const double Umax = f.parabolic ? 1.5 * f.U : f.U;
    p.bc[(int)FluidSide::Left] = {FluidBcKind::Inlet, Umax * nx, Umax * ny, f.parabolic};
    p.bc[(int)FluidSide::Right] = {FluidBcKind::Outlet};
    p.bc[(int)FluidSide::Inner] = {FluidBcKind::Wall};
    p.hydroG = 0;
    if (f.section == 0) {
        // mặt bằng: hai bờ
        const FluidBcKind wall = f.slipWalls ? FluidBcKind::Symmetry : FluidBcKind::Wall;
        p.bc[(int)FluidSide::Bottom] = {wall};
        p.bc[(int)FluidSide::Top] = {wall};
    } else {
        // mặt cắt đứng: trọng lực theo −y. Cạnh bên thấp hơn (theo cao độ trung bình) là mặt đất / đáy, cạnh kia là mặt
        // thoáng / không khí / trần; dòng chảy được phép nghiêng (kênh dốc, mái dốc). Dòng chảy thẳng đứng (ống đứng):
        // hai cạnh bên đều là thành.
        const double yl0 = -std::sin(p.flowAngle) * c[0][0] + std::cos(p.flowAngle) * c[0][1]; // y' của cạnh A–(A+Ln)
        const double yl1 = -std::sin(p.flowAngle) * c[1][0] + std::cos(p.flowAngle) * c[1][1]; // y' của cạnh B–(B+Ln)
        const double mA = 0.5 * (c[0][1] + c[3][1]), mB = 0.5 * (c[1][1] + c[2][1]);
        const double span = std::max({std::abs(c[2][1] - c[0][1]), std::abs(c[2][0] - c[0][0]), W, f.length});
        const bool vertical = std::abs(mA - mB) < 1e-6 * span;
        const bool groundIsB = mB < mA;
        // cạnh nào là "dưới" (y' nhỏ) trong hệ tọa độ dòng chảy
        const bool bLow = yl1 < yl0;
        const FluidSide sideB = bLow ? FluidSide::Bottom : FluidSide::Top, sideA = bLow ? FluidSide::Top : FluidSide::Bottom;
        const FluidSide ground = groundIsB ? sideB : sideA, top = groundIsB ? sideA : sideB;
        // đầu A/B của cửa vào nằm ở đáy, và góc hạ lưu của cạnh trên
            const double yTopOut = groundIsB ? c[3][1] : c[2][1];
        if (vertical) {
            if (f.freeSurface) return "Mặt thoáng tự do: dòng chảy không được thẳng đứng.";
            p.bc[(int)FluidSide::Bottom] = {FluidBcKind::Wall};
            p.bc[(int)FluidSide::Top] = {FluidBcKind::Wall};
        } else {
            p.bc[(int)ground] = {FluidBcKind::Wall};
            p.bc[(int)top] = {f.topWall ? FluidBcKind::Wall : FluidBcKind::Symmetry};
        }
        if (f.freeSurface) {
            // hai pha nước – không khí: biên trên thông khí, nước vào dưới mực nước ở cửa vào
            if (f.medium == 1 || f.rho < 100) return "Mặt thoáng tự do: chất lỏng phải là nước, dầu... (nặng hơn không khí).";
            if (f.waterDepth <= 0 || f.waterDepth >= W) return "Mặt thoáng tự do: độ sâu nước ở cửa vào phải > 0 và < chiều cao AB.";
            p.vof = true;
            p.rho2 = 1.204;
            p.mu2 = 1.825e-5;
            p.gravity = 9.81;
            p.inletLevel = f.waterLevelY(f.waterDepth);
            p.initLevel = f.waterLevelY(f.initDepth >= 0 ? f.initDepth : f.waterDepth);
            p.initGate = f.gate;
            p.bc[(int)top] = {FluidBcKind::Opening};
            p.hydroG = 0;
            p.Uref = f.U > 0 ? f.U : std::sqrt(9.81 * f.waterDepth);
        } else if (f.hydrostatic) {
            p.hydroG = 9.81;
            // mặt thoáng tại góc hạ lưu của cạnh trên (mặt thoáng nằm ngang khi dòng chảy nằm ngang)
            p.yFree = vertical ? std::max({c[0][1], c[1][1], c[2][1], c[3][1]}) : yTopOut;
        }
    }
    p.transient = f.transient || p.vof;
    p.dt = f.dt;
    p.endTime = f.endTime;
    if (!p.vof) p.Uref = f.U;
    p.Lref = 0;
    const double Re = f.rho * std::max(f.U, 1e-12) * W / f.mu;
    p.turbulence = f.turbulence == 2 || (f.turbulence == 0 && Re > 2300) ? TurbModel::KEpsilon : TurbModel::Laminar;
    p.turbIntensity = f.turbIntensity;
    if (p.turbulence == TurbModel::KEpsilon && note) *note += " Chảy rối: mô hình k-ε, hàm tường.";
    if (p.turbulence == TurbModel::KEpsilon && !p.transient) {
        // RANS dừng với sơ đồ bậc hai quanh vật cản tù vẫn dao động (xoáy bong); bậc nhất hội tụ ổn định
        p.scheme = ConvScheme::Upwind;
        if (note) *note += " Bài toán dừng: sơ đồ ngược dòng bậc nhất.";
    }
    if (!p.transient && p.turbulence == TurbModel::Laminar && Re > 2000) {
        // dòng chảy tầng ở Re lớn không có nghiệm dừng, lặp SIMPLE dừng phân kỳ. Tính theo thời gian (Δt theo số
        // Courant ≈ 2 với ô sát vật cản) qua khoảng 3 lần thời gian chất lỏng chảy hết miền.
        p.transient = true;
        p.timeScheme = TimeScheme::BDF2;
        p.dt = 2 * ds.h / std::max(f.U, 1e-12);
        p.endTime = 3 * f.length / std::max(f.U, 1e-12);
        p.outerIters = std::max(p.outerIters, 10);
        char t[200];
        std::snprintf(t, sizeof t, " Re lớn, dòng chảy tầng không có nghiệm dừng: tính theo thời gian, Δt = %.3g s, đến t = %.3g s.",
                      p.dt, p.endTime);
        if (note) *note += t;
    }
    if (p.transient) {
        p.alphaU = std::max(p.alphaU, 0.8);
        p.alphaP = std::max(p.alphaP, 0.5);
    }
    return {};
}

// ============================================================================ tạo lưới kênh
namespace {

// Các vị trí 0 = t0 < t1 < ... < tn = len: bước đầu h0, tăng theo cấp số nhân q đến hmax rồi giữ nguyên
std::vector<double> graded(double len, double h0, double hmax, double q = 1.1) {
    std::vector<double> d;
    double s = 0, h = h0;
    hmax = std::max(hmax, h0);
    while (s < len) {
        d.push_back(h);
        s += h;
        h = std::min(h * q, hmax);
    }
    // bỏ đoạn cuối nếu quá ngắn rồi co giãn cho khớp đúng len
    if (d.size() > 1 && s - len > 0.5 * d.back()) {
        s -= d.back();
        d.pop_back();
    }
    std::vector<double> t{0};
    for (double v : d) t.push_back(t.back() + v * len / s);
    t.back() = len;
    return t;
}

} // namespace

int meshChannel(Model& m, const ChannelSpec& s, int mat) {
    if (s.L <= 0 || s.H <= 0 || s.h <= 0) return 0;
    const double X0 = s.x0, X1 = s.x0 + s.L, Y0 = s.y0, Y1 = s.y0 + s.H;
    const double hFar = std::max(s.hFar, s.h);
    std::vector<double> xs, ys;
    // hộp bao vật cản [bx0, bx1]×[by0, by1]; lưới đều n ô mỗi cạnh trong hộp
    double bx0 = 0, bx1 = 0, by0 = 0, by1 = 0, r = 0;
    int nbx = 0, nby = 0;
    auto appendUniform = [](std::vector<double>& v, double a, double b, int n) {
        for (int i = (v.empty() ? 0 : 1); i <= n; ++i) v.push_back(a + (b - a) * i / n);
    };
    if (s.obstacle == 1) {
        r = 0.5 * s.D;
        const double a = std::min({s.D, s.xc - X0, X1 - s.xc, s.yc - Y0, Y1 - s.yc});
        if (r <= 0 || a < 1.2 * r) return 0;
        bx0 = s.xc - a, bx1 = s.xc + a, by0 = s.yc - a, by1 = s.yc + a;
        nbx = nby = std::max(8, 2 * (int)std::ceil(a / s.h)); // số chẵn: có nút tại đỉnh và đáy trụ
    } else if (s.obstacle == 2) {
        bx0 = std::max(X0, s.xc - 0.5 * s.D), bx1 = std::min(X1, s.xc + 0.5 * s.D);
        by0 = std::max(Y0, s.yc - 0.5 * s.B), by1 = std::min(Y1, s.yc + 0.5 * s.B);
        if (bx1 - bx0 < 1e-9 || by1 - by0 < 1e-9) return 0;
        if (bx1 - bx0 >= s.L - 1e-9 && by1 - by0 >= s.H - 1e-9) return 0;
        nbx = std::max(2, (int)std::ceil((bx1 - bx0) / s.h - 1e-9));
        nby = std::max(2, (int)std::ceil((by1 - by0) / s.h - 1e-9));
    }
    if (s.obstacle == 0) {
        appendUniform(xs, X0, X1, std::max(1, (int)std::ceil(s.L / s.h - 1e-9)));
        appendUniform(ys, Y0, Y1, std::max(1, (int)std::ceil(s.H / s.h - 1e-9)));
    } else {
        // từ vật cản ra xa: bước h ở sát hộp, giãn dần tới hFar
        auto axis = [&](std::vector<double>& v, double A, double B, double b0, double b1, int nb) {
            const double hb = (b1 - b0) / nb;
            if (b0 > A + 1e-12) {
                const std::vector<double> g = graded(b0 - A, hb, hFar);
                for (int i = (int)g.size() - 1; i >= 0; --i) v.push_back(b0 - g[i]);
            }
            appendUniform(v, b0, b1, nb);
            if (B > b1 + 1e-12) {
                const std::vector<double> g = graded(B - b1, hb, hFar);
                for (size_t i = 1; i < g.size(); ++i) v.push_back(b1 + g[i]);
            }
        };
        axis(xs, X0, X1, bx0, bx1, nbx);
        axis(ys, Y0, Y1, by0, by1, nby);
    }
    const int NX = (int)xs.size(), NY = (int)ys.size();
    std::vector<int> id((size_t)NX * NY, -1);
    auto node = [&](int i, int j) {
        int& k = id[(size_t)j * NX + i];
        if (k < 0) {
            k = (int)m.nodes.size();
            m.nodes.push_back({xs[i], ys[j]});
        }
        return k;
    };
    int count = 0;
    auto cell = [&](int a, int b, int c, int d, bool alt) {
        Element e;
        e.mat = mat;
        e.t = 1.0;
        e.incompatible = false;
        e.n.fill(-1);
        if (!s.triangles) {
            e.type = ElemType::Quad4;
            e.n[0] = a, e.n[1] = b, e.n[2] = c, e.n[3] = d;
            m.elements.push_back(e);
            ++count;
            return;
        }
        e.type = ElemType::Tri3;
        if (alt) {
            e.n[0] = a, e.n[1] = b, e.n[2] = c;
            m.elements.push_back(e);
            e.n[0] = a, e.n[1] = c, e.n[2] = d;
        } else {
            e.n[0] = a, e.n[1] = b, e.n[2] = d;
            m.elements.push_back(e);
            e.n[0] = b, e.n[1] = c, e.n[2] = d;
        }
        m.elements.push_back(e);
        count += 2;
    };
    auto inBox = [&](double x, double y) { return s.obstacle != 0 && x > bx0 && x < bx1 && y > by0 && y < by1; };
    for (int j = 0; j + 1 < NY; ++j)
        for (int i = 0; i + 1 < NX; ++i) {
            if (inBox(0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1]))) continue;
            cell(node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1), (i + j) % 2 == 0);
        }
    if (s.obstacle == 1) {
        // lưới O quanh trụ: vòng ngoài là chu vi hộp (dùng chung nút với lưới chữ nhật)
        int ib0 = 0, jb0 = 0;
        while (ib0 < NX && xs[ib0] < bx0 - 1e-12) ++ib0;
        while (jb0 < NY && ys[jb0] < by0 - 1e-12) ++jb0;
        const int n = nbx, nt = 4 * n;
        std::vector<int> outer(nt);
        for (int k = 0; k < n; ++k) {
            outer[k] = node(ib0 + k, jb0);
            outer[n + k] = node(ib0 + n, jb0 + k);
            outer[2 * n + k] = node(ib0 + n - k, jb0 + n);
            outer[3 * n + k] = node(ib0, jb0 + n - k);
        }
        const double a = 0.5 * (bx1 - bx0);
        const double hc = 2 * kPi * r / nt, t1 = 0.5 * hc, q = 1.12;
        const int nr = std::max(3, (int)std::ceil(std::log(1 + (q - 1) * (a - r) / t1) / std::log(q)));
        std::vector<std::vector<int>> ring(nr + 1, std::vector<int>(nt));
        ring[nr] = outer;
        for (int l = 0; l < nr; ++l) {
            const double sl = (std::pow(q, l) - 1) / (std::pow(q, nr) - 1);
            for (int k = 0; k < nt; ++k) {
                const double th = 1.25 * kPi + 2 * kPi * k / nt; // góc dưới trái 225°, ngược chiều kim đồng hồ
                const double cxp = s.xc + r * std::cos(th), cyp = s.yc + r * std::sin(th);
                const Node& B = m.nodes[outer[k]];
                ring[l][k] = (int)m.nodes.size();
                m.nodes.push_back({cxp + sl * (B.x - cxp), cyp + sl * (B.y - cyp)});
            }
        }
        for (int l = 0; l < nr; ++l)
            for (int k = 0; k < nt; ++k) {
                const int k1 = (k + 1) % nt;
                cell(ring[l][k], ring[l + 1][k], ring[l + 1][k1], ring[l][k1], (k + l) % 2 == 0);
            }
    }
    return count;
}

} // namespace fem
