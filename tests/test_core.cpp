// Kiểm chứng lõi tính toán với lời giải giải tích.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>

#include "core/Analysis.h"
#include "core/Design.h"
#include "core/Dxf.h"
#include "core/Dynamics.h"
#include "core/Report.h"
#include "core/Eigen.h"
#include "core/Examples.h"
#include "core/Fluid.h"
#include "core/TriMesher.h"
#include "core/ForceMethod.h"
#include "core/Kinematics.h"
#include "core/SectionGeom.h"
#include "core/Library.h"
#include "core/Mesher.h"
#include "core/ModelIO.h"
#include "core/PlaneElements.h"
#include "core/PostProcess.h"

using namespace fem;

static int g_fail = 0, g_pass = 0;

static void check(bool ok, const char* what, double got, double expect) {
    if (ok) {
        ++g_pass;
    } else {
        ++g_fail;
        std::printf("  FAIL %s: got %.8g, expected %.8g\n", what, got, expect);
    }
}

static void near(const char* what, double got, double expect, double relTol = 1e-6, double absTol = 1e-9) {
    check(std::abs(got - expect) <= relTol * std::abs(expect) + absTol, what, got, expect);
}

static AnalysisResult run(const Model& m) {
    AnalysisResult r = analyze(m);
    if (!r.ok) {
        ++g_fail;
        std::printf("  FAIL analysis: %s\n", r.error.c_str());
    }
    return r;
}

// Nội suy giá trị trên thanh tại x (lấy điểm mẫu gần nhất phía trái)
static double at(const std::vector<double>& xs, const std::vector<double>& v, double x) {
    for (size_t k = 0; k + 1 < xs.size(); ++k)
        if (x >= xs[k] - 1e-12 && x <= xs[k + 1] + 1e-12) {
            const double h = xs[k + 1] - xs[k];
            if (h <= 0) return v[k];
            const double t = (x - xs[k]) / h;
            return v[k] + t * (v[k + 1] - v[k]);
        }
    return v.back();
}

static Model beamModel(double L, int nEl, double E, double A, double I, double h = 0.3) {
    Model m;
    Material mat;
    mat.E = E;
    mat.nu = 0.3;
    mat.alpha = 1e-5;
    mat.gamma = 25;
    m.materials.push_back(mat);
    Section s;
    s.A = A;
    s.I = I;
    s.h = h;
    m.sections.push_back(s);
    for (int i = 0; i <= nEl; ++i) m.nodes.push_back({L * i / nEl, 0});
    for (int i = 0; i < nEl; ++i) {
        Element e;
        e.type = ElemType::Frame;
        e.n = {i, i + 1, -1, -1};
        m.elements.push_back(e);
    }
    m.cases.push_back({"TH1"});
    return m;
}

static void fix(Model& m, int node, bool x, bool y, bool r) { m.supportFor(node).fix = {x, y, r}; }

static ElemLoad dload(int e, double q1, double q2, LoadDir d = LoadDir::GlobalY) {
    ElemLoad l;
    l.elem = e;
    l.q1 = q1;
    l.q2 = q2;
    l.dir = d;
    return l;
}

static void testSimpleBeam() {
    std::puts("Dầm đơn giản, tải phân bố đều");
    const double L = 6, q = 10, E = 2.1e8, I = 7080e-8;
    Model m = beamModel(L, 1, E, 46.5e-4, I);
    fix(m, 0, true, true, false);
    fix(m, 1, false, true, false);
    m.cases[0].elem.push_back(dload(0, -q, -q));
    auto r = run(m);
    if (!r.ok) return;
    const auto& b = r.sets[0].bars[0];
    near("Mmax = qL²/8", at(b.x, b.M, L / 2), q * L * L / 8);
    near("Q(0) = qL/2", b.Q.front(), q * L / 2);
    near("Q(L) = −qL/2", b.Q.back(), -q * L / 2);
    near("M(0) = 0", b.M.front(), 0, 0, 1e-9);
    near("Ry trái", r.sets[0].R[1], q * L / 2);
    near("Độ võng giữa nhịp 5qL⁴/384EI", at(b.x, b.v, L / 2), -5 * q * std::pow(L, 4) / (384 * E * I), 1e-6);
    near("Góc xoay đầu trái qL³/24EI", r.sets[0].U[2], -q * L * L * L / (24 * E * I));
}

static void testCantilever() {
    std::puts("Dầm công-xôn, lực tập trung đầu tự do (4 phần tử)");
    const double L = 3, P = 10, E = 3e7, I = 0.3 * std::pow(0.5, 3) / 12;
    Model m = beamModel(L, 4, E, 0.15, I);
    fix(m, 0, true, true, true);
    m.cases[0].nodal.push_back({4, {0, -P, 0}});
    auto r = run(m);
    if (!r.ok) return;
    near("Độ võng đầu PL³/3EI", r.sets[0].U[3 * 4 + 1], -P * L * L * L / (3 * E * I));
    near("M ngàm = −PL", r.sets[0].bars[0].M.front(), -P * L);
    near("Q = +P", r.sets[0].bars[2].Q[5], P);
    near("Phản lực mômen", r.sets[0].R[2], P * L);
}

static void testFixedFixed() {
    std::puts("Dầm hai đầu ngàm, tải phân bố đều");
    const double L = 8, q = 12;
    Model m = beamModel(L, 1, 3e7, 0.18, 0.0054);
    fix(m, 0, true, true, true);
    fix(m, 1, true, true, true);
    m.cases[0].elem.push_back(dload(0, -q, -q));
    auto r = run(m);
    if (!r.ok) return;
    const auto& b = r.sets[0].bars[0];
    near("M đầu = −qL²/12", b.M.front(), -q * L * L / 12);
    near("M giữa = qL²/24", at(b.x, b.M, L / 2), q * L * L / 24);
    near("M cuối = −qL²/12", b.M.back(), -q * L * L / 12);
}

static void testPointAndMoment() {
    std::puts("Dầm đơn giản: lực tập trung, mômen tập trung, tải tam giác");
    const double L = 10;
    {
        Model m = beamModel(L, 1, 2e8, 0.01, 1e-4);
        fix(m, 0, true, true, false);
        fix(m, 1, false, true, false);
        ElemLoad p;
        p.elem = 0;
        p.kind = ElemLoadKind::Point;
        p.q1 = -20;
        p.a = 4;
        m.cases[0].elem.push_back(p);
        auto r = run(m);
        if (r.ok) {
            const auto& b = r.sets[0].bars[0];
            near("M dưới lực = Pab/L", at(b.x, b.M, 4), 20 * 4 * 6 / L);
            const double E = 2e8, I = 1e-4, a = 4, bb = 6;
            near("Độ võng dưới lực Pa²b²/3EIL", at(b.x, b.v, 4), -20 * a * a * bb * bb / (3 * E * I * L), 1e-5);
        }
    }
    {
        Model m = beamModel(L, 1, 2e8, 0.01, 1e-4);
        fix(m, 0, true, true, false);
        fix(m, 1, false, true, false);
        ElemLoad p;
        p.elem = 0;
        p.kind = ElemLoadKind::Moment;
        p.q1 = 30;
        p.a = L / 2;
        m.cases[0].elem.push_back(p);
        auto r = run(m);
        if (r.ok) {
            const auto& b = r.sets[0].bars[0];
            size_t k = 0;
            while (b.x[k] < L / 2 - 1e-9) ++k;
            near("M trái điểm đặt = M0/2", b.M[k], 15);
            near("M phải điểm đặt = −M0/2", b.M[k + 1], -15);
            near("Q = M0/L", b.Q[0], 3);
        }
    }
    {
        Model m = beamModel(L, 1, 2e8, 0.01, 1e-4);
        fix(m, 0, true, true, false);
        fix(m, 1, false, true, false);
        m.cases[0].elem.push_back(dload(0, 0, -9));
        auto r = run(m);
        if (r.ok) {
            const auto& b = r.sets[0].bars[0];
            double mx = 0;
            for (double v : b.M) mx = std::max(mx, v);
            near("Tải tam giác: Mmax = qL²/(9√3)", mx, 9 * L * L / (9 * std::sqrt(3.0)), 2e-3);
            near("Phản lực trái qL/6", r.sets[0].R[1], 9 * L / 6);
        }
    }
}

static void testInclinedCantilever() {
    std::puts("Công-xôn xiên 30°, lực vuông góc trục");
    const double L = 4, P = 5, E = 2e8, I = 2e-5, ang = 30 * M_PI / 180;
    Model m = beamModel(L, 2, E, 0.01, I);
    for (int i = 0; i < 3; ++i) m.nodes[i] = {L * i / 2 * std::cos(ang), L * i / 2 * std::sin(ang)};
    fix(m, 0, true, true, true);
    m.cases[0].nodal.push_back({2, {P * std::sin(ang), -P * std::cos(ang), 0}});
    auto r = run(m);
    if (!r.ok) return;
    const double d = std::hypot(r.sets[0].U[6], r.sets[0].U[7]);
    near("Độ võng PL³/3EI", d, P * L * L * L / (3 * E * I), 1e-4);
    near("M ngàm = −PL", r.sets[0].bars[0].M.front(), -P * L);
}

static void testExamples() {
    std::puts("Các ví dụ mẫu: cân bằng tổng lực và kiểm tra riêng");
    for (const auto& ex : examples()) {
        Model m = ex.build();
        auto r = analyze(m);
        if (!r.ok) {
            ++g_fail;
            std::printf("  FAIL %s: %s\n", ex.name, r.error.c_str());
            continue;
        }
        for (const auto& w : r.warnings) std::printf("  [%s] cảnh báo: %s\n", ex.name, w.c_str());
        for (const auto& s : r.sets) {
            const double scale = 1 + std::abs(s.sumFx) + std::abs(s.sumFy);
            const std::string w = std::string(ex.name) + " / " + s.name;
            // chuyển vị cưỡng bức/nhiệt cho tổng tải = 0, phản lực tự cân bằng
            // đối xứng trục: lực hướng kính cân bằng bởi ứng suất vòng, không qua phản lực
            if (m.planeMode != PlaneMode::Axisymmetric)
                check(std::abs(s.sumFx + s.sumRx) < 1e-6 * scale, (w + ": ΣX").c_str(), s.sumRx, -s.sumFx);
            check(std::abs(s.sumFy + s.sumRy) < 1e-6 * scale, (w + ": ΣY").c_str(), s.sumRy, -s.sumFy);
            // kiểm tra nội lực đầu thanh tính lại từ phương trình cân bằng
            for (size_t e = 0; e < m.elements.size(); ++e) {
                if (!m.elements[e].isBar()) continue;
                const auto& b = s.bars[e];
                const auto& f = b.endForces;
                const double sc = 1 + std::abs(f[1]) + std::abs(f[4]) + std::abs(f[2]) + std::abs(f[5]);
                if (std::abs(b.M.back() - f[5]) > 1e-6 * sc || std::abs(b.Q.back() + f[4]) > 1e-6 * sc ||
                    std::abs(b.N.back() - f[3]) > 1e-6 * sc) {
                    ++g_fail;
                    std::printf("  FAIL %s pt %zu: cân bằng thanh M %.6g/%.6g Q %.6g/%.6g\n", w.c_str(), e + 1,
                                b.M.back(), f[5], b.Q.back(), -f[4]);
                } else ++g_pass;
            }
        }
        if (std::string(ex.name) == "Dầm đơn giản") {
            const auto& b = r.sets[0].bars[0];
            near("Ví dụ dầm đơn giản Mmax", at(b.x, b.M, 3), 45);
        }
        if (std::string(ex.name) == "Dầm tĩnh định nhiều nhịp") {
            near("Mômen tại khớp 1 = 0", r.sets[0].bars[1].M.back(), 0, 0, 1e-8);
            near("Mômen tại khớp 2 = 0", r.sets[0].bars[3].M.back(), 0, 0, 1e-8);
        }
        if (std::string(ex.name) == "Khung ba khớp") {
            // phản lực đứng = ql/2 = 48, lực xô H = M_dầm/f = (q·8²/8)/6 với tải theo hình chiếu
            near("Khung ba khớp: phản lực đứng", r.sets[0].R[1], 48, 1e-6);
            near("Khung ba khớp: lực xô ngang", std::abs(r.sets[0].R[0]), 12 * 64 / 8.0 / 6.0, 1e-6);
        }
        if (std::string(ex.name) == "Dàn cầu") {
            // thanh cánh dưới khoang 3: mômen lấy với nút cánh trên x = 6: N = (R·6 − 40·3)/h
            const double R = 100, M = R * 6 - 40 * 3;
            near("Dàn: lực thanh cánh dưới giữa nhịp", r.sets[0].bars[2].N[0], M / 3, 1e-6);
        }
    }
}

static void testHinge() {
    std::puts("Khớp đầu thanh: dầm 2 đoạn có khớp giữa, công-xôn kép");
    // Hai công-xôn nối khớp tại giữa, lực P tại khớp: mỗi bên chịu P/2
    const double L = 6, P = 12;
    Model m = beamModel(L, 2, 2e8, 0.01, 1e-4);
    m.elements[0].hingeJ = true;
    fix(m, 0, true, true, true);
    fix(m, 2, true, true, true);
    m.cases[0].nodal.push_back({1, {0, -P, 0}});
    auto r = run(m);
    if (!r.ok) return;
    near("M ngàm trái = −P/2·L/2", r.sets[0].bars[0].M.front(), -P / 2 * L / 2);
    near("M tại khớp = 0", r.sets[0].bars[0].M.back(), 0, 0, 1e-9);
}

static void testThermalAndSettlement() {
    std::puts("Nhiệt độ và lún gối");
    const double L = 5, E = 3e7, A = 0.15, I = 0.003125, h = 0.5, alpha = 1e-5;
    {
        Model m = beamModel(L, 1, E, A, I, h);
        fix(m, 0, true, true, true);
        fix(m, 1, true, true, true);
        ElemLoad t;
        t.elem = 0;
        t.kind = ElemLoadKind::Temperature;
        t.q1 = 30;
        t.q2 = 20;
        m.cases[0].elem.push_back(t);
        auto r = run(m);
        if (r.ok) {
            near("Ngàm 2 đầu: N = −EAαt", r.sets[0].bars[0].N[3], -E * A * alpha * 30);
            near("Ngàm 2 đầu: M = −EIαΔt/h", r.sets[0].bars[0].M[7], -E * I * alpha * 20 / h);
        }
    }
    {
        Model m = beamModel(L, 1, E, A, I, h);
        fix(m, 0, true, true, true);
        ElemLoad t;
        t.elem = 0;
        t.kind = ElemLoadKind::Temperature;
        t.q2 = 20;
        m.cases[0].elem.push_back(t);
        auto r = run(m);
        if (r.ok) {
            near("Công-xôn: độ võng đầu = καL²/2", r.sets[0].U[4], alpha * 20 / h * L * L / 2);
            near("Công-xôn: M = 0", r.sets[0].bars[0].M[0], 0, 0, 1e-9);
        }
    }
    {
        Model m = beamModel(L, 1, E, A, I, h);
        fix(m, 0, true, true, true);
        fix(m, 1, true, true, true);
        m.cases[0].disp.push_back({1, UY, -0.01});
        auto r = run(m);
        if (r.ok) near("Lún Δ: M ngàm = 6EIΔ/L²", std::abs(r.sets[0].bars[0].M[0]), 6 * E * I * 0.01 / (L * L));
    }
}

static void testSelfWeight() {
    std::puts("Trọng lượng bản thân");
    const double L = 6;
    Model m = beamModel(L, 3, 3e7, 0.3 * 0.6, 0.3 * 0.216 / 12);
    fix(m, 0, true, true, false);
    fix(m, 3, false, true, false);
    m.cases[0].selfWeight = 1.0;
    auto r = run(m);
    if (!r.ok) return;
    const double q = 25 * 0.18;
    near("Mmax = γAL²/8", r.sets[0].bars[1].M[20], q * L * L / 8, 1e-9);
}

static void testWinkler() {
    std::puts("Dầm dài trên nền Winkler, lực tập trung giữa dầm");
    const double L = 30, P = 400, E = 2.75e7, b = 1, hh = 0.6, k = 20000;
    const double I = b * hh * hh * hh / 12;
    const int n = 120;
    Model m = beamModel(L, n, E, b * hh, I, hh);
    for (auto& e : m.elements) e.kf = k * b;
    fix(m, n / 2, true, false, false);
    m.cases[0].nodal.push_back({n / 2, {0, -P, 0}});
    auto r = run(m);
    if (!r.ok) return;
    const double beta = std::pow(k * b / (4 * E * I), 0.25);
    near("w_max = Pβ/2k", -r.sets[0].U[3 * (n / 2) + 1], P * beta / (2 * k * b), 2e-3);
    near("M_max = P/4β", r.sets[0].bars[n / 2].M.front(), P / (4 * beta), 5e-3);
}

static Model plateModel(ElemType type, int nx, int ny, double W, double H, bool distort) {
    Model m;
    Material mat;
    mat.E = 1e6;
    mat.nu = 0.25;
    m.materials.push_back(mat);
    meshRectangle(m, 0, 0, W, H, nx, ny, type, 0, 0.5);
    if (distort) {
        // dịch chuyển các nút bên trong để lưới méo
        for (auto& nd : m.nodes)
            if (nd.x > 1e-9 && nd.x < W - 1e-9 && nd.y > 1e-9 && nd.y < H - 1e-9) {
                nd.x += 0.17 * W / nx * std::sin(7 * nd.y + 3 * nd.x);
                nd.y += 0.15 * H / ny * std::cos(5 * nd.x - 2 * nd.y);
            }
    }
    m.cases.push_back({"TH1"});
    return m;
}

static void testPatch() {
    std::puts("Patch test: kéo đều (CST, Q4, Q4 không tương thích; lưới méo)");
    const double W = 4, H = 2, sig = 100, t = 0.5;
    for (int variant = 0; variant < 3; ++variant) {
        const ElemType type = variant == 0 ? ElemType::Tri3 : ElemType::Quad4;
        Model m = plateModel(type, 4, 3, W, H, true);
        for (auto& e : m.elements) e.incompatible = variant == 2;
        const int nx = 4, ny = 3;
        for (int j = 0; j <= ny; ++j) {
            fix(m, j * (nx + 1), true, j == 0, false);
            const double w = (j == 0 || j == ny) ? 0.5 : 1.0;
            m.cases[0].nodal.push_back({j * (nx + 1) + nx, {sig * t * H / ny * w, 0, 0}});
        }
        auto r = run(m);
        if (!r.ok) continue;
        double maxErr = 0;
        for (size_t e = 0; e < m.elements.size(); ++e)
            for (int k = 0; k < m.elements[e].nodeCount(); ++k) {
                const auto& s = r.sets[0].elemStress[e][k];
                maxErr = std::max(maxErr, std::abs(s[0] - sig));
                maxErr = std::max(maxErr, std::abs(s[1]));
                maxErr = std::max(maxErr, std::abs(s[2]));
            }
        const char* names[3] = {"CST: σx = 100 chính xác", "Q4: σx = 100 chính xác", "Q4 KTT: σx = 100 chính xác"};
        near(names[variant], sig + maxErr, sig, 1e-7);
    }
}

static void testPlaneCantilever() {
    std::puts("Công-xôn ứng suất phẳng so với lý thuyết dầm (Timoshenko)");
    const double L = 10, H = 1, t = 1, P = 1, E = 1e6, nu = 0.25;
    const double I = t * H * H * H / 12, G = E / (2 * (1 + nu));
    const double exact = P * L * L * L / (3 * E * I) + P * L / (5.0 / 6.0 * G * H * t);
    const int nx = 20, ny = 4;
    for (int inc = 0; inc < 2; ++inc) {
        Model m = plateModel(ElemType::Quad4, nx, ny, L, H, false);
        m.materials[0].nu = nu;
        for (auto& e : m.elements) {
            e.incompatible = inc;
            e.t = t;
        }
        for (int j = 0; j <= ny; ++j) {
            fix(m, j * (nx + 1), true, j == ny / 2, false);
            // phân bố lực cắt parabol ở đầu tự do
            const double y = H * j / ny - H / 2;
            double w = 1.5 * P / H * (1 - 4 * y * y / (H * H)) * H / ny;
            if (j == 0 || j == ny) w = 0;
            m.cases[0].nodal.push_back({j * (nx + 1) + nx, {0, -w, 0}});
        }
        // chuẩn hóa tổng lực về P
        double sum = 0;
        for (auto& l : m.cases[0].nodal) sum += -l.F[1];
        for (auto& l : m.cases[0].nodal) l.F[1] *= P / sum;
        auto r = run(m);
        if (!r.ok) continue;
        const double tip = -r.sets[0].U[3 * ((ny / 2) * (nx + 1) + nx) + 1];
        if (inc) near("Q4 không tương thích: độ võng đầu (sai số < 3%)", tip, exact, 0.03);
        else std::printf("  (Q4 thường: %.4f / lý thuyết %.4f — Q4 thường cứng hơn khi uốn)\n", tip, exact);
    }
}

static void testIO() {
    std::puts("Ghi/đọc file .fem");
    for (const auto& ex : examples()) {
        Model m = ex.build();
        std::string text = modelToText(m), err;
        Model m2;
        if (!modelFromText(text, m2, err)) {
            ++g_fail;
            std::printf("  FAIL đọc lại %s: %s\n", ex.name, err.c_str());
            continue;
        }
        auto a = analyze(m), b = analyze(m2);
        double d = 0;
        for (size_t s = 0; s < a.sets.size() && s < b.sets.size(); ++s)
            for (size_t i = 0; i < a.sets[s].U.size(); ++i) d = std::max(d, std::abs(a.sets[s].U[i] - b.sets[s].U[i]));
        check(a.ok && b.ok && a.sets.size() == b.sets.size() && d < 1e-12, ex.name, d, 0);
        check(modelToText(m2) == text, (std::string(ex.name) + " (văn bản)").c_str(), 0, 0);
    }
}

static void testMechanism() {
    std::puts("Phát hiện hệ biến hình");
    Model m = beamModel(6, 2, 2e8, 0.01, 1e-4);
    m.elements[0].hingeJ = true;
    fix(m, 0, true, true, false);
    fix(m, 2, false, true, false);
    m.cases[0].nodal.push_back({1, {0, -10, 0}});
    auto r = analyze(m);
    check(!r.ok && r.error.find("biến hình") != std::string::npos, "dầm 3 khớp thẳng hàng là hệ biến hình", 0, 0);
}

static void testSubdivide() {
    std::puts("Chia nhỏ thanh giữ nguyên kết quả");
    Model m = beamModel(8, 1, 3e7, 0.18, 0.0054);
    fix(m, 0, true, true, true);
    fix(m, 1, false, true, false);
    m.cases[0].elem.push_back(dload(0, -5, -15));
    ElemLoad p;
    p.elem = 0;
    p.kind = ElemLoadKind::Point;
    p.q1 = -20;
    p.a = 5.5;
    m.cases[0].elem.push_back(p);
    auto r1 = run(m);
    subdivideBar(m, 0, 4);
    auto r2 = run(m);
    if (r1.ok && r2.ok) {
        near("Phản lực đầu ngàm không đổi", r2.sets[0].R[2], r1.sets[0].R[2], 1e-9);
        near("Độ võng tại x=4 khớp", r2.sets[0].U[3 * 3 + 1], at(r1.sets[0].bars[0].x, r1.sets[0].bars[0].v, 4), 1e-4);
    }
}

static void testInfluenceLines() {
    std::puts("Đường ảnh hưởng");
    const double L = 6;
    {
        Model m = beamModel(L, 2, 2e8, 0.01, 1e-4);
        fix(m, 0, true, true, false);
        fix(m, 2, false, true, false);
        IlTarget t;
        t.kind = IlKind::Reaction;
        t.node = 0;
        t.dof = UY;
        auto il = influenceLine(m, {0, 1}, t, 10);
        check(il.ok, "đah tính được", 0, 0);
        if (il.ok) {
            near("đah RA: tung độ tại A = 1", il.pts.front().value, 1);
            near("đah RA: tung độ tại B = 0", il.pts.back().value, 0, 0, 1e-12);
            near("đah RA: diện tích = L/2", il.areaPos, L / 2);
        }
        t.kind = IlKind::Force;
        t.elem = 0;
        t.x = L / 2; // cuối thanh 1 = giữa nhịp
        t.qty = BarQty::M;
        il = influenceLine(m, {0, 1}, t, 10);
        if (il.ok) {
            near("đah M giữa nhịp: max = L/4", il.vmax, L / 4);
            near("đah M giữa nhịp: diện tích = L²/8", il.areaPos, L * L / 8);
        }
        t.x = 1.0;
        t.qty = BarQty::Q;
        il = influenceLine(m, {0, 1}, t, 30);
        if (il.ok) {
            near("đah Q: max = (L − x)/L", il.vmax, (L - 1) / L);
            near("đah Q: min = −x/L", il.vmin, -1 / L, 1e-6, 1e-9);
            near("đah Q: tổng diện tích = (L − 2x)/2", il.areaPos + il.areaNeg, (L - 2) / 2, 1e-6);
        }
    }
    {
        // dầm liên tục 2 nhịp bằng nhau: đah M gối giữa, min = −l·(2/3)/(4√3) tại a = l/√3
        const double l = 6;
        Model m = beamModel(2 * l, 2, 2e8, 0.01, 1e-4);
        fix(m, 0, true, true, false);
        fix(m, 1, false, true, false);
        fix(m, 2, false, true, false);
        IlTarget t;
        t.kind = IlKind::Force;
        t.elem = 0;
        t.x = l;
        t.qty = BarQty::M;
        auto il = influenceLine(m, {0, 1}, t, 60);
        if (il.ok) {
            near("đah M gối giữa dầm liên tục: min", il.vmin, -l * (2.0 / 3.0) / (4 * std::sqrt(3.0)), 2e-3);
            near("đah M gối giữa: diện tích = −l²/8 (M_B = −ql²/8)", il.areaNeg, -l * l / 8, 1e-3);
        }
    }
}

static void testEnvelope() {
    std::puts("Biểu đồ bao nội lực");
    Model m;
    for (const auto& ex : examples())
        if (std::string(ex.name) == "Dầm liên tục 3 nhịp") m = ex.build();
    auto r = run(m);
    if (!r.ok) return;
    std::vector<int> sets;
    for (size_t i = 0; i < r.sets.size(); ++i) sets.push_back((int)i);
    auto env = envelope(m, r, sets);
    bool ok = true;
    for (size_t e = 0; e < m.elements.size(); ++e)
        for (size_t s = 0; s < r.sets.size(); ++s) {
            const auto& b = r.sets[s].bars[e];
            for (size_t k = 0; k < b.x.size(); ++k) {
                const double v = b.M[k];
                const double mx = *std::max_element(env[e].max[2].begin(), env[e].max[2].end());
                const double mn = *std::min_element(env[e].min[2].begin(), env[e].min[2].end());
                ok &= v <= mx + 1e-9 && v >= mn - 1e-9;
            }
        }
    check(ok, "bao chứa mọi giá trị M", 0, 0);
    // tại gối giữa thứ nhất: Mmin = M của tổ hợp (âm lớn nhất)
    const double mEnv = env[0].min[2].back();
    double mSets = 1e300;
    for (const auto& s : r.sets) mSets = std::min(mSets, s.bars[0].M.back());
    near("Mmin tại gối = min các tổ hợp", mEnv, mSets);
}

static void testEdgeLoad() {
    std::puts("Tải trên cạnh phần tử phẳng");
    Model m;
    for (const auto& ex : examples())
        if (std::string(ex.name) == "Đập trọng lực (biến dạng phẳng)") m = ex.build();
    // áp lực thủy tĩnh p = γ(H − y) trên mặt thượng lưu x = 0, so với ví dụ (lực nút lập tay)
    std::vector<int> face;
    for (size_t i = 0; i < m.nodes.size(); ++i)
        if (std::abs(m.nodes[i].x) < 1e-9) face.push_back((int)i);
    m.cases.push_back(LoadCase("cạnh"));
    EdgeLoadSpec spec;
    spec.pressure = true;
    spec.p0 = 10 * 11;
    spec.ky = -10;
    const int ne = applyEdgeLoad(m, (int)m.cases.size() - 1, face, spec);
    check(ne == 12, "số cạnh thượng lưu", ne, 12);
    double fx1 = 0, fx2 = 0;
    for (auto& l : m.cases[1].nodal) fx1 += l.F[0];
    for (auto& l : m.cases.back().nodal) fx2 += l.F[0];
    near("Tổng áp lực nước = γH²/2", fx2, 10 * 11 * 11 / 2.0, 1e-9);
    near("Trùng với lực nút của ví dụ", fx2, fx1, 1e-9);
}

static Model columnModel(double L, double E, double A, double I) {
    Model m = beamModel(L, 1, E, A, I);
    m.nodes[1] = {0, L}; // cột thẳng đứng
    m.cases[0].nodal.push_back({1, {0, -1, 0}});
    return m;
}

static void testEigenSolver() {
    std::puts("Bộ giải trị riêng đối xứng");
    const int n = 6;
    std::vector<double> A(n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) A[i * n + j] = A[j * n + i] = std::sin(1.3 * i + 0.7 * j) + (i == j ? 3.0 * i : 0.0);
    std::vector<double> A0 = A, d, V;
    symmetricEigen(n, A, d, V);
    double res = 0, orth = 0;
    for (int k = 0; k < n; ++k)
        for (int i = 0; i < n; ++i) {
            double av = 0;
            for (int j = 0; j < n; ++j) av += A0[i * n + j] * V[j * n + k];
            res = std::max(res, std::abs(av - d[k] * V[i * n + k]));
            for (int l = 0; l < n; ++l) {
                double dot = 0;
                for (int j = 0; j < n; ++j) dot += V[j * n + k] * V[j * n + l];
                orth = std::max(orth, std::abs(dot - (k == l ? 1.0 : 0.0)));
            }
        }
    check(res < 1e-10 && orth < 1e-10, "A·v = λ·v, VᵀV = I", res + orth, 0);
    bool sorted = true;
    for (int k = 1; k < n; ++k) sorted &= d[k] >= d[k - 1];
    check(sorted, "trị riêng tăng dần", 0, 0);
}

static void testBuckling() {
    std::puts("Ổn định: lực tới hạn Euler");
    const double L = 5, E = 2.1e8, A = 46.5e-4, I = 7080e-8, EI = E * I, pi = M_PI;
    EigenOptions o;
    o.kind = EigenKind::Buckling;
    o.modes = 3;
    {
        Model m = columnModel(L, E, A, I);
        fix(m, 0, true, true, false);
        fix(m, 1, true, false, false);
        auto r = eigenAnalysis(m, o);
        check(r.ok, r.error.c_str(), 0, 0);
        if (r.ok) {
            near("Hai đầu khớp: P = π²EI/L²", r.modes[0].lambda, pi * pi * EI / (L * L), 1e-4);
            near("Dạng 2: 4π²EI/L²", r.modes[1].lambda, 4 * pi * pi * EI / (L * L), 1e-3);
            near("μ = 1", effectiveLengthFactor(m, r, 0), 1.0, 1e-4);
        }
    }
    {
        Model m = columnModel(L, E, A, I);
        fix(m, 0, true, true, true);
        auto r = eigenAnalysis(m, o);
        if (r.ok) {
            near("Công-xôn: P = π²EI/4L²", r.modes[0].lambda, pi * pi * EI / (4 * L * L), 1e-4);
            near("Công-xôn: μ = 2", effectiveLengthFactor(m, r, 0), 2.0, 1e-4);
        }
    }
    {
        Model m = columnModel(L, E, A, I);
        fix(m, 0, true, true, true);
        fix(m, 1, true, false, false);
        auto r = eigenAnalysis(m, o);
        if (r.ok) near("Ngàm – khớp: P = 20,19·EI/L² (sai số rời rạc < 0,05%)", r.modes[0].lambda, 20.1907 * EI / (L * L), 5e-4);
    }
    {
        // gối ngàm cả hai đầu nhưng thanh có khớp hai đầu → tương đương hai đầu khớp
        Model m = columnModel(L, E, A, I);
        fix(m, 0, true, true, true);
        fix(m, 1, true, false, true);
        m.elements[0].hingeI = m.elements[0].hingeJ = true;
        auto r = eigenAnalysis(m, o);
        if (r.ok) near("Khớp đầu thanh: P = π²EI/L²", r.modes[0].lambda, pi * pi * EI / (L * L), 1e-4);
    }
    {
        Model m = columnModel(L, E, A, I);
        fix(m, 0, true, true, false);
        fix(m, 1, true, false, false);
        m.cases[0].nodal[0].F[1] = +1; // kéo
        auto r = eigenAnalysis(m, o);
        check(!r.ok, "thanh chịu kéo không mất ổn định", 0, 0);
    }
}

static void testVibration() {
    std::puts("Dao động riêng");
    const double L = 6, E = 3e7, b = 0.3, h = 0.5, A = b * h, I = b * h * h * h / 12, gamma = 25;
    const double rhoA = gamma / 9.81 * A, EI = E * I, pi = M_PI;
    EigenOptions o;
    o.kind = EigenKind::Vibration;
    o.modes = 4;
    o.subdiv = 12;
    {
        Model m = beamModel(L, 1, E, A, I);
        fix(m, 0, true, true, false);
        fix(m, 1, false, true, false);
        auto r = eigenAnalysis(m, o);
        check(r.ok, r.error.c_str(), 0, 0);
        if (r.ok) {
            const double w1 = pi * pi * std::sqrt(EI / (rhoA * std::pow(L, 4)));
            near("Dầm đơn giản: ω1 = π²√(EI/mL⁴)", r.modes[0].omega, w1, 1e-4);
            near("Dầm đơn giản: ω2 = 4ω1", r.modes[1].omega, 4 * w1, 1e-3);
            near("Chu kỳ T = 2π/ω", r.modes[0].period, 2 * pi / w1, 1e-4);
        }
    }
    {
        Model m = beamModel(L, 1, E, A, I);
        fix(m, 0, true, true, true);
        auto r = eigenAnalysis(m, o);
        if (r.ok) near("Công-xôn: ω1 = 1,8751²√(EI/mL⁴)", r.modes[0].omega,
                        1.875104 * 1.875104 * std::sqrt(EI / (rhoA * std::pow(L, 4))), 1e-4);
    }
    {
        // cột không trọng lượng mang khối lượng tập trung W/g ở đỉnh: ω = √(3EI/(mL³))
        const double H = 4, W = 200;
        Model m = columnModel(H, E, A, I);
        m.cases[0].nodal[0].F[1] = -W;
        fix(m, 0, true, true, true);
        o.selfMassFactor = 0;
        o.massSet = 0;
        auto r = eigenAnalysis(m, o);
        if (r.ok) near("Hệ 1 bậc tự do: ω = √(3EI/mH³)", r.modes[0].omega, std::sqrt(3 * EI / (W / 9.81 * H * H * H)), 1e-4);
    }
}

static void testComboGenerator() {
    std::puts("Sinh tổ hợp theo TCVN 2737:1995");
    Model m = beamModel(6, 1, 2e8, 0.01, 1e-4);
    fix(m, 0, true, true, false);
    fix(m, 1, false, true, false);
    m.cases = {LoadCase("TT"), LoadCase("HT"), LoadCase("Gió trái"), LoadCase("Gió phải")};
    m.cases[0].gammaF = 1.1;
    for (int i = 1; i < 4; ++i) m.cases[i].kind = 1;
    m.cases[1].gammaF = 1.2;
    m.cases[2].group = m.cases[3].group = 1;
    m.cases[1].nodal.push_back({1, {0, -10, 0}});
    m.combos.push_back(Combo("THủ công", {{0, 1.0}}));
    const int n = generateCombos(m);
    check(n == 5, "3 tổ hợp CB1 + 2 tổ hợp CB2 (gió trái/phải loại trừ)", n, 5);
    check(m.combos.size() == 6 && !m.combos[0].generated, "giữ tổ hợp thủ công", (double)m.combos.size(), 6);
    if (m.combos.size() == 6) {
        const Combo& c = m.combos[4]; // CB2 đầu tiên: TT + 0,9(HT + Gió trái)
        check(c.terms.size() == 3, "CB2 có 3 thành phần", (double)c.terms.size(), 3);
        near("γf·ψ2 của HT = 1,2·0,9", c.terms[1].second, 1.08);
        near("γf của TT = 1,1", c.terms[0].second, 1.1);
    }
    generateCombos(m);
    check(m.combos.size() == 6, "sinh lại không nhân đôi", (double)m.combos.size(), 6);
    std::string txt = modelToText(m), err;
    Model m2;
    check(modelFromText(txt, m2, err) && m2.combos.size() == 6 && m2.combos[1].generated && m2.cases[2].group == 1 &&
              m2.cases[1].gammaF == 1.2,
          "ghi/đọc phân loại tải và tổ hợp tự sinh", 0, 0);
}

static void testSplitMerge() {
    std::puts("Cắt thanh tại điểm, gộp nút trùng");
    // dầm 8 m, tải hình thang + lực tập trung: cắt tại x = 3 không đổi kết quả
    Model m = beamModel(8, 1, 3e7, 0.18, 0.0054);
    fix(m, 0, true, true, true);
    fix(m, 1, false, true, false);
    m.cases[0].elem.push_back(dload(0, -5, -15));
    ElemLoad p;
    p.elem = 0;
    p.kind = ElemLoadKind::Point;
    p.q1 = -20;
    p.a = 5.5;
    m.cases[0].elem.push_back(p);
    auto r1 = run(m);
    const int nn = splitBarAt(m, 0, 3.0 / 8);
    auto r2 = run(m);
    check(nn == 2 && m.elements.size() == 2, "cắt thành 2 thanh", (double)m.elements.size(), 2);
    if (r1.ok && r2.ok) {
        near("phản lực ngàm không đổi", r2.sets[0].R[2], r1.sets[0].R[2], 1e-9);
        near("M tại x = 3", r2.sets[0].bars[0].M.back(), barValueAt(r1.sets[0].bars[0], BarQty::M, 3.0), 1e-9);
    }
    // hai khung giống nhau đặt chồng: gộp nút → mô hình như một khung
    Model k;
    for (const auto& ex : examples())
        if (std::string(ex.name) == "Khung 1 tầng 1 nhịp") k = ex.build();
    auto rk = run(k);
    Model k2 = k;
    const int base = (int)k2.nodes.size();
    for (int i = 0; i < base; ++i) k2.nodes.push_back(k.nodes[i]);
    for (const auto& el : k.elements) {
        Element c = el;
        for (int q = 0; q < c.nodeCount(); ++q) c.n[q] += base;
        k2.elements.push_back(c);
    }
    const int merged = mergeNodes(k2, 1e-9);
    check(merged == base && k2.nodes.size() == k.nodes.size() && k2.elements.size() == k.elements.size(),
          "gộp nút và bỏ thanh trùng", merged, base);
    auto rk2 = run(k2);
    if (rk.ok && rk2.ok) near("kết quả sau khi gộp", rk2.sets[0].U[4], rk.sets[0].U[4], 1e-12, 1e-15);
}

static void testLackOfFit() {
    std::puts("Sai số chế tạo ΔL");
    const double L = 4, E = 2e8, A = 0.01;
    Model m = beamModel(L, 1, E, A, 1e-4);
    fix(m, 0, true, true, true);
    fix(m, 1, true, true, true);
    ElemLoad l;
    l.elem = 0;
    l.kind = ElemLoadKind::LackOfFit;
    l.q1 = 0.002;
    m.cases[0].elem.push_back(l);
    auto r = run(m);
    if (r.ok) near("ngàm 2 đầu: N = −EA·ΔL/L", r.sets[0].bars[0].N[3], -E * A * 0.002 / L);
    // tĩnh định: không sinh nội lực, đầu tự do dịch chuyển ΔL
    Model c = beamModel(L, 1, E, A, 1e-4);
    fix(c, 0, true, true, true);
    c.cases[0].elem.push_back(l);
    auto rc = run(c);
    if (rc.ok) {
        near("công-xôn: N = 0", rc.sets[0].bars[0].N[0], 0, 0, 1e-9);
        near("công-xôn: dịch chuyển đầu = ΔL", rc.sets[0].U[3], 0.002, 1e-9);
    }
    // chia thanh giữ tổng ΔL
    subdivideBar(m, 0, 4);
    auto r2 = run(m);
    if (r2.ok) near("chia 4 đoạn: N không đổi", r2.sets[0].bars[2].N[0], -E * A * 0.002 / L, 1e-9);
}

static void testPDelta() {
    std::puts("Tính theo sơ đồ biến dạng (P-Δ)");
    // cột công-xôn chịu P dọc trục và H ngang: δ = H/(k³EI)(tan kL − kL), M_chân = HL + Pδ
    const double L = 5, E = 2.1e8, I = 7080e-8, A = 46.5e-4, EI = E * I, P = 1000, H = 10;
    Model m = columnModel(L, E, A, I);
    m.cases[0].nodal[0].F = {H, -P, 0};
    fix(m, 0, true, true, true);
    subdivideBar(m, 0, 10);
    auto lin = run(m);
    m.secondOrder = true;
    auto r = run(m);
    if (!r.ok || !lin.ok) return;
    const double k = std::sqrt(P / EI), exact = H / (k * k * k * EI) * (std::tan(k * L) - k * L);
    const double tip = r.sets[0].U[3 * 1];
    near("độ võng đầu cột bậc hai", tip, exact, 2e-3);
    check(tip > lin.sets[0].U[3] * 1.05, "lớn hơn tuyến tính", tip, lin.sets[0].U[3]);
    near("mômen chân cột = HL + Pδ", std::abs(r.sets[0].R[2]), H * L + P * tip, 2e-3);
    near("mômen chân cột trên biểu đồ", std::abs(r.sets[0].bars[0].M.front()), H * L + P * tip, 2e-3);
    // vượt tải trọng tới hạn → báo lỗi
    m.cases[0].nodal[0].F[1] = -1.2 * M_PI * M_PI * EI / (4 * L * L);
    auto rf = analyze(m);
    check(!rf.ok, "vượt Pcr: báo mất ổn định", 0, 0);
}

static void testFoundations() {
    std::puts("Nền Pasternak và nền chỉ chịu nén");
    const double L = 12, E = 2.75e7, b = 1, h = 0.6, I = b * h * h * h / 12, P = 300;
    const int n = 48;
    Model m = beamModel(L, n, E, b * h, I, h);
    for (auto& e : m.elements) e.kf = 20000;
    fix(m, n / 2, true, false, false);
    m.cases[0].nodal.push_back({n / 2, {0, -P, 0}});
    auto w = run(m);
    for (auto& e : m.elements) e.kG = 5000;
    auto pa = run(m);
    if (w.ok && pa.ok) {
        near("Pasternak: tổng phản lực nền = P", pa.sets[0].foundRy, P, 1e-6);
        check(std::abs(pa.sets[0].U[3 * (n / 2) + 1]) < std::abs(w.sets[0].U[3 * (n / 2) + 1]),
              "lớp trượt làm giảm độ lún", pa.sets[0].U[3 * (n / 2) + 1], w.sets[0].U[3 * (n / 2) + 1]);
        double maxErr = 0;
        for (int e = 0; e < n; ++e) {
            const auto& br = pa.sets[0].bars[e];
            const auto& f = br.endForces;
            maxErr = std::max(maxErr, std::abs(br.M.back() - f[5]));
            if (e + 1 < n && e + 1 != n / 2) maxErr = std::max(maxErr, std::abs(br.Q.back() - pa.sets[0].bars[e + 1].Q.front()));
        }
        check(maxErr < 1e-6 * P * L, "Pasternak: M khớp đầu thanh, Q liên tục giữa các đoạn", maxErr, 0);
    }
    // nền chỉ chịu nén: lực ở đầu dầm làm đầu kia vểnh lên
    Model c = beamModel(L, n, E, b * h, I, h);
    for (auto& e : c.elements) {
        e.kf = 20000;
        e.kfCompOnly = true;
    }
    fix(c, 0, true, false, false);
    c.cases[0].nodal.push_back({0, {0, -P, 0}});
    c.cases[0].nodal.push_back({0, {0, 0, -400}});
    auto rc = run(c);
    if (rc.ok) {
        near("chỉ chịu nén: tổng phản lực nền = P", rc.sets[0].foundRy, P, 1e-6);
        bool tension = false;
        for (int e = 0; e < n; ++e) {
            double fx, fy;
            barFoundationForce(barProps(c, e), rc.sets[0].bars[e], fx, fy);
            if (fy < -1e-9) tension = true;
        }
        bool lifted = false;
        for (const auto& wmsg : rc.warnings) lifted |= wmsg.find("tách khỏi nền") != std::string::npos;
        check(lifted, "có đoạn dầm tách khỏi nền", 0, 0);
        (void)tension;
    }
}

static void testKinematics() {
    std::puts("Phân tích cấu tạo hình học");
    auto ex = [](const char* name) {
        for (const auto& e : examples())
            if (std::string(e.name) == name) return e.build();
        return Model{};
    };
    struct C {
        const char* name;
        int n;
    } cases[] = {{"Dầm đơn giản", 0},      {"Dầm tĩnh định nhiều nhịp", 0}, {"Khung ba khớp", 0},
                 {"Dàn cầu", 0},           {"Khung 1 tầng 1 nhịp", 3},      {"Dầm liên tục 3 nhịp", 2},
                 {"Khung 2 tầng 2 nhịp", 12}};
    for (const auto& cs : cases) {
        const auto k = kinematicAnalysis(ex(cs.name));
        check(k.applicable && k.mechanisms == 0 && k.indeterminacy == cs.n && k.W == -cs.n,
              (std::string(cs.name) + ": bậc siêu tĩnh").c_str(), k.indeterminacy, cs.n);
    }
    // dầm hai gối có khớp giữa nhịp: biến hình, W = 1
    Model m = beamModel(6, 2, 2e8, 0.01, 1e-4);
    m.elements[0].hingeJ = true;
    fix(m, 0, true, true, false);
    fix(m, 2, false, true, false);
    auto k = kinematicAnalysis(m);
    check(k.W == 1 && k.mechanisms == 1 && !k.mechanismU.empty(), "dầm có khớp: biến hình, W = 1", k.W, 1);
    // ba gối di động song song: W = 0 nhưng biến hình tức thời
    Model t = beamModel(6, 2, 2e8, 0.01, 1e-4);
    for (int i = 0; i < 3; ++i) fix(t, i, false, true, false);
    auto kt = kinematicAnalysis(t);
    check(kt.W == 0 && kt.mechanisms == 1, "3 gối di động song song: tức thời biến hình", kt.mechanisms, 1);
}

static void testSectionGeom() {
    std::puts("Đặc trưng hình học tiết diện");
    {
        const double b = 0.2, h = 0.4;
        SectionShape sh;
        sh.kind = ShapeKind::Rect;
        sh.p = {b, h, 0, 0, 0, 0};
        const auto c = shapeContours(sh);
        const auto p = sectionProps(c);
        near("CN: A", p.A, b * h);
        near("CN: Jx = bh³/12", p.Jx, b * h * h * h / 12);
        near("CN: Jy = hb³/12", p.Jy, h * b * b * b / 12);
        near("CN: W = bh²/6", p.Wtop, b * h * h / 6);
        double ky = 0, kx = 0;
        for (const auto& v : p.kern) {
            ky = std::max(ky, v[1] - p.yc);
            kx = std::max(kx, v[0] - p.xc);
        }
        near("CN: lõi tiết diện h/6", ky, h / 6, 1e-9);
        near("CN: lõi tiết diện b/6", kx, b / 6, 1e-9);
        // τmax = 1,5 Q/A tại trọng tâm
        const double S = staticMomentAbove(c, p.yc, p.yc), w = widthAt(c, p.yc);
        near("CN: τmax = 1,5Q/A", 1.0 * S / (p.Jx * w), 1.5 / (b * h), 1e-9);
    }
    {
        // chữ T: cánh 200×20, sườn 20×180 (mm)
        const double bf = 0.2, tf = 0.02, bw = 0.02, h = 0.2;
        SectionShape sh;
        sh.kind = ShapeKind::Tee;
        sh.p = {bf, tf, bw, h, 0, 0};
        const auto c = shapeContours(sh);
        const auto p = sectionProps(c);
        const double A1 = bf * tf, A2 = bw * (h - tf), y1 = h - tf / 2, y2 = (h - tf) / 2, A = A1 + A2;
        const double yc = (A1 * y1 + A2 * y2) / A;
        const double J = bf * tf * tf * tf / 12 + A1 * (y1 - yc) * (y1 - yc) + bw * std::pow(h - tf, 3) / 12 +
                         A2 * (y2 - yc) * (y2 - yc);
        near("T: A", p.A, A);
        near("T: trọng tâm", p.yc, yc);
        near("T: Jx (định lý trục song song)", p.Jx, J);
        near("T: y thớ trên", p.yTop, h - yc);
        near("T: Jxy = 0 (đối xứng)", p.Jxy, 0, 0, 1e-15);
        const Section sec = sectionFromShape("T", sh);
        near("T: tiết diện mô hình yBot", sec.yBot, yc);
    }
    {
        SectionShape sh;
        sh.kind = ShapeKind::Pipe;
        sh.p = {0.3, 0.01, 0, 0, 0, 0};
        const auto p = sectionProps(shapeContours(sh));
        const double D = 0.3, d = 0.28;
        near("Ống: A (đa giác 180 cạnh)", p.A, M_PI * (D * D - d * d) / 4, 1e-3);
        near("Ống: J", p.Jx, M_PI * (std::pow(D, 4) - std::pow(d, 4)) / 64, 1e-3);
        SectionShape cs;
        cs.kind = ShapeKind::Circle;
        cs.p = {0.4, 0, 0, 0, 0, 0};
        const auto pc = sectionProps(shapeContours(cs));
        double r = 0;
        for (const auto& v : pc.kern) r = std::max(r, std::hypot(v[0] - pc.xc, v[1] - pc.yc));
        near("Tròn: lõi tiết diện bán kính D/8", r, 0.4 / 8, 2e-3);
    }
    {
        // thép góc đều cạnh: trục chính nghiêng 45°
        SectionShape sh;
        sh.kind = ShapeKind::Angle;
        sh.p = {0.1, 0.1, 0.01, 0, 0, 0};
        const auto p = sectionProps(shapeContours(sh));
        near("L đều cạnh: |α| = 45°", std::abs(p.alpha), 45, 1e-6);
        check(p.Jxy < 0, "L: Jxy < 0 (phần lớn diện tích ở góc phần tư II, IV quanh trọng tâm)", p.Jxy, 0);
        near("L: J1 + J2 = Jx + Jy", p.J1 + p.J2, p.Jx + p.Jy, 1e-12);
    }
    {
        // hộp = chữ nhật ngoài − chữ nhật trong
        SectionShape sh;
        sh.kind = ShapeKind::Box;
        sh.p = {0.3, 0.5, 0.02, 0, 0, 0};
        const auto c = shapeContours(sh);
        const auto p = sectionProps(c);
        near("Hộp: Jx", p.Jx, (0.3 * std::pow(0.5, 3) - 0.26 * std::pow(0.46, 3)) / 12);
        near("Hộp: bề rộng tại thành = 2t", widthAt(c, 0.25), 0.04, 1e-12);
    }
    {
        // ghi/đọc file giữ hình dạng
        Model m;
        SectionShape sh;
        sh.kind = ShapeKind::Tee;
        sh.p = {0.2, 0.02, 0.02, 0.2, 0, 0};
        m.sections.push_back(sectionFromShape("T", sh));
        SectionShape poly;
        poly.kind = ShapeKind::Polygon;
        poly.contours = {{{0, 0}, {0.1, 0}, {0.05, 0.1}}};
        m.sections.push_back(sectionFromShape("tam giác", poly));
        std::string err;
        Model r;
        const bool ok = modelFromText(modelToText(m), r, err);
        check(ok && r.sections.size() == 2 && r.sections[0].shape.kind == ShapeKind::Tee &&
                  std::abs(r.sections[0].yTop - m.sections[0].yTop) < 1e-12 && r.sections[1].shape.contours.size() == 1 &&
                  r.sections[1].shape.contours[0].size() == 3,
              "ghi/đọc tiết diện có hình dạng", 0, 0);
        near("tam giác: J = bh³/36", sectionProps(r.sections[1].shape.contours).Jx, 0.1 * 1e-3 / 36, 1e-9);
    }
}

static void testForceMethod() {
    std::puts("Phương pháp lực, Maxwell–Mohr, phương pháp chuyển vị dạng ma trận");
    auto ex = [](const char* name) {
        for (const auto& e : examples())
            if (std::string(e.name) == name) return e.build();
        return Model{};
    };
    {
        const double L = 8, q = 12;
        Model m = beamModel(L, 1, 3e7, 0.18, 0.0054);
        fix(m, 0, true, true, true);
        fix(m, 1, true, true, true);
        m.cases[0].elem.push_back(dload(0, -q, -q));
        const auto f = forceMethod(m, 0);
        check(f.ok && f.n == 3, "dầm 2 đầu ngàm: n = 3", f.n, 3);
        if (f.ok) {
            near("M đầu = −qL²/12", f.final.bars[0].M.front(), -q * L * L / 12, 1e-6);
            near("M giữa = qL²/24", barValueAt(f.final.bars[0], BarQty::M, L / 2), q * L * L / 24, 1e-6);
        }
    }
    struct C {
        const char* name;
        int set, n;
    } cases[] = {{"Khung 1 tầng 1 nhịp", 0, 3}, {"Khung 1 tầng 1 nhịp", 1, 3}, {"Khung 2 tầng 2 nhịp", 0, 12},
                 {"Dầm liên tục 3 nhịp", 6, 2}, {"Dầm tĩnh định nhiều nhịp", 0, 0}, {"Khung ba khớp", 0, 0}, {"Dàn cầu", 0, 0}};
    for (const auto& c : cases) {
        const Model m = ex(c.name);
        const auto f = forceMethod(m, c.set);
        double mmax = 1e-9;
        if (f.ok)
            for (const auto& b : f.final.bars)
                for (double v : b.M) mmax = std::max(mmax, std::abs(v));
        check(f.ok && f.n == c.n && f.maxDiffM < 1e-4 * mmax + 1e-6,
              (std::string(c.name) + ": phương pháp lực khớp PTHH").c_str(), f.ok ? f.maxDiffM : -1, 0);
    }
    {
        // dàn siêu tĩnh: thêm một thanh xiên
        Model m = ex("Dàn cầu");
        Element d = m.elements[0];
        d.n = {0, 7, -1, -1}; // nút 1 (0,0) – nút 8 (3,3)? nối thêm thanh chéo khoang đầu
        d.n = {1, 6, -1, -1};
        m.elements.push_back(d);
        const auto f = forceMethod(m, 0);
        const auto fe = analyze(m);
        check(f.ok && f.n == 1, "dàn thêm 1 thanh: n = 1", f.n, 1);
        if (f.ok && fe.ok) {
            double d0 = 0;
            for (size_t e = 0; e < m.elements.size(); ++e)
                d0 = std::max(d0, std::abs(f.final.bars[e].N[0] - fe.sets[0].bars[e].N[0]));
            check(d0 < 1e-4, "dàn siêu tĩnh: lực dọc khớp PTHH", d0, 0);
        }
    }
    {
        const double L = 6, q = 10, E = 2.1e8, I = 7080e-8;
        Model m = beamModel(L, 2, E, 46.5e-4, I);
        fix(m, 0, true, true, false);
        fix(m, 2, false, true, false);
        m.cases[0].elem.push_back(dload(0, -q, -q));
        m.cases[0].elem.push_back(dload(1, -q, -q));
        const auto f = forceMethod(m, 0);
        const auto mo = mohrDisplacement(m, f, 1, UY);
        check(mo.ok, "Maxwell–Mohr tính được", 0, 0);
        if (mo.ok) near("Maxwell–Mohr: độ võng giữa nhịp 5qL⁴/384EI", mo.value, -5 * q * std::pow(L, 4) / (384 * E * I), 2e-4);
        const Model pf = ex("Khung 1 tầng 1 nhịp");
        const auto fp = forceMethod(pf, 1);
        const auto mp = mohrDisplacement(pf, fp, 1, UX);
        if (mp.ok) near("Maxwell–Mohr: chuyển vị ngang đỉnh khung = PTHH", mp.value, mp.fem, 2e-4);
        const auto mm = matrixMethod(pf, 0);
        const auto fe = analyze(pf);
        check(mm.ok && fe.ok, "phương pháp chuyển vị dạng ma trận", 0, 0);
        if (mm.ok && fe.ok) {
            double d = 0;
            int k = 0;
            for (int i = 0; i < (int)fe.dofState.size(); ++i)
                if (fe.dofState[i] == DofFree) d = std::max(d, std::abs(mm.Z[k++] - fe.sets[0].U[i]));
            check(d < 1e-12, "Z = K⁻¹P khớp PTHH", d, 0);
        }
    }
}

static void testArchTaperTrain() {
    std::puts("Vòm, thanh vút, xếp đoàn tải");
    {
        // vòm ba khớp parabol, tải đều theo hình chiếu: trục hợp lý → M = 0 tại các nút, H = qL²/8f
        Model m;
        m.materials.push_back(Material{});
        m.sections.push_back(Section{});
        const double L = 24, f = 6, q = 10;
        const auto nodes = makeArch(m, 0, 0, L, f, 24, 0, 1, 0, 0);
        m.cases.push_back(LoadCase("q"));
        for (size_t e = 0; e < m.elements.size(); ++e) {
            ElemLoad l = dload((int)e, -q, -q);
            l.projected = true;
            m.cases[0].elem.push_back(l);
        }
        auto r = run(m);
        if (r.ok) {
            near("vòm 3 khớp: H = qL²/8f", std::abs(r.sets[0].R[0]), q * L * L / (8 * f), 1e-9);
            double mn = 0;
            for (const auto& b : r.sets[0].bars) mn = std::max({mn, std::abs(b.M.front()), std::abs(b.M.back())});
            check(mn < 1e-8 * q * L * L, "vòm parabol: M tại nút = 0 (trục hợp lý)", mn, 0);
            check(kinematicAnalysis(m).indeterminacy == 0, "vòm 3 khớp tĩnh định", 0, 0);
        }
        Model m2;
        m2.materials.push_back(Material{});
        m2.sections.push_back(Section{});
        makeArch(m2, 0, 0, L, f, 24, 1, 2, 0, 0);
        check(kinematicAnalysis(m2).indeterminacy == 3, "vòm tròn 2 ngàm: siêu tĩnh bậc 3", 0, 0);
    }
    {
        // công-xôn vút: tiết diện chữ nhật 30×60 → 30×30
        Model m = beamModel(4, 1, 3e7, 0.18, 0.0054);
        m.sections[0] = rectSection(0.3, 0.6);
        m.sections.push_back(rectSection(0.3, 0.3));
        fix(m, 0, true, true, true);
        m.cases[0].nodal.push_back({1, {0, -10, 0}});
        const int n = taperBar(m, 0, 1, 8);
        check(n == 8 && m.elements.size() == 8 && m.sections.size() == 10, "chia vút 8 đoạn", (double)m.elements.size(), 8);
        bool mono = true;
        for (int k = 0; k + 1 < 8; ++k) mono &= m.sections[m.elements[k].sec].I > m.sections[m.elements[k + 1].sec].I;
        check(mono, "I giảm dần từ ngàm ra đầu tự do", 0, 0);
        auto r = run(m);
        if (r.ok) near("M ngàm = PL (tĩnh định)", r.sets[0].bars[0].M.front(), -40, 1e-9);
    }
    {
        // dầm đơn giản L = 20: đah M giữa nhịp (tam giác, đỉnh L/4)
        const double L = 20;
        Model m = beamModel(L, 2, 2e8, 0.01, 1e-4);
        fix(m, 0, true, true, false);
        fix(m, 2, false, true, false);
        IlTarget t;
        t.kind = IlKind::Force;
        t.elem = 0;
        t.x = L / 2;
        t.qty = BarQty::M;
        const auto il = influenceLine(m, {0, 1}, t, 20);
        LoadTrain one{"P", {{100, 0}}, 0};
        near("1 trục P: Smax = PL/4", moveTrain(il, one).Smax, 100 * L / 4, 1e-9);
        LoadTrain two{"2P", {{100, 0}, {100, 4}}, 0};
        near("2 trục cách 4 m: Smax = P(L/4 + (L/2 − 4)/2)", moveTrain(il, two).Smax, 100 * (L / 4 + (L / 2 - 4) / 2), 1e-9);
        LoadTrain lane{"làn", {{0, 0}}, 10};
        near("tải làn q: S = qΩ = qL²/8", moveTrain(il, lane).Smax, 10 * L * L / 8, 1e-9);
    }
}

static void testDynamics() {
    std::puts("Động lực học: phổ phản ứng, cưỡng bức điều hòa, Newmark");
    {
        // phổ thiết kế nền C: liên tục tại TB, TC, TD
        SpectrumParams sp;
        sp.agR = 0.1;
        sp.q = 1.5;
        sp.soil = 2;
        const double ag = 0.1 * 9.81, S = 1.15;
        near("Sd(TB) = ag·S·2,5/q", spectrumSd(sp, 0.2), ag * S * 2.5 / 1.5, 1e-12);
        near("Sd(0) = ag·S·2/3", spectrumSd(sp, 0), ag * S * 2.0 / 3.0, 1e-12);
        near("Sd(TC) liên tục", spectrumSd(sp, 0.6), ag * S * 2.5 / 1.5, 1e-12);
        near("Sd(1,2 s) = ag·S·2,5/q·TC/T", spectrumSd(sp, 1.2), ag * S * 2.5 / 1.5 * 0.6 / 1.2, 1e-12);
        near("Sd(3 s) ≥ β·ag", spectrumSd(sp, 3.0), std::max(ag * S * 2.5 / 1.5 * 0.6 * 2.0 / 9.0, 0.2 * ag), 1e-12);
    }
    // cột không trọng lượng mang khối lượng W/g ở đỉnh (hệ 1 bậc tự do theo phương ngang)
    const double H = 4, E = 3e7, b = 0.3, h = 0.3, W = 200, EI = E * b * h * h * h / 12;
    Model m = columnModel(H, E, b * h, b * h * h * h / 12);
    m.cases[0].nodal[0].F = {0, -W, 0};
    m.cases.push_back(LoadCase("Ngang"));
    m.cases[1].nodal.push_back({1, {10, 0, 0}});
    fix(m, 0, true, true, true);
    EigenOptions eo;
    eo.selfMassFactor = 0;
    eo.massSet = 0;
    eo.modes = 4;
    eo.subdiv = 4;
    const double mass = W / 9.81, w = std::sqrt(3 * EI / (mass * H * H * H)), T = 2 * M_PI / w;
    {
        SpectrumParams sp;
        sp.agR = 0.15;
        sp.q = 1.5;
        sp.soil = 3;
        const auto r = spectrumAnalysis(m, eo, sp);
        check(r.ok, r.error.c_str(), 0, 0);
        if (r.ok) {
            near("chu kỳ T1", r.modes[0].T, T, 1e-4);
            near("khối lượng hữu hiệu dạng 1 = 100%", r.modes[0].ratio, 1.0, 1e-6);
            near("lực cắt đáy = m·Sd(T1)", r.baseShear, mass * spectrumSd(sp, T), 1e-4);
            near("chuyển vị đỉnh = Sd/ω²", r.combined.U[3], spectrumSd(sp, T) / (w * w), 1e-4);
            near("mômen chân cột = V·H", r.combined.bars[0].M.front(), mass * spectrumSd(sp, T) * H, 1e-4);
        }
    }
    {
        const auto r = harmonicAnalysis(m, eo, 1, 0.5 * w);
        check(r.ok, r.error.c_str(), 0, 0);
        if (r.ok) near("hệ số động μ = 1/(1 − (θ/ω)²) = 4/3", r.dynFactor, 4.0 / 3.0, 1e-4);
    }
    {
        TimeHistoryParams p;
        p.loadSet = 1;
        p.fn = TimeFn::Step;
        p.zeta = 0;
        p.dt = T / 200;
        p.tEnd = 1.2 * T;
        const auto r = timeHistory(m, eo, p);
        check(r.ok, r.error.c_str(), 0, 0);
        if (r.ok) {
            const double ust = 10 * H * H * H / (3 * EI);
            const auto h1 = historyOf(r, 1, UX);
            double mx = 0;
            for (double v : h1) mx = std::max(mx, v);
            near("tải bậc không cản: umax = 2·u_tĩnh", mx, 2 * ust, 2e-3);
            // chu kỳ: lần trở về 0 đầu tiên sau đỉnh ≈ T
            size_t k = h1.size() / 2;
            while (k + 1 < h1.size() && h1[k + 1] > 1e-3 * ust) ++k;
            near("chu kỳ dao động tự do ≈ T", r.t[k], T, 2e-2);
            const auto snap = timeSnapshot(m, r, (int)(h1.size() / 2));
            check(snap.U.size() == 6, "nội lực tại thời điểm t", 0, 0);
        }
    }
}

static void testDesign() {
    std::puts("Kiểm tra thiết kế: cấu kiện, thép TCVN 5575, BTCT TCVN 5574, độ võng, cọc");
    auto ex = [](const char* name) {
        for (const auto& e : examples())
            if (std::string(e.name) == name) return e.build();
        return Model{};
    };
    check(findMembers(ex("Khung 1 tầng 1 nhịp")).size() == 3, "khung cổng: 3 cấu kiện", 0, 0);
    check(findMembers(ex("Dầm liên tục 3 nhịp")).size() == 3, "dầm liên tục: 3 nhịp", 0, 0);
    check(findMembers(ex("Khung 2 tầng 2 nhịp")).size() == 10, "khung 2 tầng 2 nhịp: 10 cấu kiện", 0, 0);
    {
        Model m = beamModel(6, 4, 3e7, 0.18, 0.0054);
        const auto mb = findMembers(m);
        check(mb.size() == 1 && std::abs(mb[0].L - 6) < 1e-12, "4 phần tử thẳng hàng không gối giữa → 1 cấu kiện", (double)mb.size(), 1);
    }
    {
        const double fE = 230.0 / 2.1e5;
        auto p1 = [&](double lb) { return 1 - (0.073 - 5.53 * fE) * lb * std::sqrt(lb); };
        near("φ(λ̄ = 2)", steelPhi(2.0, fE), p1(2.0), 1e-12);
        check(std::abs(steelPhi(2.5 - 1e-9, fE) - steelPhi(2.5 + 1e-9, fE)) < 0.01, "φ liên tục tại λ̄ = 2,5", 0, 0);
        check(std::abs(steelPhi(4.5 - 1e-9, fE) - steelPhi(4.5 + 1e-9, fE)) < 0.01, "φ liên tục tại λ̄ = 4,5", 0, 0);
        check(steelPhi(1, fE) > steelPhi(3, fE) && steelPhi(3, fE) > steelPhi(6, fE), "φ giảm khi λ̄ tăng", 0, 0);
    }
    {
        // cột thép I20 hai đầu khớp, nén đúng tâm 300 kN, CCT38
        Model m = columnModel(4, 2.1e8, 26.8e-4, 1840e-8);
        m.materials[0] = toMaterial(materialPresets()[1]);
        m.cases[0].nodal[0].F = {0, -300, 0};
        fix(m, 0, true, true, false);
        fix(m, 1, true, false, false);
        const auto r = analyze(m);
        const auto mb = findMembers(m);
        SteelParams sp;
        const auto c = checkSteel(m, r, {0}, mb, sp);
        check(c.size() == 1, "kiểm tra 1 cấu kiện thép", 0, 0);
        if (c.size() == 1) {
            const double f = 230e3, i = std::sqrt(1840e-8 / 26.8e-4), lb = 4 / i * std::sqrt(f / 2.1e8);
            const double phi = steelPhi(lb, f / 2.1e8);
            near("hệ số sử dụng ổn định = N/(φ·A·f)", c[0].stability, 300 / (phi * 26.8e-4 * f), 1e-9);
            near("độ mảnh λ = l0/i", c[0].lambda, 4 / i, 1e-12);
        }
    }
    {
        // dầm BTCT 25×50, B20, CB400: M = 120 kNm
        const Material c = toMaterial(materialPresets()[4]);
        RcParams p;
        const auto d = rcBeamFlexure(120, 0.25, 0.5, c, p);
        const double h0 = 0.46, am = 120 / (11.5e3 * 0.25 * h0 * h0), xi = 1 - std::sqrt(1 - 2 * am);
        near("dầm: αm", d.alphaM, am, 1e-12);
        near("dầm: As = Rb·b·ξ·h0/Rs", d.As, 11.5e3 * 0.25 * xi * h0 / 350e3 * 1e4, 1e-9);
        near("ξR = 0,8/(1 + Rs/700)", rcXiR(350, 2e5), 0.8 / (1 + 350.0 / 700.0), 1e-12);
        const auto big = rcBeamFlexure(400, 0.25, 0.5, c, p);
        check(big.overReinforced && big.AsComp > 0, "M lớn: cần cốt chịu nén", big.AsComp, 0);
        // cột 30×40, N = 1200 kN, M = 150 kNm, l0 = 3,6 m: kiểm tra lại hai phương trình cân bằng
        const auto col = rcColumnSym(1200, 150, 0.3, 0.4, 3.6, c, p);
        check(col.As > 0 && col.eta >= 1, "cột: As > 0, η ≥ 1", col.As, 0);
        const double Rb = 11.5e3, Rs = 350e3, b = 0.3, hc0 = 0.36, a = 0.04, xiR = rcXiR(350, 2e5);
        const double e = col.eta * col.e0 + (hc0 - a) / 2, As = col.As * 1e-4, x = col.x;
        const double sigma = x / hc0 <= xiR ? Rs : std::clamp(Rs * (2 * (1 - x / hc0) / (1 - xiR) - 1), -Rs, Rs);
        near("cột: ΣM = N·e", Rb * b * x * (hc0 - x / 2) + Rs * As * (hc0 - a), 1200 * e, 1e-6);
        near("cột: ΣN = N", Rb * b * x + Rs * As - sigma * As, 1200, 1e-6);
        const auto sh = rcShear(150, 0.25, 0.5, c, p);
        check(sh.ok && sh.spacing > 0 && sh.spacing <= 0.23, "cốt đai: bước đai cần thiết", sh.spacing, 0);
    }
    {
        const double L = 6, q = 10, E = 2.1e8, I = 7080e-8;
        Model m = beamModel(L, 4, E, 46.5e-4, I);
        fix(m, 0, true, true, false);
        fix(m, 4, false, true, false);
        for (int e = 0; e < 4; ++e) m.cases[0].elem.push_back(dload(e, -q, -q));
        const auto r = analyze(m);
        const auto d = checkDeflection(m, r.sets[0], findMembers(m), 250);
        check(d.size() == 1, "độ võng: 1 nhịp", 0, 0);
        if (d.size() == 1) {
            near("f = 5qL⁴/384EI", d[0].f, 5 * q * std::pow(L, 4) / (384 * E * I), 1e-6);
            near("vị trí f max giữa nhịp", d[0].xAt, L / 2, 1e-9);
        }
    }
    {
        Model m;
        PileParams pp;
        makePile(m, pp);
        const auto r = analyze(m);
        check(r.ok, r.error.c_str(), 0, 0);
        if (r.ok) {
            near("cọc: phản lực nền ngang = −H", r.sets[0].foundRx, -pp.H, 1e-6);
            check(r.sets[0].U[0] > 0, "đầu cọc chuyển vị theo chiều H", r.sets[0].U[0], 0);
        }
        near("bề rộng quy ước bp = 1,5d + 0,5 (d < 0,8)", pileConventionalWidth(0.6), 1.4, 1e-12);
    }
}

static void testDxfReport() {
    std::puts("DXF và báo cáo HTML");
    auto ex = [](const char* name) {
        for (const auto& e : examples())
            if (std::string(e.name) == name) return e.build();
        return Model{};
    };
    {
        const Model pf = ex("Khung 2 tầng 2 nhịp");
        const auto r = analyze(pf);
        DxfExportOptions o;
        const std::string dxf = exportDxf(pf, &r.sets[0], o);
        check(dxf.find("ENTITIES") != std::string::npos && dxf.find("FEM_BIEUDO") != std::string::npos, "xuất DXF có biểu đồ", 0, 0);
        Model m;
        DxfImportOptions io;
        io.layer = "FEM_THANH";
        const auto ir = importDxf(dxf, m, io);
        check(ir.ok && m.elements.size() == pf.elements.size() && m.nodes.size() == pf.nodes.size(),
              "xuất rồi nhập lại lớp FEM_THANH: giữ số thanh, số nút", (double)m.elements.size(), (double)pf.elements.size());
    }
    {
        // hai đoạn thẳng chữ T (mm) và một polyline kín
        const std::string dxf =
            "0\nSECTION\n2\nENTITIES\n"
            "0\nLINE\n8\nKC\n10\n0\n20\n0\n11\n6000\n21\n0\n"
            "0\nLINE\n8\nKC\n10\n3000\n20\n0\n11\n3000\n21\n3000\n"
            "0\nLWPOLYLINE\n8\nKHAC\n90\n4\n70\n1\n10\n0\n20\n5000\n10\n1000\n20\n5000\n10\n1000\n20\n6000\n10\n0\n20\n6000\n"
            "0\nENDSEC\n0\nEOF\n";
        Model m;
        DxfImportOptions io;
        io.unit = 0.001;
        io.layer = "KC";
        auto ir = importDxf(dxf, m, io);
        check(ir.ok && m.elements.size() == 3 && m.nodes.size() == 4 && ir.splits == 1, "chữ T: dầm tự cắt tại chân cột (3 thanh, 4 nút)",
              (double)m.elements.size(), 3);
        check(std::abs(m.nodes[1].x - 6.0) < 1e-12, "đổi đơn vị mm → m", m.nodes[1].x, 6);
        Model m2;
        io.layer = "KHAC";
        ir = importDxf(dxf, m2, io);
        check(ir.ok && m2.elements.size() == 4, "polyline kín: 4 thanh", (double)m2.elements.size(), 4);
        check(dxfLayers(dxf).size() == 2, "liệt kê lớp", (double)dxfLayers(dxf).size(), 2);
    }
    {
        const Model m = ex("Khung 1 tầng 1 nhịp");
        const auto r = analyze(m);
        const std::string h = htmlReport(m, r);
        check(h.find("<svg") != std::string::npos && h.find("Biểu đồ mômen") != std::string::npos && h.size() > 5000,
              "báo cáo HTML có hình biểu đồ", (double)h.size(), 0);
        const Model d = ex("Dầm-tường (ứng suất phẳng)");
        const std::string h2 = htmlReport(d, analyze(d));
        check(h2.find("Ứng suất tại nút") != std::string::npos, "báo cáo HTML phần tử phẳng", 0, 0);
    }
}

// Nút nằm trên biên đa giác lồi cx, cy (4 góc)
static bool onPatchBoundary(const Node& p, const double cx[4], const double cy[4]) {
    for (int k = 0; k < 4; ++k) {
        const double ax = cx[k], ay = cy[k], bx = cx[(k + 1) % 4], by = cy[(k + 1) % 4];
        const double L = std::hypot(bx - ax, by - ay);
        if (std::abs((bx - ax) * (p.y - ay) - (by - ay) * (p.x - ax)) / L < 1e-9) return true;
    }
    return false;
}

static void testQuadraticPatch() {
    std::puts("Patch test phần tử bậc hai T6, Q8 (chuyển vị tuyến tính trên biên, lưới méo)");
    const double cx[4] = {0, 4, 3.6, 0.3}, cy[4] = {0, 0.4, 2.7, 2.0};
    const double E = 2e8, nu = 0.3, a = 1e-3, b = 0.5e-3, c = -0.3e-3, d = 2e-3; // ux = a x + b y, uy = c x + d y
    for (int mode = 0; mode < 2; ++mode)
        for (ElemType type : {ElemType::Tri6, ElemType::Quad8}) {
            Model m;
            Material mat;
            mat.E = E;
            mat.nu = nu;
            m.materials.push_back(mat);
            m.planeMode = mode == 0 ? PlaneMode::Stress : PlaneMode::Strain;
            meshPatch(m, cx, cy, 3, 2, type, 0, 0.4);
            m.cases.push_back({"TH1"});
            for (size_t i = 0; i < m.nodes.size(); ++i) {
                const Node& p = m.nodes[i];
                if (!onPatchBoundary(p, cx, cy)) continue;
                fix(m, (int)i, true, true, false);
                m.cases[0].disp.push_back({(int)i, UX, a * p.x + b * p.y});
                m.cases[0].disp.push_back({(int)i, UY, c * p.x + d * p.y});
            }
            auto r = run(m);
            if (!r.ok) continue;
            double du = 0;
            for (size_t i = 0; i < m.nodes.size(); ++i) {
                const Node& p = m.nodes[i];
                du = std::max(du, std::abs(r.sets[0].U[3 * i] - (a * p.x + b * p.y)));
                du = std::max(du, std::abs(r.sets[0].U[3 * i + 1] - (c * p.x + d * p.y)));
            }
            double sx, sy, txy, sz = 0;
            if (mode == 0) {
                const double k = E / (1 - nu * nu);
                sx = k * (a + nu * d);
                sy = k * (d + nu * a);
            } else {
                const double k = E / ((1 + nu) * (1 - 2 * nu));
                sx = k * ((1 - nu) * a + nu * d);
                sy = k * ((1 - nu) * d + nu * a);
                sz = nu * (sx + sy);
            }
            txy = E / (2 * (1 + nu)) * (b + c);
            double ds = 0;
            for (size_t e = 0; e < m.elements.size(); ++e)
                for (int k = 0; k < m.elements[e].nodeCount(); ++k) {
                    const auto& st = r.sets[0].elemStress[e][k];
                    ds = std::max({ds, std::abs(st[0] - sx), std::abs(st[1] - sy), std::abs(st[2] - txy), std::abs(st[3] - sz)});
                }
            const std::string nm = std::string(type == ElemType::Tri6 ? "T6" : "Q8") + (mode == 0 ? " ƯSP" : " BDP");
            check(du < 1e-12, (nm + ": chuyển vị nút trong chính xác").c_str(), du, 0);
            check(ds < 1e-6 * std::abs(sy), (nm + ": ứng suất không đổi chính xác").c_str(), ds, 0);
        }
}

static void testQuadraticCantilever() {
    std::puts("Công-xôn T6, Q8 so với lý thuyết dầm (Timoshenko), tải cắt đều ở đầu tự do");
    const double L = 10, H = 1, t = 1, P = 1, E = 1e6, nu = 0.25;
    const double I = t * H * H * H / 12, G = E / (2 * (1 + nu));
    const double exact = P * L * L * L / (3 * E * I) + P * L / (5.0 / 6.0 * G * H * t);
    for (ElemType type : {ElemType::Quad8, ElemType::Tri6}) {
        Model m;
        Material mat;
        mat.E = E;
        mat.nu = nu;
        m.materials.push_back(mat);
        meshRectangle(m, 0, 0, L, H, 10, 2, type, 0, t);
        m.cases.push_back({"TH1"});
        std::vector<int> tip;
        int mid = -1;
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            const Node& p = m.nodes[i];
            if (std::abs(p.x) < 1e-9) fix(m, (int)i, true, std::abs(p.y - H / 2) < 1e-9, false);
            if (std::abs(p.x - L) < 1e-9) {
                tip.push_back((int)i);
                if (std::abs(p.y - H / 2) < 1e-9) mid = (int)i;
            }
        }
        EdgeLoadSpec es;
        es.qy = -P / H;
        applyEdgeLoad(m, 0, tip, es);
        auto r = run(m);
        if (!r.ok || mid < 0) continue;
        const double v = -r.sets[0].U[3 * mid + 1];
        near(type == ElemType::Quad8 ? "Q8 (10×2): độ võng đầu (sai số < 2%)" : "T6 (10×2): độ võng đầu (sai số < 3%)", v,
             exact, type == ElemType::Quad8 ? 0.02 : 0.03);
    }
}

static void testAxisymmetric() {
    std::puts("Đối xứng trục: trụ dày chịu áp lực trong (Lamé, biến dạng phẳng εz = 0)");
    const double a = 1, b = 2, p = 1000, E = 2e7, nu = 0.3;
    const double sthA = p * (b * b + a * a) / (b * b - a * a);
    const double uA = (1 + nu) / E * p * a * a / (b * b - a * a) * ((1 - 2 * nu) * a + b * b / a);
    for (ElemType type : {ElemType::Quad8, ElemType::Tri6, ElemType::Quad4}) {
        Model m;
        Material mat;
        mat.E = E;
        mat.nu = nu;
        m.materials.push_back(mat);
        m.planeMode = PlaneMode::Axisymmetric;
        meshRectangle(m, a, 0, b - a, 0.1, type == ElemType::Quad4 ? 40 : 10, 1, type, 0, 1, false);
        m.cases.push_back({"TH1"});
        std::vector<int> inner;
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            fix(m, (int)i, false, true, false);
            if (std::abs(m.nodes[i].x - a) < 1e-9) inner.push_back((int)i);
        }
        EdgeLoadSpec es;
        es.pressure = true;
        es.p0 = p;
        es.clipNegative = false;
        applyEdgeLoad(m, 0, inner, es);
        double F = 0;
        for (auto& l : m.cases[0].nodal) F += l.F[0];
        const char* nm = type == ElemType::Quad8 ? "Q8" : type == ElemType::Tri6 ? "T6" : "Q4";
        near((std::string(nm) + ": tổng lực áp lực = p·2πa·h").c_str(), F, p * 2 * M_PI * a * 0.1, 1e-12);
        auto r = run(m);
        if (!r.ok) continue;
        const int n0 = inner[0];
        near((std::string(nm) + ": u(a) theo Lamé").c_str(), r.sets[0].U[3 * n0], uA, type == ElemType::Quad4 ? 2e-3 : 1e-4);
        near((std::string(nm) + ": σθ(a) = p(b²+a²)/(b²−a²)").c_str(), r.sets[0].nodeStress[n0][3], sthA,
             type == ElemType::Quad4 ? 0.03 : 0.01);
        near((std::string(nm) + ": σr(a) = −p").c_str(), r.sets[0].nodeStress[n0][0], -p, 0.0, type == ElemType::Quad4 ? 0.08 * p : 0.03 * p);
    }
    // patch test đối xứng trục (có nút trên trục): ur = α r, uz = β z ⇒ εr = εθ = α, εz = β
    for (ElemType type : {ElemType::Tri3, ElemType::Quad4, ElemType::Tri6, ElemType::Quad8}) {
        Model m;
        Material mat;
        mat.E = E;
        mat.nu = nu;
        m.materials.push_back(mat);
        m.planeMode = PlaneMode::Axisymmetric;
        const double cx[4] = {0, 2, 1.7, 0}, cy[4] = {0, 0.3, 1.8, 1.5};
        meshPatch(m, cx, cy, 3, 3, type, 0, 1, false);
        m.cases.push_back({"TH1"});
        const double al = 1e-3, be = -2e-3;
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            const Node& q = m.nodes[i];
            if (!onPatchBoundary(q, cx, cy)) continue;
            fix(m, (int)i, true, true, false);
            m.cases[0].disp.push_back({(int)i, UX, al * q.x});
            m.cases[0].disp.push_back({(int)i, UY, be * q.y});
        }
        auto r = run(m);
        if (!r.ok) continue;
        const double lam = E * nu / ((1 + nu) * (1 - 2 * nu)), mu = E / (2 * (1 + nu));
        const double tr = 2 * al + be, sr = lam * tr + 2 * mu * al, sz = lam * tr + 2 * mu * be;
        double ds = 0;
        for (size_t e = 0; e < m.elements.size(); ++e)
            for (int k = 0; k < m.elements[e].nodeCount(); ++k) {
                const auto& st = r.sets[0].elemStress[e][k];
                ds = std::max({ds, std::abs(st[0] - sr), std::abs(st[1] - sz), std::abs(st[2]), std::abs(st[3] - sr)});
            }
        const std::string nm = type == ElemType::Tri3 ? "T3" : type == ElemType::Quad4 ? "Q4" : type == ElemType::Tri6 ? "T6" : "Q8";
        check(ds < 1e-6 * std::abs(sz), (nm + ": patch test đối xứng trục").c_str(), ds, 0);
    }
}

static void testQuadraticMisc() {
    std::puts("Phần tử bậc hai: tải cạnh, ghi/đọc file, dao động riêng");
    {
        Model m;
        m.materials.push_back(Material());
        meshRectangle(m, 0, 0, 3, 1, 3, 1, ElemType::Quad8, 0, 0.5);
        m.cases.push_back({"TH1"});
        std::vector<int> top;
        for (size_t i = 0; i < m.nodes.size(); ++i)
            if (std::abs(m.nodes[i].y - 1) < 1e-9) top.push_back((int)i);
        EdgeLoadSpec es;
        es.qy = -6;
        const int ne = applyEdgeLoad(m, 0, top, es);
        check(ne == 3, "Q8: số cạnh mép trên", ne, 3);
        double sum = 0, corner = 0, midside = 0;
        for (auto& l : m.cases[0].nodal) {
            sum += l.F[1];
            if (std::abs(m.nodes[l.node].x - 0.5) < 1e-9) midside = l.F[1];
            if (std::abs(m.nodes[l.node].x) < 1e-9) corner = l.F[1];
        }
        near("Q8: tổng tải = q·L", sum, -18, 1e-12);
        near("Q8: nút giữa cạnh nhận 2/3 qa", midside, -6 * 2.0 / 3, 1e-12);
        near("Q8: nút góc biên nhận 1/6 qa", corner, -6 / 6.0, 1e-12);
    }
    {
        Model m;
        m.title = "bậc hai";
        m.materials.push_back(Material());
        m.planeMode = PlaneMode::Axisymmetric;
        meshRectangle(m, 0, 0, 2, 1, 2, 2, ElemType::Quad8, 0, 1);
        meshRectangle(m, 0, 1, 2, 1, 2, 1, ElemType::Tri6, 0, 1);
        m.cases.push_back({"TH1"});
        for (size_t i = 0; i < m.nodes.size(); ++i)
            if (m.nodes[i].y < 1e-9) fix(m, (int)i, false, true, false);
        m.cases[0].selfWeight = true;
        const std::string text = modelToText(m);
        Model m2;
        std::string err;
        const bool ok = modelFromText(text, m2, err);
        check(ok && m2.planeMode == PlaneMode::Axisymmetric && m2.elements.size() == m.elements.size() &&
                  m2.elements.back().type == ElemType::Tri6 && m2.elements[0].type == ElemType::Quad8,
              "ghi/đọc TRI6, QUAD8, PLANE AXISYM", 0, 0);
        check(ok && modelToText(m2) == text, "văn bản ghi lại trùng khớp", 0, 0);
        auto r1 = run(m);
        auto r2 = run(m2);
        double d = 0;
        if (r1.ok && r2.ok)
            for (size_t i = 0; i < r1.sets[0].U.size(); ++i) d = std::max(d, std::abs(r1.sets[0].U[i] - r2.sets[0].U[i]));
        check(r1.ok && r2.ok && d < 1e-15, "kết quả sau khi đọc lại", d, 0);
        // trọng lượng bản thân tổng = γ·V (thể tích tròn xoay π r² h)
        if (r1.ok) {
            double ry = 0;
            for (size_t i = 0; i < m.nodes.size(); ++i) ry += r1.sets[0].R[3 * i + 1];
            near("Đối xứng trục: phản lực = γ·π r²·h", ry, m.materials[0].gamma * M_PI * 4 * 2, 1e-9);
        }
    }
    {
        // công-xôn mảnh Q8: tần số riêng thứ nhất so với dầm Euler–Bernoulli
        const double L = 10, H = 0.5, E = 3e7, gamma = 25;
        Model m;
        Material mat;
        mat.E = E;
        mat.nu = 0.2;
        mat.gamma = gamma;
        m.materials.push_back(mat);
        meshRectangle(m, 0, 0, L, H, 20, 2, ElemType::Quad8, 0, 1);
        m.cases.push_back({"TH1"});
        for (size_t i = 0; i < m.nodes.size(); ++i)
            if (m.nodes[i].x < 1e-9) fix(m, (int)i, true, true, false);
        EigenOptions o;
        o.kind = EigenKind::Vibration;
        o.modes = 2;
        auto r = eigenAnalysis(m, o);
        check(r.ok, r.error.c_str(), 0, 0);
        if (r.ok) {
            const double EI = E * H * H * H / 12, rhoA = gamma / 9.81 * H;
            near("Q8: ω1 công-xôn (sai số < 2%)", r.modes[0].omega, 1.875104 * 1.875104 * std::sqrt(EI / (rhoA * std::pow(L, 4))), 0.02);
        }
    }
}

// ---------------------------------------------------------------- dòng chảy chất lỏng (FVM)
static Model fluidDomain(const ChannelSpec& s) {
    Model m;
    m.materials.push_back(Material{});
    m.cases.push_back(LoadCase("TH1"));
    meshChannel(m, s, 0);
    return m;
}

static bool solveSteady(FluidSolver& sv, int maxIt = 6000) {
    for (int i = 0; i < maxIt && !sv.converged(); ++i)
        if (!sv.step()) return false;
    return sv.converged();
}

// u nhỏ nhất trên đường thẳng đứng x = 0,5 của hộp kín
static double cavityUmin(const FluidSolver& sv) {
    double umin = 1, u, v, p;
    for (int i = 1; i < 400; ++i)
        if (sv.probe(0.5, i / 400.0, u, v, p)) umin = std::min(umin, u);
    return umin;
}

static void testFluid() {
    std::printf("Dòng chảy chất lỏng FVM: Poiseuille, hộp kín, Schäfer–Turek\n");
    // ---- lưới: nhận diện biên, phát hiện phần tử chồng lấn
    {
        ChannelSpec s;
        s.obstacle = 0;
        s.L = 1;
        s.H = 0.5;
        s.h = 0.1;
        Model m = fluidDomain(s);
        FluidMesh M;
        check(buildFluidMesh(m, M).empty(), "dựng lưới FVM", 0, 0);
        near("số ô", M.nCells, 50);
        near("mặt biên trái", M.sideFaces[(int)FluidSide::Left], 5);
        near("mặt biên dưới", M.sideFaces[(int)FluidSide::Bottom], 10);
        near("không có biên trong", M.sideFaces[(int)FluidSide::Inner], 0);
        near("số mặt trong", M.nInternal, 4 * 10 + 9 * 5);
        double V = 0;
        for (double v : M.vol) V += v;
        near("tổng diện tích ô", V, 0.5);
        m.elements.push_back(m.elements[0]);
        check(!buildFluidMesh(m, M).empty(), "phát hiện phần tử trùng nhau", 0, 0);
    }
    FluidParams pois;
    pois.rho = 1;
    pois.mu = 1;
    pois.tol = 1e-8;
    pois.bc[(int)FluidSide::Left] = {FluidBcKind::Inlet, 1, 0, true};
    pois.bc[(int)FluidSide::Right] = {FluidBcKind::Outlet};
    // ---- dòng Poiseuille: dp/dx = −8μ·Umax/H², hội tụ bậc hai trên tứ giác và tam giác lệch
    for (int tri = 0; tri < 2; ++tri) {
        double err[2];
        for (int r = 0; r < 2; ++r) {
            ChannelSpec s;
            s.obstacle = 0;
            s.L = 1;
            s.H = 0.1;
            s.h = 0.1 / (8 << r);
            s.triangles = tri;
            FluidSolver sv;
            check(sv.init(fluidDomain(s), pois), "khởi tạo Poiseuille", 0, 0);
            check(solveSteady(sv), "Poiseuille hội tụ", sv.iteration(), 0);
            double u, v, p1, p2;
            sv.probe(0.4, 0.05, u, v, p1);
            sv.probe(0.8, 0.05, u, v, p2);
            err[r] = std::abs((p2 - p1) / 0.4 + 800);
            check(sv.massImbalance() < 1e-6, "bảo toàn khối lượng", sv.massImbalance(), 0);
        }
        check(err[1] < 0.02 * 800, tri ? "Poiseuille tam giác: dp/dx sai < 2%" : "Poiseuille tứ giác: dp/dx sai < 2%", err[1], 0);
        check(err[1] < 0.4 * err[0], "Poiseuille: sai số giảm khi chia mịn (bậc hai)", err[1] / err[0], 0.25);
    }
    // ---- hộp kín nắp trượt Re = 100: u nhỏ nhất trên đường giữa −0,2140 (Botella & Peyret; Ghia: −0,2109)
    FluidParams cav;
    cav.rho = 1;
    cav.mu = 0.01;
    cav.bc[(int)FluidSide::Top] = {FluidBcKind::MovingWall, 1, 0, false};
    ChannelSpec cs;
    cs.obstacle = 0;
    cs.L = cs.H = 1;
    cs.h = 1.0 / 40;
    double uminRef = 0; // nghiệm dừng hội tụ chặt, dùng so sánh với bài toán không dừng
    for (int tri = 0; tri < 2; ++tri) {
        cs.triangles = tri;
        FluidSolver sv;
        sv.init(fluidDomain(cs), cav);
        check(solveSteady(sv), "hộp kín hội tụ", sv.iteration(), 0);
        const double umin = cavityUmin(sv);
        near(tri ? "hộp kín (tam giác): u nhỏ nhất" : "hộp kín (tứ giác): u nhỏ nhất", umin, -0.2140, 0.02);
    }
    cs.triangles = false;
    // nghiệm hội tụ không phụ thuộc hệ số giảm dư (hiệu chỉnh Majumdar)
    {
        FluidParams c2 = cav;
        c2.alphaU = 0.5;
        c2.alphaP = 0.2;
        c2.tol = 1e-7;
        FluidParams c3 = cav;
        c3.tol = 1e-7;
        FluidSolver a, b;
        a.init(fluidDomain(cs), c2);
        b.init(fluidDomain(cs), c3);
        solveSteady(a, 20000);
        solveSteady(b, 20000);
        near("nghiệm không phụ thuộc αu, αp", cavityUmin(a), cavityUmin(b), 1e-4);
        uminRef = cavityUmin(b);
    }
    // bài toán không dừng tiến tới trạng thái dừng: trùng nghiệm dừng, không phụ thuộc Δt (hiệu chỉnh Choi)
    for (TimeScheme ts : {TimeScheme::Euler, TimeScheme::BDF2}) {
        FluidParams ct = cav;
        ct.transient = true;
        ct.timeScheme = ts;
        ct.dt = 0.5;
        ct.outerIters = 30;
        ct.outerTol = 1e-7;
        ct.alphaU = 0.8;
        ct.alphaP = 0.5;
        FluidSolver sv;
        sv.init(fluidDomain(cs), ct);
        for (int i = 0; i < 80; ++i) sv.step();
        near(ts == TimeScheme::Euler ? "không dừng (Euler) → nghiệm dừng" : "không dừng (BDF2) → nghiệm dừng", cavityUmin(sv),
             uminRef, 5e-4);
    }
    // ---- Schäfer–Turek 2D-1 (Re = 20): Cd = 5,5795, Δp = 0,11752
    {
        ChannelSpec s; // mặc định: kênh 2,2×0,41, trụ D = 0,1 tại (0,2; 0,2)
        FluidParams p;
        p.rho = 1;
        p.mu = 1e-3;
        p.bc[(int)FluidSide::Left] = {FluidBcKind::Inlet, 0.3, 0, true};
        p.bc[(int)FluidSide::Right] = {FluidBcKind::Outlet};
        FluidSolver sv;
        check(sv.init(fluidDomain(s), p), "khởi tạo trụ tròn", 0, 0);
        near("Re = 20", sv.reynolds(), 20, 1e-9);
        check(solveSteady(sv), "trụ tròn hội tụ", sv.iteration(), 0);
        const FluidForces F = sv.forces();
        near("Schäfer–Turek: Cd", F.Cd, 5.5795, 0.02);
        double u, v, pa, pb;
        sv.probe(0.15, 0.2, u, v, pa);
        sv.probe(0.25, 0.2, u, v, pb);
        near("Schäfer–Turek: Δp trước – sau trụ", pa - pb, 0.11752, 0.02);
        const auto cp = sv.innerCp();
        check(!cp.empty() && std::abs(cp.front()[0]) < 180, "phân bố Cp quanh trụ", (double)cp.size(), 80);
        double cpMax = -1e9;
        for (const auto& c : cp) cpMax = std::max(cpMax, c[1]);
        check(cpMax > 1.0, "Cp lớn nhất tại điểm đón dòng > 1 (Re nhỏ)", cpMax, 1.0);
        check(!sv.streamlines(sv.defaultSeeds(10)).empty(), "đường dòng", 0, 0);
    }
    // ---- kiểm tra dữ liệu vào
    {
        FluidParams p = pois;
        p.bc[(int)FluidSide::Right] = {FluidBcKind::Wall};
        ChannelSpec s;
        s.obstacle = 0;
        s.L = 1;
        s.H = 0.1;
        s.h = 0.025;
        FluidSolver sv;
        check(!sv.init(fluidDomain(s), p), "báo lỗi: có cửa vào mà không có cửa ra", 0, 0);
        Model empty;
        check(!sv.init(empty, pois), "báo lỗi: chưa có phần tử phẳng", 0, 0);
    }
}

// Miền chảy khai báo (cạnh cửa vào AB, chiều dài L) + vật cản vẽ bằng thanh khép kín; quay toàn bộ một góc ang
static Model flowModel(double ang, bool obstacle = true) {
    Model m;
    m.materials.push_back(Material{});
    m.sections.push_back(Section{});
    m.cases.push_back(LoadCase("TH1"));
    const double c = std::cos(ang), s = std::sin(ang);
    auto R = [&](double x, double y) { return std::array<double, 2>{c * x - s * y, s * x + c * y}; };
    FlowDomain& f = m.flow;
    f.active = true;
    const auto A = R(0, 0), B = R(0, 0.41);
    f.ax = A[0], f.ay = A[1], f.bx = B[0], f.by = B[1];
    f.side = -1; // bên phải của A→B (hướng +x khi chưa quay)
    f.length = 2.2;
    f.U = 0.2;
    f.parabolic = true;
    f.rho = 1;
    f.mu = 1e-3;
    f.medium = 4;
    f.h = 0.006;
    f.hFar = 0.03;
    if (obstacle) {
        const int n = 64, n0 = (int)m.nodes.size();
        for (int i = 0; i < n; ++i) {
            const auto P = R(0.2 + 0.05 * std::cos(2 * M_PI * i / n), 0.2 + 0.05 * std::sin(2 * M_PI * i / n));
            m.nodes.push_back({P[0], P[1]});
        }
        for (int i = 0; i < n; ++i) {
            Element e;
            e.type = ElemType::Frame;
            e.n = {n0 + i, n0 + (i + 1) % n, -1, -1};
            m.elements.push_back(e);
        }
    }
    return m;
}

static void testFlowDomain() {
    std::printf("Dòng chảy khai báo: lưới tam giác Delaunay, vật cản vẽ bằng thanh\n");
    // ---- lưới tam giác: diện tích đúng, góc nhỏ nhất ≥ 20°, vật cản lõm
    {
        DomainSpec ds;
        ds.outer = {{0, 0}, {1, 0}, {1, 0.5}, {0, 0.5}};
        ds.holes = {{{0.3, 0.2}, {0.5, 0.25}, {0.3, 0.3}, {0.35, 0.25}}}; // mũi tên lõm
        ds.h = 0.01;
        ds.hFar = 0.05;
        Model m;
        std::string err;
        const int n = meshDomain(m, ds, 0, &err);
        check(n > 0, ("lưới tam giác: " + err).c_str(), n, 0);
        double A = 0, minAng = 180;
        for (const auto& e : m.elements) {
            const Node &a = m.nodes[e.n[0]], &b = m.nodes[e.n[1]], &c = m.nodes[e.n[2]];
            A += 0.5 * ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
            const double L[3] = {std::hypot(b.x - c.x, b.y - c.y), std::hypot(a.x - c.x, a.y - c.y), std::hypot(a.x - b.x, a.y - b.y)};
            for (int k = 0; k < 3; ++k) {
                const double x = L[k], y = L[(k + 1) % 3], z = L[(k + 2) % 3];
                minAng = std::min(minAng, std::acos(std::clamp((y * y + z * z - x * x) / (2 * y * z), -1.0, 1.0)) * 180 / M_PI);
            }
        }
        near("lưới tam giác: tổng diện tích = miền − vật cản", A, 0.5 - std::abs(polygonArea(ds.holes[0])), 1e-9);
        check(minAng >= 20, "lưới tam giác: góc nhỏ nhất ≥ 20°", minAng, 20);
        // vật cản cắt biên phải: chỉ phần nằm trong miền bị khoét
        ds.holes.push_back({{0.9, 0.2}, {1.2, 0.2}, {1.2, 0.3}, {0.9, 0.3}});
        Model m2;
        check(meshDomain(m2, ds, 0, &err) > 0, ("vật cản cắt biên miền: " + err).c_str(), 0, 0);
        double A2 = 0;
        for (const auto& e : m2.elements) {
            const Node &a = m2.nodes[e.n[0]], &b = m2.nodes[e.n[1]], &c = m2.nodes[e.n[2]];
            A2 += 0.5 * ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
        }
        near("vật cản cắt biên: diện tích", A2, 0.5 - std::abs(polygonArea(ds.holes[0])) - 0.1 * 0.1, 1e-9);
    }
    // ---- vòng thanh khép kín
    {
        Model m = flowModel(0);
        Element e;
        e.type = ElemType::Frame;
        m.nodes.push_back({1, 0.1});
        m.nodes.push_back({1.2, 0.1});
        e.n = {(int)m.nodes.size() - 2, (int)m.nodes.size() - 1, -1, -1};
        m.elements.push_back(e); // một thanh hở: bỏ qua
        std::vector<int> bars;
        for (int i = 0; i < (int)m.elements.size(); ++i) bars.push_back(i);
        int ignored = 0;
        const auto loops = closedLoops(m, bars, &ignored);
        check(loops.size() == 1 && loops[0].size() == 64 && ignored == 1, "tìm vòng thanh khép kín", (double)loops.size(), 1);
    }
    // ---- trụ tròn vẽ bằng 64 thanh trong miền chảy (Schäfer–Turek 2D-1); quay 30° cho cùng kết quả
    double cd0 = 0;
    for (double ang : {0.0, M_PI / 6}) {
        Model m = flowModel(ang), mesh;
        FluidParams p;
        std::string note;
        const std::string err = buildFlowProblem(m, mesh, p, &note);
        check(err.empty(), ("lập bài toán dòng chảy: " + err).c_str(), 0, 0);
        FluidSolver sv;
        check(sv.init(mesh, p), "khởi tạo dòng chảy khai báo", 0, 0);
        near("Re = 20", sv.reynolds(), 20, 1e-6);
        near("biên trong = 64 cạnh vật cản (chia nhỏ)", sv.mesh().sideFaces[(int)FluidSide::Inner] >= 64, 1);
        check(solveSteady(sv), "dòng chảy khai báo hội tụ", sv.iteration(), 0);
        const FluidForces F = sv.forces();
        if (ang == 0) {
            cd0 = F.Cd;
            near("trụ vẽ bằng thanh: Cd (Schäfer–Turek 5,58)", F.Cd, 5.5795, 0.04);
        } else
            near("miền chảy xoay 30°: Cd không đổi", F.Cd, cd0, 0.01);
    }
    // ---- ghi / đọc khai báo dòng chảy
    {
        Model m = flowModel(0.3, false), r;
        m.flow.transient = true;
        m.flow.slipWalls = true;
        std::string err;
        check(modelFromText(modelToText(m), r, err), err.c_str(), 0, 0);
        check(r.flow.active && r.flow.transient && r.flow.slipWalls && r.flow.parabolic && r.flow.side == -1, "ghi/đọc FLOW", 0, 0);
        near("ghi/đọc FLOW: xB", r.flow.bx, m.flow.bx, 1e-10);
        near("ghi/đọc FLOW: U", r.flow.U, 0.2, 1e-12);
    }
    // ---- vật cản cắt biên miền chảy (tam giác cắm qua đáy kênh)
    {
        Model m = flowModel(0, false), mesh;
        for (auto [x, y] : {std::pair{0.5, -0.05}, {0.7, -0.05}, {0.6, 0.1}}) m.nodes.push_back({x, y});
        for (int i = 0; i < 3; ++i) {
            Element e;
            e.type = ElemType::Frame;
            e.n = {i, (i + 1) % 3, -1, -1};
            m.elements.push_back(e);
        }
        FluidParams p;
        check(buildFlowProblem(m, mesh, p).empty(), "vật cản cắt biên miền chảy: phần trong miền được khoét", 0, 0);
    }
}

// Mặt cắt đứng: miền rộng (cao) H, dài L, chảy theo +x; vật cản chữ nhật [x0, x1]×[y0, y1] vẽ bằng 4 thanh
static Model sectionModel(int section, double U, double x0, double x1, double y0, double y1) {
    Model m;
    m.materials.push_back(Material{});
    m.sections.push_back(Section{});
    m.cases.push_back(LoadCase("TH1"));
    FlowDomain& f = m.flow;
    f.active = true;
    f.section = section;
    f.ax = 0, f.ay = 0, f.bx = 0, f.by = 2; // cửa vào thẳng đứng từ đáy (y = 0) lên mặt thoáng (y = 2)
    f.side = -1;
    f.length = 6;
    f.U = U;
    f.medium = 0;
    f.rho = 1000;
    f.mu = 1;
    f.h = 0.05;
    f.hFar = 0.2;
    const double P[4][2] = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
    for (auto& q : P) m.nodes.push_back({q[0], q[1]});
    for (int i = 0; i < 4; ++i) {
        Element e;
        e.type = ElemType::Frame;
        e.n = {i, (i + 1) % 4, -1, -1};
        m.elements.push_back(e);
    }
    return m;
}

static void testFlowSections() {
    std::printf("Dòng chảy: mặt bằng / mặt cắt đứng, vật cản chạm biên, áp suất thủy tĩnh\n");
    const double g = 9.81, rho = 1000;
    auto solve = [](const Model& m, FluidSolver& sv, const char* what) {
        Model mesh;
        FluidParams p;
        const std::string err = buildFlowProblem(m, mesh, p);
        check(err.empty(), (std::string(what) + ": " + err).c_str(), 0, 0);
        if (!err.empty()) return false;
        check(sv.init(mesh, p), what, 0, 0);
        return solveSteady(sv, 4000);
    };
    // ---- vật cản ngập hoàn toàn, nước gần như đứng yên: lực = lực đẩy Acsimet ρ·g·V hướng lên
    {
        FluidSolver sv;
        solve(sectionModel(1, 1e-6, 2, 3, 0.5, 1.0), sv, "vật cản ngập trong nước tĩnh");
        const FluidForces F = sv.forces();
        near("nước tĩnh: lực đẩy Acsimet Fy = ρ·g·V", F.Fy, rho * g * 0.5, 2e-3);
        near("nước tĩnh: Fx = 0", F.Fx, 0, 0, 1e-3 * rho * g * 0.5);
    }
    // ---- công trình đặt trên đáy (cạnh đáy trùng mặt đất): chỉ có áp lực nước phía trên đè xuống, ngang cân bằng
    for (double yb : {0.0, -0.3}) { // trùng đáy, hoặc vẽ cắm xuống dưới đáy
        FluidSolver sv;
        solve(sectionModel(1, 1e-6, 2, 3, yb, 1.2), sv, yb == 0 ? "công trình đặt trên đáy" : "công trình cắm xuống dưới đáy");
        const FluidForces F = sv.forces();
        near("trên đáy: Fy = −ρ·g·(mặt thoáng − đỉnh)·rộng", F.Fy, -rho * g * (2 - 1.2) * 1, 2e-3);
        near("trên đáy: Fx = 0", F.Fx, 0, 0, 1e-3 * rho * g);
        check(sv.mesh().sideFaces[(int)FluidSide::Inner] > 0, "trên đáy: có biên trong (mặt công trình)", 0, 0);
        near("trên đáy: L tham chiếu = chiều cao công trình", sv.Lref(), 1.2, 1e-9);
        double area = 0;
        for (double v : sv.mesh().vol) area += v;
        near("trên đáy: diện tích miền = miền − phần công trình trên đáy", area, 2 * 6 - 1 * 1.2, 1e-9);
    }
    // ---- có dòng chảy: mặt cắt đứng (đáy cứng, mặt thoáng trượt) — lực cản dương, công trình bị đẩy theo dòng
    {
        FluidSolver sv;
        check(solve(sectionModel(1, 0.05, 2, 3, 0, 1.2), sv, "dòng chảy qua công trình trên đáy"), "hội tụ", 0, 0);
        check(sv.forces().Fd > 0, "lực cản theo hướng dòng chảy > 0", sv.forces().Fd, 0);
        check(sv.params().bc[(int)FluidSide::Bottom].kind == FluidBcKind::Wall &&
                  sv.params().bc[(int)FluidSide::Top].kind == FluidBcKind::Symmetry,
              "mặt cắt đứng: đáy cứng, mặt thoáng trượt", 0, 0);
    }
    // ---- mặt bằng: mỏ hàn chạm bờ, hai bờ cứng; không có thủy tĩnh
    {
        FluidSolver sv;
        check(solve(sectionModel(0, 0.05, 2, 2.4, -0.2, 0.8), sv, "mỏ hàn chạm bờ (mặt bằng)"), "hội tụ", 0, 0);
        check(sv.params().hydroG == 0 && sv.params().bc[(int)FluidSide::Top].kind == FluidBcKind::Wall, "mặt bằng: hai bờ cứng, không thủy tĩnh", 0, 0);
    }
    // ---- vật cản chắn ngang toàn bộ dòng chảy; mặt cắt đứng với dòng chảy thẳng đứng
    {
        Model mesh;
        FluidParams p;
        check(!buildFlowProblem(sectionModel(1, 1, 2, 3, -0.5, 2.5), mesh, p).empty(), "báo lỗi: vật cản chắn kín dòng chảy", 0, 0);
        Model m = sectionModel(1, 1, 2, 3, 0.5, 1);
        m.flow.bx = 2, m.flow.by = 0; // cửa vào nằm ngang → dòng chảy thẳng đứng (ống đứng): hai biên bên là thành
        m.nodes.clear();
        m.elements.clear();
        check(buildFlowProblem(m, mesh, p).empty() && p.bc[(int)FluidSide::Bottom].kind == FluidBcKind::Wall &&
                  p.bc[(int)FluidSide::Top].kind == FluidBcKind::Wall,
              "mặt cắt đứng, dòng chảy thẳng đứng: hai biên bên là thành", 0, 0);
    }
    // ---- ghi / đọc
    {
        Model m = sectionModel(1, 1, 2, 3, 0, 1), r;
        m.flow.topWall = true;
        m.flow.hydrostatic = false;
        std::string err;
        check(modelFromText(modelToText(m), r, err) && r.flow.section == 1 && r.flow.topWall && !r.flow.hydrostatic,
              "ghi/đọc loại mặt cắt", 0, 0);
    }
}

// ---------------------------------------------------------------- chảy rối k-ε, mặt cắt nghiêng, mặt thoáng VOF
static void testTurbulenceFreeSurface() {
    std::printf("Chảy rối k-ε, mặt cắt đứng nghiêng, mặt thoáng tự do (VOF)\n");
    // ---- dòng chảy rối phát triển đầy đủ trong kênh phẳng 2h = 0,2 m, Re = 10⁵: hệ số ma sát theo Dean C_f = 0,073·Re^−0,25
    {
        Model m;
        m.materials.push_back(Material{});
        ChannelSpec s;
        s.obstacle = 0;
        s.L = 6;
        s.H = 0.2;
        s.h = 0.02;
        meshChannel(m, s, 0);
        FluidParams p;
        p.rho = 1000;
        p.mu = 1e-3;
        p.turbulence = TurbModel::KEpsilon;
        p.scheme = ConvScheme::Upwind;
        p.tol = 1e-5;
        p.bc[(int)FluidSide::Left] = {FluidBcKind::Inlet, 0.5, 0, false};
        p.bc[(int)FluidSide::Right] = {FluidBcKind::Outlet};
        FluidSolver sv;
        check(sv.init(m, p), "khởi tạo kênh chảy rối", 0, 0);
        check(solveSteady(sv), "kênh chảy rối hội tụ", sv.iteration(), 0);
        double u, v, p1, p2, uc;
        sv.probe(4.2, 0.1, u, v, p1);
        sv.probe(5.4, 0.1, uc, v, p2);
        const double tau = -(p2 - p1) / 1.2 * 0.1, Cf = tau / (0.5 * 1000 * 0.25);
        near("k-ε: hệ số ma sát kênh (Dean)", Cf, 0.073 * std::pow(1e5, -0.25), 0.08);
        near("k-ε: vận tốc tâm / trung bình (≈ 1,12)", uc / 0.5, 1.12, 0.05);
        double mutMax = 0;
        for (double x : sv.turbVisc()) mutMax = std::max(mutMax, x);
        check(mutMax > 10 * p.mu, "k-ε: độ nhớt rối ≫ độ nhớt phân tử", mutMax / p.mu, 10);
    }
    // ---- vật cản trong dòng nước Re ≈ 10⁶ (trước đây chảy tầng phân kỳ): k-ε dừng hội tụ, Cd của trụ vuông ≈ 2
    {
        Model m = sectionModel(0, 1, 4, 5, 2.5, 3.5), mesh;
        m.flow.by = 6;
        m.flow.length = 16;
        m.flow.mu = 1.002e-3;
        m.flow.h = 0;
        m.flow.hFar = 0;
        FluidParams p;
        std::string note;
        check(buildFlowProblem(m, mesh, p, &note).empty(), "lập bài toán Re lớn", 0, 0);
        check(p.turbulence == TurbModel::KEpsilon && !p.transient, "Re lớn: tự chọn k-ε, bài toán dừng", 0, 0);
        FluidSolver sv;
        sv.init(mesh, p);
        check(solveSteady(sv, 3000), "k-ε quanh trụ vuông hội tụ", sv.iteration(), 0);
        const double Cd = sv.forces().Cd;
        check(Cd > 1.7 && Cd < 2.7, "k-ε: Cd trụ vuông (thực nghiệm 2,1–2,2)", Cd, 2.15);
    }
    // ---- mặt cắt đứng nghiêng (kênh dốc 10°): đáy là cạnh thấp hơn, mặt thoáng trượt
    {
        Model m = sectionModel(1, 0.05, 2, 3, 0, 1), mesh;
        const double a = 10 * M_PI / 180, c = std::cos(a), s = std::sin(a);
        auto rot = [&](double& x, double& y) {
            const double X = c * x + s * y, Y = -s * x + c * y; // quay −10°: dòng chảy đi xuống dốc
            x = X, y = Y;
        };
        rot(m.flow.ax, m.flow.ay);
        rot(m.flow.bx, m.flow.by);
        for (auto& n : m.nodes) rot(n.x, n.y);
        FluidParams p;
        check(buildFlowProblem(m, mesh, p).empty(), "mặt cắt đứng nghiêng 10°", 0, 0);
        check(p.bc[(int)FluidSide::Bottom].kind == FluidBcKind::Wall && p.bc[(int)FluidSide::Top].kind == FluidBcKind::Symmetry,
              "nghiêng: đáy cứng, mặt thoáng trượt", 0, 0);
        double cc[4][2];
        m.flow.corners(cc);
        near("nghiêng: mặt thoáng thủy tĩnh tại góc hạ lưu cạnh trên", p.yFree, cc[2][1], 1e-12);
        FluidSolver sv;
        check(sv.init(mesh, p) && solveSteady(sv, 3000), "nghiêng: hội tụ", 0, 0);
    }
    // ---- VOF: nước tĩnh trong bể 1×1 m, mực 0,5 m — không có dòng giả, áp suất đáy ρ·g·h, thể tích bảo toàn
    FluidParams vp;
    vp.rho = 998.2;
    vp.mu = 1.002e-3;
    vp.vof = true;
    vp.dt = 0.01;
    vp.bc[(int)FluidSide::Left] = {FluidBcKind::Wall};
    vp.bc[(int)FluidSide::Right] = {FluidBcKind::Wall};
    vp.bc[(int)FluidSide::Bottom] = {FluidBcKind::Wall};
    vp.bc[(int)FluidSide::Top] = {FluidBcKind::Opening};
    {
        Model m;
        m.materials.push_back(Material{});
        ChannelSpec s;
        s.obstacle = 0;
        s.L = s.H = 1;
        s.h = 0.04;
        meshChannel(m, s, 0);
        FluidParams p = vp;
        p.initLevel = 0.5;
        FluidSolver sv;
        check(sv.init(m, p), "khởi tạo bể nước tĩnh", 0, 0);
        const double V0 = sv.waterVolume();
        near("VOF: thể tích nước ban đầu", V0, 0.5, 1e-9);
        while (sv.time() < 1.0)
            if (!sv.step()) break;
        double umax = 0;
        for (int c = 0; c < sv.mesh().nCells; ++c) umax = std::max(umax, std::hypot(sv.u()[c], sv.v()[c]));
        check(umax < 0.02, "VOF nước tĩnh: dòng giả < 2 cm/s", umax, 0.02);
        double u, v, pb;
        sv.probe(0.5, 0.02, u, v, pb);
        near("VOF nước tĩnh: áp suất đáy ρ·g·h", pb, 998.2 * 9.81 * 0.48, 0.01);
        near("VOF: bảo toàn thể tích nước", sv.waterVolume(), V0, 1e-6);
    }
    // ---- VOF: bể 1×1 m, mực 0,5 m, đáy thấm q = 0,05 m/s — nước thoát q·L·t, tổng thể tích bảo toàn
    {
        Model m;
        m.materials.push_back(Material{});
        ChannelSpec s;
        s.obstacle = 0;
        s.L = s.H = 1;
        s.h = 0.04;
        meshChannel(m, s, 0);
        FluidParams p = vp;
        p.initLevel = 0.5;
        p.bc[(int)FluidSide::Bottom].drain = 0.05;
        FluidSolver sv;
        check(sv.init(m, p), "khởi tạo bể đáy thấm", 0, 0);
        const double V0 = sv.waterVolume();
        while (sv.time() < 1.0 - 1e-9)
            if (!sv.step()) break;
        check(!sv.diverged(), "đáy thấm: không phân kỳ", 0, 0);
        near("đáy thấm: thể tích còn lại + đã thấm = ban đầu", sv.waterVolume() + sv.drainedVolume(), V0, 1e-6);
        near("đáy thấm: lượng thấm q·L·t", sv.drainedVolume(), 0.05 * 1 * sv.time(), 0.03);
        double umax = 0;
        for (int c = 0; c < sv.mesh().nCells; ++c)
            if (sv.alpha()[c] > 0.5) umax = std::max(umax, std::hypot(sv.u()[c], sv.v()[c]));
        check(umax < 0.1, "đáy thấm: nước hạ đều, không xoáy giả", umax, 0.05);
    }
    // ---- VOF: đáy thoát tự do — cột nước không có lực đỡ rơi tự do: lượng thoát ½·g·t²·L; bể khô: không khí đứng yên
    {
        Model m;
        m.materials.push_back(Material{});
        ChannelSpec s;
        s.obstacle = 0;
        s.L = s.H = 1;
        s.h = 0.04;
        meshChannel(m, s, 0);
        FluidParams p = vp;
        p.initLevel = 0.5;
        p.bc[(int)FluidSide::Bottom] = {FluidBcKind::FreeDrain};
        FluidSolver sv;
        check(sv.init(m, p), "khởi tạo bể đáy thoát tự do", 0, 0);
        const double V0 = sv.waterVolume();
        while (sv.time() < 0.2 - 1e-9)
            if (!sv.step()) break;
        check(!sv.diverged(), "đáy thoát tự do: không phân kỳ", 0, 0);
        near("đáy thoát tự do: thể tích còn lại + đã thoát = ban đầu", sv.waterVolume() + sv.drainedVolume(), V0, 1e-6);
        near("đáy thoát tự do: cột nước rơi tự do ½·g·t²", sv.drainedVolume(), 0.5 * 9.81 * sv.time() * sv.time(), 0.1);

        FluidParams pa = p;
        pa.initLevel = 0; // bể khô: áp suất khí quyển ở đáy khớp cột không khí, không có dòng giả
        FluidSolver sa;
        check(sa.init(m, pa), "khởi tạo bể khô đáy thoát tự do", 0, 0);
        while (sa.time() < 0.5 - 1e-9)
            if (!sa.step()) break;
        double umax = 0;
        for (int c = 0; c < sa.mesh().nCells; ++c) umax = std::max(umax, std::hypot(sa.u()[c], sa.v()[c]));
        check(umax < 0.01, "đáy thoát tự do, bể khô: không khí đứng yên", umax, 0.01);
    }
    // ---- vỡ đập Martin & Moyce (1952): cột nước a × 2a đổ ra lòng bể; mép sóng Z = x/a theo T = t·√(2g/a).
    // Thí nghiệm Z(2) ≈ 2,33; các lời giải VOF công bố ≈ 2,5–2,7 (thí nghiệm trễ do nhấc cửa).
    {
        const double a = 0.146;
        Model m;
        m.materials.push_back(Material{});
        ChannelSpec s;
        s.obstacle = 0;
        s.L = 4 * a;
        s.H = 2.4 * a;
        s.h = a / 16;
        meshChannel(m, s, 0);
        FluidParams p = vp;
        p.initLevel = 2 * a;
        p.initGate = a;
        FluidSolver sv;
        sv.init(m, p);
        const double V0 = sv.waterVolume();
        double Z1 = 0, Z2 = 0;
        while (Z2 == 0 && sv.step()) {
            const double T = sv.time() * std::sqrt(2 * 9.81 / a);
            if (T < 1 && T + sv.dt() * std::sqrt(2 * 9.81 / a) < 1) continue;
            double xf = 0;
            for (int c = 0; c < sv.mesh().nCells; ++c)
                if (sv.mesh().cy[c] < a / 8 && sv.alpha()[c] > 0.5) xf = std::max(xf, sv.mesh().cx[c]);
            if (T >= 1 && Z1 == 0) Z1 = xf / a;
            if (T >= 2) Z2 = xf / a;
        }
        check(Z2 > 2.3 && Z2 < 2.95, "vỡ đập: mép sóng Z tại T = 2", Z2, 2.6);
        check(Z1 > 1.2 && Z1 < Z2, "vỡ đập: mép sóng tiến dần", Z1, 1.5);
        near("vỡ đập: bảo toàn thể tích nước", sv.waterVolume(), V0, 1e-6);
    }
    // ---- mặt thoáng tự do khai báo trong mô hình: nước chảy qua vật cản ngập, biên trên thông khí
    {
        Model m = sectionModel(1, 0.3, 2, 2.3, 0, 0.3), mesh;
        m.flow.freeSurface = true;
        m.flow.waterDepth = 1.0;
        m.flow.h = 0.08;
        m.flow.hFar = 0.15;
        FluidParams p;
        check(buildFlowProblem(m, mesh, p).empty(), "lập bài toán mặt thoáng tự do", 0, 0);
        check(p.vof && p.transient && p.bc[(int)FluidSide::Top].kind == FluidBcKind::Opening, "mặt thoáng: VOF, biên trên thông khí", 0, 0);
        near("mặt thoáng: cao độ mặt nước cửa vào", p.inletLevel, 1.0, 1e-12);
        FluidSolver sv;
        check(sv.init(mesh, p), "khởi tạo mặt thoáng", 0, 0);
        bool ok = true;
        for (int i = 0; i < 30 && ok; ++i) ok = sv.step();
        check(ok && !sv.diverged(), "mặt thoáng: chạy 30 bước không phân kỳ", sv.iteration(), 30);
        Model bad = sectionModel(1, 0.3, 2, 2.3, 0, 0.3);
        bad.flow.freeSurface = true;
        bad.flow.waterDepth = 3; // sâu hơn chiều cao AB
        check(!buildFlowProblem(bad, mesh, p).empty(), "báo lỗi: độ sâu nước vượt chiều cao miền", 0, 0);
        // đáy thấm khai báo trong mô hình: điều kiện thấm đặt trên biên dưới (mặt đất), các biên khác không thấm
        m.flow.drainGround = true;
        m.flow.drainRate = 0.02;
        FluidParams pd;
        check(buildFlowProblem(m, mesh, pd).empty() && pd.bc[(int)FluidSide::Bottom].kind == FluidBcKind::Wall &&
                  pd.bc[(int)FluidSide::Bottom].drain == 0.02 && pd.bc[(int)FluidSide::Inner].drain == 0 &&
                  pd.bc[(int)FluidSide::Top].drain == 0,
              "đáy thấm: điều kiện thấm trên mặt đất", 0, 0);
        FluidSolver sd;
        check(sd.init(mesh, pd), "khởi tạo mặt thoáng đáy thấm", 0, 0);
        ok = true;
        for (int i = 0; i < 30 && ok; ++i) ok = sd.step();
        check(ok && !sd.diverged() && sd.drainedVolume() > 0, "đáy thấm: chạy 30 bước, nước thoát qua đáy", sd.drainedVolume(), 0);
        // đáy thoát tự do: biên dưới thành FreeDrain; không có mặt thoáng tự do thì báo lỗi
        m.flow.drainMode = 1;
        FluidParams pf;
        check(buildFlowProblem(m, mesh, pf).empty() && pf.bc[(int)FluidSide::Bottom].kind == FluidBcKind::FreeDrain &&
                  pf.bc[(int)FluidSide::Top].kind == FluidBcKind::Opening,
              "đáy thoát tự do: biên dưới thoát tự do, biên trên thông khí", 0, 0);
        m.flow.freeSurface = false;
        check(!buildFlowProblem(m, mesh, pf).empty(), "báo lỗi: đáy thoát tự do khi không có mặt thoáng tự do", 0, 0);
    }
    // ---- ghi / đọc
    {
        Model m = sectionModel(1, 1, 2, 3, 0, 1), r;
        m.flow.turbulence = 2;
        m.flow.freeSurface = true;
        m.flow.waterDepth = 1.2;
        m.flow.gate = 0.7;
        m.flow.drainGround = true;
        m.flow.drainRate = 0.003;
        m.flow.drainMode = 1;
        std::string err;
        check(modelFromText(modelToText(m), r, err) && r.flow.turbulence == 2 && r.flow.freeSurface &&
                  std::abs(r.flow.waterDepth - 1.2) < 1e-12 && std::abs(r.flow.gate - 0.7) < 1e-12 && r.flow.drainGround &&
                  std::abs(r.flow.drainRate - 0.003) < 1e-12 && r.flow.drainMode == 1,
              "ghi/đọc mô hình rối, mặt thoáng", 0, 0);
    }
}

int main() {
    testFluid();
    testFlowDomain();
    testFlowSections();
    testTurbulenceFreeSurface();
    testDxfReport();
    testDesign();
    testDynamics();
    testArchTaperTrain();
    testForceMethod();
    testSectionGeom();
    testLackOfFit();
    testPDelta();
    testFoundations();
    testKinematics();
    testSplitMerge();
    testEigenSolver();
    testBuckling();
    testVibration();
    testComboGenerator();
    testInfluenceLines();
    testEnvelope();
    testEdgeLoad();
    testSimpleBeam();
    testCantilever();
    testFixedFixed();
    testPointAndMoment();
    testInclinedCantilever();
    testHinge();
    testThermalAndSettlement();
    testSelfWeight();
    testWinkler();
    testPatch();
    testPlaneCantilever();
    testQuadraticPatch();
    testQuadraticCantilever();
    testAxisymmetric();
    testQuadraticMisc();
    testMechanism();
    testSubdivide();
    testExamples();
    testIO();
    std::printf("\n%d kiểm tra đạt, %d lỗi\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
