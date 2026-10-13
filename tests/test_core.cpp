// Kiểm chứng lõi tính toán với lời giải giải tích.
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>

#include "core/Analysis.h"
#include "core/Examples.h"
#include "core/Library.h"
#include "core/Mesher.h"
#include "core/ModelIO.h"
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

int main() {
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
    testMechanism();
    testSubdivide();
    testExamples();
    testIO();
    std::printf("\n%d kiểm tra đạt, %d lỗi\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
