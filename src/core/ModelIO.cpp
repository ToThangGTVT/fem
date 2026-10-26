#include "core/ModelIO.h"

#include "core/SectionGeom.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>

namespace fem {

namespace {

std::string quote(const std::string& s) {
    std::string r = "\"";
    for (char c : s) r += (c == '"') ? '\'' : c;
    return r + "\"";
}

std::string num(double v) {
    char buf[40];
    std::snprintf(buf, sizeof buf, "%.12g", v);
    return buf;
}

const char* dirToken(LoadDir d) {
    switch (d) {
    case LoadDir::LocalX: return "LX";
    case LoadDir::LocalY: return "LY";
    case LoadDir::GlobalX: return "GX";
    case LoadDir::GlobalY: return "GY";
    }
    return "GY";
}

bool parseDir(const std::string& t, LoadDir& d) {
    if (t == "LX") d = LoadDir::LocalX;
    else if (t == "LY") d = LoadDir::LocalY;
    else if (t == "GX") d = LoadDir::GlobalX;
    else if (t == "GY") d = LoadDir::GlobalY;
    else return false;
    return true;
}

const char* dofToken(int d) { return d == UX ? "UX" : d == UY ? "UY" : "RZ"; }

std::vector<std::string> tokenize(const std::string& line) {
    std::vector<std::string> out;
    size_t i = 0;
    while (i < line.size()) {
        while (i < line.size() && std::isspace((unsigned char)line[i])) ++i;
        if (i >= line.size() || line[i] == '#') break;
        if (line[i] == '"') {
            size_t j = line.find('"', i + 1);
            if (j == std::string::npos) j = line.size();
            out.push_back(line.substr(i + 1, j - i - 1));
            i = j + 1;
        } else {
            size_t j = i;
            while (j < line.size() && !std::isspace((unsigned char)line[j])) ++j;
            out.push_back(line.substr(i, j - i));
            i = j;
        }
    }
    return out;
}

} // namespace

std::string modelToText(const Model& m) {
    std::ostringstream o;
    o << "FEM2D 1\n";
    o << "# Đơn vị: kN, m, °C. Chỉ số đánh từ 1.\n";
    o << "TITLE " << quote(m.title) << "\n";
    o << "PLANE " << (m.planeMode == PlaneMode::Stress ? "STRESS" : m.planeMode == PlaneMode::Strain ? "STRAIN" : "AXISYM") << "\n";
    if (m.secondOrder) o << "ANALYSIS PDELTA 1\n";
    if (m.flow.active) {
        const FlowDomain& f = m.flow;
        o << "# FLOW xA yA xB yB L chiều U parabol môi_trường rho mu biên_trượt h h_xa không_dừng dt t_cuối "
             "mặt_cắt(0 bằng, 1 đứng) trần_cứng thủy_tĩnh mô_hình_rối(0 tự động, 1 tầng, 2 k-ε) cường_độ_rối "
             "mặt_thoáng_tự_do độ_sâu_cửa_vào độ_sâu_ban_đầu đoạn_vỡ_đập đáy_thấm tốc_độ_thấm kiểu_thoát_đáy(0 thấm q, 1 thoát tự do)\n";
        o << "FLOW " << num(f.ax) << " " << num(f.ay) << " " << num(f.bx) << " " << num(f.by) << " " << num(f.length) << " "
          << f.side << " " << num(f.U) << " " << (f.parabolic ? 1 : 0) << " " << f.medium << " " << num(f.rho) << " "
          << num(f.mu) << " " << (f.slipWalls ? 1 : 0) << " " << num(f.h) << " " << num(f.hFar) << " "
          << (f.transient ? 1 : 0) << " " << num(f.dt) << " " << num(f.endTime) << " " << f.section << " "
          << (f.topWall ? 1 : 0) << " " << (f.hydrostatic ? 1 : 0) << " " << f.turbulence << " " << num(f.turbIntensity)
          << " " << (f.freeSurface ? 1 : 0) << " " << num(f.waterDepth) << " " << num(f.initDepth) << " " << num(f.gate)
          << " " << (f.drainGround ? 1 : 0) << " " << num(f.drainRate) << " " << f.drainMode << "\n";
    }
    o << "\n";
    o << "# MATERIAL tên E nu alpha gamma [loại(1 thép, 2 bê tông) f/Rb Rbt (MPa)]\n";
    for (const auto& x : m.materials) {
        o << "MATERIAL " << quote(x.name) << " " << num(x.E) << " " << num(x.nu) << " " << num(x.alpha) << " "
          << num(x.gamma);
        if (x.kind != 0 || x.fd != 0) o << " " << x.kind << " " << num(x.fd) << " " << num(x.ft);
        o << "\n";
    }
    o << "# SECTION tên A I h [y_trên y_dưới] [SHAPE loại p1..p6 | POLY n x1 y1 ... (mỗi chu tuyến)]\n";
    for (const auto& x : m.sections) {
        o << "SECTION " << quote(x.name) << " " << num(x.A) << " " << num(x.I) << " " << num(x.h);
        if (x.yTop > 0 || x.yBot > 0) o << " " << num(x.yTop) << " " << num(x.yBot);
        if (x.shape.kind == ShapeKind::Polygon) {
            for (const auto& c : x.shape.contours) {
                o << " POLY " << c.size();
                for (const auto& v : c) o << " " << num(v[0]) << " " << num(v[1]);
            }
        } else if (x.shape.kind != ShapeKind::None) {
            o << " SHAPE " << (int)x.shape.kind;
            for (double v : x.shape.p) o << " " << num(v);
        }
        o << "\n";
    }
    o << "\n# NODE x y\n";
    for (const auto& n : m.nodes) o << "NODE " << num(n.x) << " " << num(n.y) << "\n";
    o << "\n# TRUSS n1 n2 mat sec | FRAME n1 n2 mat sec khớp_i khớp_j kf [kG chỉ_nén] | TRI3 n1 n2 n3 mat t | QUAD4 n1..n4 mat t im | TRI6 n1..n6 mat t | QUAD8 n1..n8 mat t\n";
    for (const auto& e : m.elements) {
        switch (e.type) {
        case ElemType::Truss:
            o << "TRUSS " << e.n[0] + 1 << " " << e.n[1] + 1 << " " << e.mat + 1 << " " << e.sec + 1 << "\n";
            break;
        case ElemType::Frame:
            o << "FRAME " << e.n[0] + 1 << " " << e.n[1] + 1 << " " << e.mat + 1 << " " << e.sec + 1 << " "
              << (int)e.hingeI << " " << (int)e.hingeJ << " " << num(e.kf);
            if (e.kG != 0 || e.kfCompOnly) o << " " << num(e.kG) << " " << (int)e.kfCompOnly;
            o << "\n";
            break;
        case ElemType::Tri3:
            o << "TRI3 " << e.n[0] + 1 << " " << e.n[1] + 1 << " " << e.n[2] + 1 << " " << e.mat + 1 << " " << num(e.t)
              << "\n";
            break;
        case ElemType::Quad4:
            o << "QUAD4 " << e.n[0] + 1 << " " << e.n[1] + 1 << " " << e.n[2] + 1 << " " << e.n[3] + 1 << " "
              << e.mat + 1 << " " << num(e.t) << " " << (int)e.incompatible << "\n";
            break;
        case ElemType::Tri6:
        case ElemType::Quad8:
            o << (e.type == ElemType::Tri6 ? "TRI6" : "QUAD8");
            for (int k = 0; k < e.nodeCount(); ++k) o << " " << e.n[k] + 1;
            o << " " << e.mat + 1 << " " << num(e.t) << "\n";
            break;
        }
    }
    o << "\n# SUPPORT nút khóa_x khóa_y khóa_φ kx ky kφ\n";
    for (const auto& s : m.supports)
        o << "SUPPORT " << s.node + 1 << " " << (int)s.fix[0] << " " << (int)s.fix[1] << " " << (int)s.fix[2] << " "
          << num(s.k[0]) << " " << num(s.k[1]) << " " << num(s.k[2]) << "\n";
    o << "\n# LOADCASE tên hệ_số_TLBT loại(0 thường xuyên, 1 tạm thời) nhóm_loại_trừ γf\n#   NLOAD nút Fx Fy Mz\n#   ELOAD pt DIST hướng q1 q2 chiếu | POINT hướng P a | MOMENT M a | TEMP t Δt | FIT ΔL\n"
         "#   DISP nút UX|UY|RZ giá_trị   (hướng: LX LY GX GY)\n";
    for (const auto& lc : m.cases) {
        o << "LOADCASE " << quote(lc.name) << " " << num(lc.selfWeight) << " " << lc.kind << " " << lc.group << " "
          << num(lc.gammaF) << "\n";
        for (const auto& l : lc.nodal)
            o << "  NLOAD " << l.node + 1 << " " << num(l.F[0]) << " " << num(l.F[1]) << " " << num(l.F[2]) << "\n";
        for (const auto& l : lc.elem) {
            o << "  ELOAD " << l.elem + 1 << " ";
            switch (l.kind) {
            case ElemLoadKind::Distributed:
                o << "DIST " << dirToken(l.dir) << " " << num(l.q1) << " " << num(l.q2) << " " << (int)l.projected;
                break;
            case ElemLoadKind::Point: o << "POINT " << dirToken(l.dir) << " " << num(l.q1) << " " << num(l.a); break;
            case ElemLoadKind::Moment: o << "MOMENT " << num(l.q1) << " " << num(l.a); break;
            case ElemLoadKind::Temperature: o << "TEMP " << num(l.q1) << " " << num(l.q2); break;
            case ElemLoadKind::LackOfFit: o << "FIT " << num(l.q1); break;
            }
            o << "\n";
        }
        for (const auto& d : lc.disp)
            o << "  DISP " << d.node + 1 << " " << dofToken(d.dof) << " " << num(d.value) << "\n";
    }
    if (!m.combos.empty()) o << "\n# COMBO tên trường_hợp hệ_số ...\n";
    for (const auto& c : m.combos) {
        o << "COMBO " << quote(c.name) << (c.generated ? " AUTO" : "");
        for (const auto& t : c.terms) o << " " << t.first + 1 << " " << num(t.second);
        o << "\n";
    }
    return o.str();
}

bool modelFromText(const std::string& text, Model& out, std::string& error) {
    Model m;
    std::istringstream in(text);
    std::string line;
    int lineNo = 0;
    LoadCase* cur = nullptr;
    auto fail = [&](const std::string& msg) {
        error = "Dòng " + std::to_string(lineNo) + ": " + msg;
        return false;
    };
    while (std::getline(in, line)) {
        ++lineNo;
        auto t = tokenize(line);
        if (t.empty()) continue;
        const std::string& k = t[0];
        auto need = [&](size_t n) { return t.size() >= n; };
        auto d = [&](size_t i, double def = 0.0) { return i < t.size() ? std::strtod(t[i].c_str(), nullptr) : def; };
        auto idx = [&](size_t i) { return i < t.size() ? std::atoi(t[i].c_str()) - 1 : -1; };
        if (k == "FEM2D") continue;
        if (k == "TITLE") m.title = t.size() > 1 ? t[1] : "";
        else if (k == "PLANE")
            m.planeMode = t.size() < 2 ? PlaneMode::Stress
                          : t[1] == "STRAIN" ? PlaneMode::Strain
                          : t[1] == "AXISYM" ? PlaneMode::Axisymmetric
                                             : PlaneMode::Stress;
        else if (k == "FLOW") {
            if (!need(8)) return fail("FLOW thiếu tham số");
            FlowDomain& f = m.flow;
            f.active = true;
            f.ax = d(1), f.ay = d(2), f.bx = d(3), f.by = d(4), f.length = d(5);
            f.side = d(6) < 0 ? -1 : 1;
            f.U = d(7);
            f.parabolic = d(8) != 0;
            f.medium = (int)d(9, 0);
            f.rho = d(10, 998.2);
            f.mu = d(11, 1.002e-3);
            f.slipWalls = d(12) != 0;
            f.h = d(13), f.hFar = d(14);
            f.transient = d(15) != 0;
            f.dt = d(16, 0.01);
            f.endTime = d(17, 5);
            f.section = d(18) != 0 ? 1 : 0;
            f.topWall = d(19) != 0;
            f.hydrostatic = d(20, 1) != 0;
            f.turbulence = std::clamp((int)d(21, 0), 0, 2);
            f.turbIntensity = d(22, 0.05);
            f.freeSurface = d(23) != 0;
            f.waterDepth = d(24);
            f.initDepth = d(25, -1);
            f.gate = d(26);
            f.drainGround = d(27) != 0;
            f.drainRate = d(28, 0.01);
            f.drainMode = d(29) != 0 ? 1 : 0;
            if (f.width() <= 0 || f.length <= 0) return fail("FLOW: kích thước miền chảy phải > 0");
        }
        else if (k == "ANALYSIS") {
            if (t.size() > 2 && t[1] == "PDELTA") m.secondOrder = d(2) != 0;
        }
        else if (k == "MATERIAL") {
            if (!need(4)) return fail("MATERIAL thiếu tham số");
            Material x;
            x.name = t[1];
            x.E = d(2);
            x.nu = d(3);
            x.alpha = d(4, 1.2e-5);
            x.gamma = d(5, 78.5);
            x.kind = (int)d(6, 0);
            x.fd = d(7, 0);
            x.ft = d(8, 0);
            m.materials.push_back(x);
        } else if (k == "SECTION") {
            if (!need(4)) return fail("SECTION thiếu tham số");
            Section x;
            x.name = t[1];
            x.A = d(2);
            x.I = d(3);
            x.h = d(4, 0.0);
            size_t i = 5;
            if (t.size() > 6 && t[5] != "SHAPE" && t[5] != "POLY") {
                x.yTop = d(5);
                x.yBot = d(6);
                i = 7;
            }
            while (i < t.size()) {
                if (t[i] == "SHAPE") {
                    x.shape.kind = (ShapeKind)(int)d(i + 1);
                    for (int k = 0; k < 6; ++k) x.shape.p[k] = d(i + 2 + k);
                    i += 8;
                } else if (t[i] == "POLY" && i + 1 < t.size()) {
                    x.shape.kind = ShapeKind::Polygon;
                    const int n = (int)d(i + 1);
                    Contour c;
                    for (int k = 0; k < n; ++k) c.push_back({d(i + 2 + 2 * k), d(i + 3 + 2 * k)});
                    x.shape.contours.push_back(c);
                    i += 2 + 2 * n;
                } else {
                    ++i;
                }
            }
            m.sections.push_back(x);
        } else if (k == "NODE") {
            if (!need(3)) return fail("NODE thiếu tọa độ");
            m.nodes.push_back({d(1), d(2)});
        } else if (k == "TRUSS" || k == "FRAME") {
            if (!need(5)) return fail(k + " thiếu tham số");
            Element e;
            e.type = k == "TRUSS" ? ElemType::Truss : ElemType::Frame;
            e.n[0] = idx(1);
            e.n[1] = idx(2);
            e.mat = idx(3);
            e.sec = idx(4);
            e.hingeI = d(5) != 0;
            e.hingeJ = d(6) != 0;
            e.kf = d(7);
            e.kG = d(8);
            e.kfCompOnly = d(9) != 0;
            m.elements.push_back(e);
        } else if (k == "TRI3") {
            if (!need(6)) return fail("TRI3 thiếu tham số");
            Element e;
            e.type = ElemType::Tri3;
            for (int i = 0; i < 3; ++i) e.n[i] = idx(1 + i);
            e.mat = idx(4);
            e.t = d(5);
            m.elements.push_back(e);
        } else if (k == "TRI6" || k == "QUAD8") {
            Element e;
            e.type = k == "TRI6" ? ElemType::Tri6 : ElemType::Quad8;
            const int nc = e.nodeCount();
            if (!need(nc + 3)) return fail(k + " thiếu tham số");
            for (int i = 0; i < nc; ++i) e.n[i] = idx(1 + i);
            e.mat = idx(nc + 1);
            e.t = d(nc + 2);
            m.elements.push_back(e);
        } else if (k == "QUAD4") {
            if (!need(7)) return fail("QUAD4 thiếu tham số");
            Element e;
            e.type = ElemType::Quad4;
            for (int i = 0; i < 4; ++i) e.n[i] = idx(1 + i);
            e.mat = idx(5);
            e.t = d(6);
            e.incompatible = d(7, 1) != 0;
            m.elements.push_back(e);
        } else if (k == "SUPPORT") {
            if (!need(5)) return fail("SUPPORT thiếu tham số");
            Support s;
            s.node = idx(1);
            for (int c = 0; c < 3; ++c) {
                s.fix[c] = d(2 + c) != 0;
                s.k[c] = d(5 + c);
            }
            m.supports.push_back(s);
        } else if (k == "LOADCASE") {
            LoadCase lc;
            lc.name = t.size() > 1 ? t[1] : "TH" + std::to_string(m.cases.size() + 1);
            lc.selfWeight = d(2);
            lc.kind = (int)d(3, 0);
            lc.group = (int)d(4, 0);
            lc.gammaF = d(5, 1.0);
            m.cases.push_back(lc);
            cur = &m.cases.back();
        } else if (k == "NLOAD" || k == "ELOAD" || k == "DISP") {
            if (!cur) return fail(k + " phải nằm sau LOADCASE");
            if (k == "NLOAD") {
                if (!need(4)) return fail("NLOAD thiếu tham số");
                cur->nodal.push_back({idx(1), {d(2), d(3), d(4)}});
            } else if (k == "DISP") {
                if (!need(4)) return fail("DISP thiếu tham số");
                PrescribedDisp p;
                p.node = idx(1);
                p.dof = t[2] == "UX" ? UX : t[2] == "UY" ? UY : RZ;
                p.value = d(3);
                cur->disp.push_back(p);
            } else {
                if (!need(4)) return fail("ELOAD thiếu tham số");
                ElemLoad l;
                l.elem = idx(1);
                const std::string& kind = t[2];
                if (kind == "DIST") {
                    l.kind = ElemLoadKind::Distributed;
                    if (!parseDir(t[3], l.dir)) return fail("hướng tải không hợp lệ: " + t[3]);
                    l.q1 = d(4);
                    l.q2 = d(5, l.q1);
                    l.projected = d(6) != 0;
                } else if (kind == "POINT") {
                    l.kind = ElemLoadKind::Point;
                    if (!parseDir(t[3], l.dir)) return fail("hướng tải không hợp lệ: " + t[3]);
                    l.q1 = d(4);
                    l.a = d(5);
                } else if (kind == "MOMENT") {
                    l.kind = ElemLoadKind::Moment;
                    l.q1 = d(3);
                    l.a = d(4);
                } else if (kind == "FIT") {
                    l.kind = ElemLoadKind::LackOfFit;
                    l.q1 = d(3);
                } else if (kind == "TEMP") {
                    l.kind = ElemLoadKind::Temperature;
                    l.q1 = d(3);
                    l.q2 = d(4);
                } else return fail("loại tải không hợp lệ: " + kind);
                cur->elem.push_back(l);
            }
        } else if (k == "COMBO") {
            Combo c;
            c.name = t.size() > 1 ? t[1] : "THop";
            size_t i0 = 2;
            if (t.size() > 2 && t[2] == "AUTO") {
                c.generated = true;
                i0 = 3;
            }
            for (size_t i = i0; i + 1 < t.size(); i += 2) c.terms.push_back({idx(i), d(i + 1)});
            m.combos.push_back(c);
        } else {
            return fail("từ khóa không hợp lệ: " + k);
        }
    }
    out = std::move(m);
    return true;
}

bool saveModel(const std::string& path, const Model& m, std::string& error) {
    std::ofstream f(path, std::ios::binary);
    if (!f) {
        error = "Không ghi được file " + path;
        return false;
    }
    f << modelToText(m);
    return true;
}

bool loadModel(const std::string& path, Model& m, std::string& error) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        error = "Không mở được file " + path;
        return false;
    }
    std::stringstream ss;
    ss << f.rdbuf();
    return modelFromText(ss.str(), m, error);
}

std::string resultReport(const Model& m, const AnalysisResult& r) {
    std::ostringstream o;
    char b[512];
    o << "KẾT QUẢ TÍNH TOÁN — " << m.title << "\n";
    o << "Đơn vị: kN, m; chuyển vị mm; góc xoay rad; ứng suất MPa\n";
    if (!r.ok) {
        o << "LỖI: " << r.error << "\n";
        return o.str();
    }
    std::snprintf(b, sizeof b, "Số nút %zu, số phần tử %zu, số phương trình %d, profile %lld\n", m.nodes.size(),
                  m.elements.size(), r.numEquations, r.profileSize);
    o << b;
    for (const auto& w : r.warnings) o << "Cảnh báo: " << w << "\n";
    for (const auto& s : r.sets) {
        o << "\n==================== " << (s.combo ? "Tổ hợp: " : "Trường hợp tải: ") << s.name << "\n";
        std::snprintf(b, sizeof b, "Tổng tải: ΣFx = %.4f, ΣFy = %.4f | Tổng phản lực: ΣRx = %.4f, ΣRy = %.4f\n",
                      s.sumFx, s.sumFy, s.sumRx, s.sumRy);
        o << b << "\nCHUYỂN VỊ NÚT\n   Nút        ux(mm)        uy(mm)       φz(rad)\n";
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            std::snprintf(b, sizeof b, "%6zu %13.5f %13.5f %13.6e\n", i + 1, s.U[3 * i] * 1e3, s.U[3 * i + 1] * 1e3,
                          s.U[3 * i + 2]);
            o << b;
        }
        o << "\nPHẢN LỰC\n   Nút         Rx(kN)         Ry(kN)       Mz(kNm)\n";
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            bool any = false;
            for (int c = 0; c < 3; ++c) any |= r.dofState[3 * i + c] == DofSupport || s.R[3 * i + c] != 0;
            if (!any) continue;
            std::snprintf(b, sizeof b, "%6zu %14.4f %14.4f %14.4f\n", i + 1, s.R[3 * i], s.R[3 * i + 1], s.R[3 * i + 2]);
            o << b;
        }
        bool header = false;
        for (size_t e = 0; e < m.elements.size(); ++e) {
            if (!m.elements[e].isBar()) continue;
            const BarResult& br = s.bars[e];
            if (br.x.empty()) continue;
            if (!header) {
                o << "\nNỘI LỰC THANH (N > 0 kéo; M > 0 căng thớ dưới)\n"
                     "   PT  Đầu       N(kN)       Q(kN)      M(kNm)   |  Mmax(kNm) tại x(m)   Mmin(kNm) tại x(m)\n";
                header = true;
            }
            const size_t n = br.x.size();
            const size_t imax = std::max_element(br.M.begin(), br.M.end()) - br.M.begin();
            const size_t imin = std::min_element(br.M.begin(), br.M.end()) - br.M.begin();
            std::snprintf(b, sizeof b, "%5zu   i  %11.4f %11.4f %11.4f   | %10.4f %7.3f   %10.4f %7.3f\n", e + 1, br.N[0],
                          br.Q[0], br.M[0], br.M[imax], br.x[imax], br.M[imin], br.x[imin]);
            o << b;
            std::snprintf(b, sizeof b, "        j  %11.4f %11.4f %11.4f\n", br.N[n - 1], br.Q[n - 1], br.M[n - 1]);
            o << b;
        }
        header = false;
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            if (s.nodeStressCount.empty() || s.nodeStressCount[i] == 0) continue;
            if (!header) {
                o << (m.planeMode == PlaneMode::Axisymmetric
                          ? "\nỨNG SUẤT TẠI NÚT (trung bình, MPa)\n   Nút        σr        σz       τrz        σθ        σ1        σ2  σ_vonMises\n"
                          : "\nỨNG SUẤT TẠI NÚT (trung bình, MPa)\n   Nút        σx        σy       τxy        σz        σ1        σ2  σ_vonMises\n");
                header = true;
            }
            const auto& v = s.nodeStress[i];
            StressInfo si = stressInfo(v[0], v[1], v[2], v[3]);
            std::snprintf(b, sizeof b, "%6zu %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %11.4f\n", i + 1, si.sx * 1e-3, si.sy * 1e-3,
                          si.txy * 1e-3, si.sz * 1e-3, si.s1 * 1e-3, si.s2 * 1e-3, si.vonMises * 1e-3);
            o << b;
        }
    }
    return o.str();
}

} // namespace fem
