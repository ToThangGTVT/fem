// Mô phỏng dòng chảy chất lỏng không nén 2D bằng phương pháp thể tích hữu hạn (FVM).
//   • Phương trình Navier–Stokes không nén: ρ(∂u/∂t + u·∇u) = −∇p + μ∇²u,  ∇·u = 0.
//   • Lưới phi cấu trúc: mỗi phần tử phẳng của mô hình (T3, Q4, T6, Q8 — chỉ lấy các nút góc) là một thể tích kiểm tra,
//     ẩn số đặt tại tâm ô (collocated). Mặt là các cạnh nối hai nút góc.
//   • Ghép áp suất – vận tốc theo SIMPLE, nội suy lưu lượng qua mặt theo Rhie–Chow (kèm hiệu chỉnh của Majumdar cho hệ số
//     giảm dư và của Choi cho số hạng thời gian để nghiệm không phụ thuộc α, Δt).
//   • Đối lưu: ngược dòng bậc nhất hoặc ngược dòng tuyến tính bậc hai (hiệu chỉnh trễ); khuếch tán: phần trực giao ẩn,
//     phần không trực giao hiệu chỉnh trễ (over-relaxed); gradient bình phương tối thiểu.
//   • Dừng (lặp SIMPLE đến hội tụ) hoặc không dừng (Euler ẩn / BDF2, vài vòng lặp ngoài mỗi bước).
//   • Chảy rối: mô hình k-ε chuẩn (Launder–Spalding) với hàm tường chuẩn, sinh rối theo Kato–Launder (tránh k quá lớn
//     tại điểm dừng trước vật cản).
//   • Mặt thoáng: hai pha nước – không khí bằng phương pháp thể tích chất lỏng (VOF), phân số thể tích α vận chuyển hiện
//     với sơ đồ nén mặt phân cách HRIC; áp suất giải dưới dạng p_rgh = p − ρ·g·x, trọng lực hướng −y.
// Ẩn số lưu lượng qua mặt là lưu lượng thể tích (m²/s cho 1 m bề dày); lưu lượng khối = ρ_mặt × lưu lượng thể tích.
// Đơn vị SI: m, s, kg/m³, Pa·s, Pa; lực trên vật cản tính cho 1 m bề dày (N/m).
#pragma once
#include <array>
#include <string>
#include <vector>

#include "core/Model.h"

namespace fem {

// Các đoạn biên được nhận diện theo hình học: bốn cạnh của hình chữ nhật bao miền tính và biên trong
// (vật cản, bậc, lỗ — mọi cạnh biên không nằm trên hình chữ nhật bao).
enum class FluidSide { Left = 0, Right, Bottom, Top, Inner, Count };
// Opening: biên thông khí (áp suất toàn phần p = 0, chất lỏng / không khí ra vào tự do) — biên trên của miền có mặt thoáng.
enum class FluidBcKind { Wall, MovingWall, Inlet, Outlet, Symmetry, Opening };

struct FluidBc {
    FluidBcKind kind = FluidBcKind::Wall;
    double u = 0, v = 0;     // vận tốc (m/s): thành di động, cửa vào (cực đại nếu phân bố parabol)
    bool parabolic = false;  // cửa vào: phân bố parabol 4s(1 − s) dọc đoạn biên
};

enum class ConvScheme { Upwind, LinearUpwind };
enum class TimeScheme { Euler, BDF2 };
enum class TurbModel { Laminar, KEpsilon };

struct FluidParams {
    double rho = 1000;    // khối lượng riêng (kg/m³)
    double mu = 1e-3;     // độ nhớt động lực (Pa·s)
    std::array<FluidBc, (int)FluidSide::Count> bc{};
    ConvScheme scheme = ConvScheme::LinearUpwind;
    double alphaU = 0.7;  // hệ số giảm dư vận tốc
    double alphaP = 0.3;  // hệ số giảm dư áp suất
    double tol = 1e-6;    // dừng: sai số dư chuẩn hóa
    bool transient = false;
    TimeScheme timeScheme = TimeScheme::BDF2;
    double dt = 0.01;     // bước thời gian (s)
    double endTime = 0;   // không dừng: thời gian kết thúc đề xuất (0: người dùng tự chọn)
    int outerIters = 20;  // số vòng lặp ngoài tối đa trong một bước thời gian
    double outerTol = 1e-4;
    double Uref = 0, Lref = 0; // vận tốc, chiều dài tham chiếu (0: tự xác định từ điều kiện biên và vật cản)
    // hướng chảy (rad): trục x' của hệ tọa độ dòng chảy. Các đoạn biên trái/phải/dưới/trên là các cạnh của hình
    // chữ nhật bao miền trong hệ x'y' (trái: cửa vào phía thượng lưu); Cd, Cl tính theo phương x', y'.
    double flowAngle = 0;
    // áp suất thủy tĩnh (mặt cắt đứng): hydroG > 0 thì áp suất toàn phần = p + ρ·g·(yFree − y), dùng cho lực lên vật cản.
    // Chất lỏng không nén một pha: trọng lực chỉ thêm phần thủy tĩnh, không làm đổi trường vận tốc.
    double hydroG = 0, yFree = 0;
    // chảy rối
    TurbModel turbulence = TurbModel::Laminar;
    double turbIntensity = 0.05; // cường độ rối ở cửa vào I: k = 1,5·(I·U)²
    double turbLength = 0;       // thang chiều dài rối ở cửa vào ℓ (m), ε = Cμ^0,75·k^1,5/ℓ; 0: 0,07 × bề rộng cửa vào
    // hai pha nước – không khí (VOF): rho, mu ở trên là của nước (pha α = 1)
    bool vof = false;
    double rho2 = 1.204, mu2 = 1.825e-5; // pha thứ hai (không khí)
    double gravity = 9.81;
    double inletLevel = 0;  // cao độ mặt nước ở cửa vào (y tổng thể): phần cửa vào cao hơn là thành
    double initLevel = 0;   // cao độ mặt nước ban đầu
    double initGate = 0;    // nước ban đầu chỉ trong đoạn dài này tính từ cửa vào theo hướng chảy (vỡ đập); 0: toàn miền
    double maxCourant = 0.5; // bước thời gian tự điều chỉnh để số Courant không vượt giá trị này (dt là bước lớn nhất)
};

// Lưới thể tích hữu hạn dựng từ các phần tử phẳng. Mặt trong đứng trước (sắp theo owner < neigh), sau đó là mặt biên.
struct FluidMesh {
    int nCells = 0, nFaces = 0, nInternal = 0;
    std::vector<int> cellElem;             // phần tử của mô hình ứng với ô
    std::vector<int> cellStart, cellNodes; // các nút góc của ô (ngược chiều kim đồng hồ), dạng CSR
    std::vector<double> px, py;            // tọa độ các nút của mô hình
    std::vector<double> cx, cy, vol;       // tâm, diện tích ô (bề dày 1 m)
    std::vector<int> owner, neigh;         // neigh = −1: mặt biên
    std::vector<int> fa, fb;               // hai nút của mặt (a → b ngược chiều kim đồng hồ quanh owner)
    std::vector<double> sx, sy;            // vectơ diện tích mặt, hướng ra khỏi owner, độ lớn = chiều dài cạnh
    std::vector<double> fx, fy;            // trung điểm mặt
    std::vector<double> w;                 // trọng số nội suy phía owner: φf = w·φP + (1 − w)·φN
    std::vector<int> side;                 // mặt biên: FluidSide
    std::vector<int> cfStart, cfFace;      // các mặt của từng ô, dạng CSR
    double x0 = 0, y0 = 0, x1 = 0, y1 = 0; // hình chữ nhật bao (tọa độ tổng thể)
    double ca = 1, sa = 0;                 // cos, sin của hướng chảy: x' = ca·x + sa·y, y' = −sa·x + ca·y
    int sideFaces[(int)FluidSide::Count] = {0, 0, 0, 0, 0};
};
// Dựng lưới từ các phần tử phẳng của mô hình. Trả về chuỗi lỗi (rỗng nếu thành công).
std::string buildFluidMesh(const Model& m, FluidMesh& mesh, double flowAngle = 0);

struct FluidForces {
    bool valid = false;   // có biên trong (vật cản)
    double Fx = 0, Fy = 0;        // tổng lực chất lỏng tác dụng lên biên trong (N/m), theo trục tổng thể
    double Fd = 0, Fl = 0;        // lực cản (theo hướng chảy) và lực nâng (vuông góc, sang trái hướng chảy)
    double Fpx = 0, Fpy = 0;      // phần do áp suất
    double Fvx = 0, Fvy = 0;      // phần do ma sát nhớt
    double Fhx = 0, Fhy = 0;      // phần do áp suất thủy tĩnh (đã gồm trong Fx, Fy; lực đẩy nổi, áp lực tĩnh)
    double Cd = 0, Cl = 0;        // hệ số cản, hệ số nâng theo Uref, Lref — chỉ phần động (không gồm thủy tĩnh)
};

struct FluidHistory {
    int iter = 0;
    double time = 0;
    double resU = 0, resV = 0, resC = 0;
    double Cd = 0, Cl = 0;
};

enum class FluidField { Speed, U, V, P, Cp, Vorticity, PTotal, Alpha, TurbK };

class FluidSolver {
public:
    // Khởi tạo từ mô hình (lưới) và thông số. Trả về false nếu lỗi, xem error().
    bool init(const Model& m, const FluidParams& p);
    const std::string& error() const { return err_; }
    const FluidMesh& mesh() const { return M_; }
    const FluidParams& params() const { return P_; }

    // Dừng: một vòng lặp SIMPLE. Không dừng: một bước thời gian. Trả về false nếu nghiệm phân kỳ.
    bool step();
    bool converged() const { return converged_; }
    bool diverged() const { return diverged_; }
    int iteration() const { return iter_; }
    double time() const { return time_; }
    const std::vector<FluidHistory>& history() const { return hist_; }

    double Uref() const { return Uref_; }
    double Lref() const { return Lref_; }
    double reynolds() const { return P_.rho * Uref_ * Lref_ / P_.mu; }
    double dt() const { return dt_; } // bước thời gian vừa dùng (VOF: tự điều chỉnh theo số Courant)
    double strouhal() const; // từ dao động của Cl theo thời gian (0 nếu chưa xác định được)

    // trường tại tâm ô
    const std::vector<double>& u() const { return u_; }
    const std::vector<double>& v() const { return v_; }
    const std::vector<double>& p() const { return p_; }        // áp suất (VOF: p_rgh = p − ρ·g·x)
    const std::vector<double>& alpha() const { return alpha_; } // VOF: phân số thể tích nước
    const std::vector<double>& turbK() const { return k_; }
    const std::vector<double>& turbEps() const { return eps_; }
    const std::vector<double>& turbVisc() const { return mut_; } // độ nhớt rối μt (Pa·s)
    double waterVolume() const;                                  // VOF: Σ α·V (m² cho 1 m bề dày)
    double pressureAt(int cell) const;                          // áp suất thực (kể cả phần thủy tĩnh) tại tâm ô
    std::vector<double> cellField(FluidField f) const;
    // nội suy về các nút của mô hình (NaN cho nút không thuộc miền); nút trên thành cứng / cửa vào lấy giá trị biên
    std::vector<double> nodeField(FluidField f, int nodeCount) const;
    // giá trị tại điểm bất kỳ (tái dựng tuyến tính trong ô chứa điểm); false nếu điểm nằm ngoài miền
    bool probe(double x, double y, double& u, double& v, double& p) const;
    int locate(double x, double y, int hint = -1) const;

    FluidForces forces() const;
    // Hệ số áp suất dọc biên trong: (góc quanh trọng tâm vật cản (độ), Cp), sắp theo góc
    std::vector<std::array<double, 2>> innerCp() const;
    // Đường dòng xuất phát từ các điểm cho trước (tích phân hai chiều)
    std::vector<std::vector<std::array<double, 2>>> streamlines(const std::vector<std::array<double, 2>>& seeds,
                                                               int maxSteps = 3000) const;
    std::vector<std::array<double, 2>> defaultSeeds(int n) const;
    double massImbalance() const; // tổng |∇·u| tương đối (kiểm tra bảo toàn khối lượng)
    double hydrostatic(double y) const { return P_.hydroG > 0 ? P_.rho * P_.hydroG * (P_.yFree - y) : 0; }

private:
    void gradient(const std::vector<double>& phi, const std::vector<double>& faceVal, std::vector<double>& gx,
                  std::vector<double>& gy) const;
    // Green–Gauss có hiệu chỉnh độ lệch (giá trị mặt φf = φ̄f + ∇φ_LSQ·r): dạng bảo toàn, dùng cho áp suất
    void gradientGG(const std::vector<double>& phi, const std::vector<double>& faceVal, std::vector<double>& gx,
                    std::vector<double>& gy) const;
    void boundaryVelocity(std::vector<double>& ub, std::vector<double>& vb) const;
    // ρ·u_f·S trên các mặt trong, u_f nội suy tuyến tính có hiệu chỉnh độ lệch (skewness) về trung điểm mặt
    void faceFlux(std::vector<double>& out) const;
    void boundaryPressure(const std::vector<double>& p, const std::vector<double>& gx, const std::vector<double>& gy,
                          std::vector<double>& pb) const;
    double inletProfile(int f) const;
    void momentum(double& resU, double& resV);
    double pressureCorrection(const std::vector<double>& majumdar);
    bool simpleIteration(double& resU, double& resV, double& resC);
    void gaussSeidel(std::vector<double>& x, const std::vector<double>& b, int sweeps) const;
    int pcg(std::vector<double>& x, const std::vector<double>& b, double relTol, int maxIt) const;
    void record(double resU, double resV, double resC);
    bool fixedPressure(int f) const;  // mặt biên có áp suất cho trước (cửa ra một pha, biên thông khí)
    double fixedPressureValue(int f) const;
    void updateProperties();          // ρ, μ của ô và của mặt theo α (VOF)
    void advectAlpha(double dt);      // vận chuyển α (HRIC, chia bước con theo số Courant)
    void boundaryAlpha(std::vector<double>& ab) const;
    void solveTurbulence();           // k-ε và hàm tường
    void initTurbulence();
    double courantRate() const;       // max_ô Σ|lưu lượng ra| / V (1/s)
    void pressureSource();            // ∇p (VOF: ∇p_rgh + g·x·∇ρ) cho phương trình động lượng

    FluidMesh M_;
    FluidParams P_;
    std::string err_;
    bool hasOutlet_ = false, converged_ = false, diverged_ = false;
    int iter_ = 0;
    double time_ = 0, Uref_ = 1, Lref_ = 1;
    double sideMin_[(int)FluidSide::Count]{}, sideMax_[(int)FluidSide::Count]{}; // phạm vi đoạn biên (cho parabol)

    std::vector<double> u_, v_, p_, F_;         // vận tốc, áp suất tại tâm ô; lưu lượng thể tích qua mặt (m²/s)
    std::vector<double> un_, vn_, Fn_;          // bước thời gian trước; Fn_ = F − ū·S (hiệu chỉnh Choi)
    std::vector<double> unn_, vnn_, Fnn_;       // hai bước trước (BDF2)
    int nOld_ = 0;                              // số bước thời gian cũ đã có
    std::vector<double> gpx_, gpy_;             // nguồn áp suất của động lượng: ∇p (VOF: ∇p_rgh + g·x·∇ρ)
    std::vector<double> pgx_, pgy_;             // ∇p (ngoại suy áp suất lên biên, giá trị tại điểm)
    std::vector<double> rhoC_, muC_, rhoF_;     // ρ, μ tại tâm ô; ρ trên mặt (lưu lượng khối = ρ_mặt·F)
    std::vector<double> alpha_, alphaF_;        // VOF: α tại ô, α trung bình trên mặt trong bước vừa qua
    std::vector<double> k_, eps_, mut_;         // k-ε: động năng rối, tốc độ tiêu tán, độ nhớt rối
    std::vector<double> kn_, en_;               // k, ε của bước thời gian trước
    std::vector<double> muWall_;                // mặt biên: độ nhớt hiệu dụng trên thành (hàm tường)
    double kIn_ = 0, epsIn_ = 0;                // k, ε ở cửa vào
    double dt_ = 0;                             // bước thời gian hiện tại
    std::vector<double> skx_, sky_;             // mặt trong: trung điểm mặt − điểm nội suy trên đoạn nối hai tâm ô
    std::vector<double> lsq_;                   // gradient bình phương tối thiểu: ma trận nghịch đảo 2×2 của từng ô
    std::vector<double> recon_;                 // dựng véc tơ ô từ thành phần pháp tuyến trên mặt: (Σ S·Sᵀ/|S|)⁻¹
    std::vector<double> diag_, upper_, lower_;  // ma trận: đường chéo, hệ số hàng owner (cột neigh), hàng neigh (cột owner)
    std::vector<double> aP_;                    // đường chéo động lượng (đã chia α) cho Rhie–Chow
    std::vector<double> at_;                    // phần đường chéo do số hạng thời gian
    double contNorm_ = 0;
    std::vector<FluidHistory> hist_;
};

// Lập bài toán từ miền chảy khai báo trong mô hình (m.flow) và các vật cản là vòng thanh khép kín vẽ bên trong miền:
// tạo lưới tam giác vào mesh (mô hình riêng, không đụng tới m), đặt tính chất chất lỏng, điều kiện biên và hướng chảy
// vào p (giữ nguyên sơ đồ đối lưu, hệ số giảm dư của p). Trả về chuỗi lỗi (rỗng nếu thành công); note: ghi chú.
std::string buildFlowProblem(const Model& m, Model& mesh, FluidParams& p, std::string* note = nullptr);

// Kênh chữ nhật [0, L]×[0, H] (góc dưới trái tại x0, y0) có thể chứa một vật cản.
struct ChannelSpec {
    double x0 = 0, y0 = 0, L = 2.2, H = 0.41;
    int obstacle = 1;              // 0: không; 1: trụ tròn; 2: hình chữ nhật (được phép chạm biên: bậc, công trình)
    double xc = 0.2, yc = 0.2;     // tâm vật cản
    double D = 0.1;                // đường kính trụ / bề rộng hình chữ nhật
    double B = 0.1;                // chiều cao hình chữ nhật
    double h = 0.01;               // kích thước ô gần vật cản (m)
    double hFar = 0.04;            // kích thước ô xa vật cản (m)
    bool triangles = false;        // chia mỗi ô tứ giác thành hai tam giác
};
// Tạo lưới phần tử phẳng (Q4 hoặc T3) cho miền chất lỏng: trụ tròn dùng lưới O quanh trụ, các vùng còn lại là lưới
// chữ nhật giãn dần theo cấp số nhân từ vật cản ra xa. Thêm vào mô hình m với vật liệu mat. Trả về số phần tử đã tạo.
int meshChannel(Model& m, const ChannelSpec& s, int mat);

const char* fluidSideName(FluidSide s);
const char* fluidBcName(FluidBcKind k);
const char* fluidFieldName(FluidField f);

} // namespace fem
