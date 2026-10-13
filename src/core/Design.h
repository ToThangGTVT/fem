// Kiểm tra thiết kế sơ bộ từ nội lực:
//   • Thép (TCVN 5575:2012): bền, ổn định khi nén (hệ số φ theo độ mảnh quy ước λ̄), cắt, độ mảnh giới hạn.
//   • Bê tông cốt thép (TCVN 5574:2018): cốt thép dọc dầm (αm, ξR), cốt đai (c = 2h0),
//     cột chữ nhật nén lệch tâm đặt thép đối xứng (kể uốn dọc η).
//   • Độ võng: f ≤ L / n.
// Đây là các công thức cơ bản của tiêu chuẩn; kết quả dùng để chọn sơ bộ, cần kiểm tra lại đầy đủ theo tiêu chuẩn.
#pragma once
#include <string>
#include <vector>

#include "core/Analysis.h"

namespace fem {

// Cấu kiện: chuỗi thanh khung thẳng hàng nối liên tiếp (không có gối, thanh khác hay khớp ở nút giữa)
struct Member {
    std::vector<int> elems;   // theo thứ tự từ nút i đến nút j của cấu kiện
    std::vector<bool> flip;   // phần tử ngược chiều với cấu kiện
    int nodeI = -1, nodeJ = -1;
    double L = 0;
    bool horizontal = false;  // gần nằm ngang (|sin| < 0,35): dầm
};
std::vector<Member> findMembers(const Model& m);

// ---------------------------------------------------------------- thép
double steelPhi(double lambdaBar, double fOverE); // hệ số uốn dọc φ (TCVN 5575:2012)

struct SteelParams {
    double gammaC = 1.0;
    double mu = 1.0;                 // hệ số chiều dài tính toán trong mặt phẳng khung (mặc định)
    std::vector<double> muMember;    // theo cấu kiện (0: dùng mu)
    double limitCompression = 120, limitTension = 400;
};
struct SteelCheck {
    int member = -1;
    double N = 0, M = 0, V = 0;      // nội lực ứng với hệ số sử dụng lớn nhất
    double lambda = 0, lambdaBar = 0, phi = 1;
    double strength = 0, stability = 0, shear = 0, slender = 0; // tỉ số sử dụng
    double util = 0;                 // lớn nhất
    std::string set;                 // bộ tải quyết định
    std::string note;
};
std::vector<SteelCheck> checkSteel(const Model& m, const AnalysisResult& r, const std::vector<int>& sets,
                                   const std::vector<Member>& members, const SteelParams& p);

// ---------------------------------------------------------------- bê tông cốt thép
struct RcParams {
    double Rs = 350, Rsc = 350, Rsw = 170;   // MPa (CB400-V, cốt đai CB240-T)
    double Es = 2e5;                          // MPa
    double a = 0.04;                          // khoảng cách trọng tâm cốt thép tới mép (m)
    double Asw = 2 * 0.503e-4;                // diện tích một lớp cốt đai (2 nhánh Φ8), m²
    double mu = 1.0;                          // hệ số chiều dài tính toán cột
    double phiL = 1.0;                        // hệ số kể tải dài hạn
    double minRatio = 0.001;                  // hàm lượng tối thiểu
};
double rcXiR(double Rs, double Es);         // ξR = 0,8 / (1 + εs,el / εb2), εb2 = 0,0035

struct RcBeamDesign {
    double As = 0, AsComp = 0;  // cm²: cốt chịu kéo, cốt chịu nén
    double M = 0;               // kNm
    double alphaM = 0, xi = 0;
    bool overReinforced = false;
};
// M > 0: căng thớ dưới; trả về cốt dưới (As) và cốt nén trên nếu cần
RcBeamDesign rcBeamFlexure(double M, double b, double h, const Material& c, const RcParams& p);

struct RcShearDesign {
    double Q = 0, Qmax = 0;     // kN; Qmax = 0,3·Rb·b·h0
    double qswReq = 0;          // kN/m
    double spacing = 0;         // m (bước đai cần thiết, 0: theo cấu tạo)
    bool ok = true;
};
RcShearDesign rcShear(double Q, double b, double h, const Material& c, const RcParams& p);

struct RcColumnDesign {
    double N = 0, M = 0;        // kN, kNm (N > 0 nén)
    double e0 = 0, eta = 1, x = 0;
    double As = 0;              // cm², mỗi phía (đối xứng)
    bool unstable = false;
};
RcColumnDesign rcColumnSym(double N, double M, double b, double h, double l0, const Material& c, const RcParams& p);

struct RcMemberResult {
    int member = -1;
    bool column = false;
    double b = 0, h = 0;
    double AsBottom = 0, AsTop = 0, AsSide = 0; // cm²
    double Mpos = 0, Mneg = 0, Nmax = 0, Qmax = 0;
    RcShearDesign shear;
    std::string note;
};
std::vector<RcMemberResult> designRc(const Model& m, const AnalysisResult& r, const std::vector<int>& sets,
                                     const std::vector<Member>& members, const RcParams& p);

// ---------------------------------------------------------------- độ võng
struct DeflectionCheck {
    int member = -1;
    double L = 0, f = 0, limit = 0, ratio = 0; // m
    double xAt = 0;
};
std::vector<DeflectionCheck> checkDeflection(const Model& m, const LoadSetResult& s, const std::vector<Member>& members,
                                             double n, bool beamsOnly = true);

// ---------------------------------------------------------------- cọc chịu tải ngang (phương pháp "m")
struct PileParams {
    double L = 20;        // chiều dài cọc trong đất, m
    double l0 = 0;        // chiều dài tự do trên mặt đất, m
    double d = 0.6;       // đường kính (hoặc cạnh) cọc, m
    bool square = false;
    double m = 5000;      // hệ số tỉ lệ m, kN/m⁴
    int segments = 40;
    double H = 100, Mh = 0; // tải ngang và mômen tại đầu cọc
    int mat = 0;
};
double pileConventionalWidth(double d);
void makePile(Model& m, const PileParams& p);

} // namespace fem
