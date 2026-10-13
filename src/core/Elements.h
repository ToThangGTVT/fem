// Ma trận độ cứng, véc tơ tải tương đương và nội lực của phần tử.
//
// Quy ước dấu nội lực thanh (theo giáo trình Cơ học kết cấu VN / Строительная механика):
//   N > 0: kéo.
//   Q > 0: làm phần tách ra quay thuận chiều kim đồng hồ.
//   M > 0: căng thớ dưới (thớ phía −y cục bộ). Biểu đồ M được vẽ về phía thớ căng.
// Trục cục bộ thanh: x từ nút i → j, y quay x ngược chiều kim đồng hồ 90°.
#pragma once
#include <array>
#include <vector>

#include "core/Model.h"

namespace fem {

using Mat6 = std::array<std::array<double, 6>, 6>;
using Vec6 = std::array<double, 6>;

struct BarProps {
    double L = 0, c = 1, s = 0;
    double EA = 0, EI = 0;
    double alpha = 0, h = 0, kf = 0;
    bool truss = false;
    bool relI = false, relJ = false; // giải phóng mômen đầu thanh
    bool hasBending = true;          // false: thanh dàn có I = 0 (chỉ chịu kéo nén)
    double kG = 0;                   // nền Pasternak (kN)
    bool second = false;             // phân tích bậc hai: thêm độ cứng hình học với lực dọc N2
    double N2 = 0;
};

// Tải trên thanh đã quy về hệ trục cục bộ.
struct BarLoads {
    double qx1 = 0, qx2 = 0, qy1 = 0, qy2 = 0; // phân bố hình thang trên toàn chiều dài
    struct Point { double a, Px, Py; };
    struct Moment { double a, M; };
    std::vector<Point> points;
    std::vector<Moment> moments;
    double T = 0;  // độ tăng nhiệt trung bình
    double dT = 0;   // t_dưới − t_trên
    double eps0 = 0; // biến dạng ban đầu do sai số chế tạo ΔL/L
    bool empty() const {
        return qx1 == 0 && qx2 == 0 && qy1 == 0 && qy2 == 0 && points.empty() && moments.empty() && T == 0 && dT == 0 &&
               eps0 == 0;
    }
};

struct BarResult {
    Vec6 endForces{}; // lực đầu thanh (nút tác dụng lên thanh), hệ cục bộ: Ni, Vi, Mi, Nj, Vj, Mj
    Vec6 uLocal{};    // chuyển vị đầu thanh cục bộ (kể cả góc xoay đầu thanh tại khớp)
    std::vector<double> x, N, Q, M, u, v; // các điểm lấy mẫu dọc thanh; u, v: chuyển vị cục bộ
};

BarProps barProps(const Model& m, int e);
void addBarLoad(const BarProps& p, const ElemLoad& l, double factor, BarLoads& out);
void addBarSelfWeight(const BarProps& p, double weightPerLength, BarLoads& out);

// đầy đủ, chưa giải phóng liên kết; gồm nền Winkler/Pasternak và độ cứng hình học khi p.second
void barLocalStiffness(const BarProps& p, Mat6& K);
void barFixedEndLoads(const BarProps& p, const BarLoads& l, Vec6& f); // tải nút tương đương, hai đầu ngàm
// Ma trận độ cứng và tải tương đương trong hệ tổng thể, đã kể khớp đầu thanh.
void barGlobal(const BarProps& p, const BarLoads& l, Mat6& Kg, Vec6& fg);
// Tính nội lực, chuyển vị dọc thanh từ chuyển vị nút tổng thể.
void barRecover(const BarProps& p, const BarLoads& l, const Vec6& uGlobal, int samples, BarResult& out);

// Ma trận dùng cho bài toán trị riêng (hệ tổng thể, đã kể khớp đầu thanh bằng phép biến đổi khử tĩnh):
//   K  — độ cứng; KG — độ cứng hình học ứng với lực dọc N (N > 0 kéo); M — khối lượng tương thích, ρA (t/m).
void barEigenMatrices(const BarProps& p, double N, double rhoA, Mat6& K, Mat6& KG, Mat6& M);

// Nội lực tại x của thanh đã tính (dùng lực đầu thanh và chuyển vị trong r)
void barForcesAt(const BarProps& p, const BarLoads& l, const BarResult& r, double x, bool after, double& N, double& Q,
                 double& M);

// Hợp lực phản lực nền Winkler tác dụng lên thanh, hệ tổng thể (Fx, Fy).
void barFoundationForce(const BarProps& p, const BarResult& r, double& Fx, double& Fy);

// Phần tử phẳng
void planeD(const Material& mat, PlaneMode mode, double D[3][3]);
bool tri3Stiffness(const double x[3], const double y[3], const double D[3][3], double t, double K[6][6]);
void tri3Stress(const double x[3], const double y[3], const double D[3][3], const double u[6], double sig[3]);
bool quad4Stiffness(const double x[4], const double y[4], const double D[3][3], double t, bool incompatible,
                    double K[8][8]);
// Ứng suất (σx, σy, τxy) tại 4 nút, ngoại suy từ điểm Gauss 2×2.
void quad4Stress(const double x[4], const double y[4], const double D[3][3], double t, bool incompatible,
                 const double u[8], double sig[4][3]);
// Lực nút tương đương của lực thể tích (bx, by) [kN/m³]
void tri3BodyForce(const double x[3], const double y[3], double t, double bx, double by, double f[6]);
void quad4BodyForce(const double x[4], const double y[4], double t, double bx, double by, double f[8]);

// Các đại lượng ứng suất dẫn xuất
struct StressInfo {
    double sx, sy, txy, sz, s1, s2, angle, vonMises;
};
StressInfo stressInfo(double sx, double sy, double txy, double sz);

} // namespace fem
