// Bài toán trị riêng: ổn định tuyến tính (tải trọng tới hạn) và dao động riêng.
//   Ổn định:   (K + λ·K_G) φ = 0, K_G ứng với lực dọc của bộ tải đã chọn; λ_cr = hệ số tải trọng tới hạn.
//   Dao động:  (K − ω²·M) φ = 0, khối lượng từ trọng lượng bản thân (γ/g) và/hoặc từ tải trọng thẳng đứng.
// Mỗi thanh được chia nhỏ bên trong để dạng riêng và trị riêng chính xác.
#pragma once
#include <string>
#include <vector>

#include "core/Model.h"

namespace fem {

enum class EigenKind { Buckling = 0, Vibration = 1, Mechanism = 2 };

struct EigenOptions {
    EigenKind kind = EigenKind::Buckling;
    int loadSet = 0;            // ổn định: bộ tải gây lực dọc (chỉ số: trường hợp tải rồi tổ hợp)
    int modes = 6;
    int subdiv = 8;             // số đoạn chia mỗi thanh
    double selfMassFactor = 1;  // dao động: hệ số khối lượng bản thân (0 = bỏ qua)
    int massSet = -1;           // dao động: bộ tải quy thành khối lượng m = |Fy|/g (−1 = không dùng)
    double massSetFactor = 1;
};

struct EigenMode {
    double lambda = 0;              // ổn định: λ_cr; dao động: ω² (rad²/s²)
    double omega = 0, freq = 0, period = 0;
    std::vector<double> U;          // dạng riêng trên mô hình đã chia nhỏ (3 bậc tự do/nút), max |u| = 1
};

struct EigenResult {
    bool ok = false;
    std::string error;
    std::vector<std::string> warnings;
    EigenKind kind = EigenKind::Buckling;
    Model sub;                              // mô hình đã chia nhỏ (nút gốc giữ nguyên chỉ số)
    std::vector<std::vector<int>> barChain; // chuỗi nút dọc theo mỗi thanh gốc (rỗng với phần tử phẳng)
    std::vector<EigenMode> modes;
    std::vector<double> memberN;            // ổn định: lực dọc trung bình của thanh gốc dưới bộ tải (kN)
    int numEquations = 0;
    double timeMs = 0;
};

EigenResult eigenAnalysis(const Model& m, const EigenOptions& opt);

// Hệ số chiều dài tính toán μ = (π/L)·√(EI/N_cr) của thanh e ở dạng mất ổn định thứ nhất (0 nếu thanh không chịu nén).
double effectiveLengthFactor(const Model& m, const EigenResult& r, int e);

// Trị riêng, véc tơ riêng của ma trận đối xứng n×n (lưu theo hàng): A = V·diag(d)·Vᵀ, d tăng dần,
// cột j của V là véc tơ riêng thứ j.
void symmetricEigen(int n, std::vector<double>& A, std::vector<double>& d, std::vector<double>& V);

} // namespace fem
