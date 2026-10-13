// Lời giải theo giáo trình Cơ học kết cấu:
//   • Phương pháp lực: hệ cơ bản, biểu đồ đơn vị M̄ᵢ, biểu đồ tải M⁰_P, hệ phương trình chính tắc δ·X + Δ_P = 0.
//   • Chuyển vị theo công thức Maxwell–Mohr (nhân biểu đồ).
//   • Phương pháp chuyển vị dạng ma trận: K·Z = P.
// Áp dụng cho hệ thanh (dàn, dầm, khung) chịu tải trọng; không xét nền đàn hồi, gối đàn hồi, nhiệt độ, lún.
#pragma once
#include <string>
#include <vector>

#include "core/Analysis.h"

namespace fem {

struct ForceMethodResult {
    bool ok = false;
    std::string error;
    std::vector<std::string> warnings;
    int n = 0;                                  // bậc siêu tĩnh
    std::vector<std::string> redundants;        // mô tả các ẩn Xᵢ
    std::vector<std::vector<double>> delta;     // δᵢₖ
    std::vector<double> deltaP;                 // Δᵢₚ
    std::vector<double> X;                      // nghiệm
    std::vector<LoadSetResult> unit;            // trạng thái đơn vị (biểu đồ M̄ᵢ, N̄ᵢ, Q̄ᵢ) trên hệ cơ bản
    LoadSetResult load;                         // trạng thái tải trên hệ cơ bản (M⁰_P)
    LoadSetResult final;                        // biểu đồ cuối cùng
    double maxDiffM = 0;                        // sai lệch |M| lớn nhất so với lời giải PTHH
    // dữ liệu nội bộ cho Maxwell–Mohr
    std::vector<int> rowOf;                     // phương trình cân bằng của từng bậc tự do (−1 nếu không có)
    std::vector<int> pivots, redundantCols;
    std::vector<std::vector<double>> Bfull;     // ma trận cân bằng (phương trình × ẩn)
    std::vector<int> barOf;                     // thanh (chỉ số phần tử) cho mỗi cột nội lực, −1: phản lực
    std::vector<int> kindOf;                    // 0: N, 1: Mi, 2: Mj, 3: phản lực
    int loadSet = 0;
    std::vector<double> qFinal;                 // nội lực (ẩn) của trạng thái cuối cùng
};

ForceMethodResult forceMethod(const Model& m, int loadSet);

struct MohrResult {
    bool ok = false;
    std::string error;
    double value = 0;                    // chuyển vị (m hoặc rad)
    double fem = 0;                      // giá trị theo PTHH để so sánh
    std::vector<std::string> terms;      // đóng góp từng thanh
    LoadSetResult unitState;             // biểu đồ do lực đơn vị trên hệ cơ bản
};
// Chuyển vị nút theo phương dof do bộ tải của fm, tính bằng nhân biểu đồ trạng thái đơn vị với biểu đồ cuối cùng.
MohrResult mohrDisplacement(const Model& m, const ForceMethodResult& fm, int node, int dof);

struct MatrixMethodResult {
    bool ok = false;
    std::string error;
    std::vector<std::string> dofLabels;     // ẩn chuyển vị Z
    std::vector<std::vector<double>> K;     // ma trận độ cứng (chỉ bậc tự do tự do)
    std::vector<double> P, Z;               // véc tơ tải, nghiệm
    struct Elem {
        int e;
        double L, c, s;
        std::vector<std::vector<double>> kLocal, kGlobal; // 6×6
        std::vector<std::string> dofs;
    };
    std::vector<Elem> elems;
};
MatrixMethodResult matrixMethod(const Model& m, int loadSet, int maxEq = 40);

} // namespace fem
