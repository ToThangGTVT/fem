// Phân tích tĩnh tuyến tính theo phương pháp chuyển vị.
#pragma once
#include <array>
#include <string>
#include <vector>

#include "core/Elements.h"
#include "core/Model.h"

namespace fem {

enum DofState : char { DofFree = 0, DofSupport = 1, DofAuto = 2 };

struct LoadSetResult {
    std::string name;
    bool combo = false;
    std::vector<double> U; // chuyển vị nút, 3 bậc tự do/nút (m, rad)
    std::vector<double> R; // phản lực (kN, kNm) tại bậc tự do bị khóa hoặc gối đàn hồi
    std::vector<BarResult> bars;                                   // theo phần tử (rỗng với phần tử phẳng)
    // phần tử phẳng: (σx, σy, τxy, σ_ngoài) tại từng nút; σ_ngoài = σz (biến dạng phẳng) hoặc σθ (đối xứng trục)
    std::vector<std::array<std::array<double, 4>, 8>> elemStress;
    std::vector<std::array<double, 4>> nodeStress;                 // trung bình tại nút
    std::vector<int> nodeStressCount;
    double sumFx = 0, sumFy = 0; // tổng tải trọng tác dụng (kiểm tra cân bằng)
    double sumRx = 0, sumRy = 0;     // tổng phản lực (kể cả phản lực nền)
    double foundRx = 0, foundRy = 0; // hợp lực phản lực nền Winkler
};

struct AnalysisResult {
    bool ok = false;
    std::string error;
    std::vector<std::string> warnings;
    int numEquations = 0;
    long long profileSize = 0;
    double timeMs = 0;
    std::vector<char> dofState; // DofState cho từng bậc tự do
    std::vector<LoadSetResult> sets; // các trường hợp tải, tiếp theo là các tổ hợp
};

AnalysisResult analyze(const Model& m, int samplesPerBar = 40);

// Sắp xếp lại nút (Reverse Cuthill–McKee) để giảm profile ma trận độ cứng.
std::vector<int> rcmOrder(const Model& m);

} // namespace fem
