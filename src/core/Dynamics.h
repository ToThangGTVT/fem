// Động lực học công trình:
//   • Phổ phản ứng theo TCVN 9386:2012 (phổ thiết kế, nền A–E, tổ hợp SRSS / CQC).
//   • Dao động cưỡng bức điều hòa (không cản): (K − θ²M)·u = F₀.
//   • Phân tích theo thời gian: Newmark gia tốc trung bình, cản Rayleigh.
// Mọi phân tích dùng mô hình đã chia nhỏ (DenseSystem); kết quả được quy về mô hình gốc để hiển thị.
#pragma once
#include <memory>
#include <string>
#include <vector>

#include "core/Analysis.h"
#include "core/Eigen.h"

namespace fem {

// ---------------------------------------------------------------- phổ phản ứng
struct SpectrumParams {
    double agR = 0.1;    // gia tốc nền tham chiếu, theo g (Phụ lục H TCVN 9386)
    double gammaI = 1.0; // hệ số tầm quan trọng
    double q = 3.9;      // hệ số ứng xử
    int soil = 2;        // 0..4: nền A, B, C, D, E
    double beta = 0.2;   // hệ số giới hạn dưới
    int direction = 0;   // 0: phương X, 1: phương Y
    bool cqc = false;    // tổ hợp CQC (mặc định SRSS)
    double zeta = 0.05;  // tỉ số cản (cho CQC)
};
void soilParams(int soil, double& S, double& TB, double& TC, double& TD);
double spectrumSd(const SpectrumParams& p, double T); // m/s²

struct ModalInfo {
    double T = 0, f = 0, omega = 0;
    double gamma = 0;     // hệ số tham gia theo phương đang xét
    double mEff = 0;      // khối lượng hữu hiệu (t)
    double ratio = 0;     // mEff / tổng khối lượng
    double Sd = 0;        // m/s²
    double baseShear = 0; // kN
};

struct SpectrumResult {
    bool ok = false;
    std::string error;
    std::vector<std::string> warnings;
    std::vector<ModalInfo> modes;
    double totalMass = 0, sumRatio = 0, baseShear = 0;
    LoadSetResult combined;            // tổ hợp (giá trị tuyệt đối) trên mô hình gốc
    std::vector<LoadSetResult> modal;  // từng dạng (có dấu)
};
SpectrumResult spectrumAnalysis(const Model& m, const EigenOptions& eo, const SpectrumParams& sp);

// ---------------------------------------------------------------- dao động cưỡng bức điều hòa
struct HarmonicResult {
    bool ok = false;
    std::string error;
    double theta = 0;            // tần số góc lực kích thích (rad/s)
    double dynFactor = 0;        // hệ số động = max|u_đ| / max|u_tĩnh|
    std::vector<double> omegas;  // vài tần số riêng thấp nhất (rad/s)
    LoadSetResult amp;           // biên độ nội lực, chuyển vị (mô hình gốc)
    LoadSetResult stat;          // lời giải tĩnh cùng tải
};
HarmonicResult harmonicAnalysis(const Model& m, const EigenOptions& eo, int loadSet, double theta);

// ---------------------------------------------------------------- theo thời gian
enum class TimeFn { Step = 0, Rect = 1, HalfSine = 2, Harmonic = 3, GroundSine = 4, GroundRecord = 5 };
struct TimeHistoryParams {
    int loadSet = 0;     // dạng phân bố tải (cho tải theo thời gian)
    TimeFn fn = TimeFn::Step;
    double amp = 1.0;    // hệ số nhân (tải) hoặc biên độ gia tốc nền (m/s²)
    double td = 0.2;     // thời gian tác dụng của xung (s)
    double freq = 2.0;   // tần số (Hz) cho tải điều hòa / gia tốc nền sin
    double dt = 0.005, tEnd = 4.0;
    double zeta = 0.05;  // tỉ số cản Rayleigh tại hai tần số riêng đầu
    int groundDir = 0;   // gia tốc nền theo X (0) hoặc Y (1)
    std::vector<std::pair<double, double>> record; // (t, a) m/s²
};
struct TimeHistoryResult {
    bool ok = false;
    std::string error;
    std::vector<std::string> warnings;
    std::vector<double> t;
    std::vector<std::vector<float>> u;   // chuyển vị các bậc tự do tự do theo từng bước
    std::vector<double> maxAbs;          // |u|max của từng bậc tự do của mô hình đã chia nhỏ (3/nút)
    std::shared_ptr<DenseSystem> sys;
    double alpha = 0, betaR = 0;         // hệ số cản Rayleigh C = αM + βK
    double omega1 = 0;
};
TimeHistoryResult timeHistory(const Model& m, const EigenOptions& eo, const TimeHistoryParams& p);
double timeFunction(const TimeHistoryParams& p, double t);
// Chuyển vị theo thời gian của một bậc tự do mô hình gốc
std::vector<double> historyOf(const TimeHistoryResult& r, int node, int dof);
// Nội lực tại bước k (lực đàn hồi K·u đặt lên mô hình đã chia nhỏ), quy về mô hình gốc
LoadSetResult timeSnapshot(const Model& m, const TimeHistoryResult& r, int k);

// Ghép kết quả trên mô hình đã chia nhỏ về mô hình gốc
LoadSetResult mergeToOriginal(const Model& orig, const DenseSystem& sys, const LoadSetResult& sub);

} // namespace fem
