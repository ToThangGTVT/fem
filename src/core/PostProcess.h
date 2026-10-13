// Hậu xử lý: biểu đồ bao nội lực, đường ảnh hưởng.
#pragma once
#include <string>
#include <vector>

#include "core/Analysis.h"

namespace fem {

enum class BarQty { N = 0, Q = 1, M = 2 };

const std::vector<double>& barValues(const BarResult& r, BarQty q);
// Giá trị tại x trên thanh (nội suy giữa các điểm mẫu). Tại bước nhảy: right = true lấy giá trị bên phải.
double barValueAt(const BarResult& r, BarQty q, double x, bool right = true);

// ---------------------------------------------------------------- biểu đồ bao
struct BarEnvelope {
    std::vector<double> x;
    std::vector<double> max[3], min[3]; // theo BarQty
};
// Bao nội lực của các thanh qua các bộ kết quả 'sets' (chỉ số trong AnalysisResult::sets).
std::vector<BarEnvelope> envelope(const Model& m, const AnalysisResult& r, const std::vector<int>& sets);

// ---------------------------------------------------------------- đường ảnh hưởng
enum class IlKind { Reaction = 0, Force = 1, Displacement = 2 };

struct IlTarget {
    IlKind kind = IlKind::Reaction;
    int node = 0;               // phản lực / chuyển vị
    int dof = UY;
    int elem = 0;               // nội lực
    double x = 0;               // vị trí tiết diện trên thanh, m (từ đầu i)
    BarQty qty = BarQty::M;
};

struct IlPoint {
    int elem;     // thanh trên đường chạy tải
    double a;     // vị trí lực P = 1 trên thanh, m
    double X, Y;  // tọa độ tổng thể của vị trí lực
    double value; // tung độ đường ảnh hưởng
};

struct InfluenceLine {
    bool ok = false;
    std::string error;
    std::vector<IlPoint> pts;      // theo thứ tự đường chạy, mỗi thanh một đoạn liên tiếp
    double areaPos = 0, areaNeg = 0; // diện tích phần dương/âm (theo hình chiếu ngang), dùng cho tải phân bố đều
    double vmax = 0, vmin = 0;
};

// Đường ảnh hưởng do lực P = 1 thẳng đứng hướng xuống di chuyển trên các thanh 'path'.
InfluenceLine influenceLine(const Model& m, const std::vector<int>& path, const IlTarget& t, int pointsPerBar = 20);

std::string ilTargetName(const Model& m, const IlTarget& t);

} // namespace fem
