// Lưới tam giác phi cấu trúc cho miền đa giác có lỗ: Delaunay ràng buộc (Bowyer–Watson, cạnh biên không bị cắt)
// và làm mịn Ruppert (chèn tâm đường tròn ngoại tiếp; tâm lấn vào cạnh biên thì chia đôi cạnh biên).
// Kích thước ô thay đổi theo khoảng cách d tới các lỗ (vật cản): h(d) = min(hFar, h + grading·d).
#pragma once
#include <array>
#include <string>
#include <vector>

#include "core/Model.h"

namespace fem {

using Polygon2 = std::vector<std::array<double, 2>>;

struct DomainSpec {
    Polygon2 outer;                   // biên ngoài
    std::vector<Polygon2> holes;      // các lỗ / vật cản, nằm trong biên ngoài, không chồng nhau
    // chỉ số nút sẵn có của mô hình ứng với các đỉnh (để lưới dùng chung nút với đường vẽ); rỗng hoặc −1: tạo nút mới
    std::vector<int> outerNodes;
    std::vector<std::vector<int>> holeNodes;
    double h = 0.01;        // kích thước ô sát vật cản (m)
    double hFar = 0.04;     // kích thước ô lớn nhất (m)
    double grading = 0.2;   // tốc độ tăng kích thước ô theo khoảng cách tới vật cản
    int maxPoints = 300000; // giới hạn số điểm (tránh chạy mãi khi góc đầu vào quá nhọn)
};

double polygonArea(const Polygon2& p); // có dấu: dương nếu ngược chiều kim đồng hồ
bool pointInPolygon(const Polygon2& p, double x, double y);

// Tạo lưới tam giác T3 (ngược chiều kim đồng hồ) vào mô hình m với vật liệu mat. Trả về số tam giác,
// hoặc 0 và ghi lý do vào err nếu dữ liệu không hợp lệ (đa giác tự cắt, lỗ nằm ngoài, lỗ chồng nhau...).
int meshDomain(Model& m, const DomainSpec& s, int mat, std::string* err = nullptr);

// Các vòng khép kín tạo bởi các thanh cho trước (mỗi nút thuộc đúng hai thanh của vòng). Trả về danh sách
// chỉ số nút theo thứ tự đi quanh vòng; ignored: số nhóm thanh liên thông không tạo thành một vòng đơn.
std::vector<std::vector<int>> closedLoops(const Model& m, const std::vector<int>& bars, int* ignored = nullptr);

} // namespace fem
