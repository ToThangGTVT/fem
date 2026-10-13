// Công cụ tạo lưới đơn giản.
#pragma once
#include "core/Model.h"

namespace fem {

// Lưới tứ giác ánh xạ song tuyến từ 4 góc (ngược chiều kim đồng hồ). type = Quad4 hoặc Tri3.
// Trả về danh sách chỉ số nút trên lưới (ny+1)×(nx+1), theo hàng từ dưới lên.
std::vector<int> meshPatch(Model& m, const double cx[4], const double cy[4], int nx, int ny, ElemType type, int mat,
                           double t, bool incompatible = true);
std::vector<int> meshRectangle(Model& m, double x0, double y0, double w, double h, int nx, int ny, ElemType type,
                               int mat, double t, bool incompatible = true);

// Chia thanh e thành n đoạn bằng nhau; tải trên thanh được phân phối lại cho các đoạn.
void subdivideBar(Model& m, int e, int n);

// Tải phân bố trên cạnh biên của phần tử phẳng, quy về lực nút tương đương.
struct EdgeLoadSpec {
    bool pressure = false;  // false: lực qx, qy (kN/m chiều dài cạnh); true: áp lực pháp tuyến
    double qx = 0, qy = 0;  // kN/m
    double p0 = 0, ky = 0;  // áp lực p(y) = p0 + ky·y (kN/m²), dương: nén vào vật thể; nhân chiều dày phần tử
    bool clipNegative = true; // bỏ phần p < 0 (ví dụ phía trên mực nước)
};
// Áp dụng cho các cạnh biên có cả hai nút nằm trong 'nodes'. Trả về số cạnh đã đặt tải.
int applyEdgeLoad(Model& m, int caseIndex, const std::vector<int>& nodes, const EdgeLoadSpec& spec);

} // namespace fem
