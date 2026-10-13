// Công cụ tạo lưới đơn giản.
#pragma once
#include "core/Model.h"

namespace fem {

// Lưới tứ giác ánh xạ song tuyến từ 4 góc (ngược chiều kim đồng hồ). type = Quad4, Tri3, Quad8 hoặc Tri6.
// Trả về danh sách chỉ số nút trên lưới (ny+1)×(nx+1), theo hàng từ dưới lên.
std::vector<int> meshPatch(Model& m, const double cx[4], const double cy[4], int nx, int ny, ElemType type, int mat,
                           double t, bool incompatible = true);
std::vector<int> meshRectangle(Model& m, double x0, double y0, double w, double h, int nx, int ny, ElemType type,
                               int mat, double t, bool incompatible = true);

// Chia thanh e thành n đoạn bằng nhau; tải trên thanh được phân phối lại cho các đoạn.
void subdivideBar(Model& m, int e, int n);

// Cắt thanh e tại vị trí t ∈ (0, 1) tính từ đầu i; tải trên thanh được chia lại. Trả về chỉ số nút mới.
int splitBarAt(Model& m, int e, double t);

// Thanh vút: chia thanh e thành n đoạn lăng trụ có tiết diện nội suy từ tiết diện hiện có (đầu i) tới secJ (đầu j).
// Tiết diện cùng dạng hình học: nội suy kích thước; ngược lại nội suy A, I, h tuyến tính. Trả về số tiết diện mới.
int taperBar(Model& m, int e, int secJ, int n);

// Vòm: nhịp L, mũi tên f, n đoạn; dạng 0: parabol, 1: cung tròn. Gối: 0 hai khớp, 1 ba khớp, 2 hai ngàm.
// Trả về chỉ số các nút từ trái sang phải; các thanh dùng vật liệu mat, tiết diện sec.
std::vector<int> makeArch(Model& m, double x0, double y0, double L, double f, int n, int shape, int supports, int mat, int sec);

// Gộp các nút trùng nhau (khoảng cách ≤ tol): nối lại phần tử, gối, tải; bỏ phần tử suy biến hoặc trùng lặp.
// Trả về số nút đã gộp.
int mergeNodes(Model& m, double tol);

// Tải phân bố trên cạnh biên của phần tử phẳng, quy về lực nút tương đương.
struct EdgeLoadSpec {
    bool pressure = false;  // false: lực qx, qy (kN/m chiều dài cạnh); true: áp lực pháp tuyến
    double qx = 0, qy = 0;  // kN/m (đối xứng trục: kN/m² trên mặt tròn xoay)
    double p0 = 0, ky = 0;  // áp lực p(y) = p0 + ky·y (kN/m²), dương: nén vào vật thể; nhân chiều dày phần tử
    bool clipNegative = true; // bỏ phần p < 0 (ví dụ phía trên mực nước)
};
// Áp dụng cho các cạnh biên có mọi nút nằm trong 'nodes'. Bài toán đối xứng trục: tải tính cho cả vành (2π r),
// không nhân chiều dày. Trả về số cạnh đã đặt tải.
int applyEdgeLoad(Model& m, int caseIndex, const std::vector<int>& nodes, const EdgeLoadSpec& spec);

} // namespace fem
