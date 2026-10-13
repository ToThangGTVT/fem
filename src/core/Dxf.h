// Nhập / xuất bản vẽ DXF (ASCII, tương thích AutoCAD R12 trở lên).
#pragma once
#include <string>

#include "core/Analysis.h"
#include "core/PostProcess.h"

namespace fem {

struct DxfExportOptions {
    bool diagram = true;       // xuất biểu đồ nội lực
    BarQty qty = BarQty::M;
    double scale = 0;          // m trên một đơn vị nội lực (0: tự động)
    bool values = true;        // ghi giá trị
    bool deformed = false;     // xuất hình biến dạng
    double deformScale = 0;    // 0: tự động
    double textHeight = 0;     // 0: tự động (theo kích thước kết cấu)
    double unit = 1.0;         // hệ số đổi đơn vị khi xuất (1: m, 1000: mm)
};
// rs có thể null (chỉ xuất sơ đồ)
std::string exportDxf(const Model& m, const LoadSetResult* rs, const DxfExportOptions& o);

struct DxfImportOptions {
    double unit = 1.0;         // hệ số đổi sang m (0,001 nếu bản vẽ theo mm)
    ElemType type = ElemType::Frame;
    std::string layer;         // chỉ nhập các đối tượng trên lớp này (rỗng: mọi lớp)
    int mat = 0, sec = 0;
};
struct DxfImportResult {
    bool ok = false;
    std::string error;
    int lines = 0, nodes = 0, elements = 0, splits = 0;
    std::vector<std::string> layers; // các lớp có trong bản vẽ
};
// Nhập các đoạn thẳng (LINE, LWPOLYLINE, POLYLINE) thành thanh; nút trùng được gộp, thanh được cắt tại nút nằm trên nó.
DxfImportResult importDxf(const std::string& text, Model& m, const DxfImportOptions& o);
// Liệt kê các lớp có đoạn thẳng trong bản vẽ
std::vector<std::string> dxfLayers(const std::string& text);

// Cắt thanh tại các nút nằm giữa thanh (dùng sau khi nhập bản vẽ). Trả về số lần cắt.
int splitBarsAtNodes(Model& m, double tol);

} // namespace fem
