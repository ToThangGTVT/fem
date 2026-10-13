// Đặc trưng hình học tiết diện bất kỳ (Sức bền vật liệu): diện tích, trọng tâm, mômen quán tính, trục chính,
// mômen kháng uốn, bán kính quán tính, lõi tiết diện; mômen tĩnh và bề rộng tại một cao độ (ứng suất tiếp).
// Hệ tọa độ tiết diện: x nằm ngang, y thẳng đứng (hướng lên, trùng trục y cục bộ của thanh).
#pragma once
#include <array>
#include <string>
#include <vector>

#include "core/Model.h"

namespace fem {

using Contour = std::vector<std::array<double, 2>>;

struct SectionProps {
    bool ok = false;
    double A = 0, xc = 0, yc = 0;           // diện tích, trọng tâm
    double Jx = 0, Jy = 0, Jxy = 0;         // mômen quán tính trung tâm (trục x ngang, y đứng)
    double J1 = 0, J2 = 0, alpha = 0;       // mômen quán tính chính, góc trục chính (độ, từ trục x)
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    double yTop = 0, yBot = 0;              // khoảng cách trọng tâm → thớ trên/dưới
    double Wtop = 0, Wbot = 0;              // mômen kháng uốn Jx/yTop, Jx/yBot
    double ix = 0, iy = 0;                  // bán kính quán tính
    Contour kern;                           // lõi tiết diện (đa giác lồi quanh trọng tâm)
};

std::vector<Contour> shapeContours(const SectionShape& s);
SectionProps sectionProps(const std::vector<Contour>& c);
// Mômen tĩnh của phần tiết diện nằm trên cao độ y đối với trục trung tâm nằm ngang (m³)
double staticMomentAbove(const std::vector<Contour>& c, double yc, double y);
// Bề rộng tiết diện tại cao độ y (m)
double widthAt(const std::vector<Contour>& c, double y);

// Tạo tiết diện từ hình dạng (tính A, I, h, yTop, yBot)
Section sectionFromShape(const std::string& name, const SectionShape& s);
const char* shapeName(ShapeKind k);
int shapeParamCount(ShapeKind k);
const char* shapeParamName(ShapeKind k, int i);

} // namespace fem
