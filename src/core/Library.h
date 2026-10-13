// Thư viện vật liệu và tiết diện theo tiêu chuẩn Việt Nam / Nga (Liên Xô).
#pragma once
#include <string>
#include <vector>

#include "core/Model.h"

namespace fem {

struct MaterialPreset {
    const char* name;
    const char* source;
    double E, nu, alpha, gamma; // kN/m², –, 1/°C, kN/m³
    int kind;                   // 1: thép, 2: bê tông
    double fd, ft;              // cường độ tính toán (MPa): thép f; bê tông Rb, Rbt
};

struct ProfilePreset {
    const char* name;
    const char* source;
    double h, b;    // mm
    double A;       // cm²
    double Ix;      // cm⁴
};

const std::vector<MaterialPreset>& materialPresets();
const std::vector<ProfilePreset>& profilePresets();

Material toMaterial(const MaterialPreset& p);
Section toSection(const ProfilePreset& p);
Section rectSection(double b, double h);           // m
Section circleSection(double d);                   // m
Section pipeSection(double D, double t);           // m
Section shapeSection(ShapeKind k, const std::array<double, 6>& p); // tiết diện T, I, hộp, [, L (m)

} // namespace fem
