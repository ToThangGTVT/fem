// Phần tử phẳng tổng quát (đẳng tham số): T3, T6, Q4, Q8; ứng suất phẳng, biến dạng phẳng, đối xứng trục.
// Bài toán đối xứng trục: x = r (bán kính, ≥ 0), y = z (trục đối xứng); độ cứng, tải tính cho cả vành tròn (2π r).
#pragma once
#include <array>
#include <vector>

#include "core/Model.h"

namespace fem {

// Hàm dạng và đạo hàm theo tọa độ tự nhiên (tam giác: ξ = L2, η = L3)
void planeShape(ElemType t, double xi, double eta, double N[8], double dxi[8], double deta[8]);
void planeNodeNatural(ElemType t, int k, double& xi, double& eta);
bool isTriangle(ElemType t);

bool planeStiffness(const Element& el, const double* x, const double* y, const Material& mat, PlaneMode mode,
                    std::vector<double>& K);
// Ứng suất tại các nút: (σx, σy, τxy, σ_ngoài)
void planeStress(const Element& el, const double* x, const double* y, const Material& mat, PlaneMode mode, const double* u,
                 std::array<std::array<double, 4>, 8>& sig);
// Lực nút tương đương của lực thể tích (kN/m³)
void planeBodyForce(const Element& el, const double* x, const double* y, PlaneMode mode, double bx, double by,
                    std::vector<double>& f);
// Khối lượng tập trung tại nút (HRZ), ρ = t/m³
void planeLumpedMass(const Element& el, const double* x, const double* y, PlaneMode mode, double rho, std::vector<double>& m);

} // namespace fem
