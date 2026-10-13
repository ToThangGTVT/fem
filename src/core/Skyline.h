// Ma trận đối xứng lưu dạng skyline (profile) và phân tích LDLᵀ (thuật toán COLSOL, Bathe).
#pragma once
#include <vector>

namespace fem {

class SkylineMatrix {
public:
    // firstRow[j] = chỉ số hàng nhỏ nhất có phần tử khác 0 trong cột j (≤ j).
    void init(const std::vector<int>& firstRow);
    int size() const { return n_; }
    long long storage() const { return (long long)a_.size(); }

    void add(int i, int j, double v); // cộng vào (i,j), tự đối xứng hóa; bỏ qua nếu nằm ngoài profile
    double get(int i, int j) const;

    // Phân tích LDLᵀ tại chỗ. Trả về -1 nếu thành công, ngược lại chỉ số phương trình có trụ ≤ 0
    // (hệ biến hình / ma trận suy biến).
    int factorize(double relTol = 1e-10);
    void solve(std::vector<double>& b) const; // giải sau khi phân tích

private:
    double& at(int i, int j) { return a_[ptr_[j] + (i - first_[j])]; }
    double at(int i, int j) const { return a_[ptr_[j] + (i - first_[j])]; }

    int n_ = 0;
    std::vector<int> first_;
    std::vector<long long> ptr_;
    std::vector<double> a_;
};

} // namespace fem
