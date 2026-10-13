#include "core/Skyline.h"

#include <algorithm>
#include <cmath>

namespace fem {

void SkylineMatrix::init(const std::vector<int>& firstRow) {
    n_ = (int)firstRow.size();
    first_ = firstRow;
    ptr_.assign(n_ + 1, 0);
    for (int j = 0; j < n_; ++j) ptr_[j + 1] = ptr_[j] + (j - first_[j] + 1);
    a_.assign((size_t)ptr_[n_], 0.0);
}

void SkylineMatrix::add(int i, int j, double v) {
    if (i > j) std::swap(i, j);
    if (i < first_[j]) return;
    at(i, j) += v;
}

double SkylineMatrix::get(int i, int j) const {
    if (i > j) std::swap(i, j);
    if (i < first_[j]) return 0.0;
    return at(i, j);
}

int SkylineMatrix::factorize(double relTol) {
    for (int j = 0; j < n_; ++j) {
        const int mj = first_[j];
        const double orig = at(j, j);
        // g_ij = a_ij − Σ l_ki g_kj
        for (int i = mj; i < j; ++i) {
            const int k0 = std::max(first_[i], mj);
            double s = 0.0;
            for (int k = k0; k < i; ++k) s += at(k, i) * at(k, j);
            at(i, j) -= s;
        }
        // l_ij = g_ij / d_i ; d_j = a_jj − Σ g_ij l_ij
        double d = orig;
        for (int i = mj; i < j; ++i) {
            const double g = at(i, j);
            const double l = g / at(i, i);
            d -= g * l;
            at(i, j) = l;
        }
        if (!(d > relTol * std::abs(orig)) || !(d > 0.0)) return j;
        at(j, j) = d;
    }
    return -1;
}

void SkylineMatrix::solve(std::vector<double>& b) const {
    for (int j = 0; j < n_; ++j) {
        double s = 0.0;
        for (int k = first_[j]; k < j; ++k) s += at(k, j) * b[k];
        b[j] -= s;
    }
    for (int j = 0; j < n_; ++j) b[j] /= at(j, j);
    for (int j = n_ - 1; j >= 0; --j) {
        const double bj = b[j];
        for (int k = first_[j]; k < j; ++k) b[k] -= at(k, j) * bj;
    }
}

} // namespace fem
