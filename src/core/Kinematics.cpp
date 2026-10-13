#include "core/Kinematics.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace fem {

namespace {

// Đưa ma trận A (rows × cols, theo hàng) về dạng bậc thang rút gọn; trả về hạng và các cột trụ.
int rref(std::vector<double>& A, int rows, int cols, std::vector<int>& pivots, double tol) {
    pivots.clear();
    int r = 0;
    for (int c = 0; c < cols && r < rows; ++c) {
        int best = r;
        for (int i = r + 1; i < rows; ++i)
            if (std::abs(A[i * cols + c]) > std::abs(A[best * cols + c])) best = i;
        if (std::abs(A[best * cols + c]) <= tol) continue;
        if (best != r)
            for (int k = 0; k < cols; ++k) std::swap(A[best * cols + k], A[r * cols + k]);
        const double pv = A[r * cols + c];
        for (int k = 0; k < cols; ++k) A[r * cols + k] /= pv;
        for (int i = 0; i < rows; ++i) {
            if (i == r) continue;
            const double f = A[i * cols + c];
            if (f == 0) continue;
            for (int k = 0; k < cols; ++k) A[i * cols + k] -= f * A[r * cols + k];
        }
        pivots.push_back(c);
        ++r;
    }
    return r;
}

} // namespace

KinematicResult kinematicAnalysis(const Model& m) {
    KinematicResult r;
    const int nn = (int)m.nodes.size();
    std::vector<int> bars;
    for (size_t e = 0; e < m.elements.size(); ++e) {
        const Element& el = m.elements[e];
        if (el.isBar()) {
            if (el.n[0] < 0 || el.n[1] < 0 || el.n[0] >= nn || el.n[1] >= nn) {
                r.note = "Thanh có nút không hợp lệ.";
                return r;
            }
            bars.push_back((int)e);
        }
    }
    if (bars.empty()) {
        r.note = "Chỉ áp dụng cho hệ thanh (dàn, khung, dầm).";
        return r;
    }
    if (bars.size() != m.elements.size()) r.note = "Phần tử phẳng được bỏ qua (chỉ xét các thanh).";
    r.applicable = true;

    // ---- đếm theo công thức Чебышёв
    std::vector<int> rigidEnds(nn, 0), hingedEnds(nn, 0);
    auto endRigid = [&](const Element& el, int k) { return el.type == ElemType::Frame && !(k == 0 ? el.hingeI : el.hingeJ); };
    for (int e : bars) {
        const Element& el = m.elements[e];
        for (int k = 0; k < 2; ++k) (endRigid(el, k) ? rigidEnds : hingedEnds)[el.n[k]]++;
    }
    r.disks = (int)bars.size();
    for (int i = 0; i < nn; ++i) {
        if (rigidEnds[i] >= 1) {
            r.rigidJoints += rigidEnds[i] - 1;
            r.hinges += hingedEnds[i];
        } else if (hingedEnds[i] >= 1) {
            r.hinges += hingedEnds[i] - 1;
        }
    }
    std::vector<int> nodeOfSupport(nn, -1);
    for (size_t i = 0; i < m.supports.size(); ++i) {
        const Support& s = m.supports[i];
        if (s.node < 0 || s.node >= nn) continue;
        nodeOfSupport[s.node] = (int)i;
        for (int c = 0; c < 3; ++c)
            if (s.fix[c] || s.k[c] > 0) ++r.supports;
    }
    r.W = 3 * r.disks - 3 * r.rigidJoints - 2 * r.hinges - r.supports;

    // ---- phương trình cân bằng nút: Fx, Fy luôn có; mômen khi nút có đầu thanh nối cứng hoặc liên kết xoay
    std::vector<int> row(3 * nn, -1);
    std::vector<char> usedNode(nn, 0);
    for (int e : bars)
        for (int k = 0; k < 2; ++k) usedNode[m.elements[e].n[k]] = 1;
    int neq = 0;
    for (int i = 0; i < nn; ++i) {
        if (!usedNode[i]) continue;
        row[3 * i] = neq++;
        row[3 * i + 1] = neq++;
        const Support* s = nodeOfSupport[i] >= 0 ? &m.supports[nodeOfSupport[i]] : nullptr;
        if (rigidEnds[i] > 0 || (s && (s->fix[2] || s->k[2] > 0))) row[3 * i + 2] = neq++;
    }
    // ---- ẩn: nội lực đầu thanh (N, Mi, Mj trừ đầu khớp) và phản lực
    std::vector<std::vector<std::pair<int, double>>> cols;
    for (int e : bars) {
        const Element& el = m.elements[e];
        const Node& a = m.nodes[el.n[0]];
        const Node& b = m.nodes[el.n[1]];
        const double L = std::hypot(b.x - a.x, b.y - a.y);
        if (L <= 0) {
            r.applicable = false;
            r.note = "Thanh có chiều dài bằng 0.";
            return r;
        }
        const double c = (b.x - a.x) / L, s = (b.y - a.y) / L;
        auto add = [&](std::vector<std::pair<int, double>>& col, int node, int dof, double v) {
            const int rr = row[3 * node + dof];
            if (rr >= 0 && v != 0) col.push_back({rr, v});
        };
        // N: thanh kéo các nút về phía nhau
        {
            std::vector<std::pair<int, double>> col;
            add(col, el.n[0], 0, c);
            add(col, el.n[0], 1, s);
            add(col, el.n[1], 0, -c);
            add(col, el.n[1], 1, -s);
            cols.push_back(col);
        }
        // mômen đầu thanh Mi (Mj) kèm lực cắt cân bằng V = M/L
        for (int k = 0; k < 2; ++k) {
            if (!endRigid(el, k)) continue;
            std::vector<std::pair<int, double>> col;
            // mômen đầu thanh M = 1 tại đầu k cần lực cắt V = 1/L ở hai đầu; lực thanh tác dụng lên nút
            // (nhân với −L để chuẩn hóa): nút i (−s, c), nút j (s, −c), mômen L tại đầu k
            add(col, el.n[k], 2, L);
            add(col, el.n[0], 0, -s);
            add(col, el.n[0], 1, c);
            add(col, el.n[1], 0, s);
            add(col, el.n[1], 1, -c);
            cols.push_back(col);
        }
    }
    for (const auto& sp : m.supports) {
        if (sp.node < 0 || sp.node >= nn || !usedNode[sp.node]) continue;
        for (int c = 0; c < 3; ++c) {
            if (!(sp.fix[c] || sp.k[c] > 0)) continue;
            const int rr = row[3 * sp.node + c];
            if (rr < 0) continue;
            cols.push_back({{rr, 1.0}});
        }
    }
    r.equations = neq;
    r.unknowns = (int)cols.size();
    if (neq > 4000 || r.unknowns > 8000) {
        r.note = "Hệ quá lớn cho phân tích đại số — chỉ dùng công thức W.";
        r.rank = -1;
        r.verdict = r.W > 0 ? "Hệ biến hình (W > 0)." : "Không kết luận được (W ≤ 0 là điều kiện cần).";
        return r;
    }
    // ma trận cân bằng B (phương trình × ẩn)
    const int R = neq, C = r.unknowns;
    std::vector<double> B((size_t)R * C, 0.0);
    for (int j = 0; j < C; ++j)
        for (const auto& [i, v] : cols[j]) B[(size_t)i * C + j] += v;
    // chuẩn hóa mỗi cột
    for (int j = 0; j < C; ++j) {
        double mx = 0;
        for (int i = 0; i < R; ++i) mx = std::max(mx, std::abs(B[(size_t)i * C + j]));
        if (mx > 0)
            for (int i = 0; i < R; ++i) B[(size_t)i * C + j] /= mx;
    }
    std::vector<double> A = B;
    std::vector<int> piv;
    r.rank = rref(A, R, C, piv, 1e-9);
    r.indeterminacy = C - r.rank;
    r.mechanisms = R - r.rank;

    if (r.mechanisms > 0) {
        r.verdict = r.W > 0 ? "Hệ biến hình: thiếu liên kết (W > 0)."
                            : "Hệ biến hình tức thời: đủ số liên kết (W ≤ 0) nhưng sắp xếp không hợp lý.";
        // dạng chuyển động: không gian rỗng của Bᵀ (chuyển vị nút không gây biến dạng)
        std::vector<double> Bt((size_t)C * R);
        for (int i = 0; i < R; ++i)
            for (int j = 0; j < C; ++j) Bt[(size_t)j * R + i] = B[(size_t)i * C + j];
        std::vector<int> pv;
        rref(Bt, C, R, pv, 1e-9);
        std::vector<char> isPivot(R, 0);
        for (int p : pv) isPivot[p] = 1;
        int freeCol = -1;
        for (int j = 0; j < R; ++j)
            if (!isPivot[j]) {
                freeCol = j;
                break;
            }
        if (freeCol >= 0) {
            std::vector<double> d(R, 0.0);
            d[freeCol] = 1;
            for (size_t k = 0; k < pv.size(); ++k) d[pv[k]] = -Bt[k * R + freeCol];
            r.mechanismU.assign(3 * nn, 0.0);
            for (int i = 0; i < 3 * nn; ++i)
                if (row[i] >= 0) r.mechanismU[i] = d[row[i]];
            // góc xoay của dạng cơ cấu trên các nút không có phương trình mômen: không hiển thị
        }
    } else if (r.indeterminacy == 0) {
        r.verdict = "Hệ bất biến hình, tĩnh định.";
    } else {
        char b[96];
        std::snprintf(b, sizeof b, "Hệ bất biến hình, siêu tĩnh bậc %d.", r.indeterminacy);
        r.verdict = b;
    }
    return r;
}

} // namespace fem
