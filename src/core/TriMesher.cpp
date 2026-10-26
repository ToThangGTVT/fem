#include "core/TriMesher.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <unordered_map>
#include <unordered_set>

namespace fem {

namespace {

struct P2 {
    double x, y;
};

double orient(const P2& a, const P2& b, const P2& c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); }

// > 0 nếu d nằm trong đường tròn ngoại tiếp tam giác abc (abc ngược chiều kim đồng hồ)
double incircle(const P2& a, const P2& b, const P2& c, const P2& d) {
    const double adx = a.x - d.x, ady = a.y - d.y, bdx = b.x - d.x, bdy = b.y - d.y, cdx = c.x - d.x, cdy = c.y - d.y;
    return (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) + (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy) +
           (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
}

void circumcenter(const P2& a, const P2& b, const P2& c, double& x, double& y) {
    const double bx = b.x - a.x, by = b.y - a.y, cx = c.x - a.x, cy = c.y - a.y;
    const double d = 2 * (bx * cy - by * cx);
    const double b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    x = a.x + (cy * b2 - by * c2) / d;
    y = a.y + (bx * c2 - cx * b2) / d;
}

bool segmentsCross(const P2& a, const P2& b, const P2& c, const P2& d) {
    const double d1 = orient(c, d, a), d2 = orient(c, d, b), d3 = orient(a, b, c), d4 = orient(a, b, d);
    return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0)) && d1 != 0 && d2 != 0 && d3 != 0 && d4 != 0;
}

double segDist(double px, double py, const P2& a, const P2& b) {
    const double dx = b.x - a.x, dy = b.y - a.y, L2 = dx * dx + dy * dy;
    double t = L2 > 0 ? ((px - a.x) * dx + (py - a.y) * dy) / L2 : 0;
    t = std::clamp(t, 0.0, 1.0);
    return std::hypot(a.x + t * dx - px, a.y + t * dy - py);
}

long long ekey(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((long long)a << 32) | (unsigned)b;
}

class Mesher {
public:
    struct Tri {
        int v[3];
        int nb[3]; // nb[i]: tam giác bên kia cạnh đối diện đỉnh v[i], tức cạnh (v[i+1], v[i+2])
        bool alive = true;
        bool domain = false;
    };
    std::vector<P2> pts;
    std::vector<Tri> tris;
    std::vector<int> freeList, vt, mark;
    int stamp = 0, last = 0;
    std::unordered_set<long long> cons; // cạnh biên (ràng buộc)
    double eps = 0;                     // ngưỡng hướng (theo bình phương kích thước miền)

    int addPoint(double x, double y) {
        pts.push_back({x, y});
        vt.push_back(-1);
        return (int)pts.size() - 1;
    }
    int newTri(int a, int b, int c) {
        int t;
        if (!freeList.empty()) {
            t = freeList.back();
            freeList.pop_back();
        } else {
            t = (int)tris.size();
            tris.emplace_back();
            mark.push_back(0);
        }
        Tri& T = tris[t];
        T.v[0] = a, T.v[1] = b, T.v[2] = c;
        T.nb[0] = T.nb[1] = T.nb[2] = -1;
        T.alive = true;
        T.domain = false;
        vt[a] = vt[b] = vt[c] = t;
        return t;
    }

    // Đi bộ tới tam giác chứa p. stopAtCons: không đi qua cạnh ràng buộc — khi bị chặn, trả −1 và cạnh chặn (ba, bb).
    // Đi bộ ngẫu nhiên có nhớ (không quay lại tam giác vừa đi qua); quá số bước thì tìm vét cạn (điểm gần như
    // thẳng hàng làm phép thử hướng không nhất quán do sai số làm tròn).
    int locate(const P2& p, int start, bool stopAtCons, int* ba = nullptr, int* bb = nullptr) {
        int t = (start >= 0 && start < (int)tris.size() && tris[start].alive) ? start : last;
        if (t < 0 || t >= (int)tris.size() || !tris[t].alive) {
            t = 0;
            while (t < (int)tris.size() && !tris[t].alive) ++t;
            if (t == (int)tris.size()) return -1;
        }
        int prev = -1;
        unsigned seed = 2463534242u;
        for (int step = 0; step < 20000; ++step) {
            const Tri& T = tris[t];
            seed ^= seed << 13, seed ^= seed >> 17, seed ^= seed << 5;
            const int r = (int)(seed % 3);
            int cross = -1;
            for (int k = 0; k < 3 && cross < 0; ++k) {
                const int i = (k + r) % 3;
                if (T.nb[i] == prev && prev >= 0) continue;
                if (orient(pts[T.v[(i + 1) % 3]], pts[T.v[(i + 2) % 3]], p) < 0) cross = i;
            }
            if (cross < 0) return t; // p nằm trong T (hoặc trên cạnh)
            const int a = T.v[(cross + 1) % 3], b = T.v[(cross + 2) % 3];
            if (stopAtCons && cons.count(ekey(a, b))) {
                if (ba) *ba = a;
                if (bb) *bb = b;
                return -1;
            }
            if (T.nb[cross] < 0) return -1;
            prev = t;
            t = T.nb[cross];
        }
        // vét cạn: tam giác có hướng nhỏ nhất (tương đối) lớn nhất
        int best = -1;
        double bestv = -1e300;
        for (int q = 0; q < (int)tris.size(); ++q) {
            if (!tris[q].alive) continue;
            const Tri& T = tris[q];
            double mn = 1e300;
            for (int i = 0; i < 3; ++i) {
                const P2 &A = pts[T.v[(i + 1) % 3]], &B = pts[T.v[(i + 2) % 3]];
                mn = std::min(mn, orient(A, B, p) / std::max(std::hypot(B.x - A.x, B.y - A.y), 1e-300));
            }
            if (mn > bestv) bestv = mn, best = q;
        }
        return best;
    }

    // Bowyer–Watson: chèn điểm pi, hốc không vượt cạnh ràng buộc (trừ allow). checkEncroach: nếu điểm lấn vào
    // đường tròn đường kính của một cạnh ràng buộc trên biên hốc thì không chèn, trả 2 và cạnh đó.
    // Trả 1 nếu chèn được, 0 nếu thất bại.
    int insert(int pi, int t0, long long allow, bool checkEncroach, std::vector<int>* created, int* ea = nullptr,
               int* eb = nullptr) {
        const P2& p = pts[pi];
        ++stamp;
        std::vector<int> cav{t0}, stack{t0};
        mark[t0] = stamp;
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            for (int i = 0; i < 3; ++i) {
                const int n = tris[t].nb[i];
                if (n < 0 || mark[n] == stamp) continue;
                const long long k = ekey(tris[t].v[(i + 1) % 3], tris[t].v[(i + 2) % 3]);
                if (k != allow && cons.count(k)) continue;
                const Tri& N = tris[n];
                if (incircle(pts[N.v[0]], pts[N.v[1]], pts[N.v[2]], p) > 0) {
                    mark[n] = stamp;
                    cav.push_back(n);
                    stack.push_back(n);
                }
            }
        }
        struct BE {
            int a, b, outer, owner;
        };
        std::vector<BE> be;
        for (int guard = 0; guard < 100; ++guard) {
            be.clear();
            bool bad = false;
            for (int t : cav) {
                if (mark[t] != stamp) continue;
                for (int i = 0; i < 3; ++i) {
                    const int n = tris[t].nb[i];
                    if (n >= 0 && mark[n] == stamp) continue;
                    const int a = tris[t].v[(i + 1) % 3], b = tris[t].v[(i + 2) % 3];
                    if (orient(pts[a], pts[b], p) <= eps) {
                        // hốc không hình sao đối với p (sai số làm tròn): bớt tam giác, hoặc thêm tam giác bên kia
                        bad = true;
                        if (t != t0) mark[t] = 0;
                        else if (n >= 0 && (ekey(a, b) == allow || !cons.count(ekey(a, b)))) {
                            mark[n] = stamp;
                            cav.push_back(n);
                        } else return 0;
                        break;
                    }
                    be.push_back({a, b, n, t});
                }
                if (bad) break;
            }
            if (!bad) break;
            if (guard == 99) return 0;
        }
        if (checkEncroach)
            for (const BE& e : be) {
                if (!cons.count(ekey(e.a, e.b))) continue;
                const P2 &A = pts[e.a], &B = pts[e.b];
                if ((A.x - p.x) * (B.x - p.x) + (A.y - p.y) * (B.y - p.y) < 0) {
                    if (ea) *ea = e.a;
                    if (eb) *eb = e.b;
                    return 2;
                }
            }
        for (int t : cav)
            if (mark[t] == stamp) {
                tris[t].alive = false;
                freeList.push_back(t);
            }
        std::unordered_map<int, int> startAt, endAt;
        std::vector<int> made;
        for (const BE& e : be) {
            const int t = newTri(pi, e.a, e.b);
            tris[t].nb[0] = e.outer;
            if (e.outer >= 0) {
                Tri& O = tris[e.outer];
                for (int j = 0; j < 3; ++j)
                    if (O.v[(j + 1) % 3] == e.b && O.v[(j + 2) % 3] == e.a) O.nb[j] = t;
            }
            startAt[e.a] = t;
            endAt[e.b] = t;
            made.push_back(t);
        }
        for (int t : made) {
            Tri& T = tris[t];
            auto f = endAt.find(T.v[1]);
            if (f != endAt.end()) T.nb[2] = f->second;
            auto g = startAt.find(T.v[2]);
            if (g != startAt.end()) T.nb[1] = g->second;
        }
        last = made.empty() ? last : made[0];
        if (created) *created = made;
        return 1;
    }

    bool edgeExists(int a, int b) {
        const int t0 = vt[a];
        if (t0 < 0) return false;
        for (int dir = 0; dir < 2; ++dir) {
            int t = t0;
            for (int guard = 0; guard < 1000 && t >= 0; ++guard) {
                const Tri& T = tris[t];
                int i = 0;
                while (i < 3 && T.v[i] != a) ++i;
                if (i == 3) return false;
                if (T.v[(i + 1) % 3] == b || T.v[(i + 2) % 3] == b) return true;
                t = dir == 0 ? T.nb[(i + 2) % 3] : T.nb[(i + 1) % 3];
                if (t == t0) return false;
            }
        }
        return false;
    }
};

} // namespace

double polygonArea(const Polygon2& p) {
    double A = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        const auto& a = p[i];
        const auto& b = p[(i + 1) % p.size()];
        A += a[0] * b[1] - b[0] * a[1];
    }
    return 0.5 * A;
}

bool pointInPolygon(const Polygon2& p, double x, double y) {
    bool in = false;
    for (size_t i = 0, j = p.size() - 1; i < p.size(); j = i++) {
        const double xi = p[i][0], yi = p[i][1], xj = p[j][0], yj = p[j][1];
        if ((yi > y) != (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) in = !in;
    }
    return in;
}

int meshDomain(Model& m, const DomainSpec& s, int mat, std::string* err) {
    auto fail = [&](const std::string& e) {
        if (err) *err = e;
        return 0;
    };
    if (s.h <= 0) return fail("Kích thước ô phải > 0.");
    // ---- kiểm tra đa giác
    std::vector<const Polygon2*> polys{&s.outer};
    for (const auto& h : s.holes) polys.push_back(&h);
    for (size_t k = 0; k < polys.size(); ++k) {
        const Polygon2& P = *polys[k];
        if (P.size() < 3 || std::abs(polygonArea(P)) < 1e-14)
            return fail(k == 0 ? "Biên ngoài không hợp lệ." : "Vật cản " + std::to_string(k) + " có diện tích bằng 0.");
    }
    auto edgesCross = [](const Polygon2& A, const Polygon2& B, bool same) {
        for (size_t i = 0; i < A.size(); ++i)
            for (size_t j = 0; j < B.size(); ++j) {
                if (same && (i == j || (i + 1) % A.size() == j || (j + 1) % A.size() == i)) continue;
                const P2 a{A[i][0], A[i][1]}, b{A[(i + 1) % A.size()][0], A[(i + 1) % A.size()][1]};
                const P2 c{B[j][0], B[j][1]}, d{B[(j + 1) % B.size()][0], B[(j + 1) % B.size()][1]};
                if (segmentsCross(a, b, c, d)) return true;
            }
        return false;
    };
    for (size_t k = 0; k < polys.size(); ++k)
        if (edgesCross(*polys[k], *polys[k], true))
            return fail(k == 0 ? "Biên ngoài tự cắt." : "Vật cản " + std::to_string(k) + " tự cắt.");
    // vật cản được phép chạm hoặc cắt biên miền (công trình đặt trên mặt đất, trụ trên đáy sông): phần nằm ngoài miền
    // bị bỏ qua; nhưng các vật cản không được chồng lên nhau
    for (size_t k = 1; k < polys.size(); ++k)
        for (size_t j = 1; j < k; ++j)
            if (edgesCross(*polys[k], *polys[j], false) || pointInPolygon(*polys[j], (*polys[k])[0][0], (*polys[k])[0][1]) ||
                pointInPolygon(*polys[k], (*polys[j])[0][0], (*polys[j])[0][1]))
                return fail("Vật cản " + std::to_string(j) + " và " + std::to_string(k) + " chồng lên nhau.");

    // ---- hàm kích thước ô
    std::vector<std::pair<P2, P2>> holeEdges;
    for (const auto& H : s.holes)
        for (size_t i = 0; i < H.size(); ++i)
            holeEdges.push_back({{H[i][0], H[i][1]}, {H[(i + 1) % H.size()][0], H[(i + 1) % H.size()][1]}});
    const double hFar = std::max(s.hFar, s.h);
    auto hAt = [&](double x, double y) {
        if (holeEdges.empty()) return s.h;
        double d = 1e300;
        for (const auto& e : holeEdges) d = std::min(d, segDist(x, y, e.first, e.second));
        return std::min(hFar, s.h + s.grading * d);
    };

    // ---- rời rạc hóa cạnh: số đoạn theo ∫ds/h
    double x0 = 1e300, y0 = 1e300, x1 = -1e300, y1 = -1e300;
    for (const auto& v : s.outer) {
        x0 = std::min(x0, v[0]), x1 = std::max(x1, v[0]);
        y0 = std::min(y0, v[1]), y1 = std::max(y1, v[1]);
    }
    const double size = std::max(x1 - x0, y1 - y0);
    Mesher M;
    M.eps = 1e-14 * size * size;
    // tam giác bao rất lớn
    const double cx = 0.5 * (x0 + x1), cy = 0.5 * (y0 + y1), R = 20 * size;
    const int s0 = M.addPoint(cx - R, cy - R), s1 = M.addPoint(cx + R, cy - R), s2 = M.addPoint(cx, cy + R);
    M.newTri(s0, s1, s2);
    std::vector<int> nodeOf; // điểm → nút sẵn có của mô hình (−1: tạo mới)
    nodeOf.assign(3, -1);
    // ---- đồ thị phẳng của các cạnh: đỉnh dùng chung theo tọa độ, cạnh biên và cạnh vật cản được chia tại giao điểm
    // và tại đỉnh nằm trên cạnh kia, các đoạn trùng nhau (vật cản tựa lên biên) chỉ giữ một
    const double tol = 1e-9 * size;
    std::vector<int> bpts;
    auto vertex = [&](double x, double y, int node) {
        for (int pi : bpts)
            if (std::abs(M.pts[pi].x - x) <= tol && std::abs(M.pts[pi].y - y) <= tol) {
                if (node >= 0 && nodeOf[pi] < 0) nodeOf[pi] = node;
                return pi;
            }
        const int pi = M.addPoint(x, y);
        nodeOf.push_back(node);
        bpts.push_back(pi);
        return pi;
    };
    struct Piece {
        int a, b;
    };
    std::vector<Piece> pieces;
    for (size_t k = 0; k < polys.size(); ++k) {
        const Polygon2& P = *polys[k];
        const std::vector<int>* ids = k == 0 ? &s.outerNodes : (k - 1 < s.holeNodes.size() ? &s.holeNodes[k - 1] : nullptr);
        for (size_t i = 0; i < P.size(); ++i) {
            const P2 a{P[i][0], P[i][1]}, b{P[(i + 1) % P.size()][0], P[(i + 1) % P.size()][1]};
            const size_t i1 = (i + 1) % P.size();
            const double dx = b.x - a.x, dy = b.y - a.y, L2 = dx * dx + dy * dy;
            std::vector<double> ts;
            // biên ngoài cắt với vật cản (vật cản không cắt nhau nên chỉ cần xét cặp biên ngoài – vật cản)
            for (size_t q = 0; q < polys.size(); ++q) {
                if (q == k || (k != 0 && q != 0)) continue;
                const Polygon2& Q = *polys[q];
                for (size_t j = 0; j < Q.size(); ++j) {
                    const P2 c{Q[j][0], Q[j][1]}, d{Q[(j + 1) % Q.size()][0], Q[(j + 1) % Q.size()][1]};
                    // đỉnh của đa giác kia nằm trên cạnh này
                    const double tc = ((c.x - a.x) * dx + (c.y - a.y) * dy) / L2;
                    if (tc > 1e-9 && tc < 1 - 1e-9 && std::abs(orient(a, b, c)) <= tol * std::sqrt(L2)) ts.push_back(tc);
                    // giao điểm thực sự
                    const double d1 = orient(c, d, a), d2 = orient(c, d, b), d3 = orient(a, b, c), d4 = orient(a, b, d);
                    const double e = tol * std::sqrt(L2);
                    if (((d1 > e && d2 < -e) || (d1 < -e && d2 > e)) && ((d3 > e && d4 < -e) || (d3 < -e && d4 > e)))
                        ts.push_back(d1 / (d1 - d2));
                }
            }
            std::sort(ts.begin(), ts.end());
            int prev = vertex(a.x, a.y, ids && i < ids->size() ? (*ids)[i] : -1);
            const int last = vertex(b.x, b.y, ids && i1 < ids->size() ? (*ids)[i1] : -1);
            for (double t : ts) {
                const int v = vertex(a.x + t * dx, a.y + t * dy, -1);
                if (v != prev) pieces.push_back({prev, v});
                prev = v;
            }
            if (last != prev) pieces.push_back({prev, last});
        }
    }
    {
        std::unordered_set<long long> seen;
        std::vector<Piece> uniq;
        for (const Piece& pc : pieces)
            if (seen.insert(ekey(pc.a, pc.b)).second) uniq.push_back(pc);
        pieces.swap(uniq);
    }
    // ---- rời rạc hóa từng đoạn: số đoạn con theo ∫ds/h
    std::vector<std::pair<int, int>> segs;
    for (const Piece& pc : pieces) {
        const P2 a = M.pts[pc.a], b = M.pts[pc.b];
        const int ns = 64;
        std::vector<double> I(ns + 1, 0);
        const double L = std::hypot(b.x - a.x, b.y - a.y);
        for (int j = 1; j <= ns; ++j) {
            const double t = (j - 0.5) / ns;
            I[j] = I[j - 1] + L / ns / hAt(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
        }
        const int n = std::max(1, (int)std::lround(I[ns]));
        int prev = pc.a;
        for (int q = 1; q < n; ++q) {
            const double target = I[ns] * q / n;
            int j = 1;
            while (j < ns && I[j] < target) ++j;
            const double t = (j - 1 + (target - I[j - 1]) / std::max(I[j] - I[j - 1], 1e-300)) / ns;
            const int v = M.addPoint(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
            nodeOf.push_back(-1);
            bpts.push_back(v);
            segs.push_back({prev, v});
            prev = v;
        }
        segs.push_back({prev, pc.b});
    }
    // ---- Delaunay của các điểm biên (thứ tự ngẫu nhiên cố định)
    std::mt19937 rng(12345);
    std::shuffle(bpts.begin(), bpts.end(), rng);
    for (int pi : bpts) {
        const int t = M.locate(M.pts[pi], M.last, false);
        if (t < 0 || M.insert(pi, t, -1, false, nullptr) != 1) return fail("Lỗi dựng lưới Delaunay (điểm biên trùng nhau?).");
    }
    // ---- khôi phục cạnh biên: cạnh chưa có thì chia đôi
    std::vector<std::pair<int, int>> queue(segs.rbegin(), segs.rend());
    int splits = 0;
    while (!queue.empty()) {
        const auto [a, b] = queue.back();
        queue.pop_back();
        if (M.edgeExists(a, b)) {
            M.cons.insert(ekey(a, b));
            continue;
        }
        if (++splits > 200000) return fail("Không khôi phục được cạnh biên.");
        const int mi = M.addPoint(0.5 * (M.pts[a].x + M.pts[b].x), 0.5 * (M.pts[a].y + M.pts[b].y));
        nodeOf.push_back(-1);
        const int t = M.locate(M.pts[mi], M.vt[a], false);
        if (t < 0 || M.insert(mi, t, -1, false, nullptr) != 1) return fail("Lỗi khôi phục cạnh biên.");
        queue.push_back({mi, b});
        queue.push_back({a, mi});
    }
    auto classify = [&](int t) {
        Mesher::Tri& T = M.tris[t];
        const double gx = (M.pts[T.v[0]].x + M.pts[T.v[1]].x + M.pts[T.v[2]].x) / 3;
        const double gy = (M.pts[T.v[0]].y + M.pts[T.v[1]].y + M.pts[T.v[2]].y) / 3;
        bool in = T.v[0] > 2 && T.v[1] > 2 && T.v[2] > 2 && pointInPolygon(s.outer, gx, gy);
        for (size_t k = 0; in && k < s.holes.size(); ++k)
            if (pointInPolygon(s.holes[k], gx, gy)) in = false;
        T.domain = in;
    };
    for (int t = 0; t < (int)M.tris.size(); ++t)
        if (M.tris[t].alive) classify(t);

    // ---- làm mịn Ruppert
    auto splitSegment = [&](int a, int b) {
        const P2 &A = M.pts[a], &B = M.pts[b];
        if (std::hypot(B.x - A.x, B.y - A.y) < 0.05 * s.h) return false;
        const int mi = M.addPoint(0.5 * (A.x + B.x), 0.5 * (A.y + B.y));
        nodeOf.push_back(-1);
        const int t = M.locate(M.pts[mi], M.vt[a], false);
        std::vector<int> made;
        const long long k = ekey(a, b);
        if (t < 0 || M.insert(mi, t, k, false, &made) != 1) {
            M.pts.pop_back();
            M.vt.pop_back();
            nodeOf.pop_back();
            return false;
        }
        M.cons.erase(k);
        M.cons.insert(ekey(a, mi));
        M.cons.insert(ekey(mi, b));
        for (int q : made) classify(q);
        return true;
    };
    const double B2 = 2.0; // (R / cạnh ngắn nhất)² > 2: góc nhỏ nhất < 20,7°
    for (int pass = 0; pass < 200; ++pass) {
        struct Bad {
            double r2;
            int t, v0, v1, v2;
        };
        std::vector<Bad> bad;
        for (int t = 0; t < (int)M.tris.size(); ++t) {
            const Mesher::Tri& T = M.tris[t];
            if (!T.alive || !T.domain) continue;
            const P2 &a = M.pts[T.v[0]], &b = M.pts[T.v[1]], &c = M.pts[T.v[2]];
            double ox, oy;
            circumcenter(a, b, c, ox, oy);
            const double r2 = (a.x - ox) * (a.x - ox) + (a.y - oy) * (a.y - oy);
            const double l2 = std::min({(a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y),
                                        (b.x - c.x) * (b.x - c.x) + (b.y - c.y) * (b.y - c.y),
                                        (c.x - a.x) * (c.x - a.x) + (c.y - a.y) * (c.y - a.y)});
            const double hc = hAt((a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3);
            const bool big = r2 > 0.36 * hc * hc; // R > 0,6·h (tam giác đều cạnh h: R = 0,577·h)
            const bool ugly = r2 > B2 * l2 && l2 > 0.04 * s.h * s.h;
            if (big || ugly) bad.push_back({r2, t, T.v[0], T.v[1], T.v[2]});
        }
        if (bad.empty()) break;
        std::sort(bad.begin(), bad.end(), [](const Bad& p, const Bad& q) { return p.r2 > q.r2; });
        int inserted = 0;
        for (const Bad& b : bad) {
            if ((int)M.pts.size() >= s.maxPoints) break;
            const Mesher::Tri& T = M.tris[b.t];
            if (!T.alive || T.v[0] != b.v0 || T.v[1] != b.v1 || T.v[2] != b.v2) continue;
            double ox, oy;
            circumcenter(M.pts[b.v0], M.pts[b.v1], M.pts[b.v2], ox, oy);
            int ba = -1, bb = -1;
            const int t = M.locate({ox, oy}, b.t, true, &ba, &bb);
            if (t < 0) {
                if (ba >= 0 && splitSegment(ba, bb)) ++inserted;
                continue;
            }
            const int pi = M.addPoint(ox, oy);
            nodeOf.push_back(-1);
            std::vector<int> made;
            const int r = M.insert(pi, t, -1, true, &made, &ba, &bb);
            if (r != 1) {
                M.pts.pop_back();
                M.vt.pop_back();
                nodeOf.pop_back();
                if (r == 2 && splitSegment(ba, bb)) ++inserted;
                continue;
            }
            for (int q : made) classify(q);
            ++inserted;
        }
        if (!inserted || (int)M.pts.size() >= s.maxPoints) break;
    }

    // ---- xuất phần tử T3
    std::vector<int> nid(M.pts.size(), -1);
    int count = 0;
    for (const auto& T : M.tris) {
        if (!T.alive || !T.domain) continue;
        Element e;
        e.type = ElemType::Tri3;
        e.mat = mat;
        e.t = 1.0;
        e.incompatible = false;
        e.n.fill(-1);
        for (int k = 0; k < 3; ++k) {
            const int v = T.v[k];
            if (nid[v] < 0) {
                if (nodeOf[v] >= 0 && nodeOf[v] < (int)m.nodes.size()) nid[v] = nodeOf[v];
                else {
                    nid[v] = (int)m.nodes.size();
                    m.nodes.push_back({M.pts[v].x, M.pts[v].y});
                }
            }
            e.n[k] = nid[v];
        }
        m.elements.push_back(e);
        ++count;
    }
    if (!count) return fail("Không tạo được phần tử nào.");
    return count;
}

std::vector<std::vector<int>> closedLoops(const Model& m, const std::vector<int>& bars, int* ignored) {
    std::map<int, std::vector<int>> adj;
    for (int e : bars) {
        if (e < 0 || e >= (int)m.elements.size() || !m.elements[e].isBar()) continue;
        const int a = m.elements[e].n[0], b = m.elements[e].n[1];
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
    std::vector<std::vector<int>> loops;
    std::unordered_set<int> seen;
    int bad = 0;
    for (const auto& [start, nbs] : adj) {
        if (seen.count(start)) continue;
        // thành phần liên thông
        std::vector<int> comp{start}, st{start};
        seen.insert(start);
        while (!st.empty()) {
            const int v = st.back();
            st.pop_back();
            for (int w : adj[v])
                if (seen.insert(w).second) {
                    comp.push_back(w);
                    st.push_back(w);
                }
        }
        bool simple = comp.size() >= 3;
        for (int v : comp) simple &= adj[v].size() == 2;
        if (!simple) {
            ++bad;
            continue;
        }
        std::vector<int> loop{start};
        int prev = start, cur = adj[start][0];
        while (cur != start && loop.size() <= comp.size()) {
            loop.push_back(cur);
            const int nx = adj[cur][0] == prev ? adj[cur][1] : adj[cur][0];
            prev = cur;
            cur = nx;
        }
        if (loop.size() == comp.size()) loops.push_back(loop);
        else ++bad;
    }
    if (ignored) *ignored = bad;
    return loops;
}

} // namespace fem
