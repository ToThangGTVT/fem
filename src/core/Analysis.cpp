#include "core/Analysis.h"

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cmath>
#include <cstdio>
#include <queue>

#include "core/PlaneElements.h"
#include "core/Skyline.h"

namespace fem {

namespace {

std::string fmt(const char* f, ...) __attribute__((format(printf, 1, 2)));
std::string fmt(const char* f, ...) {
    char buf[512];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(buf, sizeof buf, f, ap);
    va_end(ap);
    return buf;
}

struct ElemData {
    std::vector<int> dofs;
    std::vector<double> K; // ma trận độ cứng tổng thể (n×n)
    BarProps bar;
    double D[3][3];
    double x[8], y[8];
};

struct Coupling {
    int eq;    // phương trình tự do
    int dof;   // bậc tự do bị khóa
    double k;
};

std::string validate(const Model& m) {
    if (m.elements.empty()) return "Mô hình chưa có phần tử.";
    if (m.cases.empty()) return "Chưa có trường hợp tải.";
    const int nn = (int)m.nodes.size();
    for (size_t e = 0; e < m.elements.size(); ++e) {
        const Element& el = m.elements[e];
        const int id = (int)e + 1;
        for (int k = 0; k < el.nodeCount(); ++k) {
            if (el.n[k] < 0 || el.n[k] >= nn) return fmt("Phần tử %d: chỉ số nút không hợp lệ.", id);
            for (int j = 0; j < k; ++j)
                if (el.n[j] == el.n[k]) return fmt("Phần tử %d: có hai nút trùng nhau.", id);
        }
        if (el.mat < 0 || el.mat >= (int)m.materials.size()) return fmt("Phần tử %d: chưa gán vật liệu.", id);
        const Material& mat = m.materials[el.mat];
        if (mat.E <= 0) return fmt("Vật liệu \"%s\": E phải > 0.", mat.name.c_str());
        if (el.isBar()) {
            if (el.sec < 0 || el.sec >= (int)m.sections.size()) return fmt("Phần tử %d: chưa gán tiết diện.", id);
            const Section& s = m.sections[el.sec];
            if (s.A <= 0) return fmt("Tiết diện \"%s\": A phải > 0.", s.name.c_str());
            if (el.type == ElemType::Frame && s.I <= 0)
                return fmt("Tiết diện \"%s\": I phải > 0 cho thanh khung.", s.name.c_str());
            if (m.barLength((int)e) <= 1e-12) return fmt("Phần tử %d có chiều dài bằng 0.", id);
        } else {
            if (el.t <= 0) return fmt("Phần tử %d: chiều dày phải > 0.", id);
            if (mat.nu <= -1 || mat.nu >= 0.5) return fmt("Vật liệu \"%s\": hệ số Poisson phải trong (−1; 0,5).", mat.name.c_str());
            if (m.planeMode == PlaneMode::Axisymmetric)
                for (int k = 0; k < el.nodeCount(); ++k)
                    if (m.nodes[el.n[k]].x < -1e-12)
                        return fmt("Bài toán đối xứng trục: phần tử %d có nút với r = x < 0.", id);
        }
    }
    for (const auto& c : m.combos)
        for (const auto& t : c.terms)
            if (t.first < 0 || t.first >= (int)m.cases.size())
                return fmt("Tổ hợp \"%s\" tham chiếu trường hợp tải không tồn tại.", c.name.c_str());
    return {};
}

} // namespace

std::vector<int> rcmOrder(const Model& m) {
    const int nn = (int)m.nodes.size();
    std::vector<std::vector<int>> adj(nn);
    for (const auto& e : m.elements) {
        const int nc = e.nodeCount();
        for (int a = 0; a < nc; ++a)
            for (int b = 0; b < nc; ++b)
                if (a != b && e.n[a] >= 0 && e.n[b] >= 0 && e.n[a] < nn && e.n[b] < nn) adj[e.n[a]].push_back(e.n[b]);
    }
    for (auto& v : adj) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    std::vector<int> order;
    order.reserve(nn);
    std::vector<char> seen(nn, 0);
    std::vector<int> byDeg(nn);
    for (int i = 0; i < nn; ++i) byDeg[i] = i;
    std::stable_sort(byDeg.begin(), byDeg.end(), [&](int a, int b) { return adj[a].size() < adj[b].size(); });
    for (int start : byDeg) {
        if (seen[start]) continue;
        std::queue<int> q;
        q.push(start);
        seen[start] = 1;
        while (!q.empty()) {
            const int v = q.front();
            q.pop();
            order.push_back(v);
            std::vector<int> nb;
            for (int w : adj[v])
                if (!seen[w]) {
                    seen[w] = 1;
                    nb.push_back(w);
                }
            std::stable_sort(nb.begin(), nb.end(), [&](int a, int b) { return adj[a].size() < adj[b].size(); });
            for (int w : nb) q.push(w);
        }
    }
    std::reverse(order.begin(), order.end());
    return order;
}

namespace {

// Ma trận độ cứng tổng thể của phần tử (thanh: theo BarProps hiện hành)
bool elementMatrix(const Model& m, int e, ElemData& d) {
    const Element& el = m.elements[e];
    if (el.isBar()) {
        Mat6 Kg;
        Vec6 fg;
        barGlobal(d.bar, BarLoads{}, Kg, fg);
        d.K.resize(36);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) d.K[i * 6 + j] = Kg[i][j];
        return true;
    }
    return planeStiffness(el, d.x, d.y, m.materials[el.mat], m.planeMode, d.K);
}

} // namespace

AnalysisResult analyze(const Model& m, int samplesPerBar) {
    const auto t0 = std::chrono::steady_clock::now();
    AnalysisResult res;
    res.error = validate(m);
    if (!res.error.empty()) return res;

    const int nn = (int)m.nodes.size(), ne = (int)m.elements.size(), ndof = 3 * nn;
    auto warnOnce = [&](const std::string& w) {
        if (std::find(res.warnings.begin(), res.warnings.end(), w) == res.warnings.end()) res.warnings.push_back(w);
    };

    // ---- ma trận độ cứng phần tử (trạng thái ban đầu)
    std::vector<ElemData> ed(ne);
    std::vector<char> hasStiff(ndof, 0);
    bool anyCompOnly = false;
    for (int e = 0; e < ne; ++e) {
        const Element& el = m.elements[e];
        ElemData& d = ed[e];
        const int nc = el.nodeCount();
        for (int k = 0; k < nc; ++k) {
            d.x[k] = m.nodes[el.n[k]].x;
            d.y[k] = m.nodes[el.n[k]].y;
        }
        if (el.isBar()) {
            d.bar = barProps(m, e);
            for (int k = 0; k < 2; ++k)
                for (int c = 0; c < 3; ++c) d.dofs.push_back(3 * el.n[k] + c);
            anyCompOnly |= el.type == ElemType::Frame && el.kfCompOnly && (el.kf > 0 || el.kG > 0);
        } else {
            planeD(m.materials[el.mat], m.planeMode, d.D);
            for (int k = 0; k < nc; ++k)
                for (int c = 0; c < 2; ++c) d.dofs.push_back(3 * el.n[k] + c);
        }
        if (!elementMatrix(m, e, d)) {
            res.error = fmt("Phần tử %d bị suy biến (diện tích bằng 0, tứ giác lõm/xoắn%s).", e + 1,
                            m.planeMode == PlaneMode::Axisymmetric ? ", hoặc nằm ở phía r < 0" : "");
            return res;
        }
        const int n = (int)d.dofs.size();
        double kmax = 0;
        for (int i = 0; i < n; ++i) kmax = std::max(kmax, std::abs(d.K[i * n + i]));
        for (int i = 0; i < n; ++i)
            if (d.K[i * n + i] > 1e-12 * kmax) hasStiff[d.dofs[i]] = 1;
    }
    const bool nonlinear = (m.secondOrder && m.hasBars()) || anyCompOnly;

    // ---- điều kiện biên
    res.dofState.assign(ndof, DofFree);
    std::vector<double> spring(ndof, 0.0);
    for (const auto& s : m.supports) {
        if (s.node < 0 || s.node >= nn) continue;
        for (int c = 0; c < 3; ++c) {
            if (s.fix[c]) res.dofState[3 * s.node + c] = DofSupport;
            else if (s.k[c] > 0) {
                spring[3 * s.node + c] += s.k[c];
                hasStiff[3 * s.node + c] = 1;
            }
        }
    }
    std::vector<char> used(nn, 0);
    for (const auto& e : m.elements)
        for (int k = 0; k < e.nodeCount(); ++k) used[e.n[k]] = 1;
    if (m.planeMode == PlaneMode::Axisymmetric) {
        double x0, y0, x1, y1;
        m.bounds(x0, y0, x1, y1);
        const double tol = 1e-9 * std::max({1.0, x1 - x0, y1 - y0});
        for (const auto& e : m.elements)
            if (e.isPlane())
                for (int k = 0; k < e.nodeCount(); ++k)
                    if (std::abs(m.nodes[e.n[k]].x) <= tol && res.dofState[3 * e.n[k]] == DofFree)
                        res.dofState[3 * e.n[k]] = DofAuto; // trên trục đối xứng: u_r = 0
    }
    for (int d = 0; d < ndof; ++d) {
        if (res.dofState[d] != DofFree || hasStiff[d]) continue;
        res.dofState[d] = DofAuto;
        const int node = d / 3, c = d % 3;
        if (!used[node]) warnOnce(fmt("Nút %d không thuộc phần tử nào — đã bỏ qua.", node + 1));
        else if (c != RZ) warnOnce(fmt("Nút %d: bậc tự do %s không có độ cứng — đã tự động khóa.", node + 1, dofName(c)));
    }

    // ---- đánh số phương trình theo thứ tự RCM, profile
    std::vector<int> eq(ndof, -1);
    int neq = 0;
    for (int node : rcmOrder(m))
        for (int c = 0; c < 3; ++c)
            if (res.dofState[3 * node + c] == DofFree) eq[3 * node + c] = neq++;
    res.numEquations = neq;
    std::vector<int> first(neq);
    for (int i = 0; i < neq; ++i) first[i] = i;
    for (const auto& d : ed) {
        int mn = neq;
        for (int dof : d.dofs)
            if (eq[dof] >= 0) mn = std::min(mn, eq[dof]);
        for (int dof : d.dofs)
            if (eq[dof] >= 0) first[eq[dof]] = std::min(first[eq[dof]], mn);
    }

    // lắp ghép và phân tích ma trận; trả về phương trình suy biến (−1 nếu thành công)
    auto assemble = [&](const std::vector<ElemData>& eds, SkylineMatrix& K, std::vector<Coupling>& cpl) {
        K.init(first);
        cpl.clear();
        for (const auto& d : eds) {
            const int n = (int)d.dofs.size();
            for (int a = 0; a < n; ++a) {
                const int ea = eq[d.dofs[a]];
                if (ea < 0) continue;
                for (int b = 0; b < n; ++b) {
                    const double k = d.K[a * n + b];
                    if (k == 0) continue;
                    const int eb = eq[d.dofs[b]];
                    if (eb >= 0) {
                        if (ea <= eb) K.add(ea, eb, k);
                    } else if (res.dofState[d.dofs[b]] == DofSupport) {
                        cpl.push_back({ea, d.dofs[b], k});
                    }
                }
            }
        }
        for (int d = 0; d < ndof; ++d)
            if (eq[d] >= 0 && spring[d] > 0) K.add(eq[d], eq[d], spring[d]);
        return K.factorize();
    };
    auto mechanismError = [&](int bad, bool unstable) {
        int dof = 0;
        for (int d = 0; d < ndof; ++d)
            if (eq[d] == bad) dof = d;
        if (unstable)
            return fmt("Mất ổn định trong phân tích bậc hai (tải trọng vượt tải trọng tới hạn) — suy biến tại nút %d, %s.",
                       dof / 3 + 1, dofName(dof % 3));
        return fmt("Hệ biến hình (không ổn định hình học): ma trận độ cứng suy biến tại nút %d, bậc tự do %s. "
                   "Hãy kiểm tra liên kết gối tựa, khớp và cách nối các thanh.",
                   dof / 3 + 1, dofName(dof % 3));
    };

    SkylineMatrix Klin;
    std::vector<Coupling> cplLin;
    {
        const int bad = assemble(ed, Klin, cplLin);
        res.profileSize = Klin.storage();
        if (bad >= 0) {
            res.error = mechanismError(bad, false);
            return res;
        }
    }

    // ---- các bộ tải: trường hợp tải rồi tổ hợp
    struct SetDef {
        std::string name;
        bool combo;
        std::vector<std::pair<int, double>> terms;
    };
    std::vector<SetDef> defs;
    for (size_t i = 0; i < m.cases.size(); ++i) defs.push_back({m.cases[i].name, false, {{(int)i, 1.0}}});
    for (const auto& c : m.combos) defs.push_back({c.name, true, c.terms});

    // giải một bộ tải với ma trận K (đã phân tích) và trạng thái phần tử eds
    auto solveSet = [&](const SetDef& def, const std::vector<ElemData>& eds, const SkylineMatrix& K,
                        const std::vector<Coupling>& cpl, LoadSetResult& out) {
        out = LoadSetResult{};
        out.name = def.name;
        out.combo = def.combo;
        std::vector<double> Fdirect(ndof, 0.0), F(ndof, 0.0), up(ndof, 0.0);
        std::vector<BarLoads> bl(ne);
        double sw = 0;
        for (const auto& [ci, f] : def.terms) {
            const LoadCase& lc = m.cases[ci];
            sw += f * lc.selfWeight;
            for (const auto& l : lc.nodal) {
                if (l.node < 0 || l.node >= nn) continue;
                for (int c = 0; c < 3; ++c) Fdirect[3 * l.node + c] += f * l.F[c];
            }
            for (const auto& l : lc.elem) {
                if (l.elem < 0 || l.elem >= ne) continue;
                if (!m.elements[l.elem].isBar()) {
                    warnOnce(fmt("Tải trên phần tử phẳng %d chưa được hỗ trợ — đã bỏ qua (hãy dùng tải nút).", l.elem + 1));
                    continue;
                }
                addBarLoad(eds[l.elem].bar, l, f, bl[l.elem]);
            }
            for (const auto& d : lc.disp) {
                if (d.node < 0 || d.node >= nn || d.dof < 0 || d.dof > 2) continue;
                const int dof = 3 * d.node + d.dof;
                if (res.dofState[dof] != DofSupport) {
                    warnOnce(fmt("Chuyển vị cưỡng bức tại nút %d (%s) bỏ qua vì bậc tự do này không bị khóa.", d.node + 1,
                                 dofName(d.dof)));
                    continue;
                }
                up[dof] += f * d.value;
            }
        }
        std::vector<std::vector<double>> feq(ne);
        for (int e = 0; e < ne; ++e) {
            const Element& el = m.elements[e];
            const Material& mat = m.materials[el.mat];
            if (el.isBar()) {
                if (sw != 0) addBarSelfWeight(eds[e].bar, sw * mat.gamma * m.sections[el.sec].A, bl[e]);
                if (bl[e].empty()) continue;
                Mat6 Kg;
                Vec6 fg;
                barGlobal(eds[e].bar, bl[e], Kg, fg);
                feq[e].assign(fg.begin(), fg.end());
            } else if (sw != 0) {
                planeBodyForce(el, eds[e].x, eds[e].y, m.planeMode, 0, -sw * mat.gamma, feq[e]);
            }
            for (size_t k = 0; k < feq[e].size(); ++k) F[eds[e].dofs[k]] += feq[e][k];
        }
        for (int d = 0; d < ndof; ++d) {
            F[d] += Fdirect[d];
            if (d % 3 == UX) out.sumFx += F[d];
            if (d % 3 == UY) out.sumFy += F[d];
        }
        std::vector<double> rhs(neq, 0.0);
        for (int d = 0; d < ndof; ++d) {
            if (eq[d] >= 0) rhs[eq[d]] = F[d];
            else if (res.dofState[d] == DofAuto && std::abs(F[d]) > 1e-12)
                warnOnce(fmt("Tải tại nút %d theo %s bị bỏ qua vì bậc tự do này không có độ cứng.", d / 3 + 1, dofName(d % 3)));
        }
        for (const auto& c : cpl) rhs[c.eq] -= c.k * up[c.dof];
        K.solve(rhs);
        out.U.assign(ndof, 0.0);
        for (int d = 0; d < ndof; ++d) {
            if (eq[d] >= 0) out.U[d] = rhs[eq[d]];
            else if (res.dofState[d] == DofSupport) out.U[d] = up[d];
        }

        // ---- nội lực, ứng suất, phản lực
        std::vector<double> nodeForce(ndof, 0.0);
        out.bars.resize(ne);
        out.elemStress.assign(ne, {});
        out.nodeStress.assign(nn, {0, 0, 0, 0});
        out.nodeStressCount.assign(nn, 0);
        for (int e = 0; e < ne; ++e) {
            const Element& el = m.elements[e];
            const ElemData& d = eds[e];
            const int n = (int)d.dofs.size();
            std::vector<double> ue(n);
            for (int k = 0; k < n; ++k) ue[k] = out.U[d.dofs[k]];
            for (int i = 0; i < n; ++i) {
                double s = feq[e].empty() ? 0.0 : -feq[e][i];
                for (int j = 0; j < n; ++j) s += d.K[i * n + j] * ue[j];
                nodeForce[d.dofs[i]] += s;
            }
            if (el.isBar()) {
                Vec6 ug;
                for (int k = 0; k < 6; ++k) ug[k] = ue[k];
                barRecover(d.bar, bl[e], ug, samplesPerBar, out.bars[e]);
                double fx, fy;
                barFoundationForce(d.bar, out.bars[e], fx, fy);
                out.foundRx += fx;
                out.foundRy += fy;
            } else {
                std::array<std::array<double, 4>, 8> sig;
                planeStress(el, d.x, d.y, m.materials[el.mat], m.planeMode, ue.data(), sig);
                for (int k = 0; k < el.nodeCount(); ++k) {
                    out.elemStress[e][k] = sig[k];
                    auto& ns = out.nodeStress[el.n[k]];
                    for (int q = 0; q < 4; ++q) ns[q] += sig[k][q];
                    out.nodeStressCount[el.n[k]]++;
                }
            }
        }
        for (int i = 0; i < nn; ++i)
            if (out.nodeStressCount[i] > 0)
                for (auto& v : out.nodeStress[i]) v /= out.nodeStressCount[i];
        out.R.assign(ndof, 0.0);
        for (int d = 0; d < ndof; ++d) {
            if (res.dofState[d] == DofSupport) out.R[d] = nodeForce[d] - Fdirect[d];
            else if (spring[d] > 0) out.R[d] = -spring[d] * out.U[d];
            if (d % 3 == UX) out.sumRx += out.R[d];
            if (d % 3 == UY) out.sumRy += out.R[d];
        }
        out.sumRx += out.foundRx;
        out.sumRy += out.foundRy;
    };

    for (const auto& def : defs) {
        LoadSetResult out;
        if (!nonlinear) {
            solveSet(def, ed, Klin, cplLin, out);
            res.sets.push_back(std::move(out));
            continue;
        }
        // ---- lặp: phân tích bậc hai (P-Δ) và/hoặc nền chỉ chịu nén
        std::vector<ElemData> eds = ed;
        std::vector<double> N(ne, 0.0);
        std::vector<char> contact(ne, 1);
        SkylineMatrix K;
        std::vector<Coupling> cpl;
        bool converged = false;
        const int maxIt = 40;
        int it = 0;
        for (; it < maxIt; ++it) {
            for (int e = 0; e < ne; ++e) {
                const Element& el = m.elements[e];
                if (!el.isBar()) continue;
                BarProps p = ed[e].bar;
                if (el.kfCompOnly && !contact[e]) p.kf = p.kG = 0;
                if (m.secondOrder && it > 0) {
                    p.second = true;
                    p.N2 = N[e];
                }
                eds[e].bar = p;
                elementMatrix(m, e, eds[e]);
            }
            const int bad = it == 0 && !anyCompOnly ? (K = Klin, cpl = cplLin, -1) : assemble(eds, K, cpl);
            if (bad >= 0) {
                res.error = mechanismError(bad, m.secondOrder && it > 0);
                res.sets.clear();
                return res;
            }
            solveSet(def, eds, K, cpl, out);
            // cập nhật lực dọc và vùng tiếp xúc nền
            double dN = 0, Nmax = 1e-9;
            bool contactChanged = false;
            for (int e = 0; e < ne; ++e) {
                const Element& el = m.elements[e];
                if (!el.isBar()) continue;
                const BarResult& b = out.bars[e];
                if (b.N.empty()) continue;
                double nAvg = 0, vAvg = 0;
                for (double v : b.N) nAvg += v;
                for (double v : b.v) vAvg += v;
                nAvg /= b.N.size();
                vAvg /= b.v.size();
                dN = std::max(dN, std::abs(nAvg - N[e]));
                Nmax = std::max(Nmax, std::abs(nAvg));
                N[e] = nAvg;
                if (el.kfCompOnly && (el.kf > 0 || el.kG > 0)) {
                    const char c = vAvg <= 0; // thanh ép xuống nền (phía −y cục bộ)
                    if (c != contact[e]) contactChanged = true;
                    contact[e] = c;
                }
            }
            const bool nOk = !m.secondOrder || (it > 0 && dN <= 1e-7 * Nmax);
            if (nOk && !contactChanged && (it > 0 || !m.secondOrder)) {
                converged = true;
                break;
            }
        }
        if (!converged)
            warnOnce(fmt("Bộ tải \"%s\": vòng lặp phi tuyến chưa hội tụ sau %d bước.", def.name.c_str(), maxIt));
        if (anyCompOnly) {
            int lost = 0;
            for (int e = 0; e < ne; ++e)
                if (m.elements[e].kfCompOnly && !contact[e]) ++lost;
            if (lost) warnOnce(fmt("Bộ tải \"%s\": %d đoạn thanh tách khỏi nền (nền chỉ chịu nén).", def.name.c_str(), lost));
        }
        res.sets.push_back(std::move(out));
    }

    res.ok = true;
    res.timeMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return res;
}

} // namespace fem
