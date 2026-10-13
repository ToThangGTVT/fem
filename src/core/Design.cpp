#include "core/Design.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "core/Library.h"
#include "core/SectionGeom.h"

namespace fem {

namespace {
constexpr double kPi = 3.14159265358979323846;
}

// ============================================================================ cấu kiện
std::vector<Member> findMembers(const Model& m) {
    const int nn = (int)m.nodes.size(), ne = (int)m.elements.size();
    std::vector<std::vector<int>> at(nn); // các thanh nối vào nút
    for (int e = 0; e < ne; ++e)
        if (m.elements[e].isBar())
            for (int k = 0; k < 2; ++k) at[m.elements[e].n[k]].push_back(e);
    auto dir = [&](int e, int from, double& dx, double& dy) {
        const Element& el = m.elements[e];
        const Node& a = m.nodes[from];
        const Node& b = m.nodes[el.n[0] == from ? el.n[1] : el.n[0]];
        const double L = std::hypot(b.x - a.x, b.y - a.y);
        dx = (b.x - a.x) / L;
        dy = (b.y - a.y) / L;
    };
    auto hingeAt = [&](int e, int node) {
        const Element& el = m.elements[e];
        return el.type != ElemType::Frame || (el.n[0] == node ? el.hingeI : el.hingeJ);
    };
    // nút giữa của cấu kiện: đúng 2 thanh khung thẳng hàng, không gối, không khớp
    auto passThrough = [&](int node, int e, int& next) {
        if (at[node].size() != 2) return false;
        const Support* s = m.supportAt(node);
        if (s && !s->empty()) return false;
        next = at[node][0] == e ? at[node][1] : at[node][0];
        if (m.elements[next].type != ElemType::Frame || m.elements[e].type != ElemType::Frame) return false;
        if (hingeAt(e, node) || hingeAt(next, node)) return false;
        if (m.elements[e].mat != m.elements[next].mat) return false;
        double ax, ay, bx, by;
        dir(e, node, ax, ay);
        dir(next, node, bx, by);
        return ax * bx + ay * by < -0.9999; // ngược hướng = thẳng hàng
    };
    std::vector<char> used(ne, 0);
    std::vector<Member> out;
    for (int e0 = 0; e0 < ne; ++e0) {
        if (used[e0] || !m.elements[e0].isBar()) continue;
        used[e0] = 1;
        std::vector<int> chain = {e0};
        int head = m.elements[e0].n[0], tail = m.elements[e0].n[1];
        // mở rộng về phía đầu và phía cuối
        for (int side = 0; side < 2; ++side) {
            int node = side == 0 ? head : tail, cur = e0, nx;
            while (passThrough(node, cur, nx) && !used[nx]) {
                used[nx] = 1;
                if (side == 0) chain.insert(chain.begin(), nx);
                else chain.push_back(nx);
                const Element& el = m.elements[nx];
                node = el.n[0] == node ? el.n[1] : el.n[0];
                cur = nx;
            }
            (side == 0 ? head : tail) = node;
        }
        Member mb;
        mb.elems = chain;
        mb.nodeI = head;
        mb.nodeJ = tail;
        int node = head;
        for (int e : chain) {
            const Element& el = m.elements[e];
            mb.flip.push_back(el.n[0] != node);
            node = el.n[0] == node ? el.n[1] : el.n[0];
            mb.L += m.barLength(e);
        }
        const Node& a = m.nodes[head];
        const Node& b = m.nodes[tail];
        mb.horizontal = std::abs(b.y - a.y) < 0.35 * std::hypot(b.x - a.x, b.y - a.y);
        out.push_back(mb);
    }
    return out;
}

// ============================================================================ thép
double steelPhi(double lb, double fE) {
    double phi;
    if (lb <= 2.5) phi = 1 - (0.073 - 5.53 * fE) * lb * std::sqrt(lb);
    else if (lb <= 4.5) phi = 1.47 - 13.0 * fE - (0.371 - 27.3 * fE) * lb + (0.0275 - 5.53 * fE) * lb * lb;
    else phi = 332.0 / (lb * lb * (51 - lb));
    return std::clamp(phi, 1e-3, 1.0);
}

std::vector<SteelCheck> checkSteel(const Model& m, const AnalysisResult& r, const std::vector<int>& sets,
                                   const std::vector<Member>& members, const SteelParams& p) {
    std::vector<SteelCheck> out;
    for (size_t k = 0; k < members.size(); ++k) {
        const Member& mb = members[k];
        const Element& el = m.elements[mb.elems.front()];
        const Material& mat = m.materials[el.mat];
        if (mat.kind != 1) continue;
        SteelCheck c;
        c.member = (int)k;
        const Section& sec = m.sections[el.sec];
        if (mat.fd <= 0) {
            c.note = "Chưa có cường độ tính toán f của vật liệu.";
            out.push_back(c);
            continue;
        }
        const double f = mat.fd * 1e3 * p.gammaC; // kN/m²
        const double A = sec.A, I = sec.I;
        const double W = I > 0 ? I / std::max(sec.top(), sec.bot()) : 0;
        const double i = I > 0 ? std::sqrt(I / A) : 0;
        const double mu = k < p.muMember.size() && p.muMember[k] > 0 ? p.muMember[k] : p.mu;
        c.lambda = i > 0 ? mu * mb.L / i : 0;
        c.lambdaBar = c.lambda * std::sqrt(f / mat.E);
        c.phi = c.lambda > 0 ? steelPhi(c.lambdaBar, f / mat.E) : 1.0;
        // ứng suất tiếp tại trục trung hòa (cần hình dạng tiết diện)
        double tauFactor = 0; // τ = V·tauFactor
        const auto cs = shapeContours(sec.shape);
        if (!cs.empty() && I > 0) {
            const auto pr = sectionProps(cs);
            const double b = widthAt(cs, pr.yc);
            if (b > 0) tauFactor = staticMomentAbove(cs, pr.yc, pr.yc) / (I * b);
        }
        if (tauFactor == 0 && el.type == ElemType::Frame) c.note = "Không kiểm tra cắt (tiết diện không có hình dạng).";
        if (I <= 0) c.note = "Thanh không có I: chỉ kiểm tra bền kéo/nén.";
        bool compression = false;
        for (int s : sets) {
            if (s < 0 || s >= (int)r.sets.size()) continue;
            for (int e : mb.elems) {
                const BarResult& b = r.sets[s].bars[e];
                for (size_t q = 0; q < b.x.size(); ++q) {
                    const double N = b.N[q], M = el.type == ElemType::Frame ? b.M[q] : 0, V = el.type == ElemType::Frame ? b.Q[q] : 0;
                    const double st = (std::abs(N) / A + (W > 0 ? std::abs(M) / W : 0)) / f;
                    double stab = 0;
                    if (N < 0) {
                        compression = true;
                        stab = (std::abs(N) / (c.phi * A) + (W > 0 ? std::abs(M) / W : 0)) / f;
                    }
                    const double sh = tauFactor > 0 ? std::abs(V) * tauFactor / (0.58 * f) : 0;
                    const double u = std::max({st, stab, sh});
                    c.strength = std::max(c.strength, st);
                    c.stability = std::max(c.stability, stab);
                    c.shear = std::max(c.shear, sh);
                    if (u >= c.util) {
                        c.util = u;
                        c.N = N;
                        c.M = M;
                        c.V = V;
                        c.set = r.sets[s].name;
                    }
                }
            }
        }
        c.slender = c.lambda > 0 ? c.lambda / (compression ? p.limitCompression : p.limitTension) : 0;
        c.util = std::max(c.util, c.slender);
        out.push_back(c);
    }
    return out;
}

// ============================================================================ bê tông cốt thép
double rcXiR(double Rs, double Es) { return 0.8 / (1 + (Rs / Es) / 0.0035); }

RcBeamDesign rcBeamFlexure(double M, double b, double h, const Material& c, const RcParams& p) {
    RcBeamDesign d;
    d.M = M;
    const double Rb = c.fd * 1e3, Rs = p.Rs * 1e3, Rsc = p.Rsc * 1e3, h0 = h - p.a;
    const double Mabs = std::abs(M);
    if (Rb <= 0 || h0 <= 0) return d;
    const double xiR = rcXiR(p.Rs, p.Es), aR = xiR * (1 - xiR / 2);
    d.alphaM = Mabs / (Rb * b * h0 * h0);
    double As, Asc = 0;
    if (d.alphaM <= aR) {
        d.xi = 1 - std::sqrt(1 - 2 * d.alphaM);
        As = Rb * b * d.xi * h0 / Rs;
    } else {
        d.overReinforced = true;
        d.xi = xiR;
        Asc = (Mabs - aR * Rb * b * h0 * h0) / (Rsc * (h0 - p.a));
        As = (xiR * Rb * b * h0 + Rsc * Asc) / Rs;
    }
    As = std::max(As, p.minRatio * b * h0);
    d.As = As * 1e4;
    d.AsComp = Asc * 1e4;
    return d;
}

RcShearDesign rcShear(double Q, double b, double h, const Material& c, const RcParams& p) {
    RcShearDesign s;
    s.Q = std::abs(Q);
    const double Rb = c.fd * 1e3, Rbt = c.ft * 1e3, Rsw = p.Rsw * 1e3, h0 = h - p.a;
    s.Qmax = 0.3 * Rb * b * h0;
    s.ok = s.Q <= s.Qmax;
    if (s.Q <= 0.5 * Rbt * b * h0 || Rbt <= 0) return s; // bê tông đủ chịu cắt: đai theo cấu tạo
    // Q ≤ Qb + Qsw với c = c0 = 2h0: Qb = 0,75·Rbt·b·h0, Qsw = 1,5·qsw·h0
    s.qswReq = std::max((s.Q - 0.75 * Rbt * b * h0) / (1.5 * h0), 0.25 * Rbt * b);
    s.spacing = std::min({Rsw * p.Asw / s.qswReq, 0.5 * h0, 0.3});
    return s;
}

RcColumnDesign rcColumnSym(double N, double M, double b, double h, double l0, const Material& c, const RcParams& p) {
    RcColumnDesign d;
    d.N = N;
    d.M = M;
    const double Rb = c.fd * 1e3, Rs = p.Rs * 1e3, Rsc = p.Rsc * 1e3, Es = p.Es * 1e3, Eb = c.E;
    const double a = p.a, h0 = h - a;
    if (Rb <= 0 || h0 <= 0) return d;
    const double minAs = p.minRatio * b * h0;
    if (N <= 1e-9) { // kéo hoặc không có lực dọc: tính như cấu kiện chịu uốn
        d.As = rcBeamFlexure(M, b, h, c, p).As;
        return d;
    }
    const double ea = std::max({l0 / 600, h / 30, 0.01});
    d.e0 = std::abs(M) / N + ea;
    const double xiR = rcXiR(p.Rs, p.Es);
    double As = 0.005 * b * h; // giả thiết ban đầu cho độ cứng
    for (int it = 0; it < 4; ++it) {
        // uốn dọc: η = 1 / (1 − N/Ncr) khi l0/h > 4
        d.eta = 1.0;
        if (l0 / h > 4) {
            const double de = std::clamp(d.e0 / h, 0.15, 1.5);
            const double kb = 0.15 / (p.phiL * (0.3 + de));
            const double I = b * h * h * h / 12, Is = 2 * As * (h / 2 - a) * (h / 2 - a);
            const double D = kb * Eb * I + 0.7 * Es * Is;
            const double Ncr = kPi * kPi * D / (l0 * l0);
            if (N >= Ncr) {
                d.unstable = true;
                d.eta = 10;
            } else {
                d.eta = 1 / (1 - N / Ncr);
            }
        }
        const double e = d.eta * d.e0 + (h0 - a) / 2;
        auto AsOf = [&](double x) { return (N * e - Rb * b * x * (h0 - x / 2)) / (Rsc * (h0 - a)); };
        auto sigma = [&](double x) {
            const double xi = x / h0;
            if (xi <= xiR) return Rs;
            return std::clamp(Rs * (2 * (1 - xi) / (1 - xiR) - 1), -Rsc, Rs);
        };
        auto resid = [&](double x) {
            const double A = AsOf(x);
            return Rb * b * x + Rsc * A - sigma(x) * A - N;
        };
        double lo = 1e-6, hi = h;
        double x;
        if (resid(lo) * resid(hi) > 0) x = resid(hi) < 0 ? hi : lo;
        else {
            for (int k = 0; k < 100; ++k) {
                const double mid = 0.5 * (lo + hi);
                if (resid(lo) * resid(mid) <= 0) hi = mid;
                else lo = mid;
            }
            x = 0.5 * (lo + hi);
        }
        d.x = x;
        const double req = AsOf(x);
        As = std::max(req, minAs);
    }
    d.As = As * 1e4;
    return d;
}

std::vector<RcMemberResult> designRc(const Model& m, const AnalysisResult& r, const std::vector<int>& sets,
                                     const std::vector<Member>& members, const RcParams& p) {
    std::vector<RcMemberResult> out;
    for (size_t k = 0; k < members.size(); ++k) {
        const Member& mb = members[k];
        const Element& el = m.elements[mb.elems.front()];
        if (el.type != ElemType::Frame) continue;
        const Material& mat = m.materials[el.mat];
        if (mat.kind != 2) continue;
        RcMemberResult res;
        res.member = (int)k;
        res.column = !mb.horizontal;
        const Section& sec = m.sections[el.sec];
        if (sec.shape.kind != ShapeKind::Rect) {
            res.note = "Chỉ tính cho tiết diện chữ nhật.";
            out.push_back(res);
            continue;
        }
        if (mat.fd <= 0 || mat.ft <= 0) {
            res.note = "Chưa có Rb, Rbt của bê tông.";
            out.push_back(res);
            continue;
        }
        res.b = sec.shape.p[0];
        res.h = sec.shape.p[1];
        double Qmax = 0;
        for (int s : sets) {
            if (s < 0 || s >= (int)r.sets.size()) continue;
            for (int e : mb.elems) {
                const BarResult& b = r.sets[s].bars[e];
                for (size_t q = 0; q < b.x.size(); ++q) {
                    const double N = -b.N[q], M = b.M[q]; // N > 0: nén
                    Qmax = std::max(Qmax, std::abs(b.Q[q]));
                    if (!res.column) {
                        res.Mpos = std::max(res.Mpos, M);
                        res.Mneg = std::min(res.Mneg, M);
                    } else {
                        res.Nmax = std::max(res.Nmax, N);
                        const RcColumnDesign cd = rcColumnSym(N, M, res.b, res.h, p.mu * mb.L, mat, p);
                        if (cd.As > res.AsSide) {
                            res.AsSide = cd.As;
                            res.Mpos = std::abs(M);
                        }
                        if (cd.unstable) res.note = "Lực nén vượt lực tới hạn Ncr — cần tăng tiết diện.";
                    }
                }
            }
        }
        res.Qmax = Qmax;
        if (!res.column) {
            const RcBeamDesign bot = rcBeamFlexure(res.Mpos, res.b, res.h, mat, p);
            const RcBeamDesign top = rcBeamFlexure(res.Mneg, res.b, res.h, mat, p);
            res.AsBottom = std::max(bot.As, top.AsComp);
            res.AsTop = std::max(top.As, bot.AsComp);
            if (bot.overReinforced || top.overReinforced) res.note = "αm > αR: cần cốt chịu nén (đã tính).";
        }
        res.shear = rcShear(Qmax, res.b, res.h, mat, p);
        if (!res.shear.ok) res.note += (res.note.empty() ? "" : " ") + std::string("Q > 0,3·Rb·b·h0: tăng tiết diện.");
        out.push_back(res);
    }
    return out;
}

// ============================================================================ độ võng
std::vector<DeflectionCheck> checkDeflection(const Model& m, const LoadSetResult& s, const std::vector<Member>& members,
                                             double n, bool beamsOnly) {
    std::vector<DeflectionCheck> out;
    for (size_t k = 0; k < members.size(); ++k) {
        const Member& mb = members[k];
        if (beamsOnly && !mb.horizontal) continue;
        if (m.elements[mb.elems.front()].type != ElemType::Frame) continue;
        const Node& A = m.nodes[mb.nodeI];
        const Node& B = m.nodes[mb.nodeJ];
        const double L = mb.L, cm = (B.x - A.x) / L, sm = (B.y - A.y) / L;
        auto trans = [&](double dx, double dy) { return -sm * dx + cm * dy; };
        const double tI = trans(s.U[3 * mb.nodeI], s.U[3 * mb.nodeI + 1]);
        const double tJ = trans(s.U[3 * mb.nodeJ], s.U[3 * mb.nodeJ + 1]);
        DeflectionCheck d;
        d.member = (int)k;
        d.L = L;
        d.limit = L / n;
        double off = 0;
        for (size_t q = 0; q < mb.elems.size(); ++q) {
            const int e = mb.elems[q];
            const Element& el = m.elements[e];
            const Node& a = m.nodes[el.n[0]];
            const Node& b = m.nodes[el.n[1]];
            const double Le = std::hypot(b.x - a.x, b.y - a.y), c = (b.x - a.x) / Le, sn = (b.y - a.y) / Le;
            const BarResult& br = s.bars[e];
            for (size_t i = 0; i < br.x.size(); ++i) {
                const double dx = c * br.u[i] - sn * br.v[i], dy = sn * br.u[i] + c * br.v[i];
                const double pos = off + (mb.flip[q] ? Le - br.x[i] : br.x[i]);
                const double f = std::abs(trans(dx, dy) - (tI + (tJ - tI) * pos / L));
                if (f > d.f) {
                    d.f = f;
                    d.xAt = pos;
                }
            }
            off += Le;
        }
        d.ratio = d.limit > 0 ? d.f / d.limit : 0;
        out.push_back(d);
    }
    return out;
}

// ============================================================================ cọc chịu tải ngang
double pileConventionalWidth(double d) { return d >= 0.8 ? d + 1.0 : 1.5 * d + 0.5; }

void makePile(Model& m, const PileParams& p) {
    if (m.materials.empty()) m.materials.push_back(toMaterial(materialPresets()[5]));
    const int mat = std::clamp(p.mat, 0, (int)m.materials.size() - 1);
    m.sections.push_back(p.square ? rectSection(p.d, p.d) : circleSection(p.d));
    const int sec = (int)m.sections.size() - 1;
    const double bp = pileConventionalWidth(p.d);
    const int nb = std::max(4, p.segments);
    const double dz = p.L / nb;
    const int na = p.l0 > 0 ? std::max(1, (int)std::ceil(p.l0 / dz)) : 0;
    std::vector<int> nodes;
    for (int i = 0; i <= na; ++i) nodes.push_back(m.addNode(0, p.l0 - p.l0 * i / std::max(na, 1)));
    if (na == 0) nodes = {m.addNode(0, 0)};
    for (int i = 1; i <= nb; ++i) nodes.push_back(m.addNode(0, -dz * i));
    for (size_t i = 0; i + 1 < nodes.size(); ++i) {
        Element el;
        el.type = ElemType::Frame;
        el.n = {nodes[i], nodes[i + 1], -1, -1};
        el.mat = mat;
        el.sec = sec;
        const double ya = m.nodes[nodes[i]].y, yb = m.nodes[nodes[i + 1]].y;
        if (yb < 0) {
            const double z = -0.5 * (std::min(0.0, ya) + yb); // độ sâu giữa đoạn
            el.kf = p.m * z * bp;                           // hệ số nền tăng tuyến tính theo độ sâu
        }
        m.elements.push_back(el);
    }
    m.supportFor(nodes.back()).fix = {false, true, false}; // mũi cọc: chặn chuyển vị đứng
    LoadCase lc("Tải ngang đầu cọc");
    lc.nodal.push_back({nodes.front(), {p.H, 0, p.Mh}});
    m.cases.push_back(lc);
}

} // namespace fem
