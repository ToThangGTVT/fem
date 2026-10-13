#include "core/SectionGeom.h"

#include <algorithm>
#include <cmath>

namespace fem {

namespace {

constexpr double kPi = 3.14159265358979323846;

Contour rect(double x0, double y0, double x1, double y1) { return {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}}; }

Contour circle(double r, int n, bool hole) {
    Contour c;
    for (int i = 0; i < n; ++i) {
        const double a = 2 * kPi * i / n;
        c.push_back({r * std::cos(a), r * std::sin(a)});
    }
    if (hole) std::reverse(c.begin(), c.end());
    return c;
}

double cross(const std::array<double, 2>& a, const std::array<double, 2>& b) { return a[0] * b[1] - b[0] * a[1]; }

double signedArea(const Contour& c) {
    double s = 0;
    for (size_t i = 0; i < c.size(); ++i) s += cross(c[i], c[(i + 1) % c.size()]);
    return s / 2;
}

// Cắt đa giác bởi nửa mặt phẳng y ≥ y0 (Sutherland–Hodgman), giữ chiều đa giác
Contour clipAbove(const Contour& c, double y0) {
    Contour out;
    const size_t n = c.size();
    for (size_t i = 0; i < n; ++i) {
        const auto& p = c[i];
        const auto& q = c[(i + 1) % n];
        const bool pin = p[1] >= y0, qin = q[1] >= y0;
        if (pin) out.push_back(p);
        if (pin != qin) {
            const double t = (y0 - p[1]) / (q[1] - p[1]);
            out.push_back({p[0] + t * (q[0] - p[0]), y0});
        }
    }
    return out;
}

Contour convexHull(std::vector<std::array<double, 2>> pts) {
    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
    if (pts.size() < 3) return pts;
    Contour h(2 * pts.size());
    size_t k = 0;
    auto turn = [](const std::array<double, 2>& o, const std::array<double, 2>& a, const std::array<double, 2>& b) {
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    };
    for (size_t i = 0; i < pts.size(); ++i) {
        while (k >= 2 && turn(h[k - 2], h[k - 1], pts[i]) <= 0) --k;
        h[k++] = pts[i];
    }
    for (size_t i = pts.size() - 1, t = k + 1; i > 0; --i) {
        while (k >= t && turn(h[k - 2], h[k - 1], pts[i - 1]) <= 0) --k;
        h[k++] = pts[i - 1];
    }
    h.resize(k - 1);
    return h; // ngược chiều kim đồng hồ
}

} // namespace

const char* shapeName(ShapeKind k) {
    switch (k) {
    case ShapeKind::None: return "Nhập trực tiếp A, I";
    case ShapeKind::Rect: return "Chữ nhật";
    case ShapeKind::Tee: return "Chữ T";
    case ShapeKind::IBeam: return "Chữ I (có thể không đối xứng)";
    case ShapeKind::Box: return "Hộp";
    case ShapeKind::Pipe: return "Ống tròn";
    case ShapeKind::Circle: return "Tròn đặc";
    case ShapeKind::Channel: return "Chữ [";
    case ShapeKind::Angle: return "Thép góc L";
    case ShapeKind::Polygon: return "Đa giác tùy ý";
    }
    return "";
}

int shapeParamCount(ShapeKind k) {
    switch (k) {
    case ShapeKind::Rect: return 2;
    case ShapeKind::Tee: return 4;
    case ShapeKind::IBeam: return 6;
    case ShapeKind::Box: return 3;
    case ShapeKind::Pipe: return 2;
    case ShapeKind::Circle: return 1;
    case ShapeKind::Channel: return 4;
    case ShapeKind::Angle: return 3;
    default: return 0;
    }
}

const char* shapeParamName(ShapeKind k, int i) {
    static const char* rect[] = {"b", "h"};
    static const char* tee[] = {"b cánh", "t cánh", "t sườn", "h"};
    static const char* ib[] = {"b cánh trên", "t cánh trên", "t bụng", "h", "b cánh dưới", "t cánh dưới"};
    static const char* box[] = {"b", "h", "t thành"};
    static const char* pipe[] = {"D", "t thành"};
    static const char* circ[] = {"D"};
    static const char* ch[] = {"b cánh", "h", "t bụng", "t cánh"};
    static const char* ang[] = {"b (cạnh ngang)", "h (cạnh đứng)", "t"};
    switch (k) {
    case ShapeKind::Rect: return rect[i];
    case ShapeKind::Tee: return tee[i];
    case ShapeKind::IBeam: return ib[i];
    case ShapeKind::Box: return box[i];
    case ShapeKind::Pipe: return pipe[i];
    case ShapeKind::Circle: return circ[i];
    case ShapeKind::Channel: return ch[i];
    case ShapeKind::Angle: return ang[i];
    default: return "";
    }
}

std::vector<Contour> shapeContours(const SectionShape& s) {
    const auto& p = s.p;
    switch (s.kind) {
    case ShapeKind::None: return {};
    case ShapeKind::Rect: return {rect(-p[0] / 2, 0, p[0] / 2, p[1])};
    case ShapeKind::Tee: {
        const double bf = p[0], tf = p[1], bw = p[2], h = p[3];
        return {{{-bw / 2, 0}, {bw / 2, 0}, {bw / 2, h - tf}, {bf / 2, h - tf}, {bf / 2, h}, {-bf / 2, h}, {-bf / 2, h - tf},
                 {-bw / 2, h - tf}}};
    }
    case ShapeKind::IBeam: {
        const double bt = p[0], tt = p[1], tw = p[2], h = p[3], bb = p[4] > 0 ? p[4] : p[0], tb = p[5] > 0 ? p[5] : p[1];
        return {{{-bb / 2, 0},
                 {bb / 2, 0},
                 {bb / 2, tb},
                 {tw / 2, tb},
                 {tw / 2, h - tt},
                 {bt / 2, h - tt},
                 {bt / 2, h},
                 {-bt / 2, h},
                 {-bt / 2, h - tt},
                 {-tw / 2, h - tt},
                 {-tw / 2, tb},
                 {-bb / 2, tb}}};
    }
    case ShapeKind::Box: {
        const double b = p[0], h = p[1], t = p[2];
        Contour in = rect(-b / 2 + t, t, b / 2 - t, h - t);
        std::reverse(in.begin(), in.end());
        return {rect(-b / 2, 0, b / 2, h), in};
    }
    case ShapeKind::Pipe: {
        auto o = circle(p[0] / 2, 180, false), i = circle(p[0] / 2 - p[1], 180, true);
        for (auto& v : o) v[1] += p[0] / 2;
        for (auto& v : i) v[1] += p[0] / 2;
        return {o, i};
    }
    case ShapeKind::Circle: {
        auto o = circle(p[0] / 2, 180, false);
        for (auto& v : o) v[1] += p[0] / 2;
        return {o};
    }
    case ShapeKind::Channel: {
        const double b = p[0], h = p[1], tw = p[2], tf = p[3];
        return {{{0, 0}, {b, 0}, {b, tf}, {tw, tf}, {tw, h - tf}, {b, h - tf}, {b, h}, {0, h}}};
    }
    case ShapeKind::Angle: {
        const double b = p[0], h = p[1], t = p[2];
        return {{{0, 0}, {b, 0}, {b, t}, {t, t}, {t, h}, {0, h}}};
    }
    case ShapeKind::Polygon: return s.contours;
    }
    return {};
}

SectionProps sectionProps(const std::vector<Contour>& cs) {
    SectionProps r;
    double A = 0, Sx = 0, Sy = 0, Ixx = 0, Iyy = 0, Ixy = 0;
    r.xmin = r.ymin = 1e300;
    r.xmax = r.ymax = -1e300;
    std::vector<std::array<double, 2>> all;
    for (const auto& c : cs) {
        const size_t n = c.size();
        if (n < 3) continue;
        for (size_t i = 0; i < n; ++i) {
            const auto& a = c[i];
            const auto& b = c[(i + 1) % n];
            const double cr = cross(a, b);
            A += cr / 2;
            Sx += (a[1] + b[1]) * cr / 6;
            Sy += (a[0] + b[0]) * cr / 6;
            Ixx += (a[1] * a[1] + a[1] * b[1] + b[1] * b[1]) * cr / 12;
            Iyy += (a[0] * a[0] + a[0] * b[0] + b[0] * b[0]) * cr / 12;
            Ixy += (a[0] * b[1] + 2 * a[0] * a[1] + 2 * b[0] * b[1] + b[0] * a[1]) * cr / 24;
            r.xmin = std::min(r.xmin, a[0]);
            r.xmax = std::max(r.xmax, a[0]);
            r.ymin = std::min(r.ymin, a[1]);
            r.ymax = std::max(r.ymax, a[1]);
            all.push_back(a);
        }
    }
    if (!(A > 0)) return r;
    r.ok = true;
    r.A = A;
    r.xc = Sy / A;
    r.yc = Sx / A;
    r.Jx = Ixx - A * r.yc * r.yc;
    r.Jy = Iyy - A * r.xc * r.xc;
    r.Jxy = Ixy - A * r.xc * r.yc;
    const double c = 0.5 * (r.Jx + r.Jy), d = std::hypot(0.5 * (r.Jx - r.Jy), r.Jxy);
    r.J1 = c + d;
    r.J2 = c - d;
    // góc trục chính (theo giáo trình: tg2α = −2Jxy / (Jx − Jy)), trục ứng với J1
    double a2 = std::atan2(-2 * r.Jxy, r.Jx - r.Jy);
    r.alpha = 0.5 * a2 * 180 / kPi;
    r.yTop = r.ymax - r.yc;
    r.yBot = r.yc - r.ymin;
    r.Wtop = r.yTop > 0 ? r.Jx / r.yTop : 0;
    r.Wbot = r.yBot > 0 ? r.Jx / r.yBot : 0;
    r.ix = std::sqrt(std::max(0.0, r.Jx) / A);
    r.iy = std::sqrt(std::max(0.0, r.Jy) / A);

    // ---- lõi tiết diện: mỗi cạnh bao lồi coi như đường trung hòa → một đỉnh lõi
    // Làm trong hệ trục chính trung tâm (u, v), u hợp với x góc α.
    const double al = 0.5 * a2, ca = std::cos(al), sa = std::sin(al);
    const double Ju = r.Jx * ca * ca + r.Jy * sa * sa - 2 * r.Jxy * sa * ca;  // ∫v² dA
    const double Jv = r.Jx * sa * sa + r.Jy * ca * ca + 2 * r.Jxy * sa * ca;  // ∫u² dA
    const double iu2 = Ju / A, iv2 = Jv / A;
    const Contour hull = convexHull(all);
    for (size_t i = 0; i < hull.size(); ++i) {
        auto toUV = [&](const std::array<double, 2>& p) {
            const double x = p[0] - r.xc, y = p[1] - r.yc;
            return std::array<double, 2>{x * ca + y * sa, -x * sa + y * ca};
        };
        const auto P = toUV(hull[i]), Q = toUV(hull[(i + 1) % hull.size()]);
        // đường thẳng qua P, Q: nu·u + nv·v = c
        const double nu = Q[1] - P[1], nv = -(Q[0] - P[0]), cc = nu * P[0] + nv * P[1];
        if (std::abs(cc) < 1e-15) continue;
        const double uk = -iv2 * nu / cc, vk = -iu2 * nv / cc;
        r.kern.push_back({r.xc + uk * ca - vk * sa, r.yc + uk * sa + vk * ca});
    }
    return r;
}

double staticMomentAbove(const std::vector<Contour>& cs, double yc, double y) {
    double S = 0;
    for (const auto& c : cs) {
        const Contour k = clipAbove(c, y);
        const size_t n = k.size();
        if (n < 3) continue;
        for (size_t i = 0; i < n; ++i) {
            const auto& a = k[i];
            const auto& b = k[(i + 1) % n];
            const double cr = cross(a, b);
            S += ((a[1] + b[1]) / 6 - yc / 2) * cr; // ∫(y − yc) dA
        }
    }
    return S;
}

double widthAt(const std::vector<Contour>& cs, double y) {
    double w = 0;
    for (const auto& c : cs) {
        std::vector<double> xs;
        const size_t n = c.size();
        for (size_t i = 0; i < n; ++i) {
            const auto& a = c[i];
            const auto& b = c[(i + 1) % n];
            if ((a[1] <= y && b[1] > y) || (b[1] <= y && a[1] > y)) xs.push_back(a[0] + (y - a[1]) / (b[1] - a[1]) * (b[0] - a[0]));
        }
        std::sort(xs.begin(), xs.end());
        double len = 0;
        for (size_t i = 0; i + 1 < xs.size(); i += 2) len += xs[i + 1] - xs[i];
        w += signedArea(c) >= 0 ? len : -len;
    }
    return std::max(0.0, w);
}

Section sectionFromShape(const std::string& name, const SectionShape& s) {
    Section sec;
    sec.name = name;
    sec.shape = s;
    const auto p = sectionProps(shapeContours(s));
    if (p.ok) {
        sec.A = p.A;
        sec.I = p.Jx;
        sec.h = p.ymax - p.ymin;
        sec.yTop = p.yTop;
        sec.yBot = p.yBot;
    }
    return sec;
}

} // namespace fem
