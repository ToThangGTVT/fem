// Bộ đọc SVG tối giản cho biểu tượng nét: phần tử <path> (M L H V C S Q T A Z, cả dạng tương đối) và <circle>.
// Đường cong, cung được chia thành đoạn thẳng một lần rồi lưu lại; khi vẽ chỉ co giãn và nối polyline.
#include "SvgIcon.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

// Biểu tượng Lucide v0.469.0 (https://lucide.dev, cũng là bộ "Lu" của react-icons).
// ISC License — Copyright (c) Cole Bemis 2013-2022 (Feather, MIT), Lucide Contributors 2022.
const char* kSvg[] = {
    // file-plus
    R"(<path d="M15 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V7Z" />
       <path d="M14 2v4a2 2 0 0 0 2 2h4" />
       <path d="M9 15h6" />
       <path d="M12 18v-6" />)",
    // folder-open
    R"(<path d="m6 14 1.5-2.9A2 2 0 0 1 9.24 10H20a2 2 0 0 1 1.94 2.5l-1.54 6a2 2 0 0 1-1.95 1.5H4a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h3.9a2 2 0 0 1 1.69.9l.81 1.2a2 2 0 0 0 1.67.9H18a2 2 0 0 1 2 2v2" />)",
    // save
    R"(<path d="M15.2 3a2 2 0 0 1 1.4.6l3.8 3.8a2 2 0 0 1 .6 1.4V19a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2z" />
       <path d="M17 21v-7a1 1 0 0 0-1-1H8a1 1 0 0 0-1 1v7" />
       <path d="M7 3v4a1 1 0 0 0 1 1h7" />)",
    // undo-2
    R"(<path d="M9 14 4 9l5-5" />
       <path d="M4 9h10.5a5.5 5.5 0 0 1 5.5 5.5a5.5 5.5 0 0 1-5.5 5.5H11" />)",
    // redo-2
    R"(<path d="m15 14 5-5-5-5" />
       <path d="M20 9H9.5A5.5 5.5 0 0 0 4 14.5A5.5 5.5 0 0 0 9.5 20H13" />)",
    // sun
    R"(<circle cx="12" cy="12" r="4" />
       <path d="M12 2v2" /> <path d="M12 20v2" />
       <path d="m4.93 4.93 1.41 1.41" /> <path d="m17.66 17.66 1.41 1.41" />
       <path d="M2 12h2" /> <path d="M20 12h2" />
       <path d="m6.34 17.66-1.41 1.41" /> <path d="m19.07 4.93-1.41 1.41" />)",
    // moon
    R"(<path d="M12 3a6 6 0 0 0 9 9 9 9 0 1 1-9-9Z" />)",
};
static_assert(sizeof(kSvg) / sizeof(kSvg[0]) == (size_t)SvgIconId::Count, "thiếu dữ liệu biểu tượng");

constexpr float kPi = 3.14159265358979f;

struct Poly {
    std::vector<ImVec2> pts;
    bool closed = false;
};
using Shape = std::vector<Poly>;

// ---------------------------------------------------------------- đọc dữ liệu đường dẫn
struct PathReader {
    const char* s;
    void skip() {
        while (*s && (std::isspace((unsigned char)*s) || *s == ',')) ++s;
    }
    bool atNumber() {
        skip();
        return *s == '-' || *s == '+' || *s == '.' || std::isdigit((unsigned char)*s);
    }
    float number() {
        skip();
        // strtof đọc "1.4.6" thành 1.4 rồi dừng ở ".6" — đúng cách viết rút gọn của SVG
        char* end = nullptr;
        const float v = std::strtof(s, &end);
        s = end;
        return v;
    }
    bool flag() { // cờ cung là một chữ số 0/1, có thể viết liền "01"
        skip();
        const bool f = *s == '1';
        if (*s) ++s;
        return f;
    }
};

void cubic(Poly& p, ImVec2 p0, ImVec2 c1, ImVec2 c2, ImVec2 p1) {
    for (int i = 1; i <= 12; ++i) {
        const float t = i / 12.0f, u = 1 - t;
        const float a = u * u * u, b = 3 * u * u * t, c = 3 * u * t * t, d = t * t * t;
        p.pts.push_back(ImVec2(a * p0.x + b * c1.x + c * c2.x + d * p1.x, a * p0.y + b * c1.y + c * c2.y + d * p1.y));
    }
}

void quad(Poly& p, ImVec2 p0, ImVec2 c, ImVec2 p1) {
    for (int i = 1; i <= 10; ++i) {
        const float t = i / 10.0f, u = 1 - t;
        p.pts.push_back(ImVec2(u * u * p0.x + 2 * u * t * c.x + t * t * p1.x, u * u * p0.y + 2 * u * t * c.y + t * t * p1.y));
    }
}

// Cung elip theo tham số điểm đầu/cuối của SVG (phụ lục F.6.5 của đặc tả SVG)
void arc(Poly& p, ImVec2 p0, float rx, float ry, float phiDeg, bool large, bool sweep, ImVec2 p1) {
    if (rx == 0 || ry == 0) {
        p.pts.push_back(p1);
        return;
    }
    rx = std::fabs(rx);
    ry = std::fabs(ry);
    const float phi = phiDeg * kPi / 180, cp = std::cos(phi), sp = std::sin(phi);
    const float dx = (p0.x - p1.x) / 2, dy = (p0.y - p1.y) / 2;
    const float x1 = cp * dx + sp * dy, y1 = -sp * dx + cp * dy;
    const float lam = x1 * x1 / (rx * rx) + y1 * y1 / (ry * ry);
    if (lam > 1) {
        rx *= std::sqrt(lam);
        ry *= std::sqrt(lam);
    }
    const float num = rx * rx * ry * ry - rx * rx * y1 * y1 - ry * ry * x1 * x1;
    const float den = rx * rx * y1 * y1 + ry * ry * x1 * x1;
    float k = den > 0 ? std::sqrt(std::fmax(0.0f, num / den)) : 0;
    if (large == sweep) k = -k;
    const float cx1 = k * rx * y1 / ry, cy1 = -k * ry * x1 / rx;
    const float cx = cp * cx1 - sp * cy1 + (p0.x + p1.x) / 2, cy = sp * cx1 + cp * cy1 + (p0.y + p1.y) / 2;
    auto ang = [](float ux, float uy, float vx, float vy) {
        return std::atan2(ux * vy - uy * vx, ux * vx + uy * vy);
    };
    const float t0 = ang(1, 0, (x1 - cx1) / rx, (y1 - cy1) / ry);
    float dt = ang((x1 - cx1) / rx, (y1 - cy1) / ry, (-x1 - cx1) / rx, (-y1 - cy1) / ry);
    if (!sweep && dt > 0) dt -= 2 * kPi;
    if (sweep && dt < 0) dt += 2 * kPi;
    const int n = std::max(2, (int)std::ceil(std::fabs(dt) / (kPi / 12)));
    for (int i = 1; i <= n; ++i) {
        const float t = t0 + dt * i / n;
        const float ex = rx * std::cos(t), ey = ry * std::sin(t);
        p.pts.push_back(ImVec2(cp * ex - sp * ey + cx, sp * ex + cp * ey + cy));
    }
}

void parsePath(const char* d, Shape& out) {
    PathReader r{d};
    Poly cur;
    ImVec2 pos(0, 0), start(0, 0), lastCtrl(0, 0);
    char cmd = 0, prev = 0;
    auto flush = [&]() {
        if (cur.pts.size() > 1) out.push_back(cur);
        cur = Poly();
    };
    while (true) {
        r.skip();
        if (!*r.s) break;
        if (std::isalpha((unsigned char)*r.s)) cmd = *r.s++;
        else if (!cmd) break;
        const bool rel = std::islower((unsigned char)cmd);
        const ImVec2 o = rel ? pos : ImVec2(0, 0);
        auto pt = [&]() {
            const float x = r.number(), y = r.number();
            return ImVec2(o.x + x, o.y + y);
        };
        switch (std::toupper((unsigned char)cmd)) {
        case 'M':
            flush();
            pos = start = pt();
            cur.pts.push_back(pos);
            cmd = rel ? 'l' : 'L'; // các cặp số tiếp theo là lệnh L
            break;
        case 'L': pos = pt(); cur.pts.push_back(pos); break;
        case 'H': pos.x = (rel ? pos.x : 0) + r.number(); cur.pts.push_back(pos); break;
        case 'V': pos.y = (rel ? pos.y : 0) + r.number(); cur.pts.push_back(pos); break;
        case 'C': {
            const ImVec2 c1 = pt(), c2 = pt(), p1 = pt();
            cubic(cur, pos, c1, c2, p1);
            lastCtrl = c2;
            pos = p1;
            break;
        }
        case 'S': {
            const bool smooth = std::toupper((unsigned char)prev) == 'C' || std::toupper((unsigned char)prev) == 'S';
            const ImVec2 c1 = smooth ? ImVec2(2 * pos.x - lastCtrl.x, 2 * pos.y - lastCtrl.y) : pos;
            const ImVec2 c2 = pt(), p1 = pt();
            cubic(cur, pos, c1, c2, p1);
            lastCtrl = c2;
            pos = p1;
            break;
        }
        case 'Q': {
            const ImVec2 c = pt(), p1 = pt();
            quad(cur, pos, c, p1);
            lastCtrl = c;
            pos = p1;
            break;
        }
        case 'T': {
            const bool smooth = std::toupper((unsigned char)prev) == 'Q' || std::toupper((unsigned char)prev) == 'T';
            const ImVec2 c = smooth ? ImVec2(2 * pos.x - lastCtrl.x, 2 * pos.y - lastCtrl.y) : pos;
            const ImVec2 p1 = pt();
            quad(cur, pos, c, p1);
            lastCtrl = c;
            pos = p1;
            break;
        }
        case 'A': {
            const float rx = r.number(), ry = r.number(), rot = r.number();
            const bool large = r.flag(), sweep = r.flag();
            const ImVec2 p1 = pt();
            arc(cur, pos, rx, ry, rot, large, sweep, p1);
            pos = p1;
            break;
        }
        case 'Z':
            cur.closed = true;
            pos = start;
            flush();
            cur.pts.push_back(pos);
            break;
        default: return; // lệnh không hỗ trợ
        }
        prev = cmd;
        // lệnh không có tham số (Z) hoặc hết số: chờ chữ cái lệnh tiếp theo
        if (std::toupper((unsigned char)cmd) == 'Z' || !r.atNumber()) cmd = 0;
    }
    flush();
}

std::string attr(const char* b, const char* e, const char* name) {
    const std::string key = std::string(" ") + name + "=\"";
    const std::string el(b, e);
    const size_t i = el.find(key);
    if (i == std::string::npos) return "";
    const size_t j = el.find('"', i + key.size());
    return el.substr(i + key.size(), j - i - key.size());
}

Shape parseSvg(const char* svg) {
    Shape out;
    for (const char* p = std::strchr(svg, '<'); p; p = std::strchr(p + 1, '<')) {
        const char* e = std::strchr(p, '>');
        if (!e) break;
        if (!std::strncmp(p, "<path", 5)) {
            parsePath(attr(p, e, "d").c_str(), out);
        } else if (!std::strncmp(p, "<circle", 7)) {
            const float cx = std::atof(attr(p, e, "cx").c_str()), cy = std::atof(attr(p, e, "cy").c_str());
            const float rr = std::atof(attr(p, e, "r").c_str());
            Poly c;
            c.closed = true;
            for (int i = 0; i < 32; ++i)
                c.pts.push_back(ImVec2(cx + rr * std::cos(i * kPi / 16), cy + rr * std::sin(i * kPi / 16)));
            out.push_back(c);
        }
    }
    return out;
}

const Shape& shape(SvgIconId id) {
    static std::array<Shape, (size_t)SvgIconId::Count> cache;
    static std::array<bool, (size_t)SvgIconId::Count> ready{};
    if (!ready[(int)id]) {
        cache[(int)id] = parseSvg(kSvg[(int)id]);
        ready[(int)id] = true;
    }
    return cache[(int)id];
}

} // namespace

void drawSvgIcon(ImDrawList* dl, SvgIconId id, ImVec2 c, float size, ImU32 col) {
    const float k = size / 24.0f, th = std::max(1.0f, 2.0f * k);
    std::vector<ImVec2> pts;
    for (const Poly& p : shape(id)) {
        pts.clear();
        for (const ImVec2& q : p.pts) {
            const ImVec2 v(c.x + (q.x - 12) * k, c.y + (q.y - 12) * k);
            if (pts.empty() || std::fabs(v.x - pts.back().x) + std::fabs(v.y - pts.back().y) > 0.05f) pts.push_back(v);
        }
        if (p.closed && pts.size() > 2 && std::fabs(pts[0].x - pts.back().x) + std::fabs(pts[0].y - pts.back().y) < 0.05f)
            pts.pop_back();
        if (pts.size() < 2) continue;
        dl->AddPolyline(pts.data(), (int)pts.size(), col, p.closed ? ImDrawFlags_Closed : ImDrawFlags_None, th);
        // nét đầu tròn (stroke-linecap: round) và mối nối tròn ở các góc gấp
        for (size_t i = 0; i < pts.size(); ++i) {
            const bool end = !p.closed && (i == 0 || i + 1 == pts.size());
            bool corner = false;
            if (!end && pts.size() > 2) {
                const ImVec2 a = pts[(i + pts.size() - 1) % pts.size()], b = pts[i], d = pts[(i + 1) % pts.size()];
                const float ux = b.x - a.x, uy = b.y - a.y, vx = d.x - b.x, vy = d.y - b.y;
                const float lu = std::hypot(ux, uy), lv = std::hypot(vx, vy);
                corner = lu > 0 && lv > 0 && (ux * vx + uy * vy) / (lu * lv) < 0.7f;
            }
            if (end || corner) dl->AddCircleFilled(pts[i], th * 0.5f, col, 12);
        }
    }
}
