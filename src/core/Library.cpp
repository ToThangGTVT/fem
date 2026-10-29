#include "core/Library.h"

#include "core/SectionGeom.h"
#include "core/I18n.h"

#include <cmath>
#include <cstdio>

namespace fem {

const std::vector<MaterialPreset>& materialPresets() {
    // E bê tông nặng theo TCVN 5574:2018 (tương đương СП 63.13330.2018), thép theo TCVN 5575:2012.
    static const std::vector<MaterialPreset> v = {
        // thép: f với chiều dày t ≤ 20 mm; bê tông: Rb, Rbt (TCVN 5574:2018 bảng 7, 8)
        {"Thép CCT34", "TCVN 5575:2012", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 210, 0},
        {"Thép CCT38", "TCVN 5575:2012", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 230, 0},
        {"Thép С245", "СП 16.13330", 2.06e8, 0.3, 1.2e-5, 78.5, 1, 240, 0},
        {"Bê tông B15", "TCVN 5574:2018", 2.40e7, 0.2, 1.0e-5, 25.0, 2, 8.5, 0.75},
        {"Bê tông B20", "TCVN 5574:2018", 2.75e7, 0.2, 1.0e-5, 25.0, 2, 11.5, 0.90},
        {"Bê tông B25", "TCVN 5574:2018", 3.00e7, 0.2, 1.0e-5, 25.0, 2, 14.5, 1.05},
        {"Bê tông B30", "TCVN 5574:2018", 3.25e7, 0.2, 1.0e-5, 25.0, 2, 17.0, 1.15},
        {"Bê tông B35", "TCVN 5574:2018", 3.45e7, 0.2, 1.0e-5, 25.0, 2, 19.5, 1.30},
        {"Bê tông B40", "TCVN 5574:2018", 3.60e7, 0.2, 1.0e-5, 25.0, 2, 22.0, 1.40},
        {"Nhôm hợp kim", "", 7.0e7, 0.33, 2.3e-5, 27.0, 0, 0, 0},
    };
    return v;
}

const std::vector<ProfilePreset>& profilePresets() {
    static const std::vector<ProfilePreset> v = {
        // Thép chữ I cán nóng, ГОСТ 8239-89
        {"I10", "ГОСТ 8239-89", 100, 55, 12.0, 198},
        {"I12", "ГОСТ 8239-89", 120, 64, 14.7, 350},
        {"I14", "ГОСТ 8239-89", 140, 73, 17.4, 572},
        {"I16", "ГОСТ 8239-89", 160, 81, 20.2, 873},
        {"I18", "ГОСТ 8239-89", 180, 90, 23.4, 1290},
        {"I20", "ГОСТ 8239-89", 200, 100, 26.8, 1840},
        {"I22", "ГОСТ 8239-89", 220, 110, 30.6, 2550},
        {"I24", "ГОСТ 8239-89", 240, 115, 34.8, 3460},
        {"I27", "ГОСТ 8239-89", 270, 125, 40.2, 5010},
        {"I30", "ГОСТ 8239-89", 300, 135, 46.5, 7080},
        {"I33", "ГОСТ 8239-89", 330, 140, 53.8, 9840},
        {"I36", "ГОСТ 8239-89", 360, 145, 61.9, 13380},
        {"I40", "ГОСТ 8239-89", 400, 155, 72.6, 19062},
        {"I45", "ГОСТ 8239-89", 450, 160, 84.7, 27696},
        {"I50", "ГОСТ 8239-89", 500, 170, 100.0, 39727},
        {"I55", "ГОСТ 8239-89", 550, 180, 118.0, 55962},
        {"I60", "ГОСТ 8239-89", 600, 190, 138.0, 76806},
        // Thép chữ [ cán nóng, ГОСТ 8240-97 (dãy У)
        {"[10", "ГОСТ 8240-97", 100, 46, 10.9, 174},
        {"[12", "ГОСТ 8240-97", 120, 52, 13.3, 304},
        {"[14", "ГОСТ 8240-97", 140, 58, 15.6, 491},
        {"[16", "ГОСТ 8240-97", 160, 64, 18.1, 747},
        {"[18", "ГОСТ 8240-97", 180, 70, 20.7, 1090},
        {"[20", "ГОСТ 8240-97", 200, 76, 23.4, 1520},
        {"[22", "ГОСТ 8240-97", 220, 82, 26.7, 2110},
        {"[24", "ГОСТ 8240-97", 240, 90, 30.6, 2900},
        {"[27", "ГОСТ 8240-97", 270, 95, 35.2, 4160},
        {"[30", "ГОСТ 8240-97", 300, 100, 40.5, 5810},
        {"[33", "ГОСТ 8240-97", 330, 105, 46.5, 7980},
        {"[36", "ГОСТ 8240-97", 360, 110, 53.4, 10820},
        {"[40", "ГОСТ 8240-97", 400, 115, 61.5, 15220},
    };
    return v;
}

Material toMaterial(const MaterialPreset& p) {
    // tên trong bảng có sẵn là khóa tra cứu (tiếng Việt); tên vật liệu tạo ra theo ngôn ngữ hiện tại
    Material m;
    m.name = tr(p.name);
    m.E = p.E;
    m.nu = p.nu;
    m.alpha = p.alpha;
    m.gamma = p.gamma;
    m.kind = p.kind;
    m.fd = p.fd;
    m.ft = p.ft;
    return m;
}

Section toSection(const ProfilePreset& p) {
    Section s;
    s.name = p.name;
    s.A = p.A * 1e-4;
    s.I = p.Ix * 1e-8;
    s.h = p.h * 1e-3;
    return s;
}

static std::string dimName(const char* prefix, std::initializer_list<double> v) {
    std::string s = prefix;
    bool first = true;
    for (double x : v) {
        char b[32];
        std::snprintf(b, sizeof b, "%s%g", first ? " " : "x", x);
        s += b;
        first = false;
    }
    return s;
}

Section rectSection(double b, double h) {
    SectionShape s;
    s.kind = ShapeKind::Rect;
    s.p = {b, h, 0, 0, 0, 0};
    Section r = sectionFromShape(dimName(tr("CN"), {b * 100, h * 100}), s);
    // giá trị chính xác (không phụ thuộc đa giác)
    r.A = b * h;
    r.I = b * h * h * h / 12;
    return r;
}

Section circleSection(double d) {
    SectionShape s;
    s.kind = ShapeKind::Circle;
    s.p = {d, 0, 0, 0, 0, 0};
    Section r = sectionFromShape(dimName(tr("Tròn D"), {d * 1000}), s);
    const double pi = 3.14159265358979323846;
    r.A = pi * d * d / 4;
    r.I = pi * std::pow(d, 4) / 64;
    r.yTop = r.yBot = d / 2;
    return r;
}

Section pipeSection(double D, double t) {
    SectionShape s;
    s.kind = ShapeKind::Pipe;
    s.p = {D, t, 0, 0, 0, 0};
    Section r = sectionFromShape(dimName(tr("Ống"), {D * 1000, t * 1000}), s);
    const double pi = 3.14159265358979323846, d = D - 2 * t;
    r.A = pi * (D * D - d * d) / 4;
    r.I = pi * (std::pow(D, 4) - std::pow(d, 4)) / 64;
    r.yTop = r.yBot = D / 2;
    return r;
}

Section shapeSection(ShapeKind k, const std::array<double, 6>& p) {
    SectionShape s;
    s.kind = k;
    s.p = p;
    std::string n;
    switch (k) {
    case ShapeKind::Tee: n = dimName("T", {p[0] * 1000, p[1] * 1000, p[2] * 1000, p[3] * 1000}); break;
    case ShapeKind::IBeam: n = dimName("I", {p[3] * 1000, p[0] * 1000, p[1] * 1000, p[2] * 1000}); break;
    case ShapeKind::Box: n = dimName(tr("Hộp"), {p[0] * 1000, p[1] * 1000, p[2] * 1000}); break;
    case ShapeKind::Channel: n = dimName("[", {p[1] * 1000, p[0] * 1000}); break;
    case ShapeKind::Angle: n = dimName("L", {p[0] * 1000, p[1] * 1000, p[2] * 1000}); break;
    default: n = shapeName(k);
    }
    return sectionFromShape(n, s);
}

} // namespace fem
