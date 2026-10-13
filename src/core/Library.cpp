#include "core/Library.h"

#include <cmath>
#include <cstdio>

namespace fem {

const std::vector<MaterialPreset>& materialPresets() {
    // E bê tông nặng theo TCVN 5574:2018 (tương đương СП 63.13330.2018), thép theo TCVN 5575:2012.
    static const std::vector<MaterialPreset> v = {
        {"Thép CCT34", "TCVN 5575:2012", 2.1e8, 0.3, 1.2e-5, 78.5},
        {"Thép CCT38", "TCVN 5575:2012", 2.1e8, 0.3, 1.2e-5, 78.5},
        {"Thép С245", "СП 16.13330", 2.06e8, 0.3, 1.2e-5, 78.5},
        {"Bê tông B15", "TCVN 5574:2018", 2.40e7, 0.2, 1.0e-5, 25.0},
        {"Bê tông B20", "TCVN 5574:2018", 2.75e7, 0.2, 1.0e-5, 25.0},
        {"Bê tông B25", "TCVN 5574:2018", 3.00e7, 0.2, 1.0e-5, 25.0},
        {"Bê tông B30", "TCVN 5574:2018", 3.25e7, 0.2, 1.0e-5, 25.0},
        {"Bê tông B35", "TCVN 5574:2018", 3.45e7, 0.2, 1.0e-5, 25.0},
        {"Bê tông B40", "TCVN 5574:2018", 3.60e7, 0.2, 1.0e-5, 25.0},
        {"Nhôm hợp kim", "", 7.0e7, 0.33, 2.3e-5, 27.0},
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
    Material m;
    m.name = p.name;
    m.E = p.E;
    m.nu = p.nu;
    m.alpha = p.alpha;
    m.gamma = p.gamma;
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

Section rectSection(double b, double h) {
    Section s;
    char buf[64];
    std::snprintf(buf, sizeof buf, "CN %gx%g", b * 100, h * 100);
    s.name = buf;
    s.A = b * h;
    s.I = b * h * h * h / 12;
    s.h = h;
    return s;
}

Section circleSection(double d) {
    Section s;
    char buf[64];
    std::snprintf(buf, sizeof buf, "Tròn D%g", d * 1000);
    s.name = buf;
    const double pi = 3.14159265358979323846;
    s.A = pi * d * d / 4;
    s.I = pi * std::pow(d, 4) / 64;
    s.h = d;
    return s;
}

Section pipeSection(double D, double t) {
    Section s;
    char buf[64];
    std::snprintf(buf, sizeof buf, "Ống %gx%g", D * 1000, t * 1000);
    s.name = buf;
    const double pi = 3.14159265358979323846, d = D - 2 * t;
    s.A = pi * (D * D - d * d) / 4;
    s.I = pi * (std::pow(D, 4) - std::pow(d, 4)) / 64;
    s.h = D;
    return s;
}

} // namespace fem
