#include "core/Codes.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>

#include "core/I18n.h"

namespace fem {

namespace {

// Thép chữ I/H cán nóng theo kích thước (mm)
struct IShape {
    const char* name;
    const char* source;
    double h, b, tw, tf, r;
};

// A (cm²), Ix (cm⁴) của tiết diện I có bo góc bán kính r giữa bụng và cánh (khớp bảng tra trong khoảng 0,2%)
ProfilePreset toProfile(const IShape& s) {
    const double pi = 3.14159265358979323846;
    const double hw = s.h - 2 * s.tf;
    const double A = 2 * s.b * s.tf + hw * s.tw + (4 - pi) * s.r * s.r;
    double I = (s.b * s.h * s.h * s.h - (s.b - s.tw) * hw * hw * hw) / 12;
    const double af = (1 - pi / 4) * s.r * s.r, y = s.h / 2 - s.tf - 0.2234 * s.r; // một góc bo: diện tích, trọng tâm
    I += 4 * (0.0075 * std::pow(s.r, 4) + af * y * y);
    return {s.name, s.source, s.h, s.b, A / 100, I / 1e4};
}

// ---------------------------------------------------------------- vật liệu
// {tên, nguồn, E (kN/m²), ν, α, γ (kN/m³), loại (1 thép, 2 bê tông), fd, ft (MPa)}

const std::vector<MaterialPreset> kTcvnMaterials = {
    // TCVN 5575:2012 bảng 5 (t ≤ 20 mm); TCVN 5574:2018 bảng 7, 8, 10
    {"Thép CCT34", "TCVN 5575:2012", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 210, 0},
    {"Thép CCT38", "TCVN 5575:2012", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 230, 0},
    {"Thép CCT42", "TCVN 5575:2012", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 245, 0},
    {"Bê tông B15", "TCVN 5574:2018", 2.40e7, 0.2, 1.0e-5, 25.0, 2, 8.5, 0.75},
    {"Bê tông B20", "TCVN 5574:2018", 2.75e7, 0.2, 1.0e-5, 25.0, 2, 11.5, 0.90},
    {"Bê tông B25", "TCVN 5574:2018", 3.00e7, 0.2, 1.0e-5, 25.0, 2, 14.5, 1.05},
    {"Bê tông B30", "TCVN 5574:2018", 3.25e7, 0.2, 1.0e-5, 25.0, 2, 17.0, 1.15},
    {"Bê tông B35", "TCVN 5574:2018", 3.45e7, 0.2, 1.0e-5, 25.0, 2, 19.5, 1.30},
    {"Bê tông B40", "TCVN 5574:2018", 3.60e7, 0.2, 1.0e-5, 25.0, 2, 22.0, 1.40},
};

const std::vector<MaterialPreset> kEcMaterials = {
    // EN 1993-1-1 bảng 3.1 (t ≤ 40 mm), γM0 = 1,0: fyd = fy
    {"Thép S235", "EN 10025-2", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 235, 0},
    {"Thép S275", "EN 10025-2", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 275, 0},
    {"Thép S355", "EN 10025-2", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 355, 0},
    {"Thép S460", "EN 10025-2", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 460, 0},
    // EN 1992-1-1 bảng 3.1: E = Ecm; fcd = fck / 1,5; fctd = fctk,0.05 / 1,5 (αcc = αct = 1)
    {"Bê tông C20/25", "EN 1992-1-1", 3.0e7, 0.2, 1.0e-5, 25.0, 2, 13.33, 1.00},
    {"Bê tông C25/30", "EN 1992-1-1", 3.1e7, 0.2, 1.0e-5, 25.0, 2, 16.67, 1.20},
    {"Bê tông C30/37", "EN 1992-1-1", 3.3e7, 0.2, 1.0e-5, 25.0, 2, 20.00, 1.33},
    {"Bê tông C35/45", "EN 1992-1-1", 3.4e7, 0.2, 1.0e-5, 25.0, 2, 23.33, 1.47},
    {"Bê tông C40/50", "EN 1992-1-1", 3.5e7, 0.2, 1.0e-5, 25.0, 2, 26.67, 1.67},
};

const std::vector<MaterialPreset> kJisMaterials = {
    // F (cường độ tiêu chuẩn, t ≤ 40 mm) theo thông báo Bộ Xây dựng Nhật số 2464; E = 205000 MPa (AIJ)
    {"Thép SS400", "JIS G 3101", 2.05e8, 0.3, 1.2e-5, 78.5, 1, 235, 0},
    {"Thép SN400B", "JIS G 3136", 2.05e8, 0.3, 1.2e-5, 78.5, 1, 235, 0},
    {"Thép SN490B", "JIS G 3136", 2.05e8, 0.3, 1.2e-5, 78.5, 1, 325, 0},
    {"Thép SM490A", "JIS G 3106", 2.05e8, 0.3, 1.2e-5, 78.5, 1, 325, 0},
    // AIJ RC: E = 3,35·10⁴·(γ/24)²·(Fc/60)^(1/3) với γ = 23 kN/m³ (bê tông thường); fd = Fc
    {"Bê tông Fc21", "AIJ RC", 2.17e7, 0.2, 1.0e-5, 24.0, 2, 21, 0},
    {"Bê tông Fc24", "AIJ RC", 2.27e7, 0.2, 1.0e-5, 24.0, 2, 24, 0},
    {"Bê tông Fc27", "AIJ RC", 2.36e7, 0.2, 1.0e-5, 24.0, 2, 27, 0},
    {"Bê tông Fc30", "AIJ RC", 2.44e7, 0.2, 1.0e-5, 24.0, 2, 30, 0},
    {"Bê tông Fc36", "AIJ RC", 2.60e7, 0.2, 1.0e-5, 24.0, 2, 36, 0},
};

const std::vector<MaterialPreset> kKdsMaterials = {
    // KS D 3503 / KS D 3515 (t ≤ 16 mm), Es = 210000 MPa (KDS 14 31); fd = fy
    {"Thép SS275", "KS D 3503", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 275, 0},
    {"Thép SM275", "KS D 3515", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 275, 0},
    {"Thép SM355", "KS D 3515", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 355, 0},
    {"Thép SM420", "KS D 3515", 2.1e8, 0.3, 1.2e-5, 78.5, 1, 420, 0},
    // KDS 14 20 10: Ec = 8500·∛fcm, fcm = fck + 4 MPa (fck ≤ 40); fd = fck
    {"Bê tông C21", "KDS 14 20 10", 2.49e7, 0.2, 1.0e-5, 24.0, 2, 21, 0},
    {"Bê tông C24", "KDS 14 20 10", 2.58e7, 0.2, 1.0e-5, 24.0, 2, 24, 0},
    {"Bê tông C27", "KDS 14 20 10", 2.67e7, 0.2, 1.0e-5, 24.0, 2, 27, 0},
    {"Bê tông C30", "KDS 14 20 10", 2.75e7, 0.2, 1.0e-5, 24.0, 2, 30, 0},
    {"Bê tông C35", "KDS 14 20 10", 2.88e7, 0.2, 1.0e-5, 24.0, 2, 35, 0},
    {"Bê tông C40", "KDS 14 20 10", 3.00e7, 0.2, 1.0e-5, 24.0, 2, 40, 0},
};

const std::vector<MaterialPreset> kSpMaterials = {
    // СП 16.13330.2017 bảng В.5 (Ry, t nhỏ nhất của mác), E = 206000 MPa; СП 63.13330.2018 bảng 6.8, 6.11
    {"Thép С235", "СП 16.13330", 2.06e8, 0.3, 1.2e-5, 78.5, 1, 230, 0},
    {"Thép С245", "СП 16.13330", 2.06e8, 0.3, 1.2e-5, 78.5, 1, 240, 0},
    {"Thép С255", "СП 16.13330", 2.06e8, 0.3, 1.2e-5, 78.5, 1, 250, 0},
    {"Thép С345", "СП 16.13330", 2.06e8, 0.3, 1.2e-5, 78.5, 1, 340, 0},
    {"Bê tông B15", "СП 63.13330", 2.40e7, 0.2, 1.0e-5, 25.0, 2, 8.5, 0.75},
    {"Bê tông B20", "СП 63.13330", 2.75e7, 0.2, 1.0e-5, 25.0, 2, 11.5, 0.90},
    {"Bê tông B25", "СП 63.13330", 3.00e7, 0.2, 1.0e-5, 25.0, 2, 14.5, 1.05},
    {"Bê tông B30", "СП 63.13330", 3.25e7, 0.2, 1.0e-5, 25.0, 2, 17.0, 1.15},
    {"Bê tông B35", "СП 63.13330", 3.45e7, 0.2, 1.0e-5, 25.0, 2, 19.5, 1.30},
    {"Bê tông B40", "СП 63.13330", 3.60e7, 0.2, 1.0e-5, 25.0, 2, 22.0, 1.40},
};

// ---------------------------------------------------------------- thép hình (h, b, tw, tf, r — mm)

const std::vector<IShape> kEuroShapes = {
    // EN 10365: IPE, HE A, HE B
    {"IPE 100", "EN 10365", 100, 55, 4.1, 5.7, 7},   {"IPE 120", "EN 10365", 120, 64, 4.4, 6.3, 7},
    {"IPE 140", "EN 10365", 140, 73, 4.7, 6.9, 7},   {"IPE 160", "EN 10365", 160, 82, 5.0, 7.4, 9},
    {"IPE 180", "EN 10365", 180, 91, 5.3, 8.0, 9},   {"IPE 200", "EN 10365", 200, 100, 5.6, 8.5, 12},
    {"IPE 220", "EN 10365", 220, 110, 5.9, 9.2, 12}, {"IPE 240", "EN 10365", 240, 120, 6.2, 9.8, 15},
    {"IPE 270", "EN 10365", 270, 135, 6.6, 10.2, 15}, {"IPE 300", "EN 10365", 300, 150, 7.1, 10.7, 15},
    {"IPE 330", "EN 10365", 330, 160, 7.5, 11.5, 18}, {"IPE 360", "EN 10365", 360, 170, 8.0, 12.7, 18},
    {"IPE 400", "EN 10365", 400, 180, 8.6, 13.5, 21}, {"IPE 450", "EN 10365", 450, 190, 9.4, 14.6, 21},
    {"IPE 500", "EN 10365", 500, 200, 10.2, 16.0, 21}, {"IPE 550", "EN 10365", 550, 210, 11.1, 17.2, 24},
    {"IPE 600", "EN 10365", 600, 220, 12.0, 19.0, 24},
    {"HE 100 A", "EN 10365", 96, 100, 5.0, 8.0, 12},   {"HE 120 A", "EN 10365", 114, 120, 5.0, 8.0, 12},
    {"HE 140 A", "EN 10365", 133, 140, 5.5, 8.5, 12},  {"HE 160 A", "EN 10365", 152, 160, 6.0, 9.0, 15},
    {"HE 180 A", "EN 10365", 171, 180, 6.0, 9.5, 15},  {"HE 200 A", "EN 10365", 190, 200, 6.5, 10.0, 18},
    {"HE 220 A", "EN 10365", 210, 220, 7.0, 11.0, 18}, {"HE 240 A", "EN 10365", 230, 240, 7.5, 12.0, 21},
    {"HE 260 A", "EN 10365", 250, 260, 7.5, 12.5, 24}, {"HE 280 A", "EN 10365", 270, 280, 8.0, 13.0, 24},
    {"HE 300 A", "EN 10365", 290, 300, 8.5, 14.0, 27}, {"HE 320 A", "EN 10365", 310, 300, 9.0, 15.5, 27},
    {"HE 340 A", "EN 10365", 330, 300, 9.5, 16.5, 27}, {"HE 360 A", "EN 10365", 350, 300, 10.0, 17.5, 27},
    {"HE 400 A", "EN 10365", 390, 300, 11.0, 19.0, 27}, {"HE 450 A", "EN 10365", 440, 300, 11.5, 21.0, 27},
    {"HE 500 A", "EN 10365", 490, 300, 12.0, 23.0, 27},
    {"HE 100 B", "EN 10365", 100, 100, 6.0, 10.0, 12}, {"HE 120 B", "EN 10365", 120, 120, 6.5, 11.0, 12},
    {"HE 140 B", "EN 10365", 140, 140, 7.0, 12.0, 12}, {"HE 160 B", "EN 10365", 160, 160, 8.0, 13.0, 15},
    {"HE 180 B", "EN 10365", 180, 180, 8.5, 14.0, 15}, {"HE 200 B", "EN 10365", 200, 200, 9.0, 15.0, 18},
    {"HE 220 B", "EN 10365", 220, 220, 9.5, 16.0, 18}, {"HE 240 B", "EN 10365", 240, 240, 10.0, 17.0, 21},
    {"HE 260 B", "EN 10365", 260, 260, 10.0, 17.5, 24}, {"HE 280 B", "EN 10365", 280, 280, 10.5, 18.0, 24},
    {"HE 300 B", "EN 10365", 300, 300, 11.0, 19.0, 27}, {"HE 320 B", "EN 10365", 320, 300, 11.5, 20.5, 27},
    {"HE 340 B", "EN 10365", 340, 300, 12.0, 21.5, 27}, {"HE 360 B", "EN 10365", 360, 300, 12.5, 22.5, 27},
    {"HE 400 B", "EN 10365", 400, 300, 13.5, 24.0, 27}, {"HE 450 B", "EN 10365", 450, 300, 14.0, 26.0, 27},
    {"HE 500 B", "EN 10365", 500, 300, 14.5, 28.0, 27},
};

// Thép H cán nóng: JIS G 3192 (dãy rộng, trung bình, hẹp); KS D 3502 cùng dãy kích thước
std::vector<IShape> hShapes(const char* source) {
    const double d[][5] = {
        {100, 100, 6, 8, 8},     {125, 125, 6.5, 9, 8},   {150, 150, 7, 10, 8},   {175, 175, 7.5, 11, 13},
        {200, 200, 8, 12, 13},   {250, 250, 9, 14, 13},   {300, 300, 10, 15, 13}, {350, 350, 12, 19, 13},
        {400, 400, 13, 21, 22},  {148, 100, 6, 9, 8},     {194, 150, 6, 9, 8},    {244, 175, 7, 11, 13},
        {294, 200, 8, 12, 13},   {340, 250, 9, 14, 13},   {390, 300, 10, 16, 13}, {440, 300, 11, 18, 13},
        {488, 300, 11, 18, 13},  {588, 300, 12, 20, 13},  {150, 75, 5, 7, 8},     {175, 90, 5, 8, 8},
        {200, 100, 5.5, 8, 8},   {250, 125, 6, 9, 8},     {300, 150, 6.5, 9, 13}, {350, 175, 7, 11, 13},
        {400, 200, 8, 13, 13},   {450, 200, 9, 14, 13},   {500, 200, 10, 16, 13}, {600, 200, 11, 17, 13},
    };
    // giữ chuỗi tên suốt chương trình (ProfilePreset giữ const char*); deque không di chuyển phần tử khi thêm vào cuối
    static std::deque<std::string> names;
    std::vector<IShape> v;
    for (const auto& r : d) {
        char b[64];
        std::snprintf(b, sizeof b, "H-%g×%g×%g×%g", r[0], r[1], r[2], r[3]);
        names.emplace_back(b);
        v.push_back({names.back().c_str(), source, r[0], r[1], r[2], r[3], r[4]});
    }
    return v;
}

std::vector<ProfilePreset> profiles(const std::vector<IShape>& shapes, bool withGost) {
    std::vector<ProfilePreset> v;
    for (const auto& s : shapes) v.push_back(toProfile(s));
    if (withGost)
        for (const auto& p : profilePresets()) v.push_back(p);
    return v;
}

// ---------------------------------------------------------------- tổ hợp tải

std::string fmtFactor(double v) {
    char b[32];
    std::snprintf(b, sizeof b, "%g", std::round(v * 1e4) / 1e4);
    std::string s = b;
    const std::string lang = languageCode();
    if (lang == "vi" || lang == "ru")
        for (char& c : s)
            if (c == '.') c = ',';
    return s;
}

struct ComboBuilder {
    Model& m;
    int maxCombos;
    int count = 0;
    bool truncated = false;

    // terms: (trường hợp tải, hệ số); hệ số 0 bị bỏ
    void add(const std::string& label, const std::vector<std::pair<int, double>>& terms) {
        if (count >= maxCombos) {
            truncated = true;
            return;
        }
        Combo c;
        c.generated = true;
        std::string name;
        for (const auto& [i, f] : terms) {
            if (f == 0) continue;
            c.terms.push_back({i, f});
            if (!name.empty()) name += " + ";
            name += (std::abs(f - 1.0) < 1e-12 ? "" : fmtFactor(f) + "·") + m.cases[i].name;
        }
        if (c.terms.empty()) return;
        c.name = label + ": " + name;
        m.combos.push_back(c);
        ++count;
    }
};

// Khóa loại trừ: tải cùng khóa ≠ "" không xuất hiện đồng thời
std::string exclusiveKey(const LoadCase& lc) {
    if (lc.group != 0) return "g" + std::to_string(lc.group);
    const int cat = loadCategory(lc);
    if (cat == CatWind) return "W";
    if (cat == CatSeismic) return "E";
    return "";
}

// Các tập con không rỗng tương thích của [v] (mỗi khóa loại trừ tối đa một tải), xếp theo số phần tử tăng dần
std::vector<std::vector<int>> compatibleSubsets(const Model& m, const std::vector<int>& v, bool& tooMany) {
    std::vector<std::vector<int>> out;
    const int n = (int)v.size();
    tooMany = n > 14;
    if (tooMany) {
        for (int i : v) out.push_back({i});
        return out;
    }
    for (unsigned mask = 1; mask < (1u << n); ++mask) {
        std::vector<int> s;
        std::vector<std::string> keys;
        bool ok = true;
        for (int k = 0; k < n && ok; ++k) {
            if (!(mask & (1u << k))) continue;
            const std::string key = exclusiveKey(m.cases[v[k]]);
            if (!key.empty()) {
                if (std::find(keys.begin(), keys.end(), key) != keys.end()) ok = false;
                keys.push_back(key);
            }
            s.push_back(v[k]);
        }
        if (ok) out.push_back(s);
    }
    std::stable_sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.size() < b.size(); });
    return out;
}

} // namespace

const char* codeId(DesignCode c) {
    static const char* ids[] = {"tcvn", "ec", "jis", "kds", "sp"};
    return ids[(int)c];
}

const char* codeName(DesignCode c) {
    switch (c) {
    case DesignCode::TCVN: return tr("TCVN (Việt Nam)");
    case DesignCode::Eurocode: return tr("Eurocode (EN)");
    case DesignCode::JIS: return tr("JIS / AIJ (Nhật Bản)");
    case DesignCode::KDS: return tr("KDS / KS (Hàn Quốc)");
    case DesignCode::SP: return tr("СП / ГОСТ (Nga)");
    }
    return "";
}

bool parseCode(const std::string& id, DesignCode& out) {
    for (int i = 0; i < kDesignCodeCount; ++i)
        if (id == codeId((DesignCode)i)) {
            out = (DesignCode)i;
            return true;
        }
    return false;
}

DesignCode defaultCodeForLanguage(const std::string& lang) {
    const std::string l = lang.substr(0, 2);
    if (l == "en") return DesignCode::Eurocode;
    if (l == "ja") return DesignCode::JIS;
    if (l == "ko") return DesignCode::KDS;
    if (l == "ru") return DesignCode::SP;
    return DesignCode::TCVN;
}

const std::vector<MaterialPreset>& codeMaterials(DesignCode c) {
    switch (c) {
    case DesignCode::Eurocode: return kEcMaterials;
    case DesignCode::JIS: return kJisMaterials;
    case DesignCode::KDS: return kKdsMaterials;
    case DesignCode::SP: return kSpMaterials;
    default: return kTcvnMaterials;
    }
}

const std::vector<ProfilePreset>& codeProfiles(DesignCode c) {
    // khởi tạo một lần cho mọi tiêu chuẩn (an toàn luồng: hShapes ghi vào bảng tên dùng chung)
    static const std::vector<ProfilePreset> all[kDesignCodeCount] = {
        profiles(hShapes("JIS G 3192"), true), // TCVN: thép H theo JIS (dùng phổ biến) và thép I, [ theo ГОСТ
        profiles(kEuroShapes, false),
        profiles(hShapes("JIS G 3192"), false),
        profiles(hShapes("KS D 3502"), false),
        profiles({}, true),
    };
    return all[(int)c];
}

double defaultGammaF(DesignCode c, int category) {
    if (c == DesignCode::TCVN) { // TCVN 2737:1995 bảng 1, 3; gió 1,2; động đất 1,0
        static const double g[] = {1.1, 1.2, 1.4, 1.2, 1.0};
        return g[std::clamp(category, 0, 4)];
    }
    if (c == DesignCode::SP) { // СП 20.13330.2016 bảng 7.1, 8.2, 10.12, 11.1.12
        static const double g[] = {1.1, 1.2, 1.4, 1.4, 1.0};
        return g[std::clamp(category, 0, 4)];
    }
    return 1.0;
}

int generateCodeCombos(Model& m, DesignCode code, int maxCombos, std::string* note) {
    m.combos.erase(std::remove_if(m.combos.begin(), m.combos.end(), [](const Combo& c) { return c.generated; }),
                   m.combos.end());
    std::vector<int> byCat[5];
    for (size_t i = 0; i < m.cases.size(); ++i) byCat[loadCategory(m.cases[i])].push_back((int)i);
    const auto& D = byCat[CatDead];
    const auto& L = byCat[CatLive];
    const auto& S = byCat[CatSnow];
    const auto& W = byCat[CatWind];
    const auto& E = byCat[CatSeismic];
    std::vector<int> Q; // tải tạm thời (không kể động đất)
    for (int c : {CatLive, CatSnow, CatWind})
        for (int i : byCat[c]) Q.push_back(i);
    std::sort(Q.begin(), Q.end());

    ComboBuilder b{m, maxCombos};
    bool tooMany = false;
    const auto subsets = compatibleSubsets(m, Q, tooMany);
    const bool useGamma = code == DesignCode::TCVN || code == DesignCode::SP;
    auto gf = [&](int i) { return useGamma ? m.cases[i].gammaF : 1.0; };
    // hệ số (trước khi nhân γf) của mỗi trường hợp tải trong tập [s]
    using Terms = std::vector<std::pair<int, double>>;
    auto terms = [&](double fD, const std::vector<std::pair<const std::vector<int>*, double>>& groups) {
        Terms t;
        for (int i : D) t.push_back({i, fD * gf(i)});
        for (const auto& [g, f] : groups)
            for (int i : *g) t.push_back({i, f * gf(i)});
        return t;
    };
    auto with = [&](Terms t, int i, double f) {
        t.push_back({i, f * gf(i)});
        return t;
    };
    // tải tiêu chuẩn (không nhân γf): tổ hợp động đất theo TCVN 9386 / EN 1998
    auto rawTerms = [&](const std::vector<std::pair<const std::vector<int>*, double>>& groups, int e) {
        Terms t;
        for (int i : D) t.push_back({i, 1.0});
        for (const auto& [g, f] : groups)
            for (int i : *g) t.push_back({i, f});
        t.push_back({e, 1.0});
        return t;
    };

    switch (code) {
    case DesignCode::TCVN: {
        const std::string cb1 = tr("Cơ bản 1"), cb2 = tr("Cơ bản 2");
        if (Q.empty() && E.empty()) b.add(cb1, terms(1.0, {}));
        for (int i : Q) b.add(cb1, with(terms(1.0, {}), i, 1.0));
        for (const auto& s : subsets)
            if (s.size() >= 2) b.add(cb2, terms(1.0, {{&s, 0.9}}));
        for (int e : E) b.add(tr("Động đất"), rawTerms({{&L, 0.3}}, e));
        break;
    }
    case DesignCode::SP: {
        const std::string os = tr("Cơ bản");
        if (Q.empty() && E.empty()) b.add(os, terms(1.0, {}));
        for (const auto& s : subsets)
            for (int i : s) { // tải chính (ψ = 1)
                if (s.size() == 1) {
                    b.add(os, with(terms(1.0, {}), i, 1.0));
                    continue;
                }
                for (int j : s) { // tải thứ hai (ψ = 0,9), còn lại ψ = 0,7
                    if (j == i) continue;
                    Terms t = with(with(terms(1.0, {}), i, 1.0), j, 0.9);
                    for (int k : s)
                        if (k != i && k != j) t = with(t, k, 0.7);
                    b.add(os, t);
                }
            }
        for (int e : E) b.add(tr("Động đất"), with(terms(0.9, {{&L, 0.8}, {&S, 0.5}}), e, 1.0 / gf(e)));
        break;
    }
    case DesignCode::Eurocode: {
        auto psi0 = [&](int i) {
            switch (loadCategory(m.cases[i])) {
            case CatSnow: return 0.5;
            case CatWind: return 0.6;
            default: return 0.7;
            }
        };
        if (Q.empty() && E.empty()) b.add("ULS", terms(1.35, {}));
        for (const auto& s : subsets)
            for (int i : s) {
                Terms t = with(terms(1.35, {}), i, 1.5);
                for (int j : s)
                    if (j != i) t = with(t, j, 1.5 * psi0(j));
                b.add("ULS", t);
            }
        for (int e : E) b.add(tr("Động đất"), with(terms(1.0, {{&L, 0.3}}), e, 1.0));
        break;
    }
    case DesignCode::KDS: {
        // nhóm hoạt tải: các tập con tương thích của hoạt tải (không có hoạt tải: một nhóm rỗng)
        bool many = false;
        auto lsets = compatibleSubsets(m, L, many);
        if (lsets.empty()) lsets.push_back({});
        b.add("KDS", terms(1.4, {}));
        for (const auto& l : lsets) b.add("KDS", terms(1.2, {{&l, 1.6}, {&S, 0.5}}));
        if (!S.empty()) {
            for (const auto& l : lsets) b.add("KDS", terms(1.2, {{&S, 1.6}, {&l, 1.0}}));
            for (int w : W) b.add("KDS", with(terms(1.2, {{&S, 1.6}}), w, 0.65));
        }
        for (int w : W)
            for (const auto& l : lsets) b.add("KDS", with(terms(1.2, {{&l, 1.0}, {&S, 0.5}}), w, 1.3));
        for (int e : E)
            for (const auto& l : lsets) b.add("KDS", with(terms(1.2, {{&l, 1.0}, {&S, 0.2}}), e, 1.0));
        for (int w : W) b.add("KDS", with(terms(0.9, {}), w, 1.3));
        for (int e : E) b.add("KDS", with(terms(0.9, {}), e, 1.0));
        break;
    }
    case DesignCode::JIS: {
        bool many = false;
        auto lsets = compatibleSubsets(m, L, many);
        if (lsets.empty()) lsets.push_back({});
        const std::string lt = tr("Dài hạn"), st = tr("Ngắn hạn");
        for (const auto& l : lsets) b.add(lt, terms(1.0, {{&l, 1.0}}));
        for (const auto& l : lsets) {
            if (!S.empty()) b.add(st, terms(1.0, {{&l, 1.0}, {&S, 1.0}}));
            for (int w : W) b.add(st, with(terms(1.0, {{&l, 1.0}}), w, 1.0));
            for (int e : E) b.add(st, with(terms(1.0, {{&l, 1.0}}), e, 1.0));
        }
        break;
    }
    }
    if (note) {
        if (tooMany) *note = tr("Quá nhiều tải tạm thời (> 14): mỗi tổ hợp chỉ xét một tải tạm thời.");
        else if (b.truncated) *note = trf("Số tổ hợp vượt giới hạn %d — đã cắt bớt.", maxCombos);
        else note->clear();
    }
    return b.count;
}

} // namespace fem
