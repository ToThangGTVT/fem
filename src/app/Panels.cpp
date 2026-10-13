// Các panel nhập liệu, thuộc tính, kết quả và hộp thoại tạo kết cấu.
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <sstream>

#include "App.h"
#include "core/Examples.h"
#include "core/Library.h"
#include "core/Mesher.h"
#include "core/PostProcess.h"
#include "imgui_stdlib.h"

using namespace fem;

const char* contourFieldName(ContourField f);

namespace {

const char* kKinds[] = {"Phân bố", "Tập trung", "Mômen", "Nhiệt độ", "Chế tạo ΔL"};
const char* kDirs[] = {"x cục bộ", "y cục bộ", "X tổng thể", "Y tổng thể"};
const char* kDofs[] = {"ux", "uy", "φz"};

void help(const char* text) {
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    if (ImGui::BeginItemTooltip()) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

std::vector<double> parseList(const std::string& s) {
    std::vector<double> v;
    std::string t = s;
    for (char& c : t)
        if (c == ',' || c == ';') c = ' ';
    std::istringstream in(t);
    std::string tok;
    while (in >> tok) {
        // hỗ trợ dạng 3*6 (3 nhịp 6 m)
        const size_t star = tok.find('*');
        if (star != std::string::npos) {
            const int n = std::atoi(tok.substr(0, star).c_str());
            const double val = std::strtod(tok.substr(star + 1).c_str(), nullptr);
            for (int i = 0; i < n; ++i) v.push_back(val);
        } else {
            v.push_back(std::strtod(tok.c_str(), nullptr));
        }
    }
    v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !(x > 0); }), v.end());
    return v;
}

} // namespace

// ============================================================================ helpers chỉnh sửa có hoàn tác
static bool editD(App& a, const char* label, double& v, const char* fmt = "%.6g") {
    const bool ch = ImGui::InputDouble(label, &v, 0, 0, fmt);
    if (ch) {
        a.dirty = true;
        a.needAnalyze = true;
    }
    a.trackEdit();
    return ch;
}

static bool editS(App& a, const char* label, std::string& s) {
    const bool ch = ImGui::InputText(label, &s);
    if (ch) a.dirty = true;
    a.trackEdit();
    return ch;
}

static bool editB(App& a, const char* label, bool& v) {
    bool t = v;
    if (ImGui::Checkbox(label, &t)) {
        a.snapshot();
        v = t;
        a.changed();
        return true;
    }
    return false;
}

static bool editCombo(App& a, const char* label, int& v, const std::vector<std::string>& items) {
    bool ch = false;
    const char* preview = (v >= 0 && v < (int)items.size()) ? items[v].c_str() : "—";
    if (ImGui::BeginCombo(label, preview)) {
        for (int i = 0; i < (int)items.size(); ++i) {
            ImGui::PushID(i);
            if (ImGui::Selectable(items[i].c_str(), i == v) && i != v) {
                a.snapshot();
                v = i;
                a.changed();
                ch = true;
            }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }
    return ch;
}

static bool editEnum(App& a, const char* label, int cur, const char* const* names, int n, int& out) {
    std::vector<std::string> items(names, names + n);
    int v = cur;
    out = cur;
    const char* preview = (v >= 0 && v < n) ? names[v] : "—";
    bool ch = false;
    if (ImGui::BeginCombo(label, preview)) {
        for (int i = 0; i < n; ++i)
            if (ImGui::Selectable(names[i], i == v) && i != v) {
                a.snapshot();
                out = i;
                ch = true;
            }
        ImGui::EndCombo();
    }
    return ch;
}

static std::vector<std::string> matNames(const Model& m) {
    std::vector<std::string> v;
    for (size_t i = 0; i < m.materials.size(); ++i) v.push_back(std::to_string(i + 1) + ". " + m.materials[i].name);
    return v;
}

static std::vector<std::string> secNames(const Model& m) {
    std::vector<std::string> v;
    for (size_t i = 0; i < m.sections.size(); ++i) v.push_back(std::to_string(i + 1) + ". " + m.sections[i].name);
    return v;
}

static ImGuiTableFlags tableFlags() {
    return ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable |
           ImGuiTableFlags_SizingStretchProp;
}

// ============================================================================ panel mô hình
void App::modelPanel() {
    if (!ImGui::Begin("Mô hình", &showModelPanel)) {
        ImGui::End();
        return;
    }
    Model& m = model;
    if (ImGui::BeginTabBar("modeltabs", ImGuiTabBarFlags_FittingPolicyScroll)) {
        // ------------------------------------------------ chung
        if (ImGui::BeginTabItem("Chung")) {
            ImGui::SetNextItemWidth(-FLT_MIN);
            editS(*this, "##title", m.title);
            int pm = (int)m.planeMode;
            const char* modes[] = {"Ứng suất phẳng (tấm mỏng)", "Biến dạng phẳng (đập, tường chắn)"};
            int out;
            if (editEnum(*this, "Bài toán phẳng", pm, modes, 2, out)) {
                m.planeMode = (PlaneMode)out;
                changed();
            }
            ImGui::SeparatorText("Lưới vẽ");
            ImGui::SetNextItemWidth(100);
            ImGui::InputDouble("Bước lưới (m)", &gridStep, 0, 0, "%.3g");
            gridStep = std::max(gridStep, 1e-3);
            ImGui::Checkbox("Bắt dính lưới", &snapGrid);
            ImGui::SeparatorText("Thống kê");
            int nb = 0, np = 0;
            for (auto& e : m.elements) (e.isBar() ? nb : np)++;
            ImGui::Text("Nút: %zu   Thanh: %d   Phần tử phẳng: %d", m.nodes.size(), nb, np);
            ImGui::Text("Gối: %zu   Trường hợp tải: %zu   Tổ hợp: %zu", m.supports.size(), m.cases.size(), m.combos.size());
            ImGui::SeparatorText("Phân tích");
            editB(*this, "Tính theo sơ đồ biến dạng (P-Δ)", m.secondOrder);
            help("Phân tích bậc hai: kể ảnh hưởng của lực dọc lên độ cứng (P-Δ, P-δ), lặp đến khi hội tụ. Nên chia cột thành vài đoạn.");
            if (result.ok) ImGui::Text("Số phương trình: %d (profile %lld)", result.numEquations, result.profileSize);
            ImGui::EndTabItem();
        }
        // ------------------------------------------------ vật liệu
        if (ImGui::BeginTabItem("Vật liệu")) {
            if (ImGui::Button("+ Thêm")) {
                snapshot();
                m.materials.push_back(Material{});
                changed();
            }
            ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::BeginCombo("##matlib", "Thêm từ thư viện...")) {
                for (const auto& p : materialPresets()) {
                    std::string s = std::string(p.name) + (p.source[0] ? std::string("  (") + p.source + ")" : "");
                    if (ImGui::Selectable(s.c_str())) {
                        snapshot();
                        m.materials.push_back(toMaterial(p));
                        changed();
                    }
                }
                ImGui::EndCombo();
            }
            int del = -1;
            if (ImGui::BeginTable("mats", 6, tableFlags(), ImVec2(0, 0))) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("Tên", 0, 2.2f);
                ImGui::TableSetupColumn("E (kN/m²)", 0, 1.4f);
                ImGui::TableSetupColumn("ν", 0, 0.7f);
                ImGui::TableSetupColumn("α (1/°C)", 0, 1.0f);
                ImGui::TableSetupColumn("γ (kN/m³)", 0, 0.9f);
                ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, 22);
                ImGui::TableHeadersRow();
                for (size_t i = 0; i < m.materials.size(); ++i) {
                    Material& x = m.materials[i];
                    ImGui::PushID((int)i);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    editS(*this, "##n", x.name);
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    editD(*this, "##E", x.E, "%.4g");
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    editD(*this, "##nu", x.nu, "%.3g");
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    editD(*this, "##a", x.alpha, "%.3g");
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    editD(*this, "##g", x.gamma, "%.4g");
                    ImGui::TableNextColumn();
                    if (ImGui::SmallButton("x")) del = (int)i;
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
            if (del >= 0) {
                bool used = false;
                for (auto& e : m.elements) used |= e.mat == del;
                if (used || m.materials.size() == 1) addLog("Không xóa được vật liệu đang được sử dụng.");
                else {
                    snapshot();
                    m.materials.erase(m.materials.begin() + del);
                    for (auto& e : m.elements)
                        if (e.mat > del) e.mat--;
                    changed();
                }
            }
            ImGui::EndTabItem();
        }
        // ------------------------------------------------ tiết diện
        if (ImGui::BeginTabItem("Tiết diện")) {
            static double b = 0.3, h = 0.5, D = 0.2, t = 0.008;
            if (ImGui::Button("Thiết kế tiết diện...")) showSectionDesigner = true;
            ImGui::SameLine();
            if (ImGui::Button("+ Thêm")) {
                snapshot();
                m.sections.push_back(Section{});
                changed();
            }
            ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::BeginCombo("##seclib", "Thép hình ГОСТ...")) {
                for (const auto& p : profilePresets()) {
                    char s[128];
                    std::snprintf(s, sizeof s, "%-5s  A=%.1f cm²  Ix=%.0f cm⁴  (%s)", p.name, p.A, p.Ix, p.source);
                    if (ImGui::Selectable(s)) {
                        snapshot();
                        m.sections.push_back(toSection(p));
                        changed();
                    }
                }
                ImGui::EndCombo();
            }
            ImGui::SetNextItemWidth(70);
            ImGui::InputDouble("b", &b, 0, 0, "%.3g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70);
            ImGui::InputDouble("h (m)", &h, 0, 0, "%.3g");
            ImGui::SameLine();
            if (ImGui::Button("+ Chữ nhật") && b > 0 && h > 0) {
                snapshot();
                m.sections.push_back(rectSection(b, h));
                changed();
            }
            ImGui::SetNextItemWidth(70);
            ImGui::InputDouble("D", &D, 0, 0, "%.3g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70);
            ImGui::InputDouble("t (m)", &t, 0, 0, "%.3g");
            ImGui::SameLine();
            if (ImGui::Button("+ Tròn đặc") && D > 0) {
                snapshot();
                m.sections.push_back(circleSection(D));
                changed();
            }
            ImGui::SameLine();
            if (ImGui::Button("+ Ống") && D > 2 * t && t > 0) {
                snapshot();
                m.sections.push_back(pipeSection(D, t));
                changed();
            }
            int del = -1;
            if (ImGui::BeginTable("secs", 5, tableFlags())) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("Tên", 0, 1.6f);
                ImGui::TableSetupColumn("A (m²)", 0, 1.1f);
                ImGui::TableSetupColumn("I (m⁴)", 0, 1.1f);
                ImGui::TableSetupColumn("h (m)", 0, 0.8f);
                ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, 22);
                ImGui::TableHeadersRow();
                for (size_t i = 0; i < m.sections.size(); ++i) {
                    Section& x = m.sections[i];
                    ImGui::PushID((int)i);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    editS(*this, "##n", x.name);
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    editD(*this, "##A", x.A, "%.4g");
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    editD(*this, "##I", x.I, "%.4g");
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    editD(*this, "##h", x.h, "%.3g");
                    ImGui::TableNextColumn();
                    if (ImGui::SmallButton("x")) del = (int)i;
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
            if (del >= 0) {
                bool used = false;
                for (auto& e : m.elements) used |= e.isBar() && e.sec == del;
                if (used || m.sections.size() == 1) addLog("Không xóa được tiết diện đang được sử dụng.");
                else {
                    snapshot();
                    m.sections.erase(m.sections.begin() + del);
                    for (auto& e : m.elements)
                        if (e.sec > del) e.sec--;
                    changed();
                }
            }
            ImGui::EndTabItem();
        }
        // ------------------------------------------------ nút
        if (ImGui::BeginTabItem("Nút")) {
            static double nx = 0, ny = 0;
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("x", &nx, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("y", &ny, 0, 0, "%.4g");
            ImGui::SameLine();
            if (ImGui::Button("+ Thêm nút")) {
                snapshot();
                m.addNode(nx, ny);
                changed();
            }
            if (ImGui::BeginTable("nodes", 3, tableFlags())) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 44);
                ImGui::TableSetupColumn("x (m)");
                ImGui::TableSetupColumn("y (m)");
                ImGui::TableHeadersRow();
                ImGuiListClipper clip;
                clip.Begin((int)m.nodes.size());
                while (clip.Step())
                    for (int i = clip.DisplayStart; i < clip.DisplayEnd; ++i) {
                        ImGui::PushID(i);
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        char id[16];
                        std::snprintf(id, sizeof id, "%d", i + 1);
                        if (ImGui::Selectable(id, selNodes.count(i) > 0, ImGuiSelectableFlags_SpanAllColumns |
                                                                            ImGuiSelectableFlags_AllowOverlap)) {
                            if (!ImGui::GetIO().KeyShift) {
                                selNodes.clear();
                                selElems.clear();
                            }
                            selNodes.insert(i);
                        }
                        ImGui::TableNextColumn();
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        editD(*this, "##x", m.nodes[i].x, "%.4g");
                        ImGui::TableNextColumn();
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        editD(*this, "##y", m.nodes[i].y, "%.4g");
                        ImGui::PopID();
                    }
                ImGui::EndTable();
            }
            ImGui::EndTabItem();
        }
        // ------------------------------------------------ phần tử
        if (ImGui::BeginTabItem("Phần tử")) {
            const auto mats = matNames(m), secs = secNames(m);
            if (ImGui::BeginTable("elems", 5, tableFlags())) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 44);
                ImGui::TableSetupColumn("Loại", 0, 1.0f);
                ImGui::TableSetupColumn("Nút", 0, 1.0f);
                ImGui::TableSetupColumn("Vật liệu", 0, 1.2f);
                ImGui::TableSetupColumn("Tiết diện / t", 0, 1.2f);
                ImGui::TableHeadersRow();
                ImGuiListClipper clip;
                clip.Begin((int)m.elements.size());
                while (clip.Step())
                    for (int i = clip.DisplayStart; i < clip.DisplayEnd; ++i) {
                        Element& e = m.elements[i];
                        ImGui::PushID(i);
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        char id[16];
                        std::snprintf(id, sizeof id, "%d", i + 1);
                        if (ImGui::Selectable(id, selElems.count(i) > 0, ImGuiSelectableFlags_SpanAllColumns |
                                                                            ImGuiSelectableFlags_AllowOverlap)) {
                            if (!ImGui::GetIO().KeyShift) {
                                selNodes.clear();
                                selElems.clear();
                            }
                            selElems.insert(i);
                        }
                        ImGui::TableNextColumn();
                        ImGui::TextUnformatted(elemTypeName(e.type));
                        ImGui::TableNextColumn();
                        std::string ns;
                        for (int k = 0; k < e.nodeCount(); ++k) ns += (k ? "-" : "") + std::to_string(e.n[k] + 1);
                        if (e.hingeI || e.hingeJ) ns += e.hingeI && e.hingeJ ? " (khớp 2 đầu)" : e.hingeI ? " (khớp i)" : " (khớp j)";
                        ImGui::TextUnformatted(ns.c_str());
                        ImGui::TableNextColumn();
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        editCombo(*this, "##m", e.mat, mats);
                        ImGui::TableNextColumn();
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        if (e.isBar()) editCombo(*this, "##s", e.sec, secs);
                        else editD(*this, "##t", e.t, "%.3g");
                        ImGui::PopID();
                    }
                ImGui::EndTable();
            }
            ImGui::EndTabItem();
        }
        // ------------------------------------------------ gối
        if (ImGui::BeginTabItem("Gối tựa")) {
            ImGui::TextWrapped("Chọn nút trên sơ đồ rồi dùng panel Thuộc tính để gán gối. Lò xo: kN/m, kNm/rad.");
            int del = -1;
            if (ImGui::BeginTable("sups", 8, tableFlags())) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("Nút", ImGuiTableColumnFlags_WidthFixed, 40);
                ImGui::TableSetupColumn("ux", ImGuiTableColumnFlags_WidthFixed, 26);
                ImGui::TableSetupColumn("uy", ImGuiTableColumnFlags_WidthFixed, 26);
                ImGui::TableSetupColumn("φ", ImGuiTableColumnFlags_WidthFixed, 26);
                ImGui::TableSetupColumn("kx");
                ImGui::TableSetupColumn("ky");
                ImGui::TableSetupColumn("kφ");
                ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, 22);
                ImGui::TableHeadersRow();
                for (size_t i = 0; i < m.supports.size(); ++i) {
                    Support& s = m.supports[i];
                    ImGui::PushID((int)i);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::Text("%d", s.node + 1);
                    for (int c = 0; c < 3; ++c) {
                        ImGui::TableNextColumn();
                        ImGui::PushID(c);
                        editB(*this, "##f", s.fix[c]);
                        ImGui::PopID();
                    }
                    for (int c = 0; c < 3; ++c) {
                        ImGui::TableNextColumn();
                        ImGui::PushID(10 + c);
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        editD(*this, "##k", s.k[c], "%.4g");
                        ImGui::PopID();
                    }
                    ImGui::TableNextColumn();
                    if (ImGui::SmallButton("x")) del = (int)i;
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
            if (del >= 0) {
                snapshot();
                m.supports.erase(m.supports.begin() + del);
                changed();
            }
            ImGui::EndTabItem();
        }
        // ------------------------------------------------ tải trọng
        if (ImGui::BeginTabItem("Tải trọng")) {
            std::vector<std::string> cn;
            for (auto& c : m.cases) cn.push_back(c.name);
            ImGui::SetNextItemWidth(-200);
            if (ImGui::BeginCombo("##case", cn[activeCase].c_str())) {
                for (int i = 0; i < (int)cn.size(); ++i) {
                    ImGui::PushID(i);
                    if (ImGui::Selectable(cn[i].c_str(), i == activeCase)) activeCase = resultSet = i;
                    ImGui::PopID();
                }
                ImGui::EndCombo();
            }
            ImGui::SameLine();
            if (ImGui::Button("+ TH")) {
                snapshot();
                m.cases.push_back(LoadCase("TH" + std::to_string(m.cases.size() + 1)));
                activeCase = resultSet = (int)m.cases.size() - 1;
                changed();
            }
            ImGui::SameLine();
            if (ImGui::Button("Nhân bản")) {
                snapshot();
                LoadCase c = m.cases[activeCase];
                c.name += " (2)";
                m.cases.push_back(c);
                activeCase = resultSet = (int)m.cases.size() - 1;
                changed();
            }
            ImGui::SameLine();
            ImGui::BeginDisabled(m.cases.size() <= 1);
            if (ImGui::Button("Xóa TH")) {
                snapshot();
                const int d = activeCase;
                m.cases.erase(m.cases.begin() + d);
                for (auto& c : m.combos) {
                    c.terms.erase(std::remove_if(c.terms.begin(), c.terms.end(), [&](auto& t) { return t.first == d; }),
                                  c.terms.end());
                    for (auto& t : c.terms)
                        if (t.first > d) t.first--;
                }
                activeCase = 0;
                changed();
            }
            ImGui::EndDisabled();
            LoadCase& lc = m.cases[activeCase];
            ImGui::SetNextItemWidth(200);
            editS(*this, "Tên", lc.name);
            ImGui::SetNextItemWidth(80);
            editD(*this, "Hệ số trọng lượng bản thân", lc.selfWeight, "%.3g");
            help("0 = không xét. 1 = kể trọng lượng bản thân (γ·A cho thanh, γ·t cho phần tử phẳng).");
            {
                const char* kinds[] = {"Thường xuyên (tĩnh tải)", "Tạm thời (hoạt tải, gió...)"};
                int out;
                ImGui::SetNextItemWidth(200);
                if (editEnum(*this, "Loại tải", lc.kind, kinds, 2, out)) {
                    lc.kind = out;
                    changed();
                }
                if (lc.kind == 1) {
                    ImGui::SetNextItemWidth(80);
                    if (ImGui::InputInt("Nhóm loại trừ", &lc.group)) {
                        lc.group = std::max(0, lc.group);
                        dirty = true;
                    }
                    trackEdit();
                    help("Các tải tạm thời cùng nhóm (≠ 0) không tác dụng đồng thời, ví dụ gió trái và gió phải cùng nhóm 1.");
                }
                ImGui::SetNextItemWidth(80);
                editD(*this, "Hệ số độ tin cậy γf", lc.gammaF, "%.3g");
                help("Dùng khi tự sinh tổ hợp. Đặt 1 nếu giá trị tải đã nhập là tải tính toán.");
            }

            ImGui::SeparatorText("Tải nút (kN, kNm)");
            int delN = -1, delE = -1, delD = -1;
            if (ImGui::BeginTable("nl", 5, tableFlags(), ImVec2(0, std::min(160.0f, 30 + 26.0f * lc.nodal.size())))) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("Nút", ImGuiTableColumnFlags_WidthFixed, 50);
                ImGui::TableSetupColumn("Fx");
                ImGui::TableSetupColumn("Fy");
                ImGui::TableSetupColumn("Mz");
                ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, 22);
                ImGui::TableHeadersRow();
                for (size_t i = 0; i < lc.nodal.size(); ++i) {
                    auto& l = lc.nodal[i];
                    ImGui::PushID((int)i);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::Text("%d", l.node + 1);
                    for (int c = 0; c < 3; ++c) {
                        ImGui::TableNextColumn();
                        ImGui::PushID(c);
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        editD(*this, "##F", l.F[c], "%.4g");
                        ImGui::PopID();
                    }
                    ImGui::TableNextColumn();
                    if (ImGui::SmallButton("x")) delN = (int)i;
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
            ImGui::SeparatorText("Tải trên thanh");
            if (ImGui::BeginTable("el", 8, tableFlags(), ImVec2(0, std::min(220.0f, 30 + 26.0f * lc.elem.size())))) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("PT", ImGuiTableColumnFlags_WidthFixed, 34);
                ImGui::TableSetupColumn("Loại", 0, 1.1f);
                ImGui::TableSetupColumn("Hướng", 0, 1.1f);
                ImGui::TableSetupColumn("q1/P/M/t", 0, 0.9f);
                ImGui::TableSetupColumn("q2/Δt", 0, 0.9f);
                ImGui::TableSetupColumn("a (m)", 0, 0.7f);
                ImGui::TableSetupColumn("h.c", ImGuiTableColumnFlags_WidthFixed, 26);
                ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, 22);
                ImGui::TableHeadersRow();
                for (size_t i = 0; i < lc.elem.size(); ++i) {
                    auto& l = lc.elem[i];
                    ImGui::PushID((int)i);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::Text("%d", l.elem + 1);
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    int out;
                    if (editEnum(*this, "##k", (int)l.kind, kKinds, 5, out)) {
                        l.kind = (ElemLoadKind)out;
                        changed();
                    }
                    ImGui::TableNextColumn();
                    if (l.kind == ElemLoadKind::Distributed || l.kind == ElemLoadKind::Point) {
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        if (editEnum(*this, "##d", (int)l.dir, kDirs, 4, out)) {
                            l.dir = (LoadDir)out;
                            changed();
                        }
                    }
                    ImGui::TableNextColumn();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    editD(*this, "##q1", l.q1, "%.4g");
                    ImGui::TableNextColumn();
                    if (l.kind == ElemLoadKind::Distributed || l.kind == ElemLoadKind::Temperature) {
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        editD(*this, "##q2", l.q2, "%.4g");
                    }
                    ImGui::TableNextColumn();
                    if (l.kind == ElemLoadKind::Point || l.kind == ElemLoadKind::Moment) {
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        editD(*this, "##a", l.a, "%.4g");
                    }
                    ImGui::TableNextColumn();
                    if (l.kind == ElemLoadKind::Distributed) editB(*this, "##p", l.projected);
                    ImGui::TableNextColumn();
                    if (ImGui::SmallButton("x")) delE = (int)i;
                    ImGui::PopID();
                }
                ImGui::EndTable();
            }
            if (!lc.disp.empty()) {
                ImGui::SeparatorText("Chuyển vị cưỡng bức (m, rad)");
                for (size_t i = 0; i < lc.disp.size(); ++i) {
                    auto& d = lc.disp[i];
                    ImGui::PushID(1000 + (int)i);
                    ImGui::Text("Nút %d  %s =", d.node + 1, kDofs[d.dof]);
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(100);
                    editD(*this, "##v", d.value, "%.4g");
                    ImGui::SameLine();
                    if (ImGui::SmallButton("x")) delD = (int)i;
                    ImGui::PopID();
                }
            }
            if (delN >= 0 || delE >= 0 || delD >= 0) {
                snapshot();
                if (delN >= 0) lc.nodal.erase(lc.nodal.begin() + delN);
                if (delE >= 0) lc.elem.erase(lc.elem.begin() + delE);
                if (delD >= 0) lc.disp.erase(lc.disp.begin() + delD);
                changed();
            }
            ImGui::TextDisabled("Thêm tải: chọn nút/thanh trên sơ đồ → panel Thuộc tính.");
            ImGui::EndTabItem();
        }
        // ------------------------------------------------ tổ hợp
        if (ImGui::BeginTabItem("Tổ hợp")) {
            if (ImGui::Button("+ Tổ hợp")) {
                snapshot();
                Combo c;
                c.name = "THợp " + std::to_string(m.combos.size() + 1);
                for (size_t i = 0; i < m.cases.size(); ++i) c.terms.push_back({(int)i, 1.0});
                m.combos.push_back(c);
                changed();
            }
            help("Tổ hợp tuyến tính các trường hợp tải với hệ số (ví dụ hệ số độ tin cậy, hệ số tổ hợp theo TCVN 2737).");
            {
                static double psi1 = 1.0, psi2 = 0.9;
                ImGui::SeparatorText("Sinh tự động theo TCVN 2737:1995");
                ImGui::SetNextItemWidth(60);
                ImGui::InputDouble("ψ1", &psi1, 0, 0, "%.2f");
                ImGui::SameLine();
                ImGui::SetNextItemWidth(60);
                ImGui::InputDouble("ψ2", &psi2, 0, 0, "%.2f");
                ImGui::SameLine();
                if (ImGui::Button("Sinh tổ hợp")) {
                    snapshot();
                    std::string note;
                    const int n = generateCombos(m, psi1, psi2, 300, &note);
                    addLog("Đã sinh " + std::to_string(n) + " tổ hợp (cơ bản 1 và 2)." + (note.empty() ? "" : " " + note));
                    changed();
                }
                help("Cơ bản 1: tải thường xuyên + một tải tạm thời (ψ1). Cơ bản 2: tải thường xuyên + ψ2 × (từ hai tải "
                     "tạm thời trở lên, mỗi nhóm loại trừ tối đa một). Mỗi tải nhân hệ số γf của nó. Phân loại tải ở thẻ "
                     "Tải trọng. Các tổ hợp tự sinh trước đó được thay thế; tổ hợp nhập tay được giữ.");
                ImGui::Separator();
            }
            std::vector<std::string> cn;
            for (auto& c : m.cases) cn.push_back(c.name);
            int delC = -1;
            for (size_t i = 0; i < m.combos.size(); ++i) {
                Combo& c = m.combos[i];
                ImGui::PushID((int)i);
                ImGui::Separator();
                ImGui::SetNextItemWidth(-60);
                editS(*this, "##name", c.name);
                if (c.generated && ImGui::IsItemHovered()) ImGui::SetTooltip("Tổ hợp tự sinh");
                ImGui::SameLine();
                if (ImGui::SmallButton("Xóa")) delC = (int)i;
                int delT = -1;
                for (size_t k = 0; k < c.terms.size(); ++k) {
                    ImGui::PushID((int)k);
                    ImGui::SetNextItemWidth(80);
                    editD(*this, "##f", c.terms[k].second, "%.3g");
                    ImGui::SameLine();
                    ImGui::Text("×");
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(-30);
                    editCombo(*this, "##c", c.terms[k].first, cn);
                    ImGui::SameLine();
                    if (ImGui::SmallButton("x")) delT = (int)k;
                    ImGui::PopID();
                }
                if (ImGui::SmallButton("+ thành phần")) {
                    snapshot();
                    c.terms.push_back({0, 1.0});
                    changed();
                }
                if (delT >= 0) {
                    snapshot();
                    c.terms.erase(c.terms.begin() + delT);
                    changed();
                }
                ImGui::PopID();
            }
            if (delC >= 0) {
                snapshot();
                m.combos.erase(m.combos.begin() + delC);
                changed();
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
}

// ============================================================================ panel thuộc tính đối tượng chọn
void App::propsPanel() {
    if (!ImGui::Begin("Thuộc tính", &showPropsPanel)) {
        ImGui::End();
        return;
    }
    Model& m = model;
    const auto mats = matNames(m), secs = secNames(m);
    LoadCase& lc = m.cases[activeCase];

    if (selNodes.empty() && selElems.empty()) {
        ImGui::TextWrapped("Chưa chọn đối tượng. Nhấp hoặc kéo khung trên sơ đồ để chọn nút / phần tử.");
        ImGui::SeparatorText("Thanh vẽ mới");
        ImGui::SetNextItemWidth(-90);
        if (ImGui::BeginCombo("Vật liệu", curMat < (int)mats.size() ? mats[curMat].c_str() : "—")) {
            for (int i = 0; i < (int)mats.size(); ++i)
                if (ImGui::Selectable(mats[i].c_str(), i == curMat)) curMat = i;
            ImGui::EndCombo();
        }
        ImGui::SetNextItemWidth(-90);
        if (ImGui::BeginCombo("Tiết diện", curSec < (int)secs.size() ? secs[curSec].c_str() : "—")) {
            for (int i = 0; i < (int)secs.size(); ++i)
                if (ImGui::Selectable(secs[i].c_str(), i == curSec)) curSec = i;
            ImGui::EndCombo();
        }
        ImGui::SeparatorText("Thao tác");
        for (const char* t : {"Gõ lệnh ở bất cứ đâu: L (thanh), DAN, PO, E, M, CO, RO, MI, DIV, NGAM, TAIPB... Gõ ? để xem tất cả.",
                              "Enter / Space / chuột phải: kết thúc hoặc lặp lại lệnh; Esc: hủy.",
                              "Tọa độ: x,y  @dx,dy  @d<góc, hoặc gõ khoảng cách theo hướng con trỏ.",
                              "Chọn: nhấp để thêm, Shift+nhấp để bỏ; kéo trái→phải: cửa sổ, phải→trái: cắt ngang.",
                              "Nhấp grip (ô vuông xanh) để kéo nút. Lăn chuột: zoom; kéo chuột giữa/phải: pan; nhấp đúp chuột giữa: vừa khung.",
                              "F3 OSNAP, F7 lưới, F8 ORTHO, F9 SNAP, F12 DYN, Ctrl+1 thuộc tính, F5 tính."}) {
            ImGui::Bullet();
            ImGui::TextWrapped("%s", t);
        }
        ImGui::End();
        return;
    }

    // ------------------------------------------------ nút
    if (!selNodes.empty()) {
        ImGui::SeparatorText(("Nút đã chọn: " + std::to_string(selNodes.size())).c_str());
        if (selNodes.size() == 1) {
            Node& n = m.nodes[*selNodes.begin()];
            ImGui::Text("Nút %d", *selNodes.begin() + 1);
            ImGui::SetNextItemWidth(90);
            editD(*this, "x (m)", n.x, "%.4g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(90);
            editD(*this, "y (m)", n.y, "%.4g");
        }
        auto setFix = [&](bool x, bool y, bool r) {
            snapshot();
            for (int n : selNodes) m.supportFor(n).fix = {x, y, r};
            m.cleanSupports();
            changed();
        };
        ImGui::TextUnformatted("Gối tựa:");
        if (ImGui::Button("Ngàm")) setFix(true, true, true);
        ImGui::SameLine();
        if (ImGui::Button("Gối cố định")) setFix(true, true, false);
        ImGui::SameLine();
        if (ImGui::Button("Di động ↕")) setFix(false, true, false);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Gối di động: chặn chuyển vị đứng");
        ImGui::SameLine();
        if (ImGui::Button("Di động ↔")) setFix(true, false, false);
        if (ImGui::Button("Ngàm trượt")) setFix(true, false, true);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Chặn ux và φ, cho phép uy");
        ImGui::SameLine();
        if (ImGui::Button("Bỏ gối")) {
            snapshot();
            for (int n : selNodes)
                if (Support* s = m.supportAt(n)) *s = Support{n};
            m.cleanSupports();
            changed();
        }
        static double kk[3] = {0, 0, 0};
        ImGui::SetNextItemWidth(220);
        ImGui::InputScalarN("kx ky kφ", ImGuiDataType_Double, kk, 3, nullptr, nullptr, "%.4g");
        ImGui::SameLine();
        if (ImGui::Button("Gán lò xo")) {
            snapshot();
            for (int n : selNodes) m.supportFor(n).k = {kk[0], kk[1], kk[2]};
            m.cleanSupports();
            changed();
        }
        help("Gối đàn hồi: kN/m (kx, ky), kNm/rad (kφ). Chỉ có tác dụng ở bậc tự do không bị khóa.");

        ImGui::SeparatorText(("Tải nút — " + lc.name).c_str());
        static double F[3] = {0, -10, 0};
        ImGui::SetNextItemWidth(220);
        ImGui::InputScalarN("Fx Fy Mz", ImGuiDataType_Double, F, 3, nullptr, nullptr, "%.4g");
        help("Fx, Fy (kN) theo trục tổng thể (Y hướng lên: tải trọng đứng xuống là Fy < 0). Mz (kNm) dương ngược chiều kim đồng hồ.");
        if (ImGui::Button("Thêm tải nút")) {
            snapshot();
            for (int n : selNodes) lc.nodal.push_back({n, {F[0], F[1], F[2]}});
            changed();
        }
        ImGui::SameLine();
        if (ImGui::Button("Xóa tải nút")) {
            snapshot();
            lc.nodal.erase(std::remove_if(lc.nodal.begin(), lc.nodal.end(),
                                          [&](const NodalLoad& l) { return selNodes.count(l.node) > 0; }),
                           lc.nodal.end());
            changed();
        }
        static int dd = 1;
        static double dv = -0.01;
        ImGui::SetNextItemWidth(60);
        ImGui::Combo("##dof", &dd, kDofs, 3);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(90);
        ImGui::InputDouble("##dv", &dv, 0, 0, "%.4g");
        ImGui::SameLine();
        if (ImGui::Button("Chuyển vị cưỡng bức")) {
            snapshot();
            for (int n : selNodes) lc.disp.push_back({n, dd, dv});
            changed();
        }
        help("Lún / chuyển vị gối (m hoặc rad). Bậc tự do phải được khóa bởi gối tựa.");

        if (m.hasPlane() && selNodes.size() >= 2) {
            ImGui::SeparatorText(("Tải trên cạnh phần tử phẳng — " + lc.name).c_str());
            static EdgeLoadSpec spec;
            static int mode = 0;
            ImGui::RadioButton("Lực qx, qy", &mode, 0);
            ImGui::SameLine();
            ImGui::RadioButton("Áp lực pháp tuyến", &mode, 1);
            if (mode == 0) {
                double q[2] = {spec.qx, spec.qy};
                ImGui::SetNextItemWidth(160);
                if (ImGui::InputScalarN("qx qy (kN/m)", ImGuiDataType_Double, q, 2, nullptr, nullptr, "%.4g")) {
                    spec.qx = q[0];
                    spec.qy = q[1];
                }
            } else {
                ImGui::SetNextItemWidth(90);
                ImGui::InputDouble("p0", &spec.p0, 0, 0, "%.4g");
                ImGui::SameLine();
                ImGui::SetNextItemWidth(90);
                ImGui::InputDouble("k (p = p0 + k·y)", &spec.ky, 0, 0, "%.4g");
                ImGui::Checkbox("Bỏ phần p < 0", &spec.clipNegative);
                help("Áp lực kN/m² (dương: ép vào vật thể), nhân chiều dày phần tử. Áp lực nước mực H: p0 = γH, k = −γ.");
            }
            if (ImGui::Button("Đặt tải lên cạnh biên giữa các nút chọn")) {
                spec.pressure = mode == 1;
                snapshot();
                const int n = applyEdgeLoad(m, activeCase, std::vector<int>(selNodes.begin(), selNodes.end()), spec);
                addLog("Đã đặt tải lên " + std::to_string(n) + " cạnh biên (quy về lực nút).");
                if (n == 0) undoStack.pop_back();
                else changed();
            }
        }
    }

    // ------------------------------------------------ phần tử
    std::vector<int> bars, planes;
    for (int e : selElems) (m.elements[e].isBar() ? bars : planes).push_back(e);
    if (!bars.empty()) {
        ImGui::SeparatorText(("Thanh đã chọn: " + std::to_string(bars.size())).c_str());
        Element& f = m.elements[bars[0]];
        if (bars.size() == 1) {
            ImGui::Text("Phần tử %d: nút %d → %d, L = %.3f m", bars[0] + 1, f.n[0] + 1, f.n[1] + 1, m.barLength(bars[0]));
        }
        int type = f.type == ElemType::Truss ? 1 : 0;
        if (ImGui::RadioButton("Thanh khung", type == 0) && type != 0) {
            snapshot();
            for (int e : bars) m.elements[e].type = ElemType::Frame;
            changed();
        }
        ImGui::SameLine();
        if (ImGui::RadioButton("Thanh dàn", type == 1) && type != 1) {
            snapshot();
            for (int e : bars) m.elements[e].type = ElemType::Truss;
            changed();
        }
        int mat = f.mat, sec = f.sec;
        ImGui::SetNextItemWidth(-90);
        if (editCombo(*this, "Vật liệu##b", mat, mats))
            for (int e : bars) m.elements[e].mat = mat;
        ImGui::SetNextItemWidth(-90);
        if (editCombo(*this, "Tiết diện##b", sec, secs))
            for (int e : bars) m.elements[e].sec = sec;
        if (f.type == ElemType::Frame) {
            bool hi = f.hingeI, hj = f.hingeJ;
            if (editB(*this, "Khớp đầu i", hi))
                for (int e : bars) m.elements[e].hingeI = hi;
            ImGui::SameLine();
            if (editB(*this, "Khớp đầu j", hj))
                for (int e : bars) m.elements[e].hingeJ = hj;
            double kf = f.kf;
            ImGui::SetNextItemWidth(110);
            if (ImGui::InputDouble("Nền Winkler k·b (kN/m²)", &kf, 0, 0, "%.4g")) {
                for (int e : bars) m.elements[e].kf = std::max(0.0, kf);
                needAnalyze = dirty = true;
            }
            trackEdit();
            help("Hệ số nền (kN/m³) nhân bề rộng đáy dầm (m). Nên chia dầm thành nhiều đoạn (≥ 10) để kết quả chính xác.");
            double kG = f.kG;
            ImGui::SetNextItemWidth(110);
            if (ImGui::InputDouble("Nền Pasternak kG·b (kN)", &kG, 0, 0, "%.4g")) {
                for (int e : bars) m.elements[e].kG = std::max(0.0, kG);
                needAnalyze = dirty = true;
            }
            trackEdit();
            help("Hệ số thứ hai của nền (mô hình hai hệ số Pasternak / Власов): lớp trượt nối các lò xo Winkler. Đặt 0 để dùng nền Winkler.");
            bool comp = f.kfCompOnly;
            if (editB(*this, "Nền chỉ chịu nén", comp))
                for (int e : bars) m.elements[e].kfCompOnly = comp;
            help("Đoạn thanh bị nâng lên khỏi nền thì nền không làm việc (phân tích lặp). Nền nằm phía −y cục bộ, tức phía dưới khi thanh vẽ từ trái sang phải.");
        }

        ImGui::SeparatorText(("Tải trên thanh — " + lc.name).c_str());
        static int kind = 0, dir = 3;
        static double q1 = -10, q2 = -10, P = -20, a = 1, ratio = 0.5, M0 = 10, t0 = 20, dt = 10, dL = 0.001;
        static bool proj = false, useRatio = true;
        ImGui::SetNextItemWidth(110);
        ImGui::Combo("Loại##ld", &kind, kKinds, 5);
        if (kind <= 1) {
            ImGui::SameLine();
            ImGui::SetNextItemWidth(110);
            ImGui::Combo("Hướng", &dir, kDirs, 4);
        }
        switch (kind) {
        case 0:
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("q đầu i", &q1, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("q đầu j (kN/m)", &q2, 0, 0, "%.4g");
            ImGui::Checkbox("Theo hình chiếu", &proj);
            help("Tải theo phương tổng thể tính trên chiều dài hình chiếu (ví dụ tải mái, tuyết).");
            break;
        case 1:
        case 2:
            ImGui::SetNextItemWidth(90);
            if (kind == 1) ImGui::InputDouble("P (kN)", &P, 0, 0, "%.4g");
            else ImGui::InputDouble("M (kNm)", &M0, 0, 0, "%.4g");
            ImGui::Checkbox("Vị trí theo tỉ lệ", &useRatio);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(80);
            if (useRatio) ImGui::InputDouble("a/L", &ratio, 0, 0, "%.3g");
            else ImGui::InputDouble("a (m)", &a, 0, 0, "%.4g");
            break;
        case 3:
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("t (°C)", &t0, 0, 0, "%.3g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("Δt = t dưới − t trên", &dt, 0, 0, "%.3g");
            break;
        case 4:
            ImGui::SetNextItemWidth(100);
            ImGui::InputDouble("ΔL (m)", &dL, 0, 0, "%.4g");
            help("Sai số chế tạo: thanh dài hơn (ΔL > 0) hoặc ngắn hơn (ΔL < 0) so với thiết kế. Hệ siêu tĩnh sẽ phát sinh nội lực.");
            break;
        }
        if (ImGui::Button("Thêm tải cho thanh chọn")) {
            snapshot();
            for (int e : bars) {
                ElemLoad l;
                l.elem = e;
                l.kind = (ElemLoadKind)kind;
                l.dir = (LoadDir)dir;
                const double L = m.barLength(e);
                switch (kind) {
                case 0:
                    l.q1 = q1;
                    l.q2 = q2;
                    l.projected = proj;
                    break;
                case 1:
                case 2:
                    l.q1 = kind == 1 ? P : M0;
                    l.a = std::clamp(useRatio ? ratio * L : a, 0.0, L);
                    break;
                case 3:
                    l.q1 = t0;
                    l.q2 = dt;
                    break;
                case 4: l.q1 = dL; break;
                }
                lc.elem.push_back(l);
            }
            changed();
        }
        ImGui::SameLine();
        if (ImGui::Button("Xóa tải trên thanh")) {
            snapshot();
            std::set<int> s(bars.begin(), bars.end());
            lc.elem.erase(std::remove_if(lc.elem.begin(), lc.elem.end(), [&](const ElemLoad& l) { return s.count(l.elem) > 0; }),
                          lc.elem.end());
            changed();
        }
        ImGui::SeparatorText("Hình học");
        static int nDiv = 4;
        ImGui::SetNextItemWidth(80);
        ImGui::InputInt("##ndiv", &nDiv);
        nDiv = std::clamp(nDiv, 2, 200);
        ImGui::SameLine();
        if (ImGui::Button("Chia thanh")) {
            snapshot();
            std::vector<int> sorted = bars;
            std::sort(sorted.rbegin(), sorted.rend());
            for (int e : sorted) subdivideBar(m, e, nDiv);
            selElems.clear();
            changed();
        }
        ImGui::SameLine();
        if (ImGui::Button("Đảo chiều")) {
            snapshot();
            for (int e : bars) {
                Element& el = m.elements[e];
                const double L = m.barLength(e);
                std::swap(el.n[0], el.n[1]);
                std::swap(el.hingeI, el.hingeJ);
                // Trục cục bộ x, y đều đổi chiều: tải theo hướng cục bộ đổi dấu, vị trí a → L − a,
                // thớ trên/dưới hoán đổi nên Δt đổi dấu.
                for (auto& c : m.cases)
                    for (auto& l : c.elem) {
                        if (l.elem != e) continue;
                        const bool local = l.dir == LoadDir::LocalX || l.dir == LoadDir::LocalY;
                        switch (l.kind) {
                        case ElemLoadKind::Distributed:
                            std::swap(l.q1, l.q2);
                            if (local) {
                                l.q1 = -l.q1;
                                l.q2 = -l.q2;
                            }
                            break;
                        case ElemLoadKind::Point:
                            l.a = L - l.a;
                            if (local) l.q1 = -l.q1;
                            break;
                        case ElemLoadKind::Moment: l.a = L - l.a; break;
                        case ElemLoadKind::Temperature: l.q2 = -l.q2; break;
                        case ElemLoadKind::LackOfFit: break; // ΔL không phụ thuộc chiều thanh
                        }
                    }
            }
            changed();
        }
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Đổi nút i ↔ j (đổi chiều trục cục bộ, tải được giữ nguyên tác dụng)");
    }
    if (!planes.empty()) {
        ImGui::SeparatorText(("Phần tử phẳng đã chọn: " + std::to_string(planes.size())).c_str());
        Element& f = m.elements[planes[0]];
        int mat = f.mat;
        ImGui::SetNextItemWidth(-90);
        if (editCombo(*this, "Vật liệu##p", mat, mats))
            for (int e : planes) m.elements[e].mat = mat;
        double t = f.t;
        ImGui::SetNextItemWidth(100);
        if (ImGui::InputDouble("Chiều dày t (m)", &t, 0, 0, "%.4g")) {
            for (int e : planes) m.elements[e].t = t;
            needAnalyze = dirty = true;
        }
        trackEdit();
        if (f.type == ElemType::Quad4) {
            bool inc = f.incompatible;
            if (editB(*this, "Mode không tương thích (QM6)", inc))
                for (int e : planes)
                    if (m.elements[e].type == ElemType::Quad4) m.elements[e].incompatible = inc;
            help("Phần tử Q4 bổ sung mode uốn không tương thích (Wilson–Taylor): chính xác hơn nhiều khi chịu uốn.");
        }
    }
    ImGui::Separator();
    if (ImGui::Button("Xóa đối tượng chọn")) deleteSelection();
    ImGui::SameLine();
    if (ImGui::Button("Bỏ chọn")) {
        selNodes.clear();
        selElems.clear();
    }
    ImGui::End();
}

// ============================================================================ panel kết quả
void App::resultPanel() {
    if (!ImGui::Begin("Kết quả", &showResultPanel)) {
        ImGui::End();
        return;
    }
    const Model& m = model;
    {
        std::vector<std::string> names;
        for (auto& c : m.cases) names.push_back("TH: " + c.name);
        for (auto& c : m.combos) names.push_back("THợp: " + c.name);
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::BeginCombo("##rs", names.empty() ? "" : names[std::min(resultSet, (int)names.size() - 1)].c_str())) {
            for (int i = 0; i < (int)names.size(); ++i)
                if (ImGui::Selectable(names[i].c_str(), i == resultSet)) {
                    resultSet = i;
                    if (i < (int)m.cases.size()) activeCase = i;
                }
            ImGui::EndCombo();
        }
    }
    const char* views[] = {"Sơ đồ", "Biến dạng", "N", "Q", "M", "σ thanh", "Phổ màu", "Đường ảnh hưởng", "Dạng riêng"};
    int v = (int)view;
    const float right = ImGui::GetCursorScreenPos().x + ImGui::GetContentRegionAvail().x;
    for (int i = 0; i < 9; ++i) {
        const float w = ImGui::GetFrameHeight() + ImGui::GetStyle().ItemInnerSpacing.x + ImGui::CalcTextSize(views[i]).x;
        if (i && ImGui::GetItemRectMax().x + ImGui::GetStyle().ItemSpacing.x + w < right) ImGui::SameLine();
        if (ImGui::RadioButton(views[i], v == i)) {
            view = (ResultView)i;
            if (view == ResultView::Mode && !(eig.ok && eigVersion == modelVersion)) showEigenDialog = true;
        }
    }
    if (view == ResultView::Contour) {
        const ContourField fields[] = {ContourField::Ux, ContourField::Uy, ContourField::Umag, ContourField::Sx,
                                       ContourField::Sy, ContourField::Txy, ContourField::S1, ContourField::S2,
                                       ContourField::VonMises};
        ImGui::SetNextItemWidth(170);
        if (ImGui::BeginCombo("Đại lượng", contourFieldName(contour))) {
            for (auto f : fields)
                if (ImGui::Selectable(contourFieldName(f), f == contour)) contour = f;
            ImGui::EndCombo();
        }
        ImGui::Checkbox("Trung bình tại nút", &contourAveraged);
        ImGui::SameLine();
        ImGui::Checkbox("Trên hình biến dạng", &contourOnDeformed);
    }
    if (view == ResultView::N || view == ResultView::Q || view == ResultView::M) {
        ImGui::Checkbox("Biểu đồ bao", &showEnvelope);
        help("Bao nội lực (max, min) qua tất cả tổ hợp; nếu chưa có tổ hợp thì qua tất cả trường hợp tải.");
    }
    if (view == ResultView::Influence) {
        ImGui::SeparatorText("Đường ảnh hưởng (P = 1 thẳng đứng, hướng xuống)");
        const std::vector<int> path = effectiveIlPath();
        ImGui::Text("Đường chạy tải: %zu thanh%s", path.size(), ilPath.empty() ? " (tự động)" : "");
        help("Tự động: các thanh nằm ngang ở cao độ thấp nhất. Có thể chọn thanh trên sơ đồ rồi nhấn \"Dùng thanh đang chọn\".");
        if (ImGui::Button("Dùng thanh đang chọn")) {
            ilPath.clear();
            for (int e : selElems)
                if (model.elements[e].isBar()) ilPath.push_back(e);
            ilDirty = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Tự động")) {
            ilPath.clear();
            ilDirty = true;
        }
        IlTarget& t = ilTarget;
        const char* kinds[] = {"Phản lực gối", "Nội lực tại tiết diện", "Chuyển vị nút"};
        int k = (int)t.kind;
        ImGui::SetNextItemWidth(200);
        if (ImGui::Combo("Đại lượng##il", &k, kinds, 3)) {
            t.kind = (IlKind)k;
            ilDirty = true;
        }
        if (t.kind == IlKind::Force) {
            int e1 = t.elem + 1;
            ImGui::SetNextItemWidth(100);
            if (ImGui::InputInt("Thanh", &e1)) {
                t.elem = std::clamp(e1 - 1, 0, std::max(0, (int)m.elements.size() - 1));
                ilDirty = true;
            }
            ImGui::SetNextItemWidth(100);
            if (ImGui::InputDouble("x từ đầu i (m)", &t.x, 0.1, 1.0, "%.3f")) ilDirty = true;
            const char* qs[] = {"N", "Q", "M"};
            int q = (int)t.qty;
            for (int i = 0; i < 3; ++i) {
                if (i) ImGui::SameLine();
                if (ImGui::RadioButton(qs[i], q == i)) {
                    t.qty = (BarQty)i;
                    ilDirty = true;
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Lấy thanh đang chọn") && !selElems.empty()) {
                t.elem = *selElems.begin();
                if (t.elem < (int)m.elements.size() && m.elements[t.elem].isBar()) t.x = m.barLength(t.elem) / 2;
                ilDirty = true;
            }
        } else {
            int n1 = t.node + 1;
            ImGui::SetNextItemWidth(100);
            if (ImGui::InputInt("Nút", &n1)) {
                t.node = std::clamp(n1 - 1, 0, std::max(0, (int)m.nodes.size() - 1));
                ilDirty = true;
            }
            const char* ds[] = {"x", "y", "φ"};
            for (int i = 0; i < 3; ++i) {
                ImGui::SameLine();
                if (ImGui::RadioButton(ds[i], t.dof == i)) {
                    t.dof = i;
                    ilDirty = true;
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Lấy nút đang chọn") && !selNodes.empty()) {
                t.node = *selNodes.begin();
                ilDirty = true;
            }
        }
        updateInfluence();
        if (il.ok) {
            ImGui::TextUnformatted(ilTargetName(m, t).c_str());
            ImGui::Text("Tung độ: max = %.4f, min = %.4f", il.vmax, il.vmin);
            ImGui::Text("Diện tích: Ω+ = %.4f, Ω− = %.4f, Ω = %.4f", il.areaPos, il.areaNeg, il.areaPos + il.areaNeg);
            ImGui::TextDisabled("Tải đều q trên toàn đường chạy: S = q·Ω; lực P đặt tại tung độ y: S = P·y");
            ImGui::TextDisabled("(Tải q, P hướng xuống mang dấu dương)");
        } else {
            ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(200, 30, 30, 255));
            ImGui::TextWrapped("%s", il.error.c_str());
            ImGui::PopStyleColor();
        }
    }
    ImGui::SetNextItemWidth(160);
    ImGui::SliderFloat("Phóng đại biến dạng", &deformFactor, 0.05f, 20.0f, "×%.2f", ImGuiSliderFlags_Logarithmic);
    ImGui::SetNextItemWidth(160);
    ImGui::SliderFloat("Tỉ lệ biểu đồ", &diagramFactor, 0.2f, 5.0f, "×%.2f", ImGuiSliderFlags_Logarithmic);
    ImGui::Checkbox("Giá trị", &showValues);
    ImGui::SameLine();
    ImGui::Checkbox("Phản lực", &showReactions);
    ImGui::SameLine();
    ImGui::Checkbox("Sơ đồ gốc", &showUndeformed);
    ImGui::SameLine();
    ImGui::Checkbox("Kèm biến dạng", &showDeformedWithDiagram);

    if (!result.ok) {
        ImGui::Separator();
        ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(200, 30, 30, 255));
        ImGui::TextWrapped("%s", result.error.empty() ? "Chưa tính." : result.error.c_str());
        ImGui::PopStyleColor();
        if (ImGui::Button("Tính (F5)")) analyzeNow();
        ImGui::End();
        return;
    }
    const LoadSetResult* rs = currentSet();
    if (!rs) {
        ImGui::End();
        return;
    }
    for (const auto& w : result.warnings) {
        ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(190, 110, 0, 255));
        ImGui::TextWrapped("⚠ %s", w.c_str());
        ImGui::PopStyleColor();
    }

    // ---- tóm tắt
    ImGui::SeparatorText("Tóm tắt");
    double umax = 0, uxm = 0, uym = 0, Mmax = 0, Nmax = 0, Qmax = 0;
    int iu = -1;
    for (size_t i = 0; i < m.nodes.size(); ++i) {
        const double u = std::hypot(rs->U[3 * i], rs->U[3 * i + 1]);
        if (u > umax) {
            umax = u;
            iu = (int)i;
        }
        if (std::abs(rs->U[3 * i]) > std::abs(uxm)) uxm = rs->U[3 * i];
        if (std::abs(rs->U[3 * i + 1]) > std::abs(uym)) uym = rs->U[3 * i + 1];
    }
    for (size_t e = 0; e < m.elements.size(); ++e)
        for (size_t k = 0; k < rs->bars[e].x.size(); ++k) {
            const auto& b = rs->bars[e];
            if (std::abs(b.M[k]) > std::abs(Mmax)) Mmax = b.M[k];
            if (std::abs(b.N[k]) > std::abs(Nmax)) Nmax = b.N[k];
            if (std::abs(b.Q[k]) > std::abs(Qmax)) Qmax = b.Q[k];
        }
    ImGui::Text("Chuyển vị lớn nhất: %.4f mm (nút %d)", umax * 1e3, iu + 1);
    ImGui::Text("ux max = %.4f mm, uy max = %.4f mm", uxm * 1e3, uym * 1e3);
    if (m.hasBars()) ImGui::Text("|M|max = %.3f kNm   |Q|max = %.3f kN   |N|max = %.3f kN", Mmax, Qmax, Nmax);
    if (m.hasPlane()) {
        double vm = 0, s1 = -DBL_MAX, s2 = DBL_MAX;
        for (size_t i = 0; i < m.nodes.size(); ++i) {
            if (rs->nodeStressCount[i] == 0) continue;
            const auto& s = rs->nodeStress[i];
            const StressInfo si = stressInfo(s[0], s[1], s[2], s[3]);
            vm = std::max(vm, si.vonMises);
            s1 = std::max(s1, si.s1);
            s2 = std::min(s2, si.s2);
        }
        ImGui::Text("σ1 max = %.3f MPa, σ2 min = %.3f MPa, σvM max = %.3f MPa", s1 * 1e-3, s2 * 1e-3, vm * 1e-3);
    }
    ImGui::Text("ΣF = (%.3f; %.3f) kN   ΣR = (%.3f; %.3f) kN", rs->sumFx, rs->sumFy, rs->sumRx, rs->sumRy);
    if (rs->foundRy != 0) ImGui::Text("Phản lực nền Winkler: %.3f kN", rs->foundRy);
    ImGui::TextDisabled("%d phương trình, %.1f ms", result.numEquations, result.timeMs);
    if (ImGui::Button("Xuất báo cáo .txt")) requestReport();

    // ---- bảng
    if (ImGui::BeginTabBar("restabs")) {
        if (ImGui::BeginTabItem("Chuyển vị")) {
            if (ImGui::BeginTable("U", 4, tableFlags())) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("Nút", ImGuiTableColumnFlags_WidthFixed, 40);
                ImGui::TableSetupColumn("ux (mm)");
                ImGui::TableSetupColumn("uy (mm)");
                ImGui::TableSetupColumn("φz (rad)");
                ImGui::TableHeadersRow();
                ImGuiListClipper clip;
                clip.Begin((int)m.nodes.size());
                while (clip.Step())
                    for (int i = clip.DisplayStart; i < clip.DisplayEnd; ++i) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        if (ImGui::Selectable(std::to_string(i + 1).c_str(), selNodes.count(i) > 0,
                                              ImGuiSelectableFlags_SpanAllColumns)) {
                            selNodes = {i};
                            selElems.clear();
                        }
                        ImGui::TableNextColumn();
                        ImGui::Text("%.4f", rs->U[3 * i] * 1e3);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.4f", rs->U[3 * i + 1] * 1e3);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.3e", rs->U[3 * i + 2]);
                    }
                ImGui::EndTable();
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Phản lực")) {
            if (ImGui::BeginTable("R", 4, tableFlags())) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("Nút", ImGuiTableColumnFlags_WidthFixed, 40);
                ImGui::TableSetupColumn("Rx (kN)");
                ImGui::TableSetupColumn("Ry (kN)");
                ImGui::TableSetupColumn("M (kNm)");
                ImGui::TableHeadersRow();
                for (size_t i = 0; i < m.nodes.size(); ++i) {
                    bool any = false;
                    for (int c = 0; c < 3; ++c) any |= result.dofState[3 * i + c] == DofSupport || rs->R[3 * i + c] != 0;
                    if (!any) continue;
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::Text("%zu", i + 1);
                    for (int c = 0; c < 3; ++c) {
                        ImGui::TableNextColumn();
                        if (result.dofState[3 * i + c] == DofSupport || rs->R[3 * i + c] != 0)
                            ImGui::Text("%.3f", rs->R[3 * i + c]);
                        else ImGui::TextDisabled("—");
                    }
                }
                ImGui::EndTable();
            }
            ImGui::EndTabItem();
        }
        if (m.hasBars() && ImGui::BeginTabItem("Nội lực thanh")) {
            ImGui::TextDisabled("N > 0 kéo; Q > 0 quay thuận kim đồng hồ; M > 0 căng thớ dưới");
            if (ImGui::BeginTable("F", 8, tableFlags())) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("PT", ImGuiTableColumnFlags_WidthFixed, 34);
                ImGui::TableSetupColumn("Ni");
                ImGui::TableSetupColumn("Qi");
                ImGui::TableSetupColumn("Mi");
                ImGui::TableSetupColumn("Nj");
                ImGui::TableSetupColumn("Qj");
                ImGui::TableSetupColumn("Mj");
                ImGui::TableSetupColumn("M nhịp");
                ImGui::TableHeadersRow();
                for (size_t e = 0; e < m.elements.size(); ++e) {
                    if (!m.elements[e].isBar()) continue;
                    const auto& b = rs->bars[e];
                    if (b.x.empty()) continue;
                    const size_t n = b.x.size() - 1;
                    size_t kx = 0;
                    for (size_t k = 0; k <= n; ++k)
                        if (std::abs(b.M[k]) > std::abs(b.M[kx])) kx = k;
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    if (ImGui::Selectable(std::to_string(e + 1).c_str(), selElems.count((int)e) > 0,
                                          ImGuiSelectableFlags_SpanAllColumns)) {
                        selElems = {(int)e};
                        selNodes.clear();
                    }
                    const double vals[7] = {b.N[0], b.Q[0], b.M[0], b.N[n], b.Q[n], b.M[n], b.M[kx]};
                    for (double x : vals) {
                        ImGui::TableNextColumn();
                        ImGui::Text("%.2f", x);
                    }
                }
                ImGui::EndTable();
            }
            ImGui::EndTabItem();
        }
        if (m.hasPlane() && ImGui::BeginTabItem("Ứng suất")) {
            if (ImGui::BeginTable("S", 6, tableFlags())) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("Nút", ImGuiTableColumnFlags_WidthFixed, 40);
                ImGui::TableSetupColumn("σx");
                ImGui::TableSetupColumn("σy");
                ImGui::TableSetupColumn("τxy");
                ImGui::TableSetupColumn("σ1 / σ2");
                ImGui::TableSetupColumn("σvM (MPa)");
                ImGui::TableHeadersRow();
                ImGuiListClipper clip;
                clip.Begin((int)m.nodes.size());
                while (clip.Step())
                    for (int i = clip.DisplayStart; i < clip.DisplayEnd; ++i) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        ImGui::Text("%d", i + 1);
                        if (rs->nodeStressCount[i] == 0) continue;
                        const auto& s = rs->nodeStress[i];
                        const StressInfo si = stressInfo(s[0], s[1], s[2], s[3]);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.3f", si.sx * 1e-3);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.3f", si.sy * 1e-3);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.3f", si.txy * 1e-3);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.3f / %.3f", si.s1 * 1e-3, si.s2 * 1e-3);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.3f", si.vonMises * 1e-3);
                    }
                ImGui::EndTable();
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
}

void App::logPanel() {
    if (!ImGui::Begin("Nhật ký", &showLogPanel)) {
        ImGui::End();
        return;
    }
    if (ImGui::SmallButton("Xóa")) log.clear();
    ImGui::SameLine();
    ImGui::Checkbox("Tự động tính khi thay đổi", &autoAnalyze);
    ImGui::BeginChild("logscroll");
    for (const auto& l : log) {
        const bool err = l.rfind("Lỗi", 0) == 0;
        if (err) ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(200, 30, 30, 255));
        ImGui::TextUnformatted(l.c_str());
        if (err) ImGui::PopStyleColor();
    }
    if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY()) ImGui::SetScrollHereY(1.0f);
    ImGui::EndChild();
    ImGui::End();
}

// ============================================================================ hộp thoại tạo kết cấu
void App::meshDialog() {
    ImGui::SetNextWindowSize(ImVec2(380, 0), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Tạo lưới phần tử phẳng", &showMeshDialog)) {
        ImGui::End();
        return;
    }
    static double x0 = 0, y0 = 0, W = 6, H = 3, t = 0.2;
    static int nx = 12, ny = 6, type = 0;
    static bool inc = true, fixBottom = false, pinCorners = true;
    static double q = -100;
    ImGui::TextWrapped("Tạo lưới tứ giác ánh xạ trên miền chữ nhật. Nút trùng với nút có sẵn sẽ được nối.");
    ImGui::InputDouble("x0 (m)", &x0);
    ImGui::InputDouble("y0 (m)", &y0);
    ImGui::InputDouble("Rộng (m)", &W);
    ImGui::InputDouble("Cao (m)", &H);
    ImGui::InputInt("Số ô theo X", &nx);
    ImGui::InputInt("Số ô theo Y", &ny);
    nx = std::clamp(nx, 1, 400);
    ny = std::clamp(ny, 1, 400);
    const char* types[] = {"Tứ giác Q4", "Tam giác CST"};
    ImGui::Combo("Loại phần tử", &type, types, 2);
    if (type == 0) ImGui::Checkbox("Mode không tương thích (QM6)", &inc);
    ImGui::InputDouble("Chiều dày t (m)", &t);
    const auto mats = matNames(model);
    ImGui::SetNextItemWidth(200);
    if (ImGui::BeginCombo("Vật liệu", curMat < (int)mats.size() ? mats[curMat].c_str() : "—")) {
        for (int i = 0; i < (int)mats.size(); ++i)
            if (ImGui::Selectable(mats[i].c_str(), i == curMat)) curMat = i;
        ImGui::EndCombo();
    }
    ImGui::SeparatorText("Tùy chọn nhanh");
    ImGui::Checkbox("Ngàm toàn bộ cạnh dưới", &fixBottom);
    ImGui::Checkbox("Gối cố định + di động tại hai góc dưới", &pinCorners);
    ImGui::InputDouble("Tải phân bố mép trên (kN/m)", &q);
    if (ImGui::Button("Tạo lưới", ImVec2(-FLT_MIN, 0)) && W > 0 && H > 0 && t > 0) {
        snapshot();
        if (model.materials.empty()) model.materials.push_back(toMaterial(materialPresets()[5]));
        std::vector<int> g = meshRectangle(model, x0, y0, W, H, nx, ny, type == 0 ? ElemType::Quad4 : ElemType::Tri3,
                                           curMat, t, inc);
        if (fixBottom)
            for (int i = 0; i <= nx; ++i) model.supportFor(g[i]).fix = {true, true, false};
        else if (pinCorners) {
            model.supportFor(g[0]).fix = {true, true, false};
            model.supportFor(g[nx]).fix = {false, true, false};
        }
        if (q != 0) {
            LoadCase& lc = model.cases[activeCase];
            for (int i = 0; i <= nx; ++i) {
                const double w = (i == 0 || i == nx) ? 0.5 : 1.0;
                lc.nodal.push_back({g[ny * (nx + 1) + i], {0, q * w * W / nx, 0}});
            }
        }
        view = ResultView::Contour;
        fitRequested = true;
        changed();
    }
    ImGui::End();
}

void App::frameGenDialog() {
    ImGui::SetNextWindowSize(ImVec2(420, 0), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Tạo khung / dầm liên tục", &showFrameGen)) {
        ImGui::End();
        return;
    }
    static std::string spans = "6 6", heights = "3.6 3.6";
    static int base = 0, colSec = 0, beamSec = 0;
    static double qBeam = -20, wind = 10;
    static bool replace = true;
    ImGui::TextWrapped("Nhập các nhịp và chiều cao tầng cách nhau bởi dấu cách (cho phép viết 3*6). Để trống chiều cao để tạo dầm liên tục.");
    ImGui::InputText("Nhịp (m)", &spans);
    ImGui::InputText("Chiều cao tầng (m)", &heights);
    const char* bases[] = {"Ngàm", "Khớp (gối cố định)"};
    ImGui::Combo("Liên kết chân cột", &base, bases, 2);
    const auto secs = secNames(model);
    ImGui::Combo("Tiết diện cột", &colSec, [](void* d, int i) { return ((std::vector<std::string>*)d)->at(i).c_str(); },
                 (void*)&secs, (int)secs.size());
    ImGui::Combo("Tiết diện dầm", &beamSec, [](void* d, int i) { return ((std::vector<std::string>*)d)->at(i).c_str(); },
                 (void*)&secs, (int)secs.size());
    ImGui::InputDouble("Tải phân bố trên dầm (kN/m)", &qBeam);
    ImGui::InputDouble("Gió tại mỗi sàn, nút trái (kN)", &wind);
    ImGui::Checkbox("Thay thế mô hình hiện tại", &replace);
    if (ImGui::Button("Tạo", ImVec2(-FLT_MIN, 0))) {
        const auto sp = parseList(spans), hs = parseList(heights);
        if (sp.empty()) addLog("Lỗi: chưa nhập nhịp.");
        else {
            snapshot();
            if (replace) {
                Model nm;
                nm.materials = model.materials;
                nm.sections = model.sections;
                nm.title = hs.empty() ? "Dầm liên tục" : "Khung phẳng";
                model = nm;
            }
            if (model.materials.empty()) model.materials.push_back(toMaterial(materialPresets()[4]));
            if (model.sections.empty()) model.sections.push_back(rectSection(0.3, 0.5));
            colSec = std::clamp(colSec, 0, (int)model.sections.size() - 1);
            beamSec = std::clamp(beamSec, 0, (int)model.sections.size() - 1);
            std::vector<double> xs = {0}, ys = {0};
            for (double s : sp) xs.push_back(xs.back() + s);
            for (double h : hs) ys.push_back(ys.back() + h);
            const int nX = (int)xs.size(), nY = (int)ys.size();
            std::vector<std::vector<int>> id(nY, std::vector<int>(nX));
            for (int j = 0; j < nY; ++j)
                for (int i = 0; i < nX; ++i) id[j][i] = model.addNode(xs[i], ys[j]);
            auto bar = [&](int a, int b, int sec) {
                Element e;
                e.type = ElemType::Frame;
                e.n = {a, b, -1, -1};
                e.mat = std::min(curMat, (int)model.materials.size() - 1);
                e.sec = sec;
                model.elements.push_back(e);
                return (int)model.elements.size() - 1;
            };
            if (model.cases.empty()) model.cases.push_back(LoadCase("Tĩnh tải"));
            activeCase = std::min(activeCase, (int)model.cases.size() - 1);
            LoadCase& lc = model.cases[activeCase];
            for (int j = 0; j + 1 < nY; ++j)
                for (int i = 0; i < nX; ++i) bar(id[j][i], id[j + 1][i], colSec);
            const int jb0 = hs.empty() ? 0 : 1;
            for (int j = jb0; j < nY; ++j)
                for (int i = 0; i + 1 < nX; ++i) {
                    const int e = bar(id[j][i], id[j][i + 1], beamSec);
                    if (qBeam != 0) {
                        ElemLoad l;
                        l.elem = e;
                        l.q1 = l.q2 = qBeam;
                        lc.elem.push_back(l);
                    }
                }
            if (hs.empty()) {
                model.supportFor(id[0][0]).fix = {true, true, false};
                for (int i = 1; i < nX; ++i) model.supportFor(id[0][i]).fix = {false, true, false};
            } else {
                for (int i = 0; i < nX; ++i) model.supportFor(id[0][i]).fix = {true, true, base == 0};
                if (wind != 0) {
                    LoadCase g("Gió");
                    for (int j = 1; j < nY; ++j) g.nodal.push_back({id[j][0], {wind, 0, 0}});
                    model.cases.push_back(g);
                }
            }
            fitRequested = true;
            view = ResultView::M;
            changed();
        }
    }
    ImGui::End();
}

void App::trussGenDialog() {
    ImGui::SetNextWindowSize(ImVec2(400, 0), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Tạo dàn", &showTrussGen)) {
        ImGui::End();
        return;
    }
    static int kind = 0, nP = 6, sec = 0;
    static double d = 3, h = 3, P = -40;
    static bool replace = true, topLoad = false;
    const char* kinds[] = {"Pratt (xiên kéo)", "Howe (xiên nén)", "Warren (tam giác)", "Dàn tam giác (mái dốc)"};
    ImGui::Combo("Kiểu dàn", &kind, kinds, 4);
    ImGui::InputInt("Số khoang (chẵn)", &nP);
    nP = std::clamp(nP + (nP % 2), 2, 200);
    ImGui::InputDouble("Chiều dài khoang d (m)", &d);
    ImGui::InputDouble("Chiều cao dàn h (m)", &h);
    const auto secs = secNames(model);
    ImGui::Combo("Tiết diện", &sec, [](void* dd, int i) { return ((std::vector<std::string>*)dd)->at(i).c_str(); },
                 (void*)&secs, (int)secs.size());
    ImGui::InputDouble("Lực tại mỗi nút (kN)", &P);
    ImGui::Checkbox("Đặt lực ở cánh trên", &topLoad);
    ImGui::Checkbox("Thay thế mô hình hiện tại", &replace);
    if (ImGui::Button("Tạo", ImVec2(-FLT_MIN, 0)) && d > 0 && h > 0) {
        snapshot();
        if (replace) {
            Model nm;
            nm.materials = model.materials;
            nm.sections = model.sections;
            nm.title = std::string("Dàn ") + kinds[kind];
            model = nm;
        }
        if (model.materials.empty()) model.materials.push_back(toMaterial(materialPresets()[1]));
        if (model.sections.empty()) {
            Section s;
            s.name = "2L75x6";
            s.A = 17.56e-4;
            s.I = 0;
            s.h = 0.075;
            model.sections.push_back(s);
        }
        const int si = std::clamp(sec, 0, (int)model.sections.size() - 1);
        auto tr = [&](int a, int b) {
            Element e;
            e.type = ElemType::Truss;
            e.n = {a, b, -1, -1};
            e.mat = std::min(curMat, (int)model.materials.size() - 1);
            e.sec = si;
            if (model.findElement(e.type, e.n) < 0) model.elements.push_back(e);
        };
        std::vector<int> bot, top;
        for (int i = 0; i <= nP; ++i) bot.push_back(model.addNode(i * d, 0));
        if (kind == 3) {
            // dàn tam giác: cánh trên dốc lên giữa nhịp
            for (int i = 0; i <= nP; ++i) {
                const double y = h * (1.0 - std::abs(i - nP / 2.0) / (nP / 2.0));
                top.push_back(i == 0 ? bot[0] : i == nP ? bot[nP] : model.addNode(i * d, y));
            }
            for (int i = 0; i < nP; ++i) {
                tr(bot[i], bot[i + 1]);
                tr(top[i], top[i + 1]);
            }
            for (int i = 1; i < nP; ++i) tr(bot[i], top[i]);
            for (int i = 1; i < nP - 1; ++i) {
                if (i < nP / 2) tr(bot[i], top[i + 1]);
                else tr(top[i], bot[i + 1]);
            }
        } else if (kind == 2) {
            for (int i = 0; i < nP; ++i) top.push_back(model.addNode((i + 0.5) * d, h));
            for (int i = 0; i < nP; ++i) {
                tr(bot[i], bot[i + 1]);
                tr(bot[i], top[i]);
                tr(top[i], bot[i + 1]);
                if (i + 1 < nP) tr(top[i], top[i + 1]);
            }
        } else {
            for (int i = 1; i < nP; ++i) top.push_back(model.addNode(i * d, h));
            for (int i = 0; i < nP; ++i) tr(bot[i], bot[i + 1]);
            for (int i = 0; i + 1 < (int)top.size(); ++i) tr(top[i], top[i + 1]);
            tr(bot[0], top[0]);
            tr(bot[nP], top.back());
            for (int i = 1; i < nP; ++i) tr(bot[i], top[i - 1]);
            for (int i = 1; i < nP - 1; ++i) {
                const bool left = i < nP / 2;
                const bool pratt = kind == 0;
                if (left == pratt) tr(top[i - 1], bot[i + 1]);
                else tr(bot[i], top[i]);
            }
        }
        model.supportFor(bot[0]).fix = {true, true, false};
        model.supportFor(bot[nP]).fix = {false, true, false};
        if (model.cases.empty()) model.cases.push_back(LoadCase("Tải nút"));
        activeCase = std::min(activeCase, (int)model.cases.size() - 1);
        LoadCase& lc = model.cases[activeCase];
        if (P != 0) {
            if (topLoad) {
                for (int n : top)
                    if (n != bot[0] && n != bot[nP]) lc.nodal.push_back({n, {0, P, 0}});
            } else {
                for (int i = 1; i < nP; ++i) lc.nodal.push_back({bot[i], {0, P, 0}});
            }
        }
        fitRequested = true;
        view = ResultView::N;
        changed();
    }
    ImGui::End();
}

void App::aboutDialog() {
    ImGui::SetNextWindowSize(ImVec2(560, 0), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Giới thiệu & quy ước", &showAbout)) {
        ImGui::End();
        return;
    }
    ImGui::TextWrapped("FEM2D — tính toán kết cấu phẳng bằng phương pháp phần tử hữu hạn (phương pháp chuyển vị).");
    ImGui::SeparatorText("Đơn vị");
    ImGui::BulletText("Lực kN, chiều dài m, mômen kNm, E kN/m², ứng suất hiển thị MPa, chuyển vị mm.");
    ImGui::SeparatorText("Hệ tọa độ & quy ước dấu");
    ImGui::BulletText("Trục X sang phải, Y hướng lên; mômen, góc xoay dương ngược chiều kim đồng hồ.");
    ImGui::BulletText("Trục cục bộ thanh: x từ nút i đến j, y quay x ngược chiều kim đồng hồ 90°.");
    ImGui::BulletText("N > 0: kéo. Q > 0: làm phần thanh quay thuận chiều kim đồng hồ.");
    ImGui::BulletText("M > 0: căng thớ dưới (phía −y cục bộ). Biểu đồ M vẽ về phía thớ căng.");
    ImGui::SeparatorText("Phần tử");
    ImGui::BulletText("Thanh dàn: chỉ chịu kéo nén. Tải trong nhịp quy về nút như dầm đơn giản.");
    ImGui::BulletText("Thanh khung Euler–Bernoulli, khớp đầu thanh, nền Winkler, tải nhiệt độ.");
    ImGui::BulletText("Tam giác biến dạng hằng (CST), tứ giác Q4 đẳng tham số (tùy chọn mode QM6).");
    ImGui::SeparatorText("Tài liệu tham khảo");
    ImGui::BulletText("Lều Thọ Trình — Cơ học kết cấu (tập 1, 2).");
    ImGui::BulletText("Дарков А.В., Шапошников Н.Н. — Строительная механика.");
    ImGui::BulletText("Смирнов А.Ф. и др. — Строительная механика. Стержневые системы.");
    ImGui::BulletText("Варвак П.М. — Метод конечных элементов; Zienkiewicz — The Finite Element Method.");
    ImGui::BulletText("TCVN 5574:2018 (bê tông), TCVN 5575:2012 (thép), ГОСТ 8239-89, ГОСТ 8240-97.");
    ImGui::End();
}

// ============================================================================ ổn định & dao động riêng
void App::eigenDialog() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x - 500, vp->WorkPos.y + 60), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(460, 520), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Ổn định & dao động riêng", &showEigenDialog)) {
        ImGui::End();
        return;
    }
    const Model& m = model;
    int kind = (int)eigOpt.kind;
    ImGui::RadioButton("Ổn định (tải trọng tới hạn)", &kind, 0);
    ImGui::SameLine();
    ImGui::RadioButton("Dao động riêng", &kind, 1);
    eigOpt.kind = (EigenKind)kind;
    std::vector<std::string> sets;
    for (auto& c : m.cases) sets.push_back("TH: " + c.name);
    for (auto& c : m.combos) sets.push_back("THợp: " + c.name);
    auto setCombo = [&](const char* label, int& v, bool allowNone) {
        const char* pv = v >= 0 && v < (int)sets.size() ? sets[v].c_str() : "— không —";
        ImGui::SetNextItemWidth(260);
        if (ImGui::BeginCombo(label, pv)) {
            if (allowNone && ImGui::Selectable("— không —", v < 0)) v = -1;
            for (int i = 0; i < (int)sets.size(); ++i)
                if (ImGui::Selectable(sets[i].c_str(), i == v)) v = i;
            ImGui::EndCombo();
        }
    };
    if (eigOpt.kind == EigenKind::Buckling) {
        setCombo("Bộ tải gây nén", eigOpt.loadSet, false);
        help("Tính lực dọc N từ bộ tải này, rồi tìm hệ số λcr sao cho λcr × tải làm kết cấu mất ổn định "
             "(phương pháp chuyển vị, ma trận độ cứng hình học).");
    } else {
        ImGui::SetNextItemWidth(100);
        ImGui::InputDouble("Hệ số khối lượng bản thân", &eigOpt.selfMassFactor, 0, 0, "%.3g");
        help("Khối lượng phân bố ρA = γA/g của thanh, ρt của phần tử phẳng (γ lấy từ vật liệu).");
        setCombo("Quy tải thành khối lượng", eigOpt.massSet, true);
        if (eigOpt.massSet >= 0) {
            ImGui::SetNextItemWidth(100);
            ImGui::InputDouble("Hệ số", &eigOpt.massSetFactor, 0, 0, "%.3g");
        }
        help("Khối lượng tập trung tại nút m = |Fy|/g từ tải trọng thẳng đứng của bộ tải (ví dụ tĩnh tải + 0,5 hoạt tải).");
    }
    ImGui::SetNextItemWidth(100);
    ImGui::InputInt("Số dạng riêng", &eigOpt.modes);
    eigOpt.modes = std::clamp(eigOpt.modes, 1, 50);
    ImGui::SetNextItemWidth(100);
    ImGui::InputInt("Số đoạn chia mỗi thanh", &eigOpt.subdiv);
    eigOpt.subdiv = std::clamp(eigOpt.subdiv, 1, 40);
    if (ImGui::Button("Tính", ImVec2(120, 0))) runEigen();
    ImGui::SameLine();
    ImGui::Checkbox("Hoạt họa dạng riêng", &animateMode);

    if (eigVersion >= 0 && eigVersion != modelVersion) {
        ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(200, 120, 0, 255));
        ImGui::TextWrapped("Mô hình đã thay đổi sau lần tính trước — hãy tính lại.");
        ImGui::PopStyleColor();
    }
    if (!eig.ok) {
        if (!eig.error.empty()) {
            ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(200, 30, 30, 255));
            ImGui::TextWrapped("%s", eig.error.c_str());
            ImGui::PopStyleColor();
        }
        ImGui::End();
        return;
    }
    if (eig.kind == EigenKind::Mechanism) {
        ImGui::TextWrapped("Đang hiển thị dạng biến hình từ phân tích cấu tạo hình học. Bấm Tính để chạy phân tích trị riêng.");
        ImGui::End();
        return;
    }
    for (const auto& w : eig.warnings) ImGui::TextDisabled("⚠ %s", w.c_str());
    ImGui::Text("%d phương trình, %.0f ms", eig.numEquations, eig.timeMs);
    const bool buck = eig.kind == EigenKind::Buckling;
    if (ImGui::BeginTable("modes", buck ? 2 : 4, tableFlags(), ImVec2(0, std::min(200.0f, 30 + 24.0f * eig.modes.size())))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Dạng", ImGuiTableColumnFlags_WidthFixed, 50);
        if (buck) ImGui::TableSetupColumn("λcr (hệ số tải tới hạn)");
        else {
            ImGui::TableSetupColumn("ω (rad/s)");
            ImGui::TableSetupColumn("f (Hz)");
            ImGui::TableSetupColumn("T (s)");
        }
        ImGui::TableHeadersRow();
        for (int i = 0; i < (int)eig.modes.size(); ++i) {
            const EigenMode& md = eig.modes[i];
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            if (ImGui::Selectable(std::to_string(i + 1).c_str(), i == eigMode, ImGuiSelectableFlags_SpanAllColumns)) {
                eigMode = i;
                view = ResultView::Mode;
            }
            ImGui::TableNextColumn();
            if (buck) ImGui::Text("%.5g", md.lambda);
            else {
                ImGui::Text("%.5g", md.omega);
                ImGui::TableNextColumn();
                ImGui::Text("%.5g", md.freq);
                ImGui::TableNextColumn();
                ImGui::Text("%.5g", md.period);
            }
        }
        ImGui::EndTable();
    }
    if (buck && eigVersion == modelVersion) {
        ImGui::SeparatorText("Thanh chịu nén ở dạng 1");
        ImGui::TextDisabled("Ncr = λcr·N; chiều dài tính toán l0 = μ·l, μ = (π/l)·√(EI/Ncr)");
        if (ImGui::BeginTable("mu", 5, tableFlags())) {
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableSetupColumn("Thanh", ImGuiTableColumnFlags_WidthFixed, 44);
            ImGui::TableSetupColumn("N (kN)");
            ImGui::TableSetupColumn("Ncr (kN)");
            ImGui::TableSetupColumn("μ");
            ImGui::TableSetupColumn("l0 (m)");
            ImGui::TableHeadersRow();
            for (size_t e = 0; e < m.elements.size() && e < eig.memberN.size(); ++e) {
                const double mu = effectiveLengthFactor(m, eig, (int)e);
                if (mu <= 0) continue;
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                if (ImGui::Selectable(std::to_string(e + 1).c_str(), selElems.count((int)e) > 0,
                                      ImGuiSelectableFlags_SpanAllColumns)) {
                    selElems = {(int)e};
                    selNodes.clear();
                }
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", eig.memberN[e]);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", -eig.memberN[e] * eig.modes[0].lambda);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", mu);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", mu * m.barLength((int)e));
            }
            ImGui::EndTable();
        }
    }
    ImGui::End();
}

// ============================================================================ cấu tạo hình học
void App::showMechanism() {
    if (kin.mechanismU.empty()) return;
    eig = EigenResult{};
    eig.ok = true;
    eig.kind = EigenKind::Mechanism;
    eig.sub = model;
    eig.barChain.assign(model.elements.size(), {});
    for (size_t e = 0; e < model.elements.size(); ++e)
        if (model.elements[e].isBar()) eig.barChain[e] = {model.elements[e].n[0], model.elements[e].n[1]};
    EigenMode md;
    md.U = kin.mechanismU;
    double big = 0;
    for (size_t i = 0; i < md.U.size(); ++i)
        if (i % 3 != 2) big = std::max(big, std::abs(md.U[i]));
    if (big > 0)
        for (double& v : md.U) v /= big;
    eig.modes.push_back(md);
    eigVersion = modelVersion;
    eigMode = 0;
    view = ResultView::Mode;
}

void App::kinematicsDialog() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x - 520, vp->WorkPos.y + 60), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(470, 0), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Phân tích cấu tạo hình học", &showKinematics)) {
        ImGui::End();
        return;
    }
    if (kinVersion != modelVersion) {
        kin = kinematicAnalysis(model);
        kinVersion = modelVersion;
    }
    const KinematicResult& k = kin;
    if (!k.applicable) {
        ImGui::TextWrapped("%s", k.note.c_str());
        ImGui::End();
        return;
    }
    if (!k.note.empty()) ImGui::TextDisabled("%s", k.note.c_str());
    ImGui::SeparatorText("Công thức Чебышёв");
    ImGui::Text("W = 3D − 3H − 2K − C₀ = 3·%d − 3·%d − 2·%d − %d = %d", k.disks, k.rigidJoints, k.hinges, k.supports, k.W);
    ImGui::TextDisabled("D: miếng cứng (thanh), H: mối hàn đơn, K: khớp đơn, C₀: liên kết tựa");
    ImGui::TextWrapped("%s", k.W > 0 ? "W > 0: hệ thiếu liên kết → biến hình."
                                     : "W ≤ 0: đủ số liên kết (điều kiện cần); cần xét thêm cách sắp xếp liên kết.");
    if (k.rank >= 0) {
        ImGui::SeparatorText("Phân tích ma trận cân bằng");
        ImGui::Text("Số phương trình cân bằng: %d", k.equations);
        ImGui::Text("Số ẩn (nội lực đầu thanh + phản lực): %d", k.unknowns);
        ImGui::Text("Hạng ma trận: %d", k.rank);
        ImGui::Text("Bậc siêu tĩnh n = %d − %d = %d", k.unknowns, k.rank, k.indeterminacy);
        ImGui::Text("Bậc tự do biến hình = %d − %d = %d", k.equations, k.rank, k.mechanisms);
    }
    ImGui::Separator();
    const bool bad = k.mechanisms > 0;
    ImGui::PushStyleColor(ImGuiCol_Text, bad ? IM_COL32(230, 70, 60, 255) : IM_COL32(60, 180, 90, 255));
    ImGui::TextWrapped("Kết luận: %s", k.verdict.c_str());
    ImGui::PopStyleColor();
    if (bad && !k.mechanismU.empty() && ImGui::Button("Xem dạng biến hình")) showMechanism();
    ImGui::End();
}

