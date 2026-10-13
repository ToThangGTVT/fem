// Lời giải theo giáo trình: phương pháp lực (từng bước), Maxwell–Mohr, phương pháp chuyển vị dạng ma trận; tạo vòm.
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>

#include "App.h"
#include "core/ForceMethod.h"
#include "core/Library.h"
#include "core/Mesher.h"

using namespace fem;

namespace {

std::string sci(double v) {
    char b[32];
    std::snprintf(b, sizeof b, "%.4e", v);
    return b;
}

ImGuiTableFlags tflags() {
    return ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollX | ImGuiTableFlags_ScrollY |
           ImGuiTableFlags_SizingFixedFit;
}

bool setCombo(const char* label, int& v, const Model& m) {
    std::vector<std::string> sets;
    for (auto& c : m.cases) sets.push_back("TH: " + c.name);
    for (auto& c : m.combos) sets.push_back("THợp: " + c.name);
    if (sets.empty()) return false;
    v = std::clamp(v, 0, (int)sets.size() - 1);
    bool ch = false;
    ImGui::SetNextItemWidth(260);
    if (ImGui::BeginCombo(label, sets[v].c_str())) {
        for (int i = 0; i < (int)sets.size(); ++i) {
            ImGui::PushID(i);
            if (ImGui::Selectable(sets[i].c_str(), i == v)) {
                v = i;
                ch = true;
            }
            ImGui::PopID();
        }
        ImGui::EndCombo();
    }
    return ch;
}

} // namespace

const LoadSetResult* App::forceMethodShown(std::string* title) {
    if (!showForceMethod || !fm.ok || fmVersion != modelVersion || fmShow < 0) return nullptr;
    const LoadSetResult* s = nullptr;
    if (fmShow < fm.n) s = &fm.unit[fmShow];
    else if (fmShow == fm.n) s = &fm.load;
    else if (fmShow == fm.n + 1) s = &fm.final;
    else if (fmShow == fm.n + 2 && mohr.ok) s = &mohr.unitState;
    if (s && title) *title = "Phương pháp lực — biểu đồ " + s->name;
    return s;
}

void App::forceMethodDialog() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x - 560, vp->WorkPos.y + 40), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(540, 620), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Phương pháp lực", &showForceMethod)) {
        ImGui::End();
        return;
    }
    setCombo("Bộ tải##fm", fmSet, model);
    ImGui::SameLine();
    if (ImGui::Button("Giải") || (fmVersion != modelVersion && fm.ok)) {
        fm = forceMethod(model, fmSet);
        fmVersion = modelVersion;
        mohr = MohrResult{};
        fmShow = fm.ok ? fm.n + 1 : -1;
        if (fm.ok) view = ResultView::M;
        else echo("Lỗi: " + fm.error);
    }
    if (!fm.ok) {
        if (!fm.error.empty()) ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.25f, 1), "%s", fm.error.c_str());
        ImGui::TextWrapped("Giải hệ thanh siêu tĩnh theo phương pháp lực: chọn hệ cơ bản, lập biểu đồ đơn vị, "
                           "tính hệ số δik, ΔiP bằng nhân biểu đồ, giải hệ phương trình chính tắc.");
        ImGui::End();
        return;
    }
    for (const auto& w : fm.warnings) ImGui::TextDisabled("⚠ %s", w.c_str());
    ImGui::SeparatorText("1. Bậc siêu tĩnh và hệ cơ bản");
    ImGui::Text("n = %d", fm.n);
    if (fm.n == 0) ImGui::TextWrapped("Hệ tĩnh định: nội lực xác định trực tiếp từ phương trình cân bằng.");
    for (const auto& s : fm.redundants) ImGui::BulletText("%s", s.c_str());

    if (fm.n > 0) {
        ImGui::SeparatorText("2. Hệ phương trình chính tắc  δ·X + Δ_P = 0");
        ImGui::TextDisabled("δik = Σ∫M̄iM̄k/EI dx + ΣN̄iN̄kL/EA;  ΔiP = Σ∫M̄iM⁰_P/EI dx + ...");
        if (fm.n <= 8) {
            for (int i = 0; i < fm.n; ++i) {
                std::string line;
                char b[64];
                for (int k = 0; k < fm.n; ++k) {
                    std::snprintf(b, sizeof b, "%s%.4e·X%d", k ? " + " : "", fm.delta[i][k], k + 1);
                    line += b;
                }
                std::snprintf(b, sizeof b, " + (%.4e) = 0", fm.deltaP[i]);
                ImGui::TextUnformatted((line + b).c_str());
            }
        } else if (ImGui::BeginTable("##delta", fm.n + 2, tflags(), ImVec2(0, 200))) {
            ImGui::TableSetupColumn("");
            for (int k = 0; k < fm.n; ++k) ImGui::TableSetupColumn(("X" + std::to_string(k + 1)).c_str());
            ImGui::TableSetupColumn("ΔiP");
            ImGui::TableHeadersRow();
            for (int i = 0; i < fm.n; ++i) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::Text("%d", i + 1);
                for (int k = 0; k < fm.n; ++k) {
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(sci(fm.delta[i][k]).c_str());
                }
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(sci(fm.deltaP[i]).c_str());
            }
            ImGui::EndTable();
        }
        ImGui::SeparatorText("3. Nghiệm");
        for (int i = 0; i < fm.n; ++i) ImGui::Text("X%d = %.4f", i + 1, fm.X[i]);
    }
    ImGui::SeparatorText("4. Biểu đồ (vẽ trên sơ đồ)");
    ImGui::TextDisabled("M = M⁰_P + Σ Xi·M̄i");
    auto pick = [&](const char* label, int id) {
        if (ImGui::RadioButton(label, fmShow == id)) {
            fmShow = id;
            if (view != ResultView::N && view != ResultView::Q && view != ResultView::M) view = ResultView::M;
        }
    };
    pick("Không (kết quả PTHH)", -1);
    for (int i = 0; i < fm.n && i < 30; ++i) {
        ImGui::SameLine();
        pick(("M̄" + std::to_string(i + 1)).c_str(), i);
    }
    pick("M⁰_P (hệ cơ bản)", fm.n);
    ImGui::SameLine();
    pick("M cuối cùng", fm.n + 1);
    ImGui::Text("Đối chiếu PTHH: sai lệch |M| lớn nhất = %.3e kNm", fm.maxDiffM);

    ImGui::SeparatorText("5. Chuyển vị theo Maxwell–Mohr");
    ImGui::SetNextItemWidth(90);
    int n1 = mohrNode + 1;
    if (ImGui::InputInt("Nút", &n1)) mohrNode = std::clamp(n1 - 1, 0, std::max(0, (int)model.nodes.size() - 1));
    const char* ds[] = {"ux", "uy", "φz"};
    for (int i = 0; i < 3; ++i) {
        ImGui::SameLine();
        if (ImGui::RadioButton(ds[i], mohrDof == i)) mohrDof = i;
    }
    ImGui::SameLine();
    if (!selNodes.empty() && ImGui::SmallButton("Nút đang chọn")) mohrNode = *selNodes.begin();
    if (ImGui::Button("Tính chuyển vị")) {
        mohr = mohrDisplacement(model, fm, mohrNode, mohrDof);
        if (mohr.ok) fmShow = fm.n + 2;
    }
    if (mohr.ok) {
        const bool rot = mohrDof == RZ;
        ImGui::Text("Δ = Σ∫M̄kM/EI dx + ΣN̄kNL/EA = %.6g %s", rot ? mohr.value : mohr.value * 1e3, rot ? "rad" : "mm");
        ImGui::Text("Theo PTHH: %.6g %s", rot ? mohr.fem : mohr.fem * 1e3, rot ? "rad" : "mm");
        if (ImGui::TreeNode("Đóng góp từng thanh")) {
            for (const auto& t : mohr.terms) ImGui::TextUnformatted(t.c_str());
            ImGui::TreePop();
        }
        pick("Biểu đồ M̄k do lực đơn vị", fm.n + 2);
    } else if (!mohr.error.empty()) {
        ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.25f, 1), "%s", mohr.error.c_str());
    }
    ImGui::Separator();
    if (ImGui::Button("Phương pháp chuyển vị dạng ma trận...")) showMatrixMethod = true;
    ImGui::End();
}

void App::matrixMethodDialog() {
    ImGui::SetNextWindowSize(ImVec2(760, 600), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Phương pháp chuyển vị (dạng ma trận)", &showMatrixMethod)) {
        ImGui::End();
        return;
    }
    static int set = 0;
    setCombo("Bộ tải##mm", set, model);
    ImGui::SameLine();
    if (ImGui::Button("Lập") || (mmVersion != modelVersion && mm.ok)) {
        mm = matrixMethod(model, set);
        mmVersion = modelVersion;
    }
    if (!mm.ok) {
        if (!mm.error.empty()) ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.25f, 1), "%s", mm.error.c_str());
        ImGui::TextWrapped("Lập ma trận độ cứng phần tử, chuyển sang hệ tổng thể, ghép thành K, giải K·Z = P.");
        ImGui::End();
        return;
    }
    ImGui::TextWrapped("Ẩn số Z (chuyển vị nút tự do): %zu. Đơn vị: kN, m, rad.", mm.Z.size());
    if (ImGui::CollapsingHeader("Ma trận độ cứng phần tử")) {
        for (const auto& el : mm.elems) {
            char h[128];
            std::snprintf(h, sizeof h, "Thanh %d: L = %.3f m, cos = %.4f, sin = %.4f", el.e + 1, el.L, el.c, el.s);
            if (!ImGui::TreeNode(h)) continue;
            for (int which = 0; which < 2; ++which) {
                ImGui::TextUnformatted(which == 0 ? "k' (hệ cục bộ: u, v, φ tại i và j)" : "k = Tᵀ·k'·T (hệ tổng thể, đã kể khớp)");
                if (ImGui::BeginTable(which == 0 ? "##kl" : "##kg", 7, tflags() & ~ImGuiTableFlags_ScrollY)) {
                    ImGui::TableSetupColumn("");
                    for (int j = 0; j < 6; ++j) ImGui::TableSetupColumn(which ? el.dofs[j].c_str() : (j < 3 ? "i" : "j"));
                    ImGui::TableHeadersRow();
                    const auto& k = which == 0 ? el.kLocal : el.kGlobal;
                    for (int i = 0; i < 6; ++i) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        ImGui::TextUnformatted(which ? el.dofs[i].c_str() : "");
                        for (int j = 0; j < 6; ++j) {
                            ImGui::TableNextColumn();
                            ImGui::TextUnformatted(sci(k[i][j]).c_str());
                        }
                    }
                    ImGui::EndTable();
                }
            }
            ImGui::TreePop();
        }
    }
    const int n = (int)mm.Z.size();
    if (ImGui::CollapsingHeader("Ma trận độ cứng K của hệ, véc tơ tải P và nghiệm Z", ImGuiTreeNodeFlags_DefaultOpen) &&
        ImGui::BeginTable("##K", n + 3, tflags(), ImVec2(0, 360))) {
        ImGui::TableSetupScrollFreeze(1, 1);
        ImGui::TableSetupColumn("");
        for (int j = 0; j < n; ++j) ImGui::TableSetupColumn(mm.dofLabels[j].c_str());
        ImGui::TableSetupColumn("P");
        ImGui::TableSetupColumn("Z");
        ImGui::TableHeadersRow();
        for (int i = 0; i < n; ++i) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(mm.dofLabels[i].c_str());
            for (int j = 0; j < n; ++j) {
                ImGui::TableNextColumn();
                if (mm.K[i][j] == 0) ImGui::TextDisabled("0");
                else ImGui::TextUnformatted(sci(mm.K[i][j]).c_str());
            }
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(sci(mm.P[i]).c_str());
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(sci(mm.Z[i]).c_str());
        }
        ImGui::EndTable();
    }
    ImGui::End();
}

void App::archDialog() {
    ImGui::SetNextWindowSize(ImVec2(400, 0), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Tạo vòm", &showArchGen)) {
        ImGui::End();
        return;
    }
    static double L = 24, f = 6, q = -10;
    static int n = 24, shape = 0, sup = 1;
    static bool replace = true;
    ImGui::InputDouble("Nhịp L (m)", &L);
    ImGui::InputDouble("Mũi tên f (m)", &f);
    ImGui::InputInt("Số đoạn", &n);
    n = std::clamp(n, 2, 400);
    const char* shapes[] = {"Parabol y = 4f·x(L − x)/L²", "Cung tròn"};
    ImGui::Combo("Trục vòm", &shape, shapes, 2);
    const char* sups[] = {"Hai khớp", "Ba khớp (khớp đỉnh)", "Hai đầu ngàm"};
    ImGui::Combo("Liên kết", &sup, sups, 3);
    ImGui::InputDouble("Tải đều theo hình chiếu (kN/m)", &q);
    ImGui::Checkbox("Thay thế mô hình hiện tại", &replace);
    ImGui::TextDisabled("Vòm parabol ba khớp chịu tải đều là trục hợp lý: M = 0.");
    if (ImGui::Button("Tạo", ImVec2(-FLT_MIN, 0)) && L > 0 && f > 0) {
        snapshot();
        if (replace) {
            Model nm;
            nm.materials = model.materials;
            nm.sections = model.sections;
            nm.title = std::string("Vòm ") + sups[sup];
            model = nm;
        }
        if (model.materials.empty()) model.materials.push_back(Material{});
        if (model.sections.empty()) model.sections.push_back(rectSection(0.4, 0.8));
        makeArch(model, 0, 0, L, f, n, shape, sup, std::min(curMat, (int)model.materials.size() - 1),
                 std::min(curSec, (int)model.sections.size() - 1));
        if (model.cases.empty()) model.cases.push_back(LoadCase("Tải đều"));
        activeCase = std::min(activeCase, (int)model.cases.size() - 1);
        if (q != 0)
            for (size_t e = 0; e < model.elements.size(); ++e) {
                ElemLoad l;
                l.elem = (int)e;
                l.q1 = l.q2 = q;
                l.projected = true;
                model.cases[activeCase].elem.push_back(l);
            }
        fitRequested = true;
        view = ResultView::M;
        changed();
    }
    ImGui::End();
}
