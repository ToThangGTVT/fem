// Kiểm tra thiết kế (thép TCVN 5575, BTCT TCVN 5574, độ võng) và tạo cọc chịu tải ngang.
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>

#include "App.h"
#include "core/Design.h"

using namespace fem;

namespace {

ImVec4 utilColor(double u) {
    if (u <= 0.7) return ImVec4(0.25f, 0.75f, 0.4f, 1);
    if (u <= 1.0) return ImVec4(0.95f, 0.7f, 0.15f, 1);
    return ImVec4(0.95f, 0.3f, 0.25f, 1);
}

} // namespace

std::vector<int> App::designSets() const {
    std::vector<int> v;
    const int nc = (int)model.cases.size(), ns = (int)result.sets.size();
    if (designSetMode == 0) { // tất cả tổ hợp (hoặc mọi trường hợp tải nếu chưa có tổ hợp)
        for (int i = ns > nc ? nc : 0; i < ns; ++i) v.push_back(i);
    } else if (designSetMode == 1) {
        for (int i = 0; i < ns; ++i) v.push_back(i);
    } else {
        v.push_back(std::clamp(resultSet, 0, std::max(0, ns - 1)));
    }
    return v;
}

void App::designDialog() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x - 640, vp->WorkPos.y + 40), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(620, 640), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Kiểm tra thiết kế", &showDesign)) {
        ImGui::End();
        util.clear();
        return;
    }
    if (!result.ok) {
        ImGui::TextWrapped("Cần có kết quả tính toán (F5).");
        ImGui::End();
        util.clear();
        return;
    }
    const std::vector<Member> members = findMembers(model);
    // khi mở cửa sổ: chọn thẻ theo vật liệu của mô hình
    static bool pick = false;
    if (ImGui::IsWindowAppearing()) pick = true;
    bool hasSteel = false, hasConcrete = false;
    for (const auto& e : model.elements)
        if (e.isBar() && e.mat >= 0 && e.mat < (int)model.materials.size()) {
            hasSteel |= model.materials[e.mat].kind == 1;
            hasConcrete |= model.materials[e.mat].kind == 2;
        }
    ImGui::TextDisabled("Kiểm tra sơ bộ theo các công thức cơ bản của tiêu chuẩn; %zu cấu kiện.", members.size());
    const char* modes[] = {"Mọi tổ hợp (hoặc mọi TH tải nếu chưa có tổ hợp)", "Mọi trường hợp tải và tổ hợp", "Bộ kết quả đang xem"};
    ImGui::SetNextItemWidth(330);
    ImGui::Combo("Bộ tải xét", &designSetMode, modes, 3);
    ImGui::Checkbox("Tô màu hệ số sử dụng trên sơ đồ", &showUtil);
    const std::vector<int> sets = designSets();
    util.assign(model.elements.size(), -1.0);
    auto selectMember = [&](int k) {
        selElems.clear();
        selNodes.clear();
        for (int e : members[k].elems) selElems.insert(e);
    };
    if (ImGui::BeginTabBar("##design")) {
        // ================================================= thép
        if (ImGui::BeginTabItem("Thép TCVN 5575", nullptr, pick && hasSteel ? ImGuiTabItemFlags_SetSelected : 0)) {
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("γc", &steelP.gammaC, 0, 0, "%.3g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("μ mặc định", &steelP.mu, 0, 0, "%.3g");
            ImGui::SameLine();
            const bool haveMu = eig.ok && eig.kind == EigenKind::Buckling && eigVersion == modelVersion;
            ImGui::BeginDisabled(!haveMu);
            ImGui::Checkbox("μ từ phân tích ổn định", &useBucklingMu);
            ImGui::EndDisabled();
            if (!haveMu && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                ImGui::SetTooltip("Chạy Phân tích → Ổn định trước để lấy μ cho từng thanh");
            steelP.muMember.assign(members.size(), 0.0);
            if (haveMu && useBucklingMu)
                for (size_t k = 0; k < members.size(); ++k) {
                    // μ của cấu kiện: μ thanh × (chiều dài thanh / chiều dài cấu kiện)
                    const int e = members[k].elems.front();
                    const double mu = effectiveLengthFactor(model, eig, e);
                    if (mu > 0) steelP.muMember[k] = mu * model.barLength(e) / members[k].L;
                }
            const auto checks = checkSteel(model, result, sets, members, steelP);
            if (checks.empty()) ImGui::TextWrapped("Không có cấu kiện thép (vật liệu loại thép). Đặt loại vật liệu trong thẻ Vật liệu.");
            else if (ImGui::BeginTable("##steel", 9, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY,
                                       ImVec2(0, 0))) {
                ImGui::TableSetupScrollFreeze(0, 1);
                for (const char* h : {"CK", "L (m)", "λ", "φ", "Bền", "Ổn định", "Cắt", "Độ mảnh", "Max"}) ImGui::TableSetupColumn(h);
                ImGui::TableHeadersRow();
                for (const auto& c : checks) {
                    const Member& mb = members[c.member];
                    for (int e : mb.elems) util[e] = c.util;
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    char id[16];
                    std::snprintf(id, sizeof id, "%d", c.member + 1);
                    if (ImGui::Selectable(id, false, ImGuiSelectableFlags_SpanAllColumns)) selectMember(c.member);
                    if (ImGui::IsItemHovered())
                        ImGui::SetTooltip("Thanh %d…; N = %.2f kN, M = %.2f kNm, V = %.2f kN (%s)\n%s", mb.elems.front() + 1, c.N, c.M,
                                          c.V, c.set.c_str(), c.note.c_str());
                    ImGui::TableNextColumn();
                    ImGui::Text("%.2f", mb.L);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.0f", c.lambda);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.3f", c.phi);
                    for (double v : {c.strength, c.stability, c.shear, c.slender}) {
                        ImGui::TableNextColumn();
                        ImGui::TextColored(utilColor(v), "%.2f", v);
                    }
                    ImGui::TableNextColumn();
                    ImGui::TextColored(utilColor(c.util), "%.2f", c.util);
                }
                ImGui::EndTable();
            }
            ImGui::EndTabItem();
        }
        // ================================================= bê tông cốt thép
        if (ImGui::BeginTabItem("BTCT TCVN 5574", nullptr, pick && !hasSteel && hasConcrete ? ImGuiTabItemFlags_SetSelected : 0)) {
            ImGui::SetNextItemWidth(70);
            ImGui::InputDouble("Rs", &rcP.Rs, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70);
            ImGui::InputDouble("Rsc", &rcP.Rsc, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70);
            ImGui::InputDouble("Rsw (MPa)", &rcP.Rsw, 0, 0, "%.4g");
            double acm = rcP.a * 100, aswcm = rcP.Asw * 1e4;
            ImGui::SetNextItemWidth(70);
            if (ImGui::InputDouble("a (cm)", &acm, 0, 0, "%.3g")) rcP.a = acm / 100;
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70);
            if (ImGui::InputDouble("Asw một lớp đai (cm²)", &aswcm, 0, 0, "%.3g")) rcP.Asw = aswcm / 1e4;
            ImGui::SameLine();
            ImGui::SetNextItemWidth(60);
            ImGui::InputDouble("μ cột", &rcP.mu, 0, 0, "%.3g");
            ImGui::TextDisabled("ξR = %.3f; cột chữ nhật đặt thép đối xứng, kể uốn dọc η; cốt đai tính với c = 2h0.",
                                rcXiR(rcP.Rs, rcP.Es));
            const auto rc = designRc(model, result, sets, members, rcP);
            if (rc.empty()) ImGui::TextWrapped("Không có cấu kiện bê tông (vật liệu loại bê tông, tiết diện chữ nhật).");
            else if (ImGui::BeginTable("##rc", 8, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY,
                                       ImVec2(0, 0))) {
                ImGui::TableSetupScrollFreeze(0, 1);
                for (const char* h : {"CK", "Loại", "b×h (cm)", "As dưới", "As trên", "As mỗi bên", "Đai s (cm)", "Ghi chú"})
                    ImGui::TableSetupColumn(h);
                ImGui::TableHeadersRow();
                for (const auto& r : rc) {
                    const Member& mb = members[r.member];
                    // hệ số sử dụng: hàm lượng cốt thép so với 3% (chỉ để tô màu)
                    const double As = r.column ? 2 * r.AsSide : r.AsBottom + r.AsTop;
                    const double rho = r.b > 0 ? As * 1e-4 / (r.b * (r.h - rcP.a)) : 0;
                    for (int e : mb.elems) util[e] = r.shear.ok ? rho / 0.03 : 1.5;
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    char id[16];
                    std::snprintf(id, sizeof id, "%d", r.member + 1);
                    if (ImGui::Selectable(id, false, ImGuiSelectableFlags_SpanAllColumns)) selectMember(r.member);
                    if (ImGui::IsItemHovered())
                        ImGui::SetTooltip("M+ = %.2f, M− = %.2f kNm, Nmax = %.2f kN, Qmax = %.2f kN\nHàm lượng %.2f%%", r.Mpos, r.Mneg,
                                          r.Nmax, r.Qmax, rho * 100);
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(r.column ? "Cột" : "Dầm");
                    ImGui::TableNextColumn();
                    ImGui::Text("%.0f×%.0f", r.b * 100, r.h * 100);
                    ImGui::TableNextColumn();
                    if (!r.column) ImGui::Text("%.2f", r.AsBottom);
                    ImGui::TableNextColumn();
                    if (!r.column) ImGui::Text("%.2f", r.AsTop);
                    ImGui::TableNextColumn();
                    if (r.column) ImGui::Text("%.2f", r.AsSide);
                    ImGui::TableNextColumn();
                    if (!r.shear.ok) ImGui::TextColored(ImVec4(0.95f, 0.3f, 0.25f, 1), "Q > 0,3Rb·b·h0");
                    else if (r.shear.spacing > 0) ImGui::Text("%.1f", r.shear.spacing * 100);
                    else ImGui::TextDisabled("cấu tạo");
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(r.note.c_str());
                }
                ImGui::EndTable();
            }
            ImGui::EndTabItem();
        }
        // ================================================= độ võng
        if (ImGui::BeginTabItem("Độ võng")) {
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("Giới hạn L /", &deflN, 0, 0, "%.4g");
            deflN = std::max(deflN, 1.0);
            ImGui::SameLine();
            ImGui::Checkbox("Chỉ dầm (cấu kiện nằm ngang)", &deflBeamsOnly);
            ImGui::TextDisabled("Độ võng so với dây cung nối hai đầu cấu kiện; dùng bộ tải tiêu chuẩn (TTGH II).");
            double worst = 0;
            int worstSet = -1;
            std::vector<DeflectionCheck> rows(members.size());
            std::vector<std::string> setOf(members.size());
            for (int s : sets) {
                if (s < 0 || s >= (int)result.sets.size()) continue;
                for (const auto& d : checkDeflection(model, result.sets[s], members, deflN, deflBeamsOnly))
                    if (d.ratio >= rows[d.member].ratio) {
                        rows[d.member] = d;
                        setOf[d.member] = result.sets[s].name;
                        if (d.ratio > worst) {
                            worst = d.ratio;
                            worstSet = s;
                        }
                    }
            }
            (void)worstSet;
            if (ImGui::BeginTable("##defl", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, 0))) {
                ImGui::TableSetupScrollFreeze(0, 1);
                for (const char* h : {"CK", "L (m)", "f (mm)", "[f] (mm)", "f/[f]", "Bộ tải"}) ImGui::TableSetupColumn(h);
                ImGui::TableHeadersRow();
                for (size_t k = 0; k < rows.size(); ++k) {
                    const auto& d = rows[k];
                    if (d.member < 0) continue;
                    for (int e : members[k].elems) util[e] = d.ratio;
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    char id[16];
                    std::snprintf(id, sizeof id, "%zu", k + 1);
                    if (ImGui::Selectable(id, false, ImGuiSelectableFlags_SpanAllColumns)) selectMember((int)k);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.2f", d.L);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.2f", d.f * 1e3);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.2f", d.limit * 1e3);
                    ImGui::TableNextColumn();
                    ImGui::TextColored(utilColor(d.ratio), "%.2f", d.ratio);
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(setOf[k].c_str());
                }
                ImGui::EndTable();
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    pick = false;
    ImGui::End();
}

void App::pileDialog() {
    ImGui::SetNextWindowSize(ImVec2(420, 0), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Cọc chịu tải ngang (phương pháp m)", &showPileGen)) {
        ImGui::End();
        return;
    }
    static PileParams p;
    static bool replace = true;
    ImGui::TextWrapped("Hệ số nền tăng tuyến tính theo độ sâu: C_z = m·z (TCVN 10304:2014, Phụ lục A); "
                       "bề rộng quy ước bp = d + 1 (d ≥ 0,8 m), 1,5d + 0,5 (d < 0,8 m).");
    ImGui::InputDouble("Chiều dài trong đất L (m)", &p.L);
    ImGui::InputDouble("Chiều dài tự do l0 (m)", &p.l0);
    ImGui::InputDouble("Đường kính / cạnh d (m)", &p.d);
    ImGui::Checkbox("Cọc vuông", &p.square);
    ImGui::InputDouble("Hệ số m (kN/m⁴)", &p.m);
    ImGui::InputInt("Số đoạn trong đất", &p.segments);
    p.segments = std::clamp(p.segments, 4, 400);
    ImGui::InputDouble("Lực ngang H (kN)", &p.H);
    ImGui::InputDouble("Mômen M đầu cọc (kNm)", &p.Mh);
    ImGui::Checkbox("Thay thế mô hình hiện tại", &replace);
    ImGui::TextDisabled("bp = %.2f m", pileConventionalWidth(p.d));
    if (ImGui::Button("Tạo cọc", ImVec2(-FLT_MIN, 0)) && p.L > 0 && p.d > 0) {
        snapshot();
        if (replace) {
            Model nm;
            nm.materials = model.materials;
            nm.title = "Cọc chịu tải ngang (phương pháp m)";
            model = nm;
        }
        p.mat = std::min(curMat, std::max(0, (int)model.materials.size() - 1));
        makePile(model, p);
        activeCase = resultSet = (int)model.cases.size() - 1;
        fitRequested = true;
        view = ResultView::M;
        changed();
    }
    ImGui::End();
}
