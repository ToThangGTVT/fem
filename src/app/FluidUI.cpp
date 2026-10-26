// Cửa sổ dòng chảy chất lỏng: miền tính (bài toán mẫu hoặc các phần tử phẳng có sẵn), tính chất chất lỏng,
// điều kiện biên, thông số bộ giải FVM (SIMPLE), chạy trực tiếp và xem lực, hệ số Cd, Cl, Cp, lịch sử hội tụ.
#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>

#include "App.h"
#include "core/Fluid.h"
#include "imgui_internal.h"

using namespace fem;

namespace {

struct Medium {
    const char* name;
    double rho, mu;
};
const Medium kMedia[] = {{"Nước 20°C", 998.2, 1.002e-3},
                         {"Không khí 20°C", 1.204, 1.825e-5},
                         {"Dầu động cơ 20°C", 880, 0.29},
                         {"Glyxerin 20°C", 1260, 1.41},
                         {"Tùy chọn", 0, 0}};
constexpr int kMediumCustom = 4;

struct Preset {
    const char* name;
    const char* note;
};
const Preset kPresets[] = {
    {"Hộp kín nắp trượt (lid-driven cavity), Re = 100",
     "Hộp 1×1 m, nắp trên trượt với u = 1 m/s. So sánh: Ghia và cộng sự (1982), u nhỏ nhất trên đường giữa ≈ −0,21."},
    {"Trụ tròn trong kênh, Re = 20 (dừng)",
     "Bài toán chuẩn Schäfer–Turek 2D-1: kênh 2,2×0,41 m, trụ D = 0,1 m. Tham chiếu Cd = 5,58, Cl = 0,0106."},
    {"Trụ tròn trong kênh, Re = 100 (dòng xoáy Kármán)",
     "Schäfer–Turek 2D-2, không dừng: xoáy bong khỏi trụ luân phiên. Tham chiếu Cd max ≈ 3,23, Cl max ≈ 1,0, St ≈ 0,30."},
    {"Trụ vuông trong kênh, Re = 100 (không dừng)",
     "Tiết diện vuông cạnh 0,1 m (dạng mặt cắt cột, trụ cầu): xoáy bong tại các góc sắc."},
    {"Bậc lùi (backward-facing step), Re = 100",
     "Kênh mở rộng đột ngột gấp đôi: vùng xoáy sau bậc và điểm dòng bám lại thành dưới."},
    {"Kênh phẳng (dòng Poiseuille)",
     "Dòng chảy tầng giữa hai bản phẳng: phân bố vận tốc parabol, dp/dx = −8μ·Umax/H²."},
};
constexpr int kPresetCount = (int)(sizeof kPresets / sizeof kPresets[0]);

} // namespace

ImU32 fluidBcColor(FluidBcKind k) {
    switch (k) {
    case FluidBcKind::Wall: return IM_COL32(150, 155, 165, 255);
    case FluidBcKind::MovingWall: return IM_COL32(255, 170, 40, 255);
    case FluidBcKind::Inlet: return IM_COL32(60, 150, 255, 255);
    case FluidBcKind::Outlet: return IM_COL32(240, 80, 70, 255);
    case FluidBcKind::Symmetry: return IM_COL32(80, 200, 120, 255);
    case FluidBcKind::Opening: return IM_COL32(150, 220, 255, 255);
    case FluidBcKind::FreeDrain: return IM_COL32(205, 150, 75, 255);
    }
    return IM_COL32_WHITE;
}

namespace {

// Đồ thị nhiều đường y(x); logY: trục log10
void plotSeries(ImDrawList* dl, ImVec2 p0, ImVec2 sz, const std::vector<double>& x, const std::vector<std::vector<double>>& ys,
                const std::vector<ImU32>& cols, const std::vector<const char*>& names, bool logY, const App& a, const char* xl) {
    dl->AddRectFilled(p0, ImVec2(p0.x + sz.x, p0.y + sz.y), a.pal.bg);
    dl->AddRect(p0, ImVec2(p0.x + sz.x, p0.y + sz.y), a.pal.textDim);
    if (x.size() < 2) return;
    double x0 = x.front(), x1 = x.back(), y0 = DBL_MAX, y1 = -DBL_MAX;
    auto tr = [&](double v) { return logY ? (v > 0 ? std::log10(v) : NAN) : v; };
    for (const auto& y : ys)
        for (double v : y) {
            const double t = tr(v);
            if (!std::isfinite(t)) continue;
            y0 = std::min(y0, t);
            y1 = std::max(y1, t);
        }
    if (y0 > y1) return;
    if (logY) {
        y0 = std::floor(y0);
        y1 = std::ceil(y1);
    }
    if (y1 - y0 < 1e-12) y1 = y0 + 1, y0 -= 1;
    if (x1 <= x0) x1 = x0 + 1;
    const float L = 48, B = 18;
    auto P = [&](double xx, double yy) {
        return ImVec2((float)(p0.x + L + (xx - x0) / (x1 - x0) * (sz.x - L - 6)),
                      (float)(p0.y + 6 + (y1 - yy) / (y1 - y0) * (sz.y - B - 12)));
    };
    if (!logY && y0 < 0 && y1 > 0) dl->AddLine(P(x0, 0), P(x1, 0), a.pal.grid);
    const size_t step = std::max<size_t>(1, x.size() / 1500);
    for (size_t s = 0; s < ys.size(); ++s) {
        std::vector<ImVec2> pts;
        for (size_t i = 0; i < x.size() && i < ys[s].size(); i += step)
            if (std::isfinite(tr(ys[s][i]))) pts.push_back(P(x[i], tr(ys[s][i])));
        dl->AddPolyline(pts.data(), (int)pts.size(), cols[s], 0, 1.5f);
    }
    char b[48];
    std::snprintf(b, sizeof b, logY ? "1e%.0f" : "%.3g", y1);
    dl->AddText(ImVec2(p0.x + 3, p0.y + 3), a.pal.text, b);
    std::snprintf(b, sizeof b, logY ? "1e%.0f" : "%.3g", y0);
    dl->AddText(ImVec2(p0.x + 3, p0.y + sz.y - B - 14), a.pal.text, b);
    std::snprintf(b, sizeof b, "%s %.4g", xl, x1);
    dl->AddText(ImVec2(p0.x + sz.x - ImGui::CalcTextSize(b).x - 4, p0.y + sz.y - B + 2), a.pal.text, b);
    float lx = p0.x + L + 4;
    for (size_t s = 0; s < names.size() && s < cols.size(); ++s) {
        dl->AddLine(ImVec2(lx, p0.y + sz.y - B + 9), ImVec2(lx + 14, p0.y + sz.y - B + 9), cols[s], 2);
        dl->AddText(ImVec2(lx + 17, p0.y + sz.y - B + 2), a.pal.text, names[s]);
        lx += 24 + ImGui::CalcTextSize(names[s]).x;
    }
}

// Ghi chú chữ mờ, tự xuống dòng theo bề rộng cửa sổ
void note(const char* fmt, ...) {
    va_list a;
    va_start(a, fmt);
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
    ImGui::TextWrappedV(fmt, a);
    ImGui::PopStyleColor();
    va_end(a);
}

double nowSec() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

} // namespace

const FluidMesh& App::currentFluidMesh() {
    if (fluidMeshVersion != modelVersion) {
        fluidMeshError = buildFluidMesh(model, fluidMesh);
        fluidMeshVersion = modelVersion;
    }
    return fluidMesh;
}

void App::fluidApplyPreset(int k) {
    ChannelSpec s;
    FluidParams p;
    p.rho = 1;
    p.mu = 1e-3;
    auto bc = [&](FluidSide sd, FluidBcKind kind, double u = 0, double v = 0, bool par = false) {
        p.bc[(int)sd] = {kind, u, v, par};
    };
    bc(FluidSide::Left, FluidBcKind::Inlet, 0.3, 0, true);
    bc(FluidSide::Right, FluidBcKind::Outlet);
    fluidEndTime = 8;
    switch (k) {
    case 0: // hộp kín nắp trượt
        s.obstacle = 0;
        s.L = s.H = 1;
        s.h = 1.0 / 48;
        p.mu = 0.01;
        bc(FluidSide::Left, FluidBcKind::Wall);
        bc(FluidSide::Right, FluidBcKind::Wall);
        bc(FluidSide::Top, FluidBcKind::MovingWall, 1, 0);
        break;
    case 1: break; // Schäfer–Turek 2D-1: thông số mặc định của ChannelSpec
    case 2:
    case 3:
        bc(FluidSide::Left, FluidBcKind::Inlet, 1.5, 0, true);
        if (k == 3) {
            s.obstacle = 2;
            s.D = s.B = 0.1;
        }
        p.transient = true;
        p.dt = 0.01;
        p.alphaU = 0.8;
        p.alphaP = 0.5;
        p.outerIters = 10;
        break;
    case 4: // bậc lùi: bậc cao 0,1 m dài 0,3 m ở góc dưới trái
        s.L = 2.5;
        s.H = 0.2;
        s.obstacle = 2;
        s.xc = 0.15;
        s.D = 0.3;
        s.yc = 0.05;
        s.B = 0.1;
        s.h = 0.01;
        s.hFar = 0.02;
        bc(FluidSide::Left, FluidBcKind::Inlet, 1.5, 0, true);
        break;
    case 5: // kênh phẳng
        s.obstacle = 0;
        s.L = 1;
        s.H = 0.1;
        s.h = 0.01;
        p.mu = 0.01;
        bc(FluidSide::Left, FluidBcKind::Inlet, 1, 0, true);
        break;
    }
    fluidPreset = std::clamp(k, 0, kPresetCount - 1);
    fluidSpec = s;
    fluidP = p;
    fluidMedium = kMediumCustom;

    snapshot();
    Model nm;
    nm.title = kPresets[std::clamp(k, 0, kPresetCount - 1)].name;
    Material mat;
    mat.name = "Chất lỏng";
    mat.E = 1;
    mat.nu = 0.3;
    mat.gamma = p.rho * 9.81e-3;
    nm.materials.push_back(mat);
    nm.cases.push_back(LoadCase("TH1"));
    if (meshChannel(nm, s, 0) <= 0) {
        undoStack.pop_back();
        addLog("Lỗi: không tạo được lưới miền tính.");
        return;
    }
    model = std::move(nm);
    activeCase = resultSet = 0;
    selNodes.clear();
    selElems.clear();
    cmd.reset();
    view = ResultView::None;
    fitRequested = true;
    changed();
    analyzeNow(); // không có gối tựa: lỗi được bỏ qua, chỉ để cập nhật trạng thái
    fluidField = (int)(p.transient ? FluidField::Vorticity : FluidField::Speed);
    if (fluidInit()) fluidRun = true;
}

bool App::fluidInit() {
    fluid = std::make_unique<FluidSolver>();
    ++fluidSerial;
    fluidVersion = modelVersion;
    fluidRun = false;
    const Model* src = &model;
    if (model.flow.active) {
        // miền chảy khai báo: lưới tam giác bám theo các vật cản vẽ bằng thanh khép kín
        std::string note;
        const std::string err = buildFlowProblem(model, fluidMeshModel, fluidP, &note);
        if (!err.empty()) {
            addLog("Lỗi: " + err);
            fluid.reset();
            return false;
        }
        addLog("Dòng chảy: " + note);
        fluidEndTime = fluidP.endTime > 0 ? fluidP.endTime : model.flow.endTime;
        if (fluidP.vof) fluidField = (int)FluidField::Alpha;
        src = &fluidMeshModel;
    }
    if (!fluid->init(*src, fluidP)) {
        addLog("Lỗi: " + fluid->error());
        fluid.reset();
        return false;
    }
    char b[200];
    std::snprintf(b, sizeof b, "Dòng chảy: %d ô, %d mặt; Re = %.4g (U = %.4g m/s, L = %.4g m).", fluid->mesh().nCells,
                  fluid->mesh().nFaces, fluid->reynolds(), fluid->Uref(), fluid->Lref());
    addLog(b);
    return true;
}

void App::editFlow() {
    flowDraft = model.flow;
    flowPopupRequest = true;
}

bool App::flowFields(FlowDomain& f, bool inPanel) {
    bool ch = false;
    // inPanel: mỗi ô sửa thẳng vào mô hình, có hoàn tác (trackEdit)
    auto after = [&]() {
        if (inPanel) trackEdit();
        ch |= ImGui::IsItemEdited();
    };
    auto instant = [&](auto&& apply) { // thay đổi tức thời (nút, ô chọn)
        if (inPanel) snapshot();
        apply();
        if (inPanel) changed();
        ch = true;
    };
    ImGui::PushID(inPanel ? "flowpanel" : "flowpopup");
    const double W = f.width();
    ImGui::TextUnformatted("Loại mặt cắt:");
    ImGui::SameLine();
    if (ImGui::RadioButton("Mặt bằng", f.section == 0)) instant([&] { f.section = 0; });
    ImGui::SameLine();
    if (ImGui::RadioButton("Mặt cắt đứng", f.section == 1)) instant([&] { f.section = 1; });
    double nx, ny;
    f.direction(nx, ny);
    if (f.section == 1 && std::abs(nx) < 0.5)
        ImGui::TextColored(ImVec4(0.95f, 0.3f, 0.25f, 1), "Mặt cắt đứng: dòng chảy phải theo phương ngang (AB thẳng đứng).");
    ImGui::Text(f.section == 0 ? "Cửa vào AB rộng %.4g m" : "Cửa vào AB: độ sâu dòng chảy %.4g m", W);
    ImGui::SetNextItemWidth(90);
    ImGui::InputDouble("Chiều dài miền chảy L (m)", &f.length, 0, 0, "%.4g");
    after();
    f.length = std::max(f.length, 1e-6);
    if (ImGui::Button("Đổi chiều chảy")) instant([&] { f.side = -f.side; });
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Chảy sang phía bên kia của cạnh cửa vào AB");

    ImGui::SeparatorText("Vận tốc");
    ImGui::SetNextItemWidth(90);
    ImGui::InputDouble("U trung bình (m/s)", &f.U, 0, 0, "%.4g");
    after();
    bool par = f.parabolic;
    if (ImGui::Checkbox("Phân bố parabol (Umax = 1,5·U)", &par)) instant([&] { f.parabolic = par; });

    ImGui::SeparatorText("Chất lỏng");
    int med = std::clamp(f.medium, 0, kMediumCustom);
    ImGui::SetNextItemWidth(170);
    if (ImGui::Combo("##medium", &med, [](void*, int i) { return kMedia[i].name; }, nullptr, 5))
        instant([&] {
            f.medium = med;
            if (med != kMediumCustom) {
                f.rho = kMedia[med].rho;
                f.mu = kMedia[med].mu;
            }
        });
    ImGui::SetNextItemWidth(80);
    ImGui::InputDouble("ρ (kg/m³)", &f.rho, 0, 0, "%.5g");
    if (ImGui::IsItemEdited()) f.medium = kMediumCustom;
    after();
    ImGui::SameLine();
    ImGui::SetNextItemWidth(80);
    ImGui::InputDouble("μ (Pa·s)", &f.mu, 0, 0, "%.4g");
    if (ImGui::IsItemEdited()) f.medium = kMediumCustom;
    after();

    if (f.section == 0) {
        ImGui::SeparatorText("Hai bờ");
        if (ImGui::RadioButton("Thành cứng (kênh, bờ kè)", !f.slipWalls)) instant([&] { f.slipWalls = false; });
        ImGui::SameLine();
        if (ImGui::RadioButton("Trượt (sông rộng, dòng tự do)", f.slipWalls)) instant([&] { f.slipWalls = true; });
    } else {
        ImGui::SeparatorText("Biên dưới / biên trên");
        note("Biên dưới: mặt đất / đáy sông (thành cứng). Công trình được phép đặt trên hoặc cắm xuống đáy.");
        bool dr = f.drainGround;
        if (ImGui::Checkbox("Đáy thấm nước (nước thoát qua mặt đất)", &dr)) instant([&] { f.drainGround = dr; });
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Nước thoát ra khỏi miền qua mặt đất / đáy với tốc độ thấm q (vận tốc tiếp tuyến trên đáy vẫn bằng 0). "
                              "Có mặt thoáng tự do: chỉ phần đáy ngập nước thoát.");
        if (f.drainGround) {
            ImGui::Indent();
            if (ImGui::RadioButton("Thấm với tốc độ q", f.drainMode == 0)) instant([&] { f.drainMode = 0; });
            ImGui::SameLine();
            ImGui::BeginDisabled(!f.freeSurface);
            if (ImGui::RadioButton("Thoát tự do (không giới hạn)", f.drainMode == 1)) instant([&] { f.drainMode = 1; });
            ImGui::EndDisabled();
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                ImGui::SetTooltip("Đáy như tấm lưới hở: nước tới đáy bao nhiêu thoát hết bấy nhiêu, rơi tự do theo trọng lực "
                                  "(không bị hút cưỡng bức). Cần mặt thoáng tự do.");
            if (f.drainMode == 1 && !f.freeSurface)
                ImGui::TextColored(ImVec4(0.95f, 0.3f, 0.25f, 1), "Thoát tự do cần Mặt thoáng tự do (VOF).");
            if (f.drainMode == 0) {
                ImGui::SetNextItemWidth(80);
                ImGui::InputDouble("Tốc độ thấm q (m/s)", &f.drainRate, 0, 0, "%.3g");
                after();
                f.drainRate = std::max(f.drainRate, 0.0);
                note("q ≈ hệ số thấm K của đất (thấm thẳng đứng): sỏi ~10⁻², cát thô ~10⁻³, cát mịn ~10⁻⁴ m/s.");
            } else
                note("Nước chạm đáy thoát ngay, không có lực đỡ: vũng nước nằm trên đáy sẽ rơi hết xuống. Chỗ đáy khô, "
                     "không khí đứng yên.");
            ImGui::Unindent();
        }
        const int top = f.freeSurface ? 1 : f.topWall ? 2 : 0;
        if (ImGui::RadioButton("Mặt thoáng phẳng (trượt)", top == 0)) instant([&] { f.freeSurface = f.topWall = false; });
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Mặt nước coi như cố định; hoặc không khí tự do phía trên khi tính gió");
        ImGui::SameLine();
        if (ImGui::RadioButton("Mặt thoáng tự do (VOF)", top == 1)) instant([&] {
                f.freeSurface = true;
                f.topWall = false;
                if (f.waterDepth <= 0 || f.waterDepth >= W) f.waterDepth = 0.6 * W;
            });
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Hai pha nước – không khí: mực nước thay đổi, sóng, nước dâng trước công trình, dòng qua đập tràn, "
                              "vỡ đập. AB gồm cả lớp không khí phía trên mặt nước.");
        ImGui::SameLine();
        if (ImGui::RadioButton("Trần cứng (cống, ống)", top == 2)) instant([&] {
                f.topWall = true;
                f.freeSurface = false;
            });
        if (f.freeSurface) {
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("Độ sâu nước ở cửa vào (m)", &f.waterDepth, 0, 0, "%.4g");
            after();
            if (f.waterDepth <= 0 || f.waterDepth >= W)
                ImGui::TextColored(ImVec4(0.95f, 0.3f, 0.25f, 1), "Độ sâu phải > 0 và < chiều cao AB (%.4g m) để có lớp không khí.", W);
            bool same = f.initDepth < 0;
            if (ImGui::Checkbox("Mực nước ban đầu bằng ở cửa vào", &same)) instant([&] { f.initDepth = same ? -1 : f.waterDepth; });
            if (!same) {
                ImGui::SetNextItemWidth(80);
                ImGui::InputDouble("Độ sâu nước ban đầu (m)", &f.initDepth, 0, 0, "%.4g");
                after();
                f.initDepth = std::max(f.initDepth, 0.0);
            }
            ImGui::SetNextItemWidth(80);
            ImGui::InputDouble("Nước ban đầu chỉ trong đoạn dài (m)", &f.gate, 0, 0, "%.4g");
            after();
            f.gate = std::max(f.gate, 0.0);
            note("0: toàn miền. > 0: cột nước sau cửa đập (vỡ đập, xả cửa cống), phần sau đó ban đầu khô.");
            if (f.medium == 1)
                ImGui::TextColored(ImVec4(0.95f, 0.3f, 0.25f, 1), "Mặt thoáng tự do cần chất lỏng nặng (nước, dầu...).");
        } else {
            bool hs = f.hydrostatic;
            if (ImGui::Checkbox("Cộng áp suất thủy tĩnh ρ·g·h (mặt thoáng tại biên trên)", &hs)) instant([&] { f.hydrostatic = hs; });
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Lực lên công trình gồm cả áp lực tĩnh và lực đẩy nổi; vận tốc không đổi. Bỏ chọn với không khí.");
        }
    }

    ImGui::SeparatorText("Mô hình dòng chảy");
    if (ImGui::RadioButton("Tự động", f.turbulence == 0)) instant([&] { f.turbulence = 0; });
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Chảy rối k-ε khi Re > 2300, chảy tầng khi Re nhỏ");
    ImGui::SameLine();
    if (ImGui::RadioButton("Chảy tầng", f.turbulence == 1)) instant([&] { f.turbulence = 1; });
    ImGui::SameLine();
    if (ImGui::RadioButton("Chảy rối k-ε", f.turbulence == 2)) instant([&] { f.turbulence = 2; });
    if (f.turbulence != 1) {
        double I = f.turbIntensity * 100;
        ImGui::SetNextItemWidth(80);
        if (ImGui::InputDouble("Cường độ rối ở cửa vào (%)", &I, 0, 0, "%.3g")) f.turbIntensity = std::clamp(I, 0.1, 50.0) / 100;
        after();
    }

    ImGui::SeparatorText("Lưới tam giác");
    ImGui::SetNextItemWidth(80);
    ImGui::InputDouble("Ô sát vật cản (m)", &f.h, 0, 0, "%.4g");
    after();
    ImGui::SameLine();
    ImGui::SetNextItemWidth(80);
    ImGui::InputDouble("Ô xa (m)", &f.hFar, 0, 0, "%.4g");
    after();
    f.h = std::max(f.h, 0.0);
    f.hFar = std::max(f.hFar, 0.0);
    note("0: tự chọn theo bề rộng dòng chảy (W/40) và kích thước vật cản.");

    ImGui::SeparatorText("Phân tích");
    const bool vof = f.section == 1 && f.freeSurface;
    ImGui::BeginDisabled(vof);
    if (ImGui::RadioButton("Dừng", !f.transient && !vof)) instant([&] { f.transient = false; });
    ImGui::SameLine();
    if (ImGui::RadioButton("Theo thời gian", f.transient || vof)) instant([&] { f.transient = true; });
    ImGui::EndDisabled();
    if (vof) note("Mặt thoáng tự do luôn tính theo thời gian; Δt tự giảm theo số Courant (Δt dưới đây là bước lớn nhất).");
    if (f.transient || vof) {
        ImGui::SetNextItemWidth(80);
        ImGui::InputDouble("Δt (s)", &f.dt, 0, 0, "%.4g");
        after();
        ImGui::SameLine();
        ImGui::SetNextItemWidth(80);
        ImGui::InputDouble("Đến t (s)", &f.endTime, 0, 0, "%.4g");
        after();
        f.dt = std::max(f.dt, 1e-6);
    }
    if (f.rho > 0 && f.mu > 0) {
        const double Re = f.rho * f.U * W / f.mu;
        note(f.section == 0 ? "Re theo bề rộng dòng chảy = ρ·U·W/μ = %.4g" : "Re theo chiều cao miền = ρ·U·H/μ = %.4g", Re);
        const bool turb = f.turbulence == 2 || (f.turbulence == 0 && Re > 2300);
        if (turb) note("→ chảy rối, mô hình k-ε với hàm tường%s.", f.transient ? "" : " (bài toán dừng dùng sơ đồ bậc nhất)");
        else if (Re > 2300)
            ImGui::TextColored(ImVec4(0.95f, 0.6f, 0.1f, 1),
                               "Re lớn mà chọn chảy tầng: chỉ cho kết quả định tính; bài toán dừng sẽ được tính theo thời gian.");
    }
    ImGui::PopID();
    return ch;
}

void App::flowPopup() {
    const char* title = "Khai báo dòng chảy";
    if (flowPopupRequest) {
        ImGui::OpenPopup(title);
        flowPopupRequest = false;
    }
    if (!ImGui::BeginPopupModal(title, nullptr, ImGuiWindowFlags_AlwaysAutoResize)) return;
    flowFields(flowDraft, false);
    ImGui::Separator();
    const bool ok = ImGui::Button("Đồng ý", ImVec2(120, 0)) || ImGui::IsKeyPressed(ImGuiKey_Enter) ||
                    ImGui::IsKeyPressed(ImGuiKey_KeypadEnter);
    ImGui::SameLine();
    const bool cancel = ImGui::Button("Hủy", ImVec2(120, 0)) || ImGui::IsKeyPressed(ImGuiKey_Escape);
    if (ok && flowDraft.U > 0 && flowDraft.rho > 0 && flowDraft.mu > 0) {
        const bool first = !model.flow.active;
        snapshot();
        model.flow = flowDraft;
        model.flow.active = true;
        if (model.cases.empty()) model.cases.push_back(LoadCase("TH1"));
        changed();
        fitRequested = first;
        if (first)
            echo("Đã khai báo dòng chảy. Vẽ vật cản bằng lệnh L (đường khép kín, gõ C để khép), rồi bấm Tính (F5).");
        else echo("Đã cập nhật dòng chảy. Bấm Tính (F5) để tính lại.");
        ImGui::CloseCurrentPopup();
    } else if (ok)
        ImGui::TextColored(ImVec4(0.95f, 0.3f, 0.25f, 1), "U, ρ, μ phải > 0.");
    if (cancel) ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
}

void App::fluidDialog() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x - 560, vp->WorkPos.y + 40), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(540, 760), ImGuiCond_FirstUseEver);
    // lần đầu mở trong phiên: ghép thành thẻ cạnh panel Kết quả để không che sơ đồ
    static bool docked = false;
    // (chờ vài khung hình: bố cục dock dựng ở khung hình đầu tiên)
    if (ImGuiWindow* res = ImGui::FindWindowByName("Kết quả"); !docked && res && res->DockId && ImGui::GetFrameCount() > 3) {
        ImGui::SetNextWindowDockID(res->DockId, ImGuiCond_Always);
        docked = true;
    }
    if (!ImGui::Begin("Dòng chảy chất lỏng (FVM)", &showFluid)) {
        ImGui::End();
        return;
    }
    const ImVec4 warn(0.95f, 0.6f, 0.1f, 1);

    // ---------------------------------------------------------------- chạy (mỗi khung hình tối đa ~30 ms)
    if (fluidRun && fluid && fluidVersion == modelVersion) {
        const double t0 = nowSec();
        const bool steady = !fluid->params().transient;
        while (nowSec() - t0 < 0.03) {
            if (steady && fluid->converged()) break;
            if (!steady && fluid->time() >= fluidEndTime - 1e-9) break;
            if (!fluid->step()) break;
        }
        if (fluid->diverged()) {
            fluidRun = false;
            addLog("Lỗi: lời giải phân kỳ — giảm hệ số giảm dư αu, αp, giảm Δt hoặc làm mịn lưới.");
        } else if (steady && fluid->converged()) {
            fluidRun = false;
            addLog("Dòng chảy: hội tụ sau " + std::to_string(fluid->iteration()) + " vòng lặp.");
        } else if (!steady && fluid->time() >= fluidEndTime - 1e-9) {
            fluidRun = false;
            addLog("Dòng chảy: đã tính đến t = " + std::to_string(fluid->time()).substr(0, 6) + " s.");
        }
    } else if (fluidVersion != modelVersion)
        fluidRun = false;

    // ---------------------------------------------------------------- miền chảy khai báo (Trang chủ › Dòng chảy)
    if (model.flow.active) {
        const FlowDomain& f = model.flow;
        ImGui::SeparatorText("Dòng chảy khai báo");
        ImGui::TextWrapped("%s: %s %.4g m, dài %.4g m; U = %.4g m/s%s; %s%s; ρ = %.5g kg/m³, μ = %.4g Pa·s; %s.",
                           f.section == 0 ? "Mặt bằng" : "Mặt cắt đứng", f.section == 0 ? "rộng" : "cao", f.width(), f.length,
                           f.U, f.parabolic ? " (parabol)" : "",
                           f.section == 0 ? (f.slipWalls ? "bờ trượt" : "bờ cứng")
                           : f.freeSurface ? "mặt thoáng tự do (VOF)"
                           : f.topWall     ? "đáy cứng, trần cứng"
                                           : "đáy cứng, mặt thoáng phẳng",
                           f.section != 1 || !f.drainGround ? ""
                           : f.drainMode == 1               ? ", đáy thoát tự do"
                           : f.drainRate > 0                ? ", đáy thấm nước"
                                                            : "",
                           f.rho, f.mu, f.turbulence == 1 ? "chảy tầng" : f.turbulence == 2 ? "chảy rối k-ε" : "mô hình tự động");
        if (ImGui::Button("Sửa khai báo...")) editFlow();
        ImGui::SameLine();
        note("Vật cản: các đường khép kín vẽ bằng lệnh L bên trong miền chảy.");
    } else {
    if (ImGui::CollapsingHeader("Miền tính và lưới", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::BeginCombo("##preset", kPresets[fluidPreset].name)) {
            for (int i = 0; i < kPresetCount; ++i)
                if (ImGui::Selectable(kPresets[i].name, i == fluidPreset)) fluidPreset = i;
            ImGui::EndCombo();
        }
        note("%s", kPresets[fluidPreset].note);
        if (ImGui::Button("Tạo bài toán mẫu và chạy")) fluidApplyPreset(fluidPreset);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Thay mô hình hiện tại bằng miền tính, lưới và điều kiện biên của bài toán mẫu");

        if (ImGui::TreeNode("Kênh chữ nhật có vật cản (tự tạo lưới)")) {
            ChannelSpec& s = fluidSpec;
            ImGui::PushItemWidth(80);
            ImGui::InputDouble("Dài L (m)", &s.L, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::InputDouble("Cao H (m)", &s.H, 0, 0, "%.4g");
            const char* obs[] = {"Không", "Trụ tròn", "Hình chữ nhật"};
            ImGui::SetNextItemWidth(140);
            ImGui::Combo("Vật cản", &s.obstacle, obs, 3);
            if (s.obstacle) {
                ImGui::InputDouble("Tâm x (m)", &s.xc, 0, 0, "%.4g");
                ImGui::SameLine();
                ImGui::InputDouble("Tâm y (m)", &s.yc, 0, 0, "%.4g");
                ImGui::InputDouble(s.obstacle == 1 ? "Đường kính D (m)" : "Rộng (m)", &s.D, 0, 0, "%.4g");
                if (s.obstacle == 2) {
                    ImGui::SameLine();
                    ImGui::InputDouble("Cao (m)", &s.B, 0, 0, "%.4g");
                }
            }
            ImGui::InputDouble("Ô gần vật cản (m)", &s.h, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::InputDouble("Ô xa (m)", &s.hFar, 0, 0, "%.4g");
            ImGui::PopItemWidth();
            ImGui::Checkbox("Lưới tam giác (chia mỗi tứ giác thành 2 tam giác)", &s.triangles);
            note("Hình chữ nhật được phép chạm biên: bậc, công trình đặt trên mặt đất.");
            if (ImGui::Button("Tạo lưới (thay mô hình hiện tại)")) {
                snapshot();
                Model nm;
                nm.title = "Miền chất lỏng";
                Material mat;
                mat.name = "Chất lỏng";
                mat.E = 1;
                nm.materials.push_back(mat);
                nm.cases.push_back(LoadCase("TH1"));
                const int n = meshChannel(nm, s, 0);
                if (n <= 0) {
                    undoStack.pop_back();
                    addLog("Lỗi: kích thước không hợp lệ hoặc vật cản không nằm trong kênh (trụ tròn cần cách biên ≥ 0,6 D).");
                } else {
                    model = std::move(nm);
                    activeCase = resultSet = 0;
                    selNodes.clear();
                    selElems.clear();
                    view = ResultView::None;
                    fitRequested = true;
                    changed();
                    addLog("Đã tạo lưới miền chất lỏng: " + std::to_string(n) + " phần tử.");
                }
            }
            ImGui::TreePop();
        }
        const FluidMesh& fm = currentFluidMesh();
        if (!fluidMeshError.empty()) ImGui::TextColored(warn, "%s", fluidMeshError.c_str());
        else
            note("Miền tính: %d ô (phần tử phẳng), %d mặt. Có thể dùng lưới bất kỳ từ Lưới phẳng / DXF, "
                                "xóa phần tử để tạo vật cản.",
                                fm.nCells, fm.nFaces);
    }

    // ---------------------------------------------------------------- chất lỏng
    if (ImGui::CollapsingHeader("Chất lỏng", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::SetNextItemWidth(170);
        if (ImGui::Combo("##medium", &fluidMedium, [](void*, int i) { return kMedia[i].name; }, nullptr, 5) &&
            fluidMedium != kMediumCustom) {
            fluidP.rho = kMedia[fluidMedium].rho;
            fluidP.mu = kMedia[fluidMedium].mu;
        }
        ImGui::SetNextItemWidth(70);
        if (ImGui::InputDouble("ρ (kg/m³)", &fluidP.rho, 0, 0, "%.5g")) fluidMedium = kMediumCustom;
        ImGui::SameLine();
        ImGui::SetNextItemWidth(70);
        if (ImGui::InputDouble("μ (Pa·s)", &fluidP.mu, 0, 0, "%.4g")) fluidMedium = kMediumCustom;
        note("Độ nhớt động học ν = μ/ρ = %.4g m²/s", fluidP.rho > 0 ? fluidP.mu / fluidP.rho : 0.0);
    }

    // ---------------------------------------------------------------- điều kiện biên
    if (ImGui::CollapsingHeader("Điều kiện biên", ImGuiTreeNodeFlags_DefaultOpen)) {
        const FluidMesh& fm = currentFluidMesh();
        // mỗi đoạn biên: tên + loại; dòng thứ hai (khi cần) nhập vận tốc
        const float nameW = ImGui::CalcTextSize("Trong (vật cản)").x + 40;
        for (int s = 0; s < (int)FluidSide::Count; ++s) {
            FluidBc& b = fluidP.bc[s];
            ImGui::PushID(s);
            const ImVec2 c = ImGui::GetCursorScreenPos();
            const float fh = ImGui::GetFrameHeight();
            ImGui::GetWindowDrawList()->AddRectFilled(ImVec2(c.x, c.y + fh * 0.3f), ImVec2(c.x + 10, c.y + fh * 0.7f), fluidBcColor(b.kind));
            ImGui::Dummy(ImVec2(14, fh));
            ImGui::SameLine();
            ImGui::AlignTextToFramePadding();
            const int nf = fm.nCells ? fm.sideFaces[s] : 0;
            if (nf) ImGui::TextUnformatted(fluidSideName((FluidSide)s));
            else ImGui::TextDisabled("%s", fluidSideName((FluidSide)s));
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("%d mặt biên", nf);
            ImGui::SameLine(nameW);
            int k = (int)b.kind;
            ImGui::SetNextItemWidth(std::max(120.0f, ImGui::GetContentRegionAvail().x));
            if (ImGui::Combo("##k", &k, [](void*, int i) { return fluidBcName((FluidBcKind)i); }, nullptr, 6)) b.kind = (FluidBcKind)k;
            if (b.kind == FluidBcKind::Inlet || b.kind == FluidBcKind::MovingWall) {
                const float indent = nameW - ImGui::GetStyle().WindowPadding.x;
                ImGui::SetCursorPosX(ImGui::GetCursorPosX() + indent);
                ImGui::SetNextItemWidth(55);
                ImGui::InputDouble("u##u", &b.u, 0, 0, "%.4g");
                ImGui::SameLine();
                ImGui::SetNextItemWidth(55);
                ImGui::InputDouble("v (m/s)##v", &b.v, 0, 0, "%.4g");
                if (b.kind == FluidBcKind::Inlet) {
                    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + indent);
                    ImGui::Checkbox("phân bố parabol", &b.parabolic);
                }
            }
            ImGui::PopID();
        }
        note("Bốn cạnh của hình chữ nhật bao miền; biên trong là mọi cạnh biên còn lại (vật cản, bậc, lỗ).");
        note("Cửa vào parabol: u, v là vận tốc lớn nhất (vận tốc trung bình = 2/3).");
    }

    } // !model.flow.active

    // ---------------------------------------------------------------- bộ giải
    if (ImGui::CollapsingHeader("Phương pháp giải")) {
        int sch = (int)fluidP.scheme;
        const char* schemes[] = {"Ngược dòng bậc nhất (ổn định)", "Ngược dòng tuyến tính bậc hai (chính xác)"};
        ImGui::SetNextItemWidth(280);
        if (ImGui::Combo("Sơ đồ đối lưu", &sch, schemes, 2)) fluidP.scheme = (ConvScheme)sch;
        if (!model.flow.active) {
            ImGui::Checkbox("Không dừng (theo thời gian)", &fluidP.transient);
            int tm = (int)fluidP.turbulence;
            const char* tms[] = {"Chảy tầng", "Chảy rối k-ε (hàm tường)"};
            ImGui::SetNextItemWidth(220);
            if (ImGui::Combo("Mô hình dòng chảy", &tm, tms, 2)) fluidP.turbulence = (TurbModel)tm;
        }
        ImGui::PushItemWidth(80);
        if (fluidP.transient && model.flow.active) {
            ImGui::InputInt("Số vòng lặp ngoài tối đa mỗi bước", &fluidP.outerIters);
            fluidP.outerIters = std::clamp(fluidP.outerIters, 1, 200);
        } else if (fluidP.transient) {
            ImGui::InputDouble("Δt (s)", &fluidP.dt, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::InputDouble("Đến t (s)", &fluidEndTime, 0, 0, "%.4g");
            int ts = (int)fluidP.timeScheme;
            const char* tss[] = {"Euler ẩn (bậc 1)", "BDF2 (bậc 2)"};
            ImGui::SetNextItemWidth(160);
            if (ImGui::Combo("Sơ đồ thời gian", &ts, tss, 2)) fluidP.timeScheme = (TimeScheme)ts;
            ImGui::InputInt("Số vòng lặp ngoài tối đa mỗi bước", &fluidP.outerIters);
            fluidP.outerIters = std::clamp(fluidP.outerIters, 1, 200);
        } else {
            ImGui::InputDouble("Sai số dừng", &fluidP.tol, 0, 0, "%.1e");
        }
        ImGui::InputDouble("αu", &fluidP.alphaU, 0, 0, "%.2f");
        ImGui::SameLine();
        ImGui::InputDouble("αp", &fluidP.alphaP, 0, 0, "%.2f");
        ImGui::SameLine();
        note("hệ số giảm dư (SIMPLE)");
        if (!model.flow.active) {
            ImGui::InputDouble("U tham chiếu", &fluidP.Uref, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::InputDouble("L tham chiếu", &fluidP.Lref, 0, 0, "%.4g");
        }
        ImGui::PopItemWidth();
        if (!model.flow.active) note("0: tự xác định (vận tốc cửa vào / thành di động; chiều cao vật cản).");
    }

    // ---------------------------------------------------------------- điều khiển
    ImGui::Separator();
    const bool stale = fluid && fluidVersion != modelVersion;
    if (ImGui::Button("Khởi tạo")) fluidInit();
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Dựng lưới FVM từ mô hình hiện tại và đặt vận tốc, áp suất ban đầu bằng 0");
    ImGui::SameLine();
    if (!fluidRun) {
        if (ImGui::Button("▶ Chạy") && (fluid && !stale ? true : fluidInit())) {
            if (fluid->params().transient && fluid->time() >= fluidEndTime - 1e-9) fluidEndTime = fluid->time() + 5;
            fluidRun = true;
        }
    } else if (ImGui::Button("■ Dừng"))
        fluidRun = false;
    ImGui::SameLine();
    ImGui::BeginDisabled(!fluid || stale || fluidRun);
    if (ImGui::Button("Một bước")) fluid->step();
    ImGui::EndDisabled();
    if (stale) ImGui::TextColored(warn, "Mô hình hoặc lưới đã thay đổi — bấm Khởi tạo.");

    if (fluid && !stale) {
        const FluidSolver& S = *fluid;
        const auto& H = S.history();
        if (S.params().transient)
            ImGui::Text("t = %.4g s, bước %d   Re = %.4g", S.time(), S.iteration(), S.reynolds());
        else
            ImGui::Text("Vòng lặp %d%s   Re = %.4g", S.iteration(), S.converged() ? " — đã hội tụ" : "", S.reynolds());
        if (!H.empty())
            note("Sai số dư: u %.2e   v %.2e   liên tục %.2e", H.back().resU, H.back().resV, H.back().resC);
        if (S.diverged()) ImGui::TextColored(ImVec4(0.95f, 0.3f, 0.25f, 1), "Phân kỳ — giảm αu, αp hoặc Δt rồi khởi tạo lại.");

        // ---- hiển thị
        ImGui::SeparatorText("Hiển thị");
        // các đại lượng có nghĩa với bài toán đang tính
        std::vector<std::pair<int, const char*>> fields = {{(int)FluidField::Speed, "|u| vận tốc"}, {(int)FluidField::U, "u"},
                                                           {(int)FluidField::V, "v"},
                                                           {(int)FluidField::P, S.params().vof ? "p áp suất (thực)" : "p áp suất"},
                                                           {(int)FluidField::Cp, "Cp"}, {(int)FluidField::Vorticity, "ω xoáy"}};
        if (S.params().hydroG > 0) fields.push_back({(int)FluidField::PTotal, "p toàn phần"}); // p + ρ·g·h
        if (S.params().vof) fields.push_back({(int)FluidField::Alpha, "α tỉ lệ nước"});
        if (S.params().turbulence != TurbModel::Laminar) fields.push_back({(int)FluidField::TurbK, "k rối"});
        bool known = false;
        for (auto& fl : fields) known |= fl.first == fluidField;
        if (!known) fluidField = (int)FluidField::Speed;
        for (size_t i = 0; i < fields.size(); ++i) {
            if (i % 3) ImGui::SameLine(); // hai hàng, vừa panel hẹp
            ImGui::RadioButton(fields[i].second, &fluidField, fields[i].first);
        }
        if (S.params().vof)
            note("Thể tích nước %.5g m² (1 m bề dày); Δt hiện tại %.3g s; đường trắng: mặt nước α = 0,5.", S.waterVolume(), S.dt());
        if (S.params().vof && S.drainedVolume() > 0) note("Nước đã thoát qua đáy: %.5g m²", S.drainedVolume());
        ImGui::Checkbox("Véc tơ vận tốc", &fluidVectors);
        ImGui::SameLine();
        ImGui::Checkbox("Đường dòng", &fluidStreams);
        ImGui::SameLine();
        ImGui::Checkbox("Lưới", &fluidMeshLines);

        // ---- lực lên vật cản
        const FluidForces F = S.forces();
        if (F.valid) {
            ImGui::SeparatorText("Lực chất lỏng lên biên trong (1 m bề dày)");
            ImGui::TextWrapped("Lực cản (theo hướng chảy) = %.5g N/m, lực nâng = %.5g N/m", F.Fd, F.Fl);
            note("Theo trục tổng thể: Fx = %.5g (áp suất %.4g, nhớt %.4g), Fy = %.5g (áp suất %.4g, nhớt %.4g) N/m", F.Fx,
                 F.Fpx, F.Fvx, F.Fy, F.Fpy, F.Fvy);
            ImGui::Text("Cd = %.4f   Cl = %.4f", F.Cd, F.Cl);
            if (S.params().transient) {
                const double st = S.strouhal();
                if (st > 0) {
                    ImGui::SameLine();
                    ImGui::Text("  St = f·L/U = %.4f", st);
                }
            }
            note("Cd = lực cản / (½ρU²L), U = %.4g m/s, L = %.4g m (bề rộng vật cản vuông góc hướng chảy); kN/m = N/m × 10⁻³",
                 S.Uref(), S.Lref());
        }

        // ---- đồ thị
        if (H.size() >= 2) {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const float w = ImGui::GetContentRegionAvail().x;
            std::vector<double> x, ru, rv, rc;
            for (const auto& h : H) {
                x.push_back(S.params().transient ? h.time : h.iter);
                ru.push_back(h.resU);
                rv.push_back(h.resV);
                rc.push_back(h.resC);
            }
            ImGui::SeparatorText("Sai số dư");
            ImVec2 p0 = ImGui::GetCursorScreenPos();
            ImGui::Dummy(ImVec2(w, 130));
            plotSeries(dl, p0, ImVec2(w, 130), x, {ru, rv, rc},
                       {IM_COL32(80, 160, 255, 255), IM_COL32(90, 210, 120, 255), IM_COL32(255, 120, 60, 255)},
                       {"u", "v", "liên tục"}, true, *this, S.params().transient ? "t (s)" : "vòng lặp");
            if (F.valid) {
                std::vector<double> cd, cl;
                for (const auto& h : H) {
                    cd.push_back(h.Cd);
                    cl.push_back(h.Cl);
                }
                ImGui::SeparatorText("Hệ số lực");
                p0 = ImGui::GetCursorScreenPos();
                ImGui::Dummy(ImVec2(w, 130));
                plotSeries(dl, p0, ImVec2(w, 130), x, {cd, cl}, {IM_COL32(240, 80, 70, 255), IM_COL32(80, 160, 255, 255)},
                           {"Cd", "Cl"}, false, *this, S.params().transient ? "t (s)" : "vòng lặp");
                const auto cp = S.innerCp();
                if (cp.size() > 2) {
                    std::vector<double> th, v;
                    for (const auto& c : cp) {
                        th.push_back(c[0]);
                        v.push_back(c[1]);
                    }
                    ImGui::SeparatorText("Hệ số áp suất Cp trên biên trong");
                    p0 = ImGui::GetCursorScreenPos();
                    ImGui::Dummy(ImVec2(w, 130));
                    plotSeries(dl, p0, ImVec2(w, 130), th, {v}, {IM_COL32(255, 170, 40, 255)}, {"Cp"}, false, *this, "θ (°)");
                    note("θ tính từ điểm đón dòng (phía −x), dương ở nửa trên.");
                }
            }
        }
    }
    ImGui::End();
    if (!showFluid) fluidRun = false;
}
