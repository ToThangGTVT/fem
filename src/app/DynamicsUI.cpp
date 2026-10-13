// Cửa sổ động lực học: phổ phản ứng TCVN 9386, dao động cưỡng bức điều hòa, phân tích theo thời gian (Newmark).
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <sstream>

#include "App.h"
#include "core/Dynamics.h"
#include "imgui_stdlib.h"

using namespace fem;

namespace {

constexpr double kPi = 3.14159265358979323846;

bool setCombo(const char* label, int& v, const Model& m) {
    std::vector<std::string> sets;
    for (auto& c : m.cases) sets.push_back("TH: " + c.name);
    for (auto& c : m.combos) sets.push_back("THợp: " + c.name);
    if (sets.empty()) return false;
    v = std::clamp(v, 0, (int)sets.size() - 1);
    bool ch = false;
    ImGui::SetNextItemWidth(240);
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

// Đồ thị đơn giản y(x)
void plot(ImDrawList* dl, ImVec2 p0, ImVec2 sz, const std::vector<double>& x, const std::vector<double>& y, ImU32 col, ImU32 axis,
          ImU32 text, const char* xl, const char* yl, double markX = NAN) {
    dl->AddRect(p0, ImVec2(p0.x + sz.x, p0.y + sz.y), axis);
    if (x.size() < 2) return;
    double x0 = x.front(), x1 = x.back(), y0 = 0, y1 = 0;
    for (double v : y) {
        y0 = std::min(y0, v);
        y1 = std::max(y1, v);
    }
    if (y1 - y0 < 1e-30) y1 = y0 + 1;
    const float L = 46, B = 18;
    auto P = [&](double xx, double yy) {
        return ImVec2((float)(p0.x + L + (xx - x0) / (x1 - x0) * (sz.x - L - 6)),
                      (float)(p0.y + 6 + (y1 - yy) / (y1 - y0) * (sz.y - B - 12)));
    };
    dl->AddLine(P(x0, 0), P(x1, 0), axis);
    std::vector<ImVec2> pts;
    const size_t step = std::max<size_t>(1, x.size() / 2000);
    for (size_t i = 0; i < x.size(); i += step) pts.push_back(P(x[i], y[i]));
    dl->AddPolyline(pts.data(), (int)pts.size(), col, 1.6f);
    char b[48];
    std::snprintf(b, sizeof b, "%.3g", y1);
    dl->AddText(ImVec2(p0.x + 3, p0.y + 3), text, b);
    std::snprintf(b, sizeof b, "%.3g", y0);
    dl->AddText(ImVec2(p0.x + 3, p0.y + sz.y - B - 14), text, b);
    std::snprintf(b, sizeof b, "%s %.3g", xl, x1);
    dl->AddText(ImVec2(p0.x + sz.x - ImGui::CalcTextSize(b).x - 4, p0.y + sz.y - B + 2), text, b);
    dl->AddText(ImVec2(p0.x + L, p0.y + sz.y - B + 2), text, yl);
    if (!std::isnan(markX)) {
        const ImVec2 a = P(markX, y0), c = P(markX, y1);
        dl->AddLine(a, c, IM_COL32(255, 160, 40, 255), 1.5f);
    }
}

} // namespace

const LoadSetResult* App::overrideShown(std::string* title) {
    if (const LoadSetResult* f = forceMethodShown(title)) return f;
    if (!showDynamics || dynVersion != modelVersion) return nullptr;
    const LoadSetResult* s = nullptr;
    switch (dynShow) {
    case 1:
        if (spec.ok) s = &spec.combined;
        break;
    case 2:
        if (spec.ok && specMode >= 0 && specMode < (int)spec.modal.size()) s = &spec.modal[specMode];
        break;
    case 3:
        if (harm.ok) s = &harm.amp;
        break;
    case 4:
        if (th.ok && !thSnap.U.empty()) s = &thSnap;
        break;
    default: break;
    }
    if (s && title) *title = "Động lực học — " + s->name;
    return s;
}

void App::dynamicsDialog() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x - 600, vp->WorkPos.y + 40), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(580, 680), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Động lực học", &showDynamics)) {
        ImGui::End();
        return;
    }
    if (dynVersion != modelVersion && (spec.ok || harm.ok || th.ok)) {
        ImGui::TextColored(ImVec4(0.9f, 0.6f, 0.1f, 1), "Mô hình đã thay đổi — hãy tính lại.");
    }
    // ---- khối lượng (chung)
    if (ImGui::CollapsingHeader("Khối lượng và rời rạc hóa", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::SetNextItemWidth(90);
        ImGui::InputDouble("Hệ số khối lượng bản thân (γ/g)", &dynOpt.selfMassFactor, 0, 0, "%.3g");
        bool useSet = dynOpt.massSet >= 0;
        if (ImGui::Checkbox("Quy tải thẳng đứng thành khối lượng:", &useSet)) dynOpt.massSet = useSet ? 0 : -1;
        if (dynOpt.massSet >= 0) {
            ImGui::SameLine();
            setCombo("##massset", dynOpt.massSet, model);
            ImGui::SetNextItemWidth(90);
            ImGui::InputDouble("Hệ số (ví dụ TT + 0,5 HT: dùng tổ hợp)", &dynOpt.massSetFactor, 0, 0, "%.3g");
        }
        ImGui::SetNextItemWidth(90);
        ImGui::InputInt("Số dạng dao động", &dynOpt.modes);
        dynOpt.modes = std::clamp(dynOpt.modes, 1, 60);
        ImGui::SetNextItemWidth(90);
        ImGui::InputInt("Số đoạn chia mỗi thanh", &dynOpt.subdiv);
        dynOpt.subdiv = std::clamp(dynOpt.subdiv, 1, 20);
    }
    auto showRadio = [&](const char* label, int id) {
        if (ImGui::RadioButton(label, dynShow == id)) {
            dynShow = id;
            if (view != ResultView::N && view != ResultView::Q && view != ResultView::M && view != ResultView::Deformed)
                view = ResultView::M;
        }
    };
    if (ImGui::BeginTabBar("##dyn")) {
        // ================================================= phổ phản ứng
        if (ImGui::BeginTabItem("Phổ phản ứng TCVN 9386")) {
            SpectrumParams& sp = specParams;
            ImGui::SetNextItemWidth(90);
            ImGui::InputDouble("agR (g)", &sp.agR, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(90);
            ImGui::InputDouble("γI", &sp.gammaI, 0, 0, "%.3g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(90);
            ImGui::InputDouble("q", &sp.q, 0, 0, "%.3g");
            const char* soils[] = {"Nền A", "Nền B", "Nền C", "Nền D", "Nền E"};
            ImGui::SetNextItemWidth(90);
            ImGui::Combo("Loại nền", &sp.soil, soils, 5);
            ImGui::SameLine();
            ImGui::RadioButton("Phương X", &sp.direction, 0);
            ImGui::SameLine();
            ImGui::RadioButton("Phương Y", &sp.direction, 1);
            ImGui::Checkbox("Tổ hợp CQC (mặc định SRSS)", &sp.cqc);
            double S, TB, TC, TD;
            soilParams(sp.soil, S, TB, TC, TD);
            ImGui::TextDisabled("S = %.2f, TB = %.2f s, TC = %.2f s, TD = %.1f s; ag = γI·agR·g = %.3f m/s²", S, TB, TC, TD,
                                sp.gammaI * sp.agR * 9.81);
            // đồ thị phổ
            {
                std::vector<double> T, Sd;
                for (int i = 0; i <= 200; ++i) {
                    T.push_back(4.0 * i / 200);
                    Sd.push_back(spectrumSd(sp, T.back()));
                }
                const ImVec2 p0 = ImGui::GetCursorScreenPos(), sz(ImGui::GetContentRegionAvail().x, 150);
                ImGui::Dummy(sz);
                ImDrawList* dl = ImGui::GetWindowDrawList();
                plot(dl, p0, sz, T, Sd, pal.posLine, pal.textDim, pal.text, "T (s)", "Sd (m/s²)");
                if (spec.ok)
                    for (const auto& mi : spec.modes) {
                        if (mi.T > 4) continue;
                        const float x = p0.x + 46 + (float)(mi.T / 4.0) * (sz.x - 52);
                        dl->AddLine(ImVec2(x, p0.y + 6), ImVec2(x, p0.y + sz.y - 18), IM_COL32(255, 160, 40, 160));
                    }
            }
            if (ImGui::Button("Tính phổ phản ứng")) {
                spec = spectrumAnalysis(model, dynOpt, sp);
                dynVersion = modelVersion;
                dynShow = spec.ok ? 1 : 0;
                if (spec.ok) {
                    view = ResultView::M;
                    char b[96];
                    std::snprintf(b, sizeof b, "Phổ phản ứng: lực cắt đáy %.2f kN, Σ khối lượng hữu hiệu %.1f%%", spec.baseShear,
                                  spec.sumRatio * 100);
                    echo(b);
                } else echo("Lỗi: " + spec.error);
            }
            if (spec.ok) {
                for (const auto& w : spec.warnings) ImGui::TextColored(ImVec4(0.9f, 0.6f, 0.1f, 1), "⚠ %s", w.c_str());
                ImGui::Text("Tổng khối lượng %.3f t; Σ khối lượng hữu hiệu %.1f%%; lực cắt đáy %.3f kN", spec.totalMass,
                            spec.sumRatio * 100, spec.baseShear);
                if (ImGui::BeginTable("##modes", 7, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY,
                                      ImVec2(0, 170))) {
                    ImGui::TableSetupScrollFreeze(0, 1);
                    for (const char* h : {"Dạng", "T (s)", "f (Hz)", "Γ", "mhh %", "Sd (m/s²)", "V (kN)"}) ImGui::TableSetupColumn(h);
                    ImGui::TableHeadersRow();
                    for (size_t i = 0; i < spec.modes.size(); ++i) {
                        const auto& mi = spec.modes[i];
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn();
                        if (ImGui::Selectable(std::to_string(i + 1).c_str(), dynShow == 2 && specMode == (int)i,
                                              ImGuiSelectableFlags_SpanAllColumns)) {
                            dynShow = 2;
                            specMode = (int)i;
                        }
                        ImGui::TableNextColumn();
                        ImGui::Text("%.4f", mi.T);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.3f", mi.f);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.3f", mi.gamma);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.1f", mi.ratio * 100);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.3f", mi.Sd);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.2f", mi.baseShear);
                    }
                    ImGui::EndTable();
                }
                showRadio("Hiển thị tổ hợp (giá trị tuyệt đối)", 1);
                ImGui::SameLine();
                showRadio("Không hiển thị", 0);
                ImGui::TextDisabled("Chọn một dòng trong bảng để xem nội lực của từng dạng.");
            } else if (!spec.error.empty()) {
                ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.25f, 1), "%s", spec.error.c_str());
            }
            ImGui::EndTabItem();
        }
        // ================================================= dao động cưỡng bức
        if (ImGui::BeginTabItem("Cưỡng bức điều hòa")) {
            ImGui::TextWrapped("Lực P(t) = P₀·sin θt, P₀ phân bố theo bộ tải chọn; tính biên độ dao động ổn định (không cản).");
            setCombo("Bộ tải P₀##h", harmSet, model);
            ImGui::SetNextItemWidth(100);
            ImGui::InputDouble("Tần số lực kích thích (Hz)", &harmFreq, 0, 0, "%.4g");
            if (ImGui::Button("Tính dao động cưỡng bức")) {
                harm = harmonicAnalysis(model, dynOpt, harmSet, 2 * kPi * harmFreq);
                dynVersion = modelVersion;
                dynShow = harm.ok ? 3 : 0;
                if (!harm.ok) echo("Lỗi: " + harm.error);
            }
            if (harm.ok) {
                ImGui::Text("θ = %.3f rad/s; hệ số động μ = max|u_đ| / max|u_tĩnh| = %.4f", harm.theta, harm.dynFactor);
                std::string ws;
                char b[48];
                for (double w : harm.omegas) {
                    std::snprintf(b, sizeof b, "%.2f  ", w);
                    ws += b;
                }
                ImGui::TextWrapped("Tần số riêng ω (rad/s): %s", ws.c_str());
                for (double w : harm.omegas)
                    if (std::abs(w - harm.theta) < 0.1 * w) {
                        ImGui::TextColored(ImVec4(0.9f, 0.4f, 0.2f, 1), "⚠ Gần cộng hưởng (θ ≈ ω = %.2f rad/s)", w);
                        break;
                    }
                showRadio("Hiển thị biên độ", 3);
                ImGui::SameLine();
                showRadio("Không hiển thị##h", 0);
            } else if (!harm.error.empty()) {
                ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.25f, 1), "%s", harm.error.c_str());
            }
            ImGui::EndTabItem();
        }
        // ================================================= theo thời gian
        if (ImGui::BeginTabItem("Theo thời gian (Newmark)")) {
            TimeHistoryParams& p = thParams;
            const char* fns[] = {"Tải bậc thang", "Xung chữ nhật", "Xung nửa sin", "Tải điều hòa", "Gia tốc nền dạng sin",
                                 "Gia tốc nền tự nhập (t, a)"};
            int fn = (int)p.fn;
            ImGui::SetNextItemWidth(220);
            if (ImGui::Combo("Dạng tác động", &fn, fns, 6)) p.fn = (TimeFn)fn;
            const bool ground = p.fn == TimeFn::GroundSine || p.fn == TimeFn::GroundRecord;
            if (!ground) setCombo("Bộ tải (phân bố)##t", p.loadSet, model);
            ImGui::SetNextItemWidth(90);
            ImGui::InputDouble(ground ? "Biên độ / hệ số gia tốc (m/s²)" : "Hệ số nhân tải", &p.amp, 0, 0, "%.4g");
            if (p.fn == TimeFn::Rect || p.fn == TimeFn::HalfSine) {
                ImGui::SetNextItemWidth(90);
                ImGui::InputDouble("Thời gian tác dụng td (s)", &p.td, 0, 0, "%.4g");
            }
            if (p.fn == TimeFn::Harmonic || p.fn == TimeFn::GroundSine) {
                ImGui::SetNextItemWidth(90);
                ImGui::InputDouble("Tần số (Hz)", &p.freq, 0, 0, "%.4g");
            }
            if (ground) {
                ImGui::RadioButton("Theo X##g", &p.groundDir, 0);
                ImGui::SameLine();
                ImGui::RadioButton("Theo Y##g", &p.groundDir, 1);
            }
            if (p.fn == TimeFn::GroundRecord) {
                ImGui::TextDisabled("Mỗi dòng: t (s), a (m/s²)");
                if (ImGui::InputTextMultiline("##rec", &thRecordText, ImVec2(-FLT_MIN, 80))) {
                    p.record.clear();
                    std::istringstream in(thRecordText);
                    std::string line;
                    while (std::getline(in, line)) {
                        for (char& c : line)
                            if (c == ',' || c == ';' || c == '\t') c = ' ';
                        std::istringstream ls(line);
                        double t, a;
                        if (ls >> t >> a) p.record.push_back({t, a});
                    }
                }
            }
            ImGui::SetNextItemWidth(90);
            ImGui::InputDouble("Δt (s)", &p.dt, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(90);
            ImGui::InputDouble("Thời gian (s)", &p.tEnd, 0, 0, "%.4g");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70);
            ImGui::InputDouble("ζ", &p.zeta, 0, 0, "%.3g");
            if (ImGui::Button("Tính theo thời gian")) {
                th = timeHistory(model, dynOpt, p);
                dynVersion = modelVersion;
                thStep = 0;
                thSnapStep = -1;
                if (th.ok) {
                    // nút có chuyển vị lớn nhất theo phương X
                    double best = -1;
                    for (size_t d = 0; d < th.maxAbs.size() && d / 3 < model.nodes.size(); ++d)
                        if (d % 3 != RZ && th.maxAbs[d] > best) {
                            best = th.maxAbs[d];
                            thNode = (int)(d / 3);
                            thDof = (int)(d % 3);
                        }
                    dynShow = 4;
                } else echo("Lỗi: " + th.error);
            }
            if (th.ok) {
                for (const auto& w : th.warnings) ImGui::TextColored(ImVec4(0.9f, 0.6f, 0.1f, 1), "⚠ %s", w.c_str());
                ImGui::TextDisabled("Cản Rayleigh: C = %.4g·M + %.4g·K; ω1 = %.3f rad/s", th.alpha, th.betaR, th.omega1);
                ImGui::SetNextItemWidth(80);
                int n1 = thNode + 1;
                if (ImGui::InputInt("Nút##th", &n1)) thNode = std::clamp(n1 - 1, 0, (int)model.nodes.size() - 1);
                const char* ds[] = {"ux", "uy", "φz"};
                for (int i = 0; i < 3; ++i) {
                    ImGui::SameLine();
                    if (ImGui::RadioButton(ds[i], thDof == i)) thDof = i;
                }
                const std::vector<double> h = historyOf(th, thNode, thDof);
                std::vector<double> hmm(h.size());
                double peak = 0;
                int kPeak = 0;
                for (size_t k = 0; k < h.size(); ++k) {
                    hmm[k] = thDof == RZ ? h[k] : h[k] * 1e3;
                    if (std::abs(h[k]) > peak) {
                        peak = std::abs(h[k]);
                        kPeak = (int)k;
                    }
                }
                const ImVec2 p0 = ImGui::GetCursorScreenPos(), sz(ImGui::GetContentRegionAvail().x, 170);
                ImGui::Dummy(sz);
                plot(ImGui::GetWindowDrawList(), p0, sz, th.t, hmm, pal.posLine, pal.textDim, pal.text, "t (s)",
                     thDof == RZ ? "φ (rad)" : "u (mm)", th.t[std::clamp(thStep, 0, (int)th.t.size() - 1)]);
                ImGui::Text("|u|max = %.4f %s tại t = %.3f s", thDof == RZ ? peak : peak * 1e3, thDof == RZ ? "rad" : "mm", th.t[kPeak]);
                ImGui::SetNextItemWidth(-140);
                ImGui::SliderInt("##step", &thStep, 0, (int)th.t.size() - 1, "bước %d");
                ImGui::SameLine();
                if (ImGui::Button("Đỉnh")) thStep = kPeak;
                ImGui::SameLine();
                ImGui::Checkbox("Chạy", &thPlay);
                if (thPlay) {
                    thStep += std::max(1, (int)th.t.size() / 300);
                    if (thStep >= (int)th.t.size()) thStep = 0;
                }
                ImGui::Text("t = %.3f s", th.t[std::clamp(thStep, 0, (int)th.t.size() - 1)]);
                if (thSnapStep != thStep && dynShow == 4) {
                    thSnap = timeSnapshot(model, th, thStep);
                    thSnapStep = thStep;
                }
                showRadio("Hiển thị nội lực tại thời điểm t", 4);
                ImGui::SameLine();
                showRadio("Không hiển thị##t", 0);
            } else if (!th.error.empty()) {
                ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.25f, 1), "%s", th.error.c_str());
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
}
