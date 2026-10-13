#include "App.h"

#include <SDL3/SDL.h>

#include <algorithm>
#include <cstdio>
#include <fstream>

#include "core/Examples.h"
#include "core/Library.h"
#include "core/ModelIO.h"
#include "imgui_internal.h"

using namespace fem;

namespace {

constexpr size_t kMaxUndo = 60;

const SDL_DialogFileFilter kFemFilter[] = {{"Mô hình FEM2D (*.fem)", "fem"}, {"Tất cả", "*"}};
const SDL_DialogFileFilter kTxtFilter[] = {{"Báo cáo (*.txt)", "txt"}};

std::string fileName(const std::string& p) {
    const size_t k = p.find_last_of("/\\");
    return k == std::string::npos ? p : p.substr(k + 1);
}

} // namespace

Palette makePalette(bool dark) {
    Palette p;
    if (!dark) {
        p.bg = IM_COL32(255, 255, 255, 255);
        p.grid = IM_COL32(232, 235, 240, 255);
        p.gridMajor = IM_COL32(210, 215, 224, 255);
        p.axis = IM_COL32(160, 168, 180, 255);
        p.bar = IM_COL32(30, 40, 60, 255);
        p.truss = IM_COL32(60, 90, 130, 255);
        p.plane = IM_COL32(190, 205, 225, 255);
        p.planeEdge = IM_COL32(110, 125, 150, 255);
        p.node = IM_COL32(20, 25, 35, 255);
        p.text = IM_COL32(20, 25, 35, 255);
        p.textDim = IM_COL32(110, 115, 125, 255);
        p.select = IM_COL32(255, 140, 0, 255);
        p.hover = IM_COL32(0, 150, 255, 255);
        p.support = IM_COL32(40, 120, 60, 255);
        p.load = IM_COL32(210, 40, 40, 255);
        p.reaction = IM_COL32(20, 140, 70, 255);
        p.deformed = IM_COL32(200, 30, 140, 255);
        p.posFill = IM_COL32(60, 120, 220, 70);
        p.posLine = IM_COL32(30, 80, 190, 255);
        p.negFill = IM_COL32(225, 70, 60, 70);
        p.negLine = IM_COL32(190, 40, 30, 255);
        p.hinge = IM_COL32(255, 255, 255, 255);
    } else {
        p.bg = IM_COL32(28, 30, 36, 255);
        p.grid = IM_COL32(40, 43, 50, 255);
        p.gridMajor = IM_COL32(54, 58, 68, 255);
        p.axis = IM_COL32(90, 96, 110, 255);
        p.bar = IM_COL32(225, 230, 240, 255);
        p.truss = IM_COL32(150, 190, 240, 255);
        p.plane = IM_COL32(60, 75, 100, 255);
        p.planeEdge = IM_COL32(120, 135, 160, 255);
        p.node = IM_COL32(235, 235, 240, 255);
        p.text = IM_COL32(235, 235, 240, 255);
        p.textDim = IM_COL32(150, 155, 165, 255);
        p.select = IM_COL32(255, 160, 40, 255);
        p.hover = IM_COL32(80, 190, 255, 255);
        p.support = IM_COL32(110, 210, 130, 255);
        p.load = IM_COL32(255, 100, 90, 255);
        p.reaction = IM_COL32(90, 220, 140, 255);
        p.deformed = IM_COL32(255, 90, 200, 255);
        p.posFill = IM_COL32(80, 140, 255, 80);
        p.posLine = IM_COL32(120, 170, 255, 255);
        p.negFill = IM_COL32(255, 90, 80, 80);
        p.negLine = IM_COL32(255, 120, 110, 255);
        p.hinge = IM_COL32(28, 30, 36, 255);
    }
    return p;
}

void App::addLog(const std::string& s) {
    log.push_back(s);
    if (log.size() > 500) log.erase(log.begin());
}

std::string App::windowTitle() const {
    std::string t = "FEM2D — ";
    t += filePath.empty() ? (model.title.empty() ? "Mô hình mới" : model.title) : fileName(filePath);
    if (dirty) t += " *";
    return t;
}

void App::clampIndices() {
    if (model.cases.empty()) model.cases.push_back(LoadCase("TH1"));
    activeCase = std::clamp(activeCase, 0, (int)model.cases.size() - 1);
    curMat = std::clamp(curMat, 0, std::max(0, (int)model.materials.size() - 1));
    curSec = std::clamp(curSec, 0, std::max(0, (int)model.sections.size() - 1));
    const int nsets = (int)(model.cases.size() + model.combos.size());
    resultSet = std::clamp(resultSet, 0, std::max(0, nsets - 1));
    for (auto it = selNodes.begin(); it != selNodes.end();)
        it = (*it >= (int)model.nodes.size()) ? selNodes.erase(it) : std::next(it);
    for (auto it = selElems.begin(); it != selElems.end();)
        it = (*it >= (int)model.elements.size()) ? selElems.erase(it) : std::next(it);
    if (pendingNode >= (int)model.nodes.size()) pendingNode = -1;
}

const LoadSetResult* App::currentSet() const {
    if (!resultValid || resultSet < 0 || resultSet >= (int)result.sets.size()) return nullptr;
    const LoadSetResult& s = result.sets[resultSet];
    // bảo vệ khi mô hình vừa thay đổi trong khung hình hiện tại
    if (s.U.size() != 3 * model.nodes.size() || s.bars.size() != model.elements.size()) return nullptr;
    return &s;
}

void App::analyzeNow() {
    clampIndices();
    const std::string prevError = result.error;
    result = analyze(model);
    resultValid = result.ok;
    ++resultVersion;
    needAnalyze = false;
    if (!result.ok && result.error != prevError) addLog("Lỗi: " + result.error);
    clampIndices();
}

void App::snapshot() {
    undoStack.push_back(model);
    if (undoStack.size() > kMaxUndo) undoStack.erase(undoStack.begin());
    redoStack.clear();
}

void App::changed() {
    dirty = true;
    ++modelVersion;
    needAnalyze = true;
    clampIndices();
}

void App::trackEdit() {
    if (ImGui::IsItemActivated()) {
        editSnapshot = model;
        editSnapshotValid = true;
    }
    if (ImGui::IsItemDeactivatedAfterEdit() && editSnapshotValid) {
        undoStack.push_back(std::move(editSnapshot));
        if (undoStack.size() > kMaxUndo) undoStack.erase(undoStack.begin());
        redoStack.clear();
        editSnapshotValid = false;
        changed();
    } else if (ImGui::IsItemDeactivated()) {
        editSnapshotValid = false;
    }
}

void App::undo() {
    if (undoStack.empty()) return;
    redoStack.push_back(std::move(model));
    model = std::move(undoStack.back());
    undoStack.pop_back();
    changed();
}

void App::redo() {
    if (redoStack.empty()) return;
    undoStack.push_back(std::move(model));
    model = std::move(redoStack.back());
    redoStack.pop_back();
    changed();
}

void App::deleteSelection() {
    if (selNodes.empty() && selElems.empty()) return;
    snapshot();
    model.removeElements(std::vector<int>(selElems.begin(), selElems.end()));
    model.removeNodes(std::vector<int>(selNodes.begin(), selNodes.end()));
    selNodes.clear();
    selElems.clear();
    pendingNode = -1;
    changed();
}

void App::newModel() {
    model = Model{};
    model.title = "Mô hình mới";
    model.materials.push_back(toMaterial(materialPresets()[0]));
    model.sections.push_back(toSection(profilePresets()[5]));
    model.cases.push_back(LoadCase("TH1"));
    filePath.clear();
    undoStack.clear();
    redoStack.clear();
    selNodes.clear();
    selElems.clear();
    dirty = false;
    activeCase = resultSet = 0;
    cam.cx = 5;
    cam.cy = 2;
    cam.scale = 60;
    fitRequested = false;
    tool = Tool::AddFrame;
    ++modelVersion;
    ilPath.clear();
    analyzeNow();
}

void App::loadExample(int i) {
    const auto& ex = examples();
    if (i < 0 || i >= (int)ex.size()) return;
    model = ex[i].build();
    filePath.clear();
    undoStack.clear();
    redoStack.clear();
    selNodes.clear();
    selElems.clear();
    dirty = false;
    activeCase = resultSet = 0;
    fitRequested = true;
    tool = Tool::Select;
    ++modelVersion;
    ilPath.clear();
    ilTarget = IlTarget{};
    view = model.hasPlane() ? ResultView::Contour : model.hasBars() && model.elements[0].type == ElemType::Truss ? ResultView::N
                                                                                                                  : ResultView::M;
    if (model.hasPlane()) contour = ContourField::Sx;
    analyzeNow();
    addLog(std::string("Đã mở ví dụ: ") + ex[i].name);
}

bool App::openFile(const std::string& path) {
    Model m;
    std::string err;
    if (!loadModel(path, m, err)) {
        addLog("Lỗi mở file: " + err);
        return false;
    }
    model = std::move(m);
    filePath = path;
    undoStack.clear();
    redoStack.clear();
    selNodes.clear();
    selElems.clear();
    dirty = false;
    activeCase = resultSet = 0;
    fitRequested = true;
    tool = Tool::Select;
    ++modelVersion;
    ilPath.clear();
    ilTarget = IlTarget{};
    analyzeNow();
    addLog("Đã mở " + path);
    return true;
}

bool App::saveFile(const std::string& path) {
    std::string err;
    if (!saveModel(path, model, err)) {
        addLog("Lỗi lưu file: " + err);
        return false;
    }
    filePath = path;
    dirty = false;
    addLog("Đã lưu " + path);
    return true;
}

void App::requestOpen() {
    SDL_ShowOpenFileDialog(
        [](void* ud, const char* const* list, int) {
            if (!list || !list[0]) return;
            App* a = (App*)ud;
            std::lock_guard<std::mutex> lk(a->dialogMutex);
            a->dialogOpenPath = list[0];
        },
        this, window, kFemFilter, 2, nullptr, false);
}

void App::requestSave(bool saveAs) {
    if (!saveAs && !filePath.empty()) {
        saveFile(filePath);
        return;
    }
    SDL_ShowSaveFileDialog(
        [](void* ud, const char* const* list, int) {
            if (!list || !list[0]) return;
            App* a = (App*)ud;
            std::lock_guard<std::mutex> lk(a->dialogMutex);
            a->dialogSavePath = list[0];
        },
        this, window, kFemFilter, 1, filePath.empty() ? "mo_hinh.fem" : filePath.c_str());
}

void App::requestReport() {
    SDL_ShowSaveFileDialog(
        [](void* ud, const char* const* list, int) {
            if (!list || !list[0]) return;
            App* a = (App*)ud;
            std::lock_guard<std::mutex> lk(a->dialogMutex);
            a->dialogReportPath = list[0];
        },
        this, window, kTxtFilter, 1, "ket_qua.txt");
}

void App::processDialogs() {
    std::string open, save, report;
    {
        std::lock_guard<std::mutex> lk(dialogMutex);
        open.swap(dialogOpenPath);
        save.swap(dialogSavePath);
        report.swap(dialogReportPath);
    }
    if (!open.empty()) openFile(open);
    if (!save.empty() && pending != Pending::None) {
        // lưu xong (từ hộp thoại "Lưu thay đổi?") thì tiếp tục thao tác đang chờ
        if (save.size() < 4 || save.substr(save.size() - 4) != ".fem") save += ".fem";
        if (saveFile(save)) runPending();
        else pending = Pending::None;
        save.clear();
    }
    if (!save.empty()) {
        if (save.size() < 4 || save.substr(save.size() - 4) != ".fem") save += ".fem";
        saveFile(save);
    }
    if (!report.empty()) {
        if (report.size() < 4 || report.substr(report.size() - 4) != ".txt") report += ".txt";
        if (needAnalyze) analyzeNow();
        std::ofstream f(report, std::ios::binary);
        if (f) {
            f << resultReport(model, result);
            addLog("Đã xuất báo cáo " + report);
        } else addLog("Không ghi được " + report);
    }
}

void App::shortcuts() {
    using namespace ImGui;
    const ImGuiInputFlags g = ImGuiInputFlags_RouteGlobal;
    if (Shortcut(ImGuiMod_Ctrl | ImGuiKey_Z, g)) undo();
    if (Shortcut(ImGuiMod_Ctrl | ImGuiKey_Y, g) || Shortcut(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_Z, g)) redo();
    if (Shortcut(ImGuiMod_Ctrl | ImGuiKey_S, g)) requestSave(false);
    if (Shortcut(ImGuiMod_Ctrl | ImGuiKey_O, g) && confirmDiscard(Pending::Open)) requestOpen();
    if (Shortcut(ImGuiMod_Ctrl | ImGuiKey_N, g) && confirmDiscard(Pending::New)) newModel();
    if (Shortcut(ImGuiKey_F5, g)) {
        analyzeNow();
        if (result.ok)
            addLog("Tính xong: " + std::to_string(result.numEquations) + " phương trình, " +
                   std::to_string((int)result.timeMs) + " ms");
    }
}

void App::frame() {
    if (needAnalyze && autoAnalyze) analyzeNow();
    shortcuts();
    mainMenu();
    const ImGuiID dock = ImGui::DockSpaceOverViewport(0, nullptr, ImGuiDockNodeFlags_None);
    dockLayout(dock);

    viewport();
    if (showModelPanel) modelPanel();
    if (showPropsPanel) propsPanel();
    if (showResultPanel) resultPanel();
    if (showLogPanel) logPanel();
    if (showMeshDialog) meshDialog();
    if (showFrameGen) frameGenDialog();
    if (showTrussGen) trussGenDialog();
    if (showAbout) aboutDialog();
    if (showDemo) ImGui::ShowDemoWindow(&showDemo);
    savePrompt();
    clampIndices();
}

void App::dockLayout(unsigned int dockId) {
    if (!resetLayout) return;
    resetLayout = false;
    ImGuiDockNode* node = ImGui::DockBuilderGetNode(dockId);
    static bool firstFrame = true;
    const bool force = !firstFrame;
    firstFrame = false;
    if (node && node->IsSplitNode() && !force) return; // đã có bố cục lưu trong file .ini
    ImGui::DockBuilderRemoveNode(dockId);
    ImGui::DockBuilderAddNode(dockId, ImGuiDockNodeFlags_DockSpace);
    ImGui::DockBuilderSetNodeSize(dockId, ImGui::GetMainViewport()->WorkSize);
    ImGuiID center = dockId;
    ImGuiID left = ImGui::DockBuilderSplitNode(center, ImGuiDir_Left, 0.24f, nullptr, &center);
    ImGuiID right = ImGui::DockBuilderSplitNode(center, ImGuiDir_Right, 0.30f, nullptr, &center);
    ImGuiID bottom = ImGui::DockBuilderSplitNode(center, ImGuiDir_Down, 0.16f, nullptr, &center);
    ImGuiID leftBottom = ImGui::DockBuilderSplitNode(left, ImGuiDir_Down, 0.45f, nullptr, &left);
    ImGui::DockBuilderDockWindow("Sơ đồ kết cấu", center);
    ImGui::DockBuilderDockWindow("Mô hình", left);
    ImGui::DockBuilderDockWindow("Thuộc tính", leftBottom);
    ImGui::DockBuilderDockWindow("Kết quả", right);
    ImGui::DockBuilderDockWindow("Nhật ký", bottom);
    ImGui::DockBuilderFinish(dockId);
}

void App::mainMenu() {
    if (!ImGui::BeginMainMenuBar()) return;
    if (ImGui::BeginMenu("Tệp")) {
        if (ImGui::MenuItem("Mới", "Ctrl+N") && confirmDiscard(Pending::New)) newModel();
        if (ImGui::MenuItem("Mở...", "Ctrl+O") && confirmDiscard(Pending::Open)) requestOpen();
        if (ImGui::MenuItem("Lưu", "Ctrl+S")) requestSave(false);
        if (ImGui::MenuItem("Lưu thành...")) requestSave(true);
        ImGui::Separator();
        if (ImGui::BeginMenu("Ví dụ mẫu")) {
            const auto& ex = examples();
            for (size_t i = 0; i < ex.size(); ++i) {
                if (ImGui::MenuItem(ex[i].name) && confirmDiscard(Pending::Example, (int)i)) loadExample((int)i);
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", ex[i].note);
            }
            ImGui::EndMenu();
        }
        ImGui::Separator();
        if (ImGui::MenuItem("Xuất báo cáo kết quả...")) requestReport();
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Sửa")) {
        if (ImGui::MenuItem("Hoàn tác", "Ctrl+Z", false, !undoStack.empty())) undo();
        if (ImGui::MenuItem("Làm lại", "Ctrl+Y", false, !redoStack.empty())) redo();
        ImGui::Separator();
        if (ImGui::MenuItem("Chọn tất cả", "Ctrl+A")) {
            for (int i = 0; i < (int)model.nodes.size(); ++i) selNodes.insert(i);
            for (int i = 0; i < (int)model.elements.size(); ++i) selElems.insert(i);
        }
        if (ImGui::MenuItem("Xóa đối tượng chọn", "Delete", false, !selNodes.empty() || !selElems.empty()))
            deleteSelection();
        if (ImGui::MenuItem("Xóa nút thừa")) {
            snapshot();
            model.removeUnusedNodes();
            changed();
        }
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Tạo")) {
        if (ImGui::MenuItem("Khung / dầm liên tục...")) showFrameGen = true;
        if (ImGui::MenuItem("Dàn...")) showTrussGen = true;
        if (ImGui::MenuItem("Lưới phần tử phẳng...")) showMeshDialog = true;
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Tính toán")) {
        if (ImGui::MenuItem("Tính ngay", "F5")) {
            analyzeNow();
            if (result.ok) addLog("Tính xong: " + std::to_string(result.numEquations) + " phương trình");
        }
        ImGui::MenuItem("Tự động tính khi thay đổi", nullptr, &autoAnalyze);
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Xem")) {
        ImGui::MenuItem("Mô hình", nullptr, &showModelPanel);
        ImGui::MenuItem("Thuộc tính", nullptr, &showPropsPanel);
        ImGui::MenuItem("Kết quả", nullptr, &showResultPanel);
        ImGui::MenuItem("Nhật ký", nullptr, &showLogPanel);
        ImGui::Separator();
        ImGui::MenuItem("Số hiệu nút", nullptr, &showNodeIds);
        ImGui::MenuItem("Số hiệu phần tử", nullptr, &showElemIds);
        ImGui::MenuItem("Tải trọng", nullptr, &showLoads);
        ImGui::MenuItem("Gối tựa", nullptr, &showSupports);
        ImGui::MenuItem("Lưới", nullptr, &showGrid);
        ImGui::MenuItem("Nền tối", nullptr, &darkCanvas);
        ImGui::Separator();
        if (ImGui::MenuItem("Vừa khung nhìn", "F")) fitRequested = true;
        if (ImGui::MenuItem("Khôi phục bố cục")) resetLayout = true;
        ImGui::MenuItem("ImGui demo", nullptr, &showDemo);
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Trợ giúp")) {
        if (ImGui::MenuItem("Giới thiệu & quy ước")) showAbout = true;
        ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
}

// ---------------------------------------------------------------- hỏi lưu thay đổi

bool App::confirmDiscard(Pending p, int example, const std::string& path) {
    pending = p;
    pendingExample = example;
    pendingPath = path;
    if (!dirty) {
        runPending();
        return p == Pending::New || p == Pending::Open; // thao tác do nơi gọi tự thực hiện
    }
    openSavePrompt = true;
    return false;
}

void App::runPending() {
    const Pending p = pending;
    pending = Pending::None;
    switch (p) {
    case Pending::None: break;
    case Pending::New:
    case Pending::Open: break; // được thực hiện tại nơi gọi khi không cần hỏi
    case Pending::OpenPath: openFile(pendingPath); break;
    case Pending::Example: loadExample(pendingExample); break;
    case Pending::Quit: quitNow = true; break;
    }
}

void App::savePrompt() {
    if (openSavePrompt) {
        ImGui::OpenPopup("Lưu thay đổi?");
        openSavePrompt = false;
    }
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (!ImGui::BeginPopupModal("Lưu thay đổi?", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) return;
    ImGui::Text("Mô hình \"%s\" có thay đổi chưa được lưu.", model.title.c_str());
    ImGui::Spacing();
    auto proceed = [&]() {
        const Pending p = pending;
        dirty = false;
        if (p == Pending::New) {
            pending = Pending::None;
            newModel();
        } else if (p == Pending::Open) {
            pending = Pending::None;
            requestOpen();
        } else {
            runPending();
        }
    };
    if (ImGui::Button("Lưu", ImVec2(110, 0))) {
        ImGui::CloseCurrentPopup();
        if (!filePath.empty()) {
            if (saveFile(filePath)) proceed();
            else pending = Pending::None;
        } else {
            requestSave(true); // tiếp tục sau khi chọn được tên file (processDialogs)
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Không lưu", ImVec2(110, 0))) {
        ImGui::CloseCurrentPopup();
        proceed();
    }
    ImGui::SameLine();
    if (ImGui::Button("Hủy", ImVec2(110, 0)) || ImGui::IsKeyPressed(ImGuiKey_Escape)) {
        pending = Pending::None;
        ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
}

// ---------------------------------------------------------------- biểu đồ bao, đường ảnh hưởng

void App::updateEnvelope() {
    if (envVersion == resultVersion) return;
    envVersion = resultVersion;
    env.clear();
    envSets.clear();
    if (!result.ok) return;
    // bao theo các tổ hợp nếu có, nếu không thì theo các trường hợp tải
    const int nc = (int)model.cases.size(), ns = (int)result.sets.size();
    for (int i = (ns > nc ? nc : 0); i < ns; ++i) envSets.push_back(i);
    env = envelope(model, result, envSets);
}

std::vector<int> App::effectiveIlPath() const {
    std::vector<int> path;
    for (int e : ilPath)
        if (e >= 0 && e < (int)model.elements.size() && model.elements[e].isBar()) path.push_back(e);
    if (path.empty()) {
        // tự động: các thanh nằm ngang ở cao độ thấp nhất (dầm, cánh dưới dàn, dầm tầng 1)
        double ymin = 1e300;
        std::vector<int> horiz;
        for (int e = 0; e < (int)model.elements.size(); ++e) {
            const Element& el = model.elements[e];
            if (!el.isBar()) continue;
            const Node& a = model.nodes[el.n[0]];
            const Node& b = model.nodes[el.n[1]];
            const double L = std::hypot(b.x - a.x, b.y - a.y);
            if (L <= 0 || std::abs(b.y - a.y) > 1e-6 * L) continue;
            horiz.push_back(e);
            ymin = std::min(ymin, a.y);
        }
        for (int e : horiz)
            if (std::abs(model.nodes[model.elements[e].n[0]].y - ymin) < 1e-9 * std::max(1.0, std::abs(ymin)))
                path.push_back(e);
    }
    std::sort(path.begin(), path.end(), [&](int a, int b) {
        auto mid = [&](int e) {
            const Element& el = model.elements[e];
            return model.nodes[el.n[0]].x + model.nodes[el.n[1]].x;
        };
        return mid(a) < mid(b);
    });
    return path;
}

void App::updateInfluence() {
    if (!ilDirty && ilVersion == modelVersion) return;
    ilDirty = false;
    ilVersion = modelVersion;
    il = influenceLine(model, effectiveIlPath(), ilTarget, 20);
}

