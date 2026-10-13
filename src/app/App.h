// Trạng thái ứng dụng và khai báo các thành phần giao diện.
#pragma once
#include <mutex>
#include <set>
#include <string>
#include <vector>

#include "core/Analysis.h"
#include "core/Model.h"
#include "core/PostProcess.h"
#include "imgui.h"

struct SDL_Window;

enum class Tool { Select, AddNode, AddFrame, AddTruss };
enum class ResultView { None, Deformed, N, Q, M, Sigma, Contour, Influence };
enum class Pending { None, New, Open, OpenPath, Example, Quit };
enum class ContourField { Ux, Uy, Umag, Sx, Sy, Txy, S1, S2, VonMises };

struct Camera {
    double cx = 5, cy = 2;   // tâm nhìn (m)
    double scale = 60;       // px / m
    ImVec2 c0{0, 0}, c1{0, 0}; // khung canvas trên màn hình
    ImVec2 toScreen(double x, double y) const {
        return ImVec2((float)((c0.x + c1.x) * 0.5 + (x - cx) * scale), (float)((c0.y + c1.y) * 0.5 - (y - cy) * scale));
    }
    void toWorld(ImVec2 p, double& x, double& y) const {
        x = cx + (p.x - (c0.x + c1.x) * 0.5) / scale;
        y = cy - (p.y - (c0.y + c1.y) * 0.5) / scale;
    }
};

struct Palette {
    ImU32 bg, grid, gridMajor, axis, bar, truss, plane, planeEdge, node, text, textDim, select, hover, support, load,
        reaction, deformed, posFill, posLine, negFill, negLine, hinge;
};

struct App {
    // ---- mô hình & lịch sử
    fem::Model model;
    std::string filePath;
    bool dirty = false;
    std::vector<fem::Model> undoStack, redoStack;
    fem::Model editSnapshot;
    bool editSnapshotValid = false;

    // ---- kết quả
    fem::AnalysisResult result;
    bool resultValid = false;
    bool needAnalyze = false;
    bool autoAnalyze = true;
    int modelVersion = 0, resultVersion = 0;
    int activeCase = 0; // trường hợp tải đang soạn
    int resultSet = 0;  // bộ kết quả đang xem (trường hợp tải hoặc tổ hợp)

    // ---- chọn
    std::set<int> selNodes, selElems;
    int hoverNode = -1, hoverElem = -1;

    // ---- nhìn & công cụ
    Camera cam;
    bool fitRequested = true;
    int fitFrames = 0; // lặp lại vài khung hình vì kích thước canvas ổn định sau khi dock
    Tool tool = Tool::Select;
    int pendingNode = -1;
    bool pendingValid = false;
    double pendingX = 0, pendingY = 0;
    double gridStep = 0.5;
    bool snapGrid = true;
    int curMat = 0, curSec = 0;

    // ---- tùy chọn hiển thị
    bool showNodeIds = false, showElemIds = false, showLoads = true, showSupports = true, showValues = true;
    bool showReactions = true, showUndeformed = true, showDeformedWithDiagram = false, showGrid = true;
    bool darkCanvas = false;
    ResultView view = ResultView::M;
    ContourField contour = ContourField::VonMises;
    bool contourAveraged = true, contourOnDeformed = false;
    float deformFactor = 1.0f, diagramFactor = 1.0f;
    double deformAuto = 1.0; // hệ số phóng đại tự động (tính khi vẽ)

    // ---- biểu đồ bao
    bool showEnvelope = false;
    std::vector<fem::BarEnvelope> env;
    std::vector<int> envSets;
    int envVersion = -1;

    // ---- đường ảnh hưởng
    fem::IlTarget ilTarget;
    std::vector<int> ilPath; // rỗng: tự động chọn
    fem::InfluenceLine il;
    int ilVersion = -1;
    bool ilDirty = true;

    // ---- hỏi lưu thay đổi
    Pending pending = Pending::None;
    int pendingExample = -1;
    std::string pendingPath;
    bool openSavePrompt = false;
    bool quitNow = false;

    // ---- cửa sổ phụ
    bool showModelPanel = true, showPropsPanel = true, showResultPanel = true, showLogPanel = true;
    bool showMeshDialog = false, showFrameGen = false, showTrussGen = false, showAbout = false, showDemo = false;
    bool resetLayout = true;

    std::vector<std::string> log;
    SDL_Window* window = nullptr;
    std::mutex dialogMutex;
    std::string dialogOpenPath, dialogSavePath, dialogReportPath;

    // ---- thao tác
    void newModel();
    void loadExample(int i);
    bool openFile(const std::string& path);
    bool saveFile(const std::string& path);
    void requestOpen();
    void requestSave(bool saveAs);
    void requestReport();
    void processDialogs();

    void snapshot();          // lưu trạng thái hiện tại vào ngăn hoàn tác
    void changed();           // mô hình đã thay đổi → tính lại
    void trackEdit();         // gọi ngay sau một ô nhập liệu để hỗ trợ hoàn tác
    void undo();
    void redo();
    void analyzeNow();
    void deleteSelection();
    void clampIndices();
    void addLog(const std::string& s);
    bool confirmDiscard(Pending p, int example = -1, const std::string& path = {});
    void runPending();
    void savePrompt();
    void updateEnvelope();
    void updateInfluence();
    std::vector<int> effectiveIlPath() const;
    std::string windowTitle() const;
    const fem::LoadSetResult* currentSet() const;

    // ---- giao diện (định nghĩa trong Panels.cpp / View.cpp)
    void frame();
    void mainMenu();
    void dockLayout(unsigned int dockId);
    void modelPanel();
    void propsPanel();
    void resultPanel();
    void logPanel();
    void meshDialog();
    void frameGenDialog();
    void trussGenDialog();
    void aboutDialog();
    void viewport();
    void shortcuts();
};

Palette makePalette(bool dark);
