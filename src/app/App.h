// Trạng thái ứng dụng và khai báo các thành phần giao diện.
#pragma once
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <vector>

#include "core/Analysis.h"
#include "core/Design.h"
#include "core/Dxf.h"
#include "core/Dynamics.h"
#include "core/Eigen.h"
#include "core/ForceMethod.h"
#include "core/Kinematics.h"
#include "core/Model.h"
#include "core/PostProcess.h"
#include "Commands.h"
#include "imgui.h"

struct SDL_Window;

enum class ResultView { None, Deformed, N, Q, M, Sigma, Contour, Influence, Mode };
enum class Pending { None, New, Open, OpenPath, Example, Quit };
enum class ContourField { Ux, Uy, Umag, Sx, Sy, Txy, S1, S2, VonMises, Sz };

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
    ImU32 axisX, axisY, crosshair, snap, grip, gripHover, gripHot, windowSel, crossingSel;
};

// Màu giao diện ribbon / thanh trạng thái theo chủ đề
struct UiColors {
    ImU32 tabRow, ribbon, panelTitle, separator, text, textDim, hover, active, accent, cmdBg, statusBg;
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
    double gridStep = 0.5;
    bool snapGrid = true;
    int curMat = 0, curSec = 0;

    // ---- tùy chọn hiển thị
    bool showNodeIds = false, showElemIds = false, showLoads = true, showSupports = true, showValues = true;
    bool showReactions = true, showUndeformed = true, showDeformedWithDiagram = false, showGrid = true;
    bool darkTheme = true; // giao diện tối (kiểu AutoCAD) / sáng
    Palette pal;
    UiColors ui;
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

    // ---- ổn định & dao động riêng
    fem::EigenOptions eigOpt;
    fem::EigenResult eig;
    int eigVersion = -1, eigMode = 0;
    bool animateMode = true;
    bool showEigenDialog = false;

    // ---- cấu tạo hình học
    fem::KinematicResult kin;
    int kinVersion = -1;
    bool showKinematics = false;
    bool showSectionDesigner = false, showSectionStress = false;

    // ---- lời giải theo giáo trình
    fem::ForceMethodResult fm;
    int fmVersion = -1, fmSet = 0, fmShow = -1;
    fem::MohrResult mohr;
    int mohrNode = 0, mohrDof = 1;
    fem::MatrixMethodResult mm;
    int mmVersion = -1;
    bool showForceMethod = false, showMatrixMethod = false, showArchGen = false;
    // xếp đoàn tải trên đường ảnh hưởng
    int trainPreset = 1;
    fem::LoadTrain customTrain{"Tự nhập", {{100, 0}, {100, 4}}, 0};
    std::string customTrainText = "100@0 100@4";
    bool trainShow = false;

    // ---- động lực học
    fem::EigenOptions dynOpt;
    fem::SpectrumParams specParams;
    fem::SpectrumResult spec;
    fem::HarmonicResult harm;
    fem::TimeHistoryParams thParams;
    fem::TimeHistoryResult th;
    fem::LoadSetResult thSnap;
    std::string thRecordText = "0 0\n0.1 1\n0.2 -1\n0.3 0.5\n0.4 0\n";
    int dynVersion = -1, dynShow = 0, specMode = 0, harmSet = 0, thStep = 0, thSnapStep = -1, thNode = 0, thDof = 0;
    double harmFreq = 2.0;
    bool showDynamics = false, thPlay = false;

    // ---- kiểm tra thiết kế
    bool showDesign = false, showPileGen = false, showUtil = true, useBucklingMu = false, deflBeamsOnly = true;
    int designSetMode = 0;
    double deflN = 250;
    fem::SteelParams steelP;
    fem::RcParams rcP;
    std::vector<double> util; // hệ số sử dụng từng phần tử để tô màu (−1: không có)

    // ---- dòng lệnh kiểu AutoCAD
    std::unique_ptr<Command> cmd;
    std::string cmdLine, lastCommand;
    std::vector<std::string> typedHistory;
    int typedPos = -1;
    bool osnap = true, ortho = false, dynInput = true, cmdExpanded = false;
    SnapPoint cursor;
    bool cursorInCanvas = false;
    double cursorWX = 0, cursorWY = 0;
    int ribbonTab = 0;
    // giá trị mặc định của các lệnh (nhớ lần nhập trước)
    double defNodal[3] = {0, -10, 0};
    double defQ1 = -10, defQ2 = -10, defP = -20, defM = 10;
    int defDir = 3, defPDir = 3, defDiv = 2;
    bool defProj = false;

    // ---- hỏi lưu thay đổi
    Pending pending = Pending::None;
    int pendingExample = -1;
    std::string pendingPath;
    bool openSavePrompt = false;
    bool quitNow = false;

    // ---- cửa sổ phụ
    bool showModelPanel = true, showPropsPanel = true, showResultPanel = true, showLogPanel = false;
    bool showMeshDialog = false, showFrameGen = false, showTrussGen = false, showAbout = false, showDemo = false;
    bool resetLayout = true;

    std::vector<std::string> log;
    SDL_Window* window = nullptr;
    std::mutex dialogMutex;
    std::string dialogOpenPath, dialogSavePath, dialogReportPath, dialogDxfIn, dialogDxfOut, dialogHtml;
    // nhập DXF
    bool showDxfImport = false, dxfReplace = false;
    std::string dxfText, dxfPath;
    std::vector<std::string> dxfLayerList;
    int dxfLayer = -1, dxfUnit = 0, dxfType = 0;

    // ---- thao tác
    void newModel();
    void loadExample(int i);
    bool openFile(const std::string& path);
    bool saveFile(const std::string& path);
    void requestOpen();
    void requestSave(bool saveAs);
    void requestReport();
    void requestDxfIn();
    void requestDxfOut();
    void requestHtml();
    void dxfImportDialog();
    void processDialogs();

    void snapshot();          // lưu trạng thái hiện tại vào ngăn hoàn tác
    void changed();           // mô hình đã thay đổi → tính lại
    void trackEdit();         // gọi ngay sau một ô nhập liệu để hỗ trợ hoàn tác
    void undo();
    void redo();
    void analyzeNow(bool explicitRun = false); // explicitRun: ghi lỗi vào lịch sử lệnh
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

    // ---- lệnh (Ribbon.cpp)
    void startCommand(const std::string& name);
    void submitCommandLine();
    void cancelCommand();
    void commandKeys();
    void echo(const std::string& s);
    std::string currentPrompt();
    int nodeForSnap(const SnapPoint& p); // nút tại điểm: dùng nút có sẵn, cắt thanh hoặc tạo nút mới
    std::set<int> selectionNodes() const;
    void applyTheme();

    // ---- giao diện (định nghĩa trong Ribbon.cpp / Panels.cpp / View.cpp)
    void ribbon();
    void statusBar();
    void commandLine();
    void appMenu();
    void frame();
    void dockLayout(unsigned int dockId);
    void modelPanel();
    void propsPanel();
    void resultPanel();
    void logPanel();
    void meshDialog();
    void frameGenDialog();
    void trussGenDialog();
    void aboutDialog();
    void eigenDialog();
    void kinematicsDialog();
    void sectionDesigner();
    void sectionStressDialog();
    void forceMethodDialog();
    void matrixMethodDialog();
    void archDialog();
    void dynamicsDialog();
    void designDialog();
    void pileDialog();
    std::vector<int> designSets() const;
    const fem::LoadSetResult* overrideShown(std::string* title = nullptr);
    const fem::LoadSetResult* forceMethodShown(std::string* title = nullptr);
    const fem::LoadTrain& currentTrain() const;
    void showMechanism();
    void runEigen();
    void viewport();
    void shortcuts();
};

Palette makePalette(bool dark);
