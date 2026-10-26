// FEM2D — chương trình phần tử hữu hạn cho kết cấu phẳng. Giao diện: SDL3 + Dear ImGui.
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <vector>
#include <string>

#include "App.h"
#include "core/ForceMethod.h"
#include "core/Model.h"
#include "imgui.h"
#include "imgui_impl_sdl3.h"
#include "imgui_impl_sdlrenderer3.h"

static void setupStyle(float scale) {
    ImGuiStyle& s = ImGui::GetStyle();
    s.FramePadding = ImVec2(6, 3);
    s.ItemSpacing = ImVec2(7, 5);
    s.WindowBorderSize = 1;
    s.ScaleAllSizes(scale);
    s.FontScaleDpi = scale;
}

static void loadFonts() {
    ImGuiIO& io = ImGui::GetIO();
    // Phông chữ có đủ dấu tiếng Việt, chữ Hy Lạp (σ, τ, φ) và chữ Nga (ГОСТ).
    const char* candidates[] = {
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
        "/Library/Fonts/Arial Unicode.ttf",
        "C:\\Windows\\Fonts\\segoeui.ttf",
        "C:\\Windows\\Fonts\\arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/TTF/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf",
        "/usr/share/fonts/noto/NotoSans-Regular.ttf",
    };
    ImFont* font = nullptr;
    for (const char* path : candidates) {
        if (FILE* f = std::fopen(path, "rb")) {
            std::fclose(f);
            font = io.Fonts->AddFontFromFileTTF(path, 15.0f);
            if (font) break;
        }
    }
    if (!font) {
        io.Fonts->AddFontDefault();
        return;
    }
    // phông dự phòng cho ký hiệu thiếu (chỉ số trên ⁴, mũi tên...)
    const char* fallbacks[] = {"/System/Library/Fonts/Supplemental/Arial Unicode.ttf", "/Library/Fonts/Arial Unicode.ttf",
                               "C:\\Windows\\Fonts\\seguisym.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"};
    for (const char* path : fallbacks) {
        if (FILE* f = std::fopen(path, "rb")) {
            std::fclose(f);
            ImFontConfig cfg;
            cfg.MergeMode = true;
            if (io.Fonts->AddFontFromFileTTF(path, 15.0f, &cfg)) break;
        }
    }
    ImGui::GetStyle().FontSizeBase = 15.0f;
}

int main(int argc, char** argv) {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        std::printf("SDL_Init lỗi: %s\n", SDL_GetError());
        return 1;
    }
    const float scale = SDL_GetDisplayContentScale(SDL_GetPrimaryDisplay());
    SDL_WindowFlags flags = SDL_WINDOW_RESIZABLE | SDL_WINDOW_HIDDEN | SDL_WINDOW_HIGH_PIXEL_DENSITY;
    SDL_Window* window = SDL_CreateWindow("FEM2D", (int)(1440 * scale), (int)(900 * scale), flags);
    if (!window) {
        std::printf("SDL_CreateWindow lỗi: %s\n", SDL_GetError());
        return 1;
    }
    SDL_Renderer* renderer = SDL_CreateRenderer(window, nullptr);
    if (!renderer) {
        std::printf("SDL_CreateRenderer lỗi: %s\n", SDL_GetError());
        return 1;
    }
    SDL_SetRenderVSync(renderer, 1);
    SDL_SetWindowPosition(window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
    SDL_ShowWindow(window);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    // không bật điều hướng bàn phím: phím gõ được chuyển vào dòng lệnh (kiểu AutoCAD)
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.IniFilename = "fem2d_layout_v2.ini"; // v2: bố cục kiểu AutoCAD
    setupStyle(scale);
    ImGui_ImplSDL3_InitForSDLRenderer(window, renderer);
    ImGui_ImplSDLRenderer3_Init(renderer);
    loadFonts();

    App app;
    app.window = window;
    app.applyTheme();
    // --shot <ảnh.bmp> <số ví dụ> <chế độ xem 0..6> [đại lượng màu 0..8] [cờ]: chụp màn hình rồi thoát
    std::string shotPath;
    int shotFrames = 0;
    // --smoke: thao tác như người dùng qua dòng lệnh và chuột, kiểm tra kết quả rồi thoát
    bool smoke = false;
    int smokeFrame = 0;
    std::vector<std::function<void()>> script;
    ImVec2 smokeMouse(-1, -1); // vị trí chuột giả lập, gửi lại mỗi khung hình
    if (argc > 2 && std::string(argv[1]) == "--shot") {
        shotPath = argv[2];
        app.loadExample(argc > 3 ? std::atoi(argv[3]) : 0);
        if (argc > 4) app.view = (ResultView)std::atoi(argv[4]);
        if (argc > 5) app.contour = (ContourField)std::atoi(argv[5]);
        if (argc > 6) {
            // cờ: e = biểu đồ bao; h = thẻ và cửa sổ Hướng dẫn; s<k> = bộ kết quả thứ k; m<thanh>:<x> = đah M tại tiết diện; q<thanh>:<x> = đah Q
            const std::string f = argv[6];
            app.showEnvelope = f.find('e') != std::string::npos;
            if (f.find('k') != std::string::npos) app.startCommand("CTHH");
            if (f.find('d') != std::string::npos) app.showSectionDesigner = true;
            if (f.find('f') != std::string::npos) {
                app.showForceMethod = true;
                app.fm = fem::forceMethod(app.model, 0);
                app.fmVersion = app.modelVersion;
                app.fmShow = 0;
                app.view = ResultView::M;
            }
            if (f.find('r') != std::string::npos) app.trainShow = true;
            if (f.find('s') != std::string::npos) app.resultSet = std::atoi(f.c_str() + f.find('s') + 1); // s<k>: bộ kết quả k
            if (f.find('g') != std::string::npos) app.showDesign = true;
            if (f.find('h') != std::string::npos) {
                app.ribbonTab = 5;
                app.openGuide(GuideTopic::Start);
                app.openGuide(GuideTopic::Influence);
            }
            if (f.find('p') != std::string::npos) {
                app.showDynamics = true;
                app.dynOpt.modes = 6;
                app.spec = fem::spectrumAnalysis(app.model, app.dynOpt, app.specParams);
                app.dynVersion = app.modelVersion;
                app.dynShow = 1;
                app.view = ResultView::M;
            }
            if (f.find('t') != std::string::npos) {
                app.selElems = {1};
                app.showSectionStress = true;
            }
            if (f[0] == 'u') {
                // u: mặt cắt đứng, nhà đặt trên mặt đất, gió không dừng
                fem::Model m;
                m.title = "Gió quanh nhà — mặt cắt đứng";
                m.materials.push_back(fem::Material{});
                m.sections.push_back(fem::Section{});
                m.cases.push_back(fem::LoadCase("TH1"));
                fem::FlowDomain& fl = m.flow;
                fl.active = true;
                fl.section = 1;
                fl.ax = 0, fl.ay = 0, fl.bx = 0, fl.by = 30;
                fl.side = -1;
                fl.length = 80;
                fl.U = 10;
                fl.medium = 1;
                fl.rho = 1.204;
                fl.mu = 1.825e-5; // không khí thật: Re ≈ 2·10⁷, mô hình rối k-ε (tự động)
                fl.hydrostatic = false;
                fl.transient = false;
                const double P[5][2] = {{20, 0}, {30, 0}, {30, 12}, {25, 16}, {20, 12}};
                for (auto& q : P) m.nodes.push_back({q[0], q[1]});
                for (int i = 0; i < 5; ++i) {
                    fem::Element e;
                    e.type = fem::ElemType::Frame;
                    e.n = {i, (i + 1) % 5, -1, -1};
                    m.elements.push_back(e);
                }
                app.model = m;
                app.changed();
                app.fitRequested = true;
                app.view = ResultView::None;
                app.startCommand("TINH");
                const int steps = f.size() > 1 ? std::atoi(f.c_str() + 1) : 600;
                for (int i = 0; app.fluid && i < steps && !app.fluid->converged(); ++i) app.fluid->step();
                app.fluidRun = false;
                app.fluidField = (int)fem::FluidField::Speed;
            }
            if (f[0] == 'y' || f[0] == 'Y' || f[0] == 'Z') {
                // y<t>: mặt cắt đứng có mặt thoáng tự do — cột nước vỡ đập đổ vào khối công trình; Y<t>: đáy thấm nước;
                // Z<t>: đáy thoát tự do
                fem::Model m;
                m.title = "Vỡ đập — mặt thoáng tự do";
                m.materials.push_back(fem::Material{});
                m.sections.push_back(fem::Section{});
                m.cases.push_back(fem::LoadCase("TH1"));
                fem::FlowDomain& fl = m.flow;
                fl.active = true;
                fl.section = 1;
                fl.ax = 0, fl.ay = 0, fl.bx = 0, fl.by = 1.2;
                fl.side = -1;
                fl.length = 4;
                fl.U = 0;
                fl.freeSurface = true;
                fl.waterDepth = 0.8;
                fl.gate = 1.2;
                fl.h = 0.03;
                fl.hFar = 0.05;
                fl.dt = 0.01;
                fl.endTime = 5;
                fl.drainGround = f[0] == 'Y' || f[0] == 'Z';
                fl.drainMode = f[0] == 'Z' ? 1 : 0;
                fl.drainRate = 0.2;
                const double P[4][2] = {{2.6, -0.05}, {2.9, -0.05}, {2.9, 0.3}, {2.6, 0.3}};
                for (auto& q : P) m.nodes.push_back({q[0], q[1]});
                for (int i = 0; i < 4; ++i) {
                    fem::Element e;
                    e.type = fem::ElemType::Frame;
                    e.n = {i, (i + 1) % 4, -1, -1};
                    m.elements.push_back(e);
                }
                app.model = m;
                app.changed();
                app.fitRequested = true;
                app.view = ResultView::None;
                app.startCommand("TINH");
                const double tEnd = f.size() > 1 ? std::atof(f.c_str() + 1) : 1.0;
                while (app.fluid && app.fluid->time() < tEnd && app.fluid->step()) {
                }
                app.fluidRun = false;
            }
            if (f[0] == 'w' || f[0] == 'W') {
                // w: trụ tròn vẽ bằng thanh trong miền chảy khai báo, tính đến hội tụ; W: chỉ khai báo, chưa tính
                fem::Model m;
                m.title = "Dòng chảy quanh trụ";
                m.materials.push_back(fem::Material{});
                m.sections.push_back(fem::Section{});
                m.cases.push_back(fem::LoadCase("TH1"));
                m.flow.active = true;
                m.flow.ax = 0, m.flow.ay = 0, m.flow.bx = 0, m.flow.by = 0.41;
                m.flow.length = 2.2;
                m.flow.U = 0.2;
                m.flow.parabolic = true;
                m.flow.medium = 4;
                m.flow.rho = 1;
                m.flow.mu = 1e-3;
                const double pi = 3.14159265358979323846;
                for (int i = 0; i < 48; ++i)
                    m.nodes.push_back({0.2 + 0.05 * std::cos(2 * pi * i / 48), 0.2 + 0.05 * std::sin(2 * pi * i / 48)});
                for (int i = 0; i < 48; ++i) {
                    fem::Element e;
                    e.type = fem::ElemType::Frame;
                    e.n = {i, (i + 1) % 48, -1, -1};
                    m.elements.push_back(e);
                }
                app.model = m;
                app.changed();
                app.fitRequested = true;
                app.view = ResultView::None;
                if (f.find('q') != std::string::npos) app.editFlow(); // mở cửa sổ khai báo
                if (f[0] == 'w') {
                    app.startCommand("TINH");
                    for (int i = 0; app.fluid && i < 3000 && !app.fluid->converged(); ++i) app.fluid->step();
                    app.fluidRun = false;
                }
            }
            if (f[0] == 'c') {
                // c<k>[:<số bước>]: bài toán chất lỏng mẫu k, chạy trước số bước (mặc định đến khi hội tụ / hết thời gian)
                app.showFluid = true;
                app.fluidApplyPreset(std::atoi(f.c_str() + 1));
                const size_t colon = f.find(':');
                const int steps = colon != std::string::npos ? std::atoi(f.c_str() + colon + 1) : 100000;
                for (int i = 0; app.fluid && i < steps; ++i) {
                    if (!app.fluid->params().transient && app.fluid->converged()) break;
                    if (app.fluid->params().transient && app.fluid->time() >= app.fluidEndTime) break;
                    if (!app.fluid->step()) break;
                }
                app.fluidRun = false;
                if (argc > 5) app.fluidField = std::atoi(argv[5]);
            }
            if (f.find('l') != std::string::npos) {
                app.darkTheme = false;
                app.applyTheme();
            }
            if (f[0] == 'b' || f[0] == 'v') {
                app.eigOpt.kind = f[0] == 'b' ? fem::EigenKind::Buckling : fem::EigenKind::Vibration;
                app.eigMode = 0;
                app.animateMode = false;
                app.runEigen();
                app.animateMode = false;
                if (f.size() > 1) app.eigMode = std::atoi(f.c_str() + 1) - 1;
                app.showEigenDialog = true;
            }
            if (f[0] == 'm' || f[0] == 'q') {
                app.ilTarget.kind = fem::IlKind::Force;
                app.ilTarget.qty = f[0] == 'm' ? fem::BarQty::M : fem::BarQty::Q;
                app.ilTarget.elem = std::atoi(f.c_str() + 1) - 1;
                app.ilTarget.x = std::atof(f.c_str() + f.find(':') + 1);
            }
        }
        io.IniFilename = nullptr;
    } else if (argc > 1 && std::string(argv[1]) == "--smoke") {
        smoke = true;
        io.IniFilename = nullptr;
        io.ConfigInputTrickleEventQueue = false;
        app.newModel();
        app.cancelCommand();
        auto wait = [&](int n) {
            for (int i = 0; i < n; ++i) script.push_back([] {});
        };
        auto at = [&app](double x, double y) { return app.cam.toScreen(x, y); };
        auto type = [&](const char* t) {
            script.push_back([&io, t] { io.AddInputCharactersUTF8(t); });
            script.push_back([&io] { io.AddKeyEvent(ImGuiKey_Enter, true); });
            script.push_back([&io] { io.AddKeyEvent(ImGuiKey_Enter, false); });
        };
        auto click = [&](double x, double y) {
            script.push_back([&smokeMouse, at, x, y] { smokeMouse = at(x, y); });
            script.push_back([&io] { io.AddMouseButtonEvent(0, true); });
            script.push_back([&io] { io.AddMouseButtonEvent(0, false); });
        };
        auto drag = [&](double x0, double y0, double x1, double y1) {
            script.push_back([&smokeMouse, at, x0, y0] { smokeMouse = at(x0, y0); });
            script.push_back([&io] { io.AddMouseButtonEvent(0, true); });
            script.push_back([&smokeMouse, at, x1, y1] { smokeMouse = at(x1, y1); });
            script.push_back([&io] { io.AddMouseButtonEvent(0, false); });
        };
        auto check = [&](const char* what, std::function<bool()> f) {
            script.push_back([what, f] { std::printf("smoke: %-58s %s\n", what, f() ? "OK" : "LỖI"); });
        };
        wait(4);
        script.push_back([&app] {
            app.fitRequested = false;
            app.fitFrames = 0;
            app.cam.cx = 3;
            app.cam.cy = 2;
            app.cam.scale = 50;
        });
        wait(2);
        type("L");
        click(0, 0);
        click(0, 4);
        click(6, 4);
        click(6, 0);
        type(""); // Enter: kết thúc lệnh
        check("vẽ khung bằng lệnh L: 4 nút, 3 thanh", [&app] { return app.model.nodes.size() == 4 && app.model.elements.size() == 3; });
        type("L");
        click(3.02, 4.03); // OSNAP trung điểm dầm
        type("@0,-3");     // tọa độ tương đối
        type("");
        check("nối vào trung điểm dầm (cắt dầm): 6 nút, 5 thanh", [&app] { return app.model.nodes.size() == 6 && app.model.elements.size() == 5; });
        drag(-0.5, -0.5, 6.5, 0.5); // chọn cửa sổ trái→phải: 2 nút chân cột
        check("chọn cửa sổ: 2 nút, 0 thanh", [&app] { return app.selNodes.size() == 2 && app.selElems.empty(); });
        type("NGAM");
        check("NGAM với tập chọn sẵn: 2 gối ngàm", [&app] { return app.model.supports.size() == 2 && app.model.supports[0].fix[2]; });
        click(1.5, 4);
        click(4.5, 4);
        type("TAIPB");
        type("-10");
        type(""); // hướng mặc định GY
        check("TAIPB -10 trên 2 nửa dầm", [&app] { return app.model.cases[0].elem.size() == 2; });
        drag(5.8, 3.7, 0.2, 4.3); // cắt ngang phải→trái: 2 nửa dầm + cột giữa + nút giữa
        check("chọn cắt ngang: 3 thanh", [&app] { return app.selElems.size() == 3; });
        script.push_back([&app] {
            app.selNodes.clear();
            app.selElems.clear();
            app.analyzeNow();
        });
        check("tính: ΣRy = 60 kN", [&app] { return app.result.ok && std::abs(app.result.sets[0].sumRy - 60) < 1e-9; });
        type("Z");
        type("E");
        script.push_back([&app] { app.undo(); });
        check("hoàn tác bỏ tải", [&app] { return app.model.cases[0].elem.empty(); });
        type("CHATLONG");
        check("CHATLONG mở cửa sổ dòng chảy", [&app] { return app.showFluid; });
        script.push_back([&app] { app.fluidApplyPreset(0); });
        wait(20);
        check("dòng chảy: chạy trực tiếp qua các khung hình", [&app] {
            return app.fluidShown() && app.fluidRun && app.fluid->iteration() > 5 && !app.fluid->diverged();
        });
        script.push_back([&app] { app.undo(); });
        check("hoàn tác: trả lại khung, ẩn kết quả dòng chảy", [&app] { return app.model.elements.size() == 5 && !app.fluidShown(); });
        // khai báo dòng chảy (Trang chủ › Dòng chảy), vẽ vật cản bằng L, bấm Tính
        script.push_back([&app] {
            app.newModel();
            app.cancelCommand();
            app.showFluid = false;
            app.fitRequested = false;
            app.fitFrames = 0;
            app.cam.cx = 3;
            app.cam.cy = 1;
            app.cam.scale = 60;
        });
        wait(2);
        type("DONGCHAY");
        click(0, 0);
        click(0, 2);
        click(4, 1); // phía +x, dài 4 m
        wait(3);
        script.push_back([&io] { io.AddKeyEvent(ImGuiKey_Enter, true); });
        script.push_back([&io] { io.AddKeyEvent(ImGuiKey_Enter, false); });
        wait(2);
        check("DONGCHAY + popup: miền chảy rộng 2 m, dài 4 m, chảy theo +x", [&app] {
            double nx, ny;
            app.model.flow.direction(nx, ny);
            return app.model.flow.active && std::abs(app.model.flow.width() - 2) < 1e-9 &&
                   std::abs(app.model.flow.length - 4) < 1e-9 && nx > 0.999;
        });
        type("L");
        click(2, 0.5);
        click(2.5, 0.5);
        click(2.5, 1.5);
        click(2, 1.5);
        type("C");
        check("vẽ vật cản khép kín bằng L ... C: 4 thanh", [&app] { return app.model.elements.size() == 4; });
        type("TINH");
        wait(20);
        check("Tính (TINH): lưới quanh vật cản, giải dòng chảy trực tiếp", [&app] {
            return app.fluidShown() && app.fluid->mesh().sideFaces[(int)fem::FluidSide::Inner] >= 4 &&
                   app.fluid->iteration() > 3 && !app.fluid->diverged();
        });
        // mặt cắt đứng: công trình đặt trên mặt đất
        script.push_back([&app] {
            app.newModel();
            app.cancelCommand();
            app.showFluid = false;
        });
        wait(2);
        type("DONGCHAY");
        type("D");
        click(0, 0);
        click(0, 2);
        click(4, 1);
        wait(3);
        script.push_back([&io] { io.AddKeyEvent(ImGuiKey_Enter, true); });
        script.push_back([&io] { io.AddKeyEvent(ImGuiKey_Enter, false); });
        wait(2);
        check("DONGCHAY D: mặt cắt đứng", [&app] { return app.model.flow.active && app.model.flow.section == 1; });
        type("L");
        click(2, 0);
        click(2.5, 0);
        click(2.5, 1);
        click(2, 1);
        type("C");
        type("TINH");
        wait(20);
        check("mặt cắt đứng: công trình trên mặt đất, có áp suất thủy tĩnh", [&app] {
            return app.fluidShown() && app.fluid->params().hydroG > 0 &&
                   app.fluid->mesh().sideFaces[(int)fem::FluidSide::Inner] > 0 && !app.fluid->diverged();
        });
        // mặt thoáng tự do (VOF) và chảy rối k-ε trên cùng mô hình
        script.push_back([&app] {
            app.snapshot();
            app.model.flow.freeSurface = true;
            app.model.flow.waterDepth = 1.2;
            app.model.flow.U = 0.5;
            app.model.flow.turbulence = 2;
            app.changed();
        });
        type("TINH");
        wait(15);
        check("mặt thoáng tự do + k-ε: VOF chạy trực tiếp, nước bảo toàn", [&app] {
            return app.fluidShown() && app.fluid->params().vof && app.fluid->params().turbulence == fem::TurbModel::KEpsilon &&
                   app.fluid->iteration() > 2 && !app.fluid->diverged() && app.fluid->waterVolume() > 0;
        });
        script.push_back([&app] {
            if (std::getenv("FEM_SMOKE_LOG"))
                for (const auto& l : app.log) std::printf("  | %s\n", l.c_str());
            app.quitNow = true;
        });
    } else if (argc > 1) {
        app.openFile(argv[1]);
    } else {
        app.loadExample(3);
    }

    bool done = false;
    while (!done) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            ImGui_ImplSDL3_ProcessEvent(&ev);
            const bool closeReq = ev.type == SDL_EVENT_QUIT ||
                                  (ev.type == SDL_EVENT_WINDOW_CLOSE_REQUESTED && ev.window.windowID == SDL_GetWindowID(window));
            if (closeReq && (smoke || !shotPath.empty() || app.confirmDiscard(Pending::Quit))) done = true;
            if (ev.type == SDL_EVENT_DROP_FILE && ev.drop.data) app.confirmDiscard(Pending::OpenPath, -1, ev.drop.data);
        }
        if (SDL_GetWindowFlags(window) & SDL_WINDOW_MINIMIZED) {
            SDL_Delay(10);
            continue;
        }
        app.processDialogs();
        if (app.quitNow) done = true;
        SDL_SetWindowTitle(window, app.windowTitle().c_str());

        // luôn nhận ký tự để gõ lệnh ở bất cứ đâu (backend ImGui tắt text input khi không có ô nhập)
        if (!SDL_TextInputActive(window)) SDL_StartTextInput(window);
        ImGui_ImplSDLRenderer3_NewFrame();
        ImGui_ImplSDL3_NewFrame();
        // giả lập sau backend để vị trí chuột thật không ghi đè sự kiện giả lập
        if (smoke && smokeMouse.x >= 0) io.AddMousePosEvent(smokeMouse.x, smokeMouse.y);
        if (smoke && smokeFrame < (int)script.size()) script[smokeFrame++]();
        if (smoke && smokeMouse.x >= 0) io.AddMousePosEvent(smokeMouse.x, smokeMouse.y);
        ImGui::NewFrame();
        app.frame();
        ImGui::Render();
        SDL_SetRenderScale(renderer, io.DisplayFramebufferScale.x, io.DisplayFramebufferScale.y);
        SDL_SetRenderDrawColor(renderer, 240, 240, 242, 255);
        SDL_RenderClear(renderer);
        ImGui_ImplSDLRenderer3_RenderDrawData(ImGui::GetDrawData(), renderer);
        if (!shotPath.empty() && ++shotFrames == 6) {
            if (SDL_Surface* surf = SDL_RenderReadPixels(renderer, nullptr)) {
                SDL_SaveBMP(surf, shotPath.c_str());
                SDL_DestroySurface(surf);
            }
            done = true;
        }
        SDL_RenderPresent(renderer);
    }

    ImGui_ImplSDLRenderer3_Shutdown();
    ImGui_ImplSDL3_Shutdown();
    ImGui::DestroyContext();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
