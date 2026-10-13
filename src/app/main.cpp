// FEM2D — chương trình phần tử hữu hạn cho kết cấu phẳng. Giao diện: SDL3 + Dear ImGui.
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#include <cstdio>
#include <cstdlib>
#include <functional>
#include <vector>
#include <string>

#include "App.h"
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
    // --shot <ảnh.bmp> <số ví dụ> <chế độ xem 0..6> [đại lượng màu 0..8]: chụp màn hình rồi thoát
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
            // cờ: e = biểu đồ bao; m<thanh>:<x> = đah M tại tiết diện; q<thanh>:<x> = đah Q
            const std::string f = argv[6];
            app.showEnvelope = f.find('e') != std::string::npos;
            if (f.find('k') != std::string::npos) app.startCommand("CTHH");
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
