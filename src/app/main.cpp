// FEM2D — chương trình phần tử hữu hạn cho kết cấu phẳng. Giao diện: SDL3 + Dear ImGui.
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#include <cstdio>
#include <cstdlib>
#include <string>

#include "App.h"
#include "core/Model.h"
#include "imgui.h"
#include "imgui_impl_sdl3.h"
#include "imgui_impl_sdlrenderer3.h"

static void setupStyle(float scale) {
    ImGui::StyleColorsLight();
    ImGuiStyle& s = ImGui::GetStyle();
    s.WindowRounding = 4;
    s.FrameRounding = 3;
    s.GrabRounding = 3;
    s.TabRounding = 3;
    s.FramePadding = ImVec2(6, 3);
    s.ItemSpacing = ImVec2(7, 5);
    s.WindowBorderSize = 1;
    s.FrameBorderSize = 1;
    ImVec4* c = s.Colors;
    c[ImGuiCol_WindowBg] = ImVec4(0.96f, 0.96f, 0.97f, 1);
    c[ImGuiCol_ChildBg] = ImVec4(1, 1, 1, 0);
    c[ImGuiCol_Header] = ImVec4(0.26f, 0.47f, 0.80f, 0.22f);
    c[ImGuiCol_HeaderHovered] = ImVec4(0.26f, 0.47f, 0.80f, 0.35f);
    c[ImGuiCol_HeaderActive] = ImVec4(0.26f, 0.47f, 0.80f, 0.50f);
    c[ImGuiCol_Button] = ImVec4(0.88f, 0.90f, 0.93f, 1);
    c[ImGuiCol_ButtonHovered] = ImVec4(0.26f, 0.47f, 0.80f, 0.35f);
    c[ImGuiCol_ButtonActive] = ImVec4(0.26f, 0.47f, 0.80f, 0.60f);
    c[ImGuiCol_TitleBgActive] = ImVec4(0.80f, 0.85f, 0.93f, 1);
    c[ImGuiCol_Border] = ImVec4(0, 0, 0, 0.18f);
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
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard | ImGuiConfigFlags_DockingEnable;
    io.IniFilename = "fem2d_layout.ini";
    setupStyle(scale);
    ImGui_ImplSDL3_InitForSDLRenderer(window, renderer);
    ImGui_ImplSDLRenderer3_Init(renderer);
    loadFonts();

    App app;
    app.window = window;
    // --shot <ảnh.bmp> <số ví dụ> <chế độ xem 0..6> [đại lượng màu 0..8]: chụp màn hình rồi thoát
    std::string shotPath;
    int shotFrames = 0;
    // --smoke: giả lập chuột vẽ khung cổng bằng công cụ "Thanh khung", gán gối, rồi in kết quả
    bool smoke = false;
    int smokeFrame = 0;
    if (argc > 2 && std::string(argv[1]) == "--shot") {
        shotPath = argv[2];
        app.loadExample(argc > 3 ? std::atoi(argv[3]) : 0);
        if (argc > 4) app.view = (ResultView)std::atoi(argv[4]);
        if (argc > 5) app.contour = (ContourField)std::atoi(argv[5]);
        if (argc > 6) {
            // cờ: e = biểu đồ bao; m<thanh>:<x> = đah M tại tiết diện; q<thanh>:<x> = đah Q
            const std::string f = argv[6];
            app.showEnvelope = f.find('e') != std::string::npos;
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
        app.newModel();
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

        if (smoke) {
            const double pts[4][2] = {{0, 0}, {0, 4}, {6, 4}, {6, 0}};
            const int k = smokeFrame / 4, ph = smokeFrame % 4;
            if (smokeFrame == 3) app.fitRequested = false;
            if (k >= 1 && k <= 4) {
                const ImVec2 p = app.cam.toScreen(pts[k - 1][0], pts[k - 1][1]);
                if (ph == 0) io.AddMousePosEvent(p.x, p.y);
                if (ph == 1) io.AddMouseButtonEvent(0, true);
                if (ph == 2) io.AddMouseButtonEvent(0, false);
            } else if (k == 5) {
                if (ph == 0) {
                    app.tool = Tool::Select;
                    const ImVec2 a = app.cam.toScreen(-0.5, -0.5);
                    io.AddMousePosEvent(a.x, a.y);
                }
                if (ph == 1) io.AddMouseButtonEvent(0, true);
                if (ph == 2) {
                    const ImVec2 b = app.cam.toScreen(6.5, 0.5); // kéo khung chọn 2 nút chân cột
                    io.AddMousePosEvent(b.x, b.y);
                }
                if (ph == 3) io.AddMouseButtonEvent(0, false);
            } else if (k == 6 && ph == 0) {
                app.snapshot();
                for (int n : app.selNodes) app.model.supportFor(n).fix = {true, true, true};
                fem::ElemLoad l;
                l.elem = 1;
                l.q1 = l.q2 = -10;
                app.model.cases[0].elem.push_back(l);
                app.changed();
            } else if (k == 7) {
                app.analyzeNow();
                std::printf("smoke: nút %zu, phần tử %zu, chọn %zu nút, gối %zu, tính %s, ΣRy = %.3f\n",
                            app.model.nodes.size(), app.model.elements.size(), app.selNodes.size(),
                            app.model.supports.size(), app.result.ok ? "OK" : app.result.error.c_str(),
                            app.result.ok ? app.result.sets[0].sumRy : 0.0);
                app.undo();
                app.analyzeNow();
                std::printf("smoke: sau hoàn tác: gối %zu, tải %zu\n", app.model.supports.size(),
                            app.model.cases[0].elem.size());
                done = true;
            }
            ++smokeFrame;
        }
        ImGui_ImplSDLRenderer3_NewFrame();
        ImGui_ImplSDL3_NewFrame();
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
