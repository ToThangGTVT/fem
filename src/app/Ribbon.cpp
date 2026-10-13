// Giao diện kiểu AutoCAD: chủ đề sáng/tối, ribbon, nút ứng dụng, thanh trạng thái, dòng lệnh và xử lý lệnh.
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>

#include "App.h"
#include "Commands.h"
#include "core/Examples.h"
#include "core/Mesher.h"
#include "imgui_internal.h"

using namespace fem;

namespace {

constexpr float kPi = 3.14159265358979f;

ImU32 rgb(int r, int g, int b, int a = 255) { return IM_COL32(r, g, b, a); }
ImVec4 v4(int r, int g, int b, float a = 1.0f) { return ImVec4(r / 255.0f, g / 255.0f, b / 255.0f, a); }

std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

void appendUtf8(std::string& s, unsigned c) {
    if (c < 0x80) s += (char)c;
    else if (c < 0x800) {
        s += (char)(0xC0 | (c >> 6));
        s += (char)(0x80 | (c & 0x3F));
    } else if (c < 0x10000) {
        s += (char)(0xE0 | (c >> 12));
        s += (char)(0x80 | ((c >> 6) & 0x3F));
        s += (char)(0x80 | (c & 0x3F));
    } else {
        s += (char)(0xF0 | (c >> 18));
        s += (char)(0x80 | ((c >> 12) & 0x3F));
        s += (char)(0x80 | ((c >> 6) & 0x3F));
        s += (char)(0x80 | (c & 0x3F));
    }
}

void popUtf8(std::string& s) {
    while (!s.empty()) {
        const unsigned char c = (unsigned char)s.back();
        s.pop_back();
        if ((c & 0xC0) != 0x80) break;
    }
}

// ============================================================================ biểu tượng vẽ bằng vector
enum class Icon {
    Line, Truss, Node, Mesh, FrameGen, TrussGen, Erase, Move, Copy, Rotate, Mirror, Divide, Merge, Fix, Pin, Roller,
    RollerX, Free, NLoad, DLoad, PLoad, MLoad, Temp, Settle, SelfW, Case, Combo, Run, Buckle, Modal, Influence, Schema,
    Deformed, DiagN, DiagQ, DiagM, Sigma, Contour, Envelope, Report, ZoomExt, Dist, Id, List, Props, New, Open, Save,
    Undo, Redo, Sun, Moon, Help
};

struct IconColors {
    ImU32 fg, blue, red, green, orange, dim;
};

IconColors iconColors(bool dark) {
    if (dark) return {rgb(225, 230, 238), rgb(80, 165, 255), rgb(255, 95, 85), rgb(90, 210, 120), rgb(255, 180, 60), rgb(140, 150, 165)};
    return {rgb(45, 52, 64), rgb(0, 110, 215), rgb(215, 50, 40), rgb(30, 150, 60), rgb(225, 130, 0), rgb(140, 148, 160)};
}

void arrowHead(ImDrawList* dl, ImVec2 tip, ImVec2 dir, float sz, ImU32 c) {
    const float L = std::sqrt(dir.x * dir.x + dir.y * dir.y);
    if (L <= 0) return;
    const ImVec2 d(dir.x / L, dir.y / L), n(-d.y, d.x);
    dl->AddTriangleFilled(tip, ImVec2(tip.x - d.x * sz + n.x * sz * 0.5f, tip.y - d.y * sz + n.y * sz * 0.5f),
                          ImVec2(tip.x - d.x * sz - n.x * sz * 0.5f, tip.y - d.y * sz - n.y * sz * 0.5f), c);
}

void drawIcon(ImDrawList* dl, Icon ic, ImVec2 c, float s, const IconColors& k) {
    const float h = s * 0.5f;
    auto P = [&](float x, float y) { return ImVec2(c.x + x * h, c.y + y * h); };
    const float th = std::max(1.2f, s / 18.0f);
    auto line = [&](float x0, float y0, float x1, float y1, ImU32 col, float t = 0) {
        dl->AddLine(P(x0, y0), P(x1, y1), col, t > 0 ? t : th);
    };
    auto dot = [&](float x, float y, ImU32 col, float r = 0.13f) { dl->AddCircleFilled(P(x, y), r * h, col); };
    auto arrow = [&](float x0, float y0, float x1, float y1, ImU32 col) {
        line(x0, y0, x1, y1, col);
        arrowHead(dl, P(x1, y1), ImVec2(x1 - x0, y1 - y0), 0.32f * h, col);
    };
    auto hatch = [&](float x0, float x1, float y, ImU32 col) {
        line(x0, y, x1, y, col);
        for (float x = x0 + 0.05f; x <= x1; x += 0.3f) line(x, y, x - 0.22f, y + 0.25f, col, th * 0.7f);
    };
    auto tri = [&](float y0, ImU32 col) {
        dl->AddTriangle(P(0, y0), P(-0.45f, y0 + 0.65f), P(0.45f, y0 + 0.65f), col, th);
    };
    auto text = [&](const char* t, float x, float y, ImU32 col, float scale) {
        ImFont* f = ImGui::GetFont();
        const float fs = s * scale;
        const ImVec2 ts = f->CalcTextSizeA(fs, FLT_MAX, 0, t);
        dl->AddText(f, fs, ImVec2(P(x, y).x - ts.x / 2, P(x, y).y - ts.y / 2), col, t);
    };
    auto diagram = [&](ImU32 col, int shape) {
        line(-0.9f, 0.1f, 0.9f, 0.1f, k.fg);
        std::vector<ImVec2> pts;
        for (int i = 0; i <= 12; ++i) {
            const float x = -0.9f + 1.8f * i / 12;
            float v = 0;
            if (shape == 0) v = -0.55f;                                  // N: hằng
            else if (shape == 1) v = (x < 0 ? -0.5f : 0.5f);             // Q: bậc
            else v = 0.75f * (1 - (x / 0.9f) * (x / 0.9f)) - 0.05f;      // M: parabol (vẽ xuống, thớ căng)
            pts.push_back(P(x, 0.1f + v));
        }
        for (size_t i = 0; i + 1 < pts.size(); ++i) {
            dl->AddQuadFilled(P(-0.9f + 1.8f * i / 12, 0.1f), pts[i], pts[i + 1], P(-0.9f + 1.8f * (i + 1) / 12, 0.1f),
                              (col & 0x00FFFFFF) | 0x50000000);
            dl->AddLine(pts[i], pts[i + 1], col, th);
        }
    };
    switch (ic) {
    case Icon::Line:
        line(-0.75f, 0.7f, 0.75f, -0.7f, k.fg);
        dot(-0.75f, 0.7f, k.blue, 0.17f);
        dot(0.75f, -0.7f, k.blue, 0.17f);
        break;
    case Icon::Truss:
        line(-0.9f, 0.6f, 0.9f, 0.6f, k.fg);
        line(-0.9f, 0.6f, 0, -0.6f, k.fg);
        line(0.9f, 0.6f, 0, -0.6f, k.fg);
        line(0, -0.6f, 0, 0.6f, k.fg);
        for (auto [x, y] : {std::pair{-0.9f, 0.6f}, {0.9f, 0.6f}, {0.f, -0.6f}, {0.f, 0.6f}}) {
            dl->AddCircleFilled(P(x, y), 0.14f * h, ImGui::GetColorU32(ImGuiCol_WindowBg));
            dl->AddCircle(P(x, y), 0.14f * h, k.blue, 0, th);
        }
        break;
    case Icon::Node:
        line(-0.7f, 0, 0.7f, 0, k.dim);
        line(0, -0.7f, 0, 0.7f, k.dim);
        dot(0, 0, k.blue, 0.25f);
        break;
    case Icon::Mesh:
        for (int i = 0; i <= 3; ++i) {
            const float t = -0.8f + 1.6f * i / 3;
            line(t, -0.8f, t, 0.8f, k.fg, th * 0.8f);
            line(-0.8f, t, 0.8f, t, k.fg, th * 0.8f);
        }
        dl->AddRectFilled(P(-0.27f, -0.27f), P(0.27f, 0.27f), (k.blue & 0x00FFFFFF) | 0x80000000);
        break;
    case Icon::FrameGen:
        line(-0.75f, 0.85f, -0.75f, -0.75f, k.fg);
        line(0.75f, 0.85f, 0.75f, -0.75f, k.fg);
        line(0, 0.85f, 0, -0.75f, k.fg);
        line(-0.75f, -0.75f, 0.75f, -0.75f, k.blue);
        line(-0.75f, 0.05f, 0.75f, 0.05f, k.blue);
        break;
    case Icon::TrussGen:
        line(-0.95f, 0.5f, 0.95f, 0.5f, k.fg);
        line(-0.5f, -0.4f, 0.5f, -0.4f, k.fg);
        line(-0.95f, 0.5f, -0.5f, -0.4f, k.fg);
        line(-0.5f, -0.4f, 0, 0.5f, k.fg);
        line(0, 0.5f, 0.5f, -0.4f, k.fg);
        line(0.5f, -0.4f, 0.95f, 0.5f, k.fg);
        break;
    case Icon::Erase:
        line(-0.8f, 0.5f, 0.5f, -0.7f, k.dim);
        line(0.15f, 0.2f, 0.8f, 0.85f, k.red, th * 1.5f);
        line(0.8f, 0.2f, 0.15f, 0.85f, k.red, th * 1.5f);
        break;
    case Icon::Move:
        arrow(0, 0, 0.85f, 0, k.fg);
        arrow(0, 0, -0.85f, 0, k.fg);
        arrow(0, 0, 0, 0.85f, k.fg);
        arrow(0, 0, 0, -0.85f, k.fg);
        dot(0, 0, k.blue, 0.18f);
        break;
    case Icon::Copy:
        dl->AddRect(P(-0.8f, -0.3f), P(0.3f, 0.8f), k.dim, 0, 0, th);
        dl->AddRect(P(-0.3f, -0.8f), P(0.8f, 0.3f), k.blue, 0, 0, th);
        break;
    case Icon::Rotate: {
        dl->PathArcTo(P(0, 0), 0.65f * h, -0.4f * kPi, 1.25f * kPi, 20);
        dl->PathStroke(k.fg, 0, th);
        const ImVec2 tip = P(0.65f * std::cos(-0.4f * kPi), 0.65f * std::sin(-0.4f * kPi));
        arrowHead(dl, tip, ImVec2(1, 0.4f), 0.38f * h, k.fg);
        dot(0, 0, k.blue, 0.16f);
        break;
    }
    case Icon::Mirror:
        for (float y = -0.9f; y < 0.9f; y += 0.35f) line(0, y, 0, y + 0.18f, k.dim);
        dl->AddTriangleFilled(P(-0.2f, -0.6f), P(-0.85f, 0.6f), P(-0.2f, 0.6f), k.blue);
        dl->AddTriangle(P(0.2f, -0.6f), P(0.85f, 0.6f), P(0.2f, 0.6f), k.fg, th);
        break;
    case Icon::Divide:
        line(-0.9f, 0.1f, 0.9f, 0.1f, k.fg);
        for (float x : {-0.9f, -0.3f, 0.3f, 0.9f}) dot(x, 0.1f, k.blue, 0.15f);
        break;
    case Icon::Merge:
        dl->AddCircle(P(-0.3f, 0), 0.45f * h, k.dim, 0, th);
        dl->AddCircle(P(0.3f, 0), 0.45f * h, k.dim, 0, th);
        dot(0, 0, k.blue, 0.2f);
        break;
    case Icon::Fix:
        line(0, -0.85f, 0, 0.25f, k.fg, th * 1.5f);
        hatch(-0.7f, 0.7f, 0.25f, k.green);
        break;
    case Icon::Pin:
        tri(-0.55f, k.green);
        hatch(-0.65f, 0.65f, 0.15f, k.green);
        break;
    case Icon::Roller:
        tri(-0.65f, k.green);
        dl->AddCircle(P(-0.25f, 0.15f), 0.12f * h, k.green, 0, th);
        dl->AddCircle(P(0.25f, 0.15f), 0.12f * h, k.green, 0, th);
        hatch(-0.65f, 0.65f, 0.32f, k.green);
        break;
    case Icon::RollerX:
        dl->AddTriangle(P(0.45f, 0), P(-0.2f, -0.45f), P(-0.2f, 0.45f), k.green, th);
        dl->AddCircle(P(-0.33f, -0.25f), 0.12f * h, k.green, 0, th);
        dl->AddCircle(P(-0.33f, 0.25f), 0.12f * h, k.green, 0, th);
        line(-0.5f, -0.7f, -0.5f, 0.7f, k.green);
        break;
    case Icon::Free:
        tri(-0.55f, k.dim);
        line(-0.6f, -0.6f, 0.6f, 0.6f, k.red, th * 1.4f);
        line(0.6f, -0.6f, -0.6f, 0.6f, k.red, th * 1.4f);
        break;
    case Icon::NLoad:
        arrow(0, -0.9f, 0, 0.35f, k.red);
        dot(0, 0.55f, k.fg, 0.2f);
        break;
    case Icon::DLoad:
        line(-0.9f, 0.7f, 0.9f, 0.7f, k.fg, th * 1.3f);
        line(-0.8f, -0.6f, 0.8f, -0.6f, k.red);
        for (float x : {-0.8f, -0.27f, 0.27f, 0.8f}) arrow(x, -0.6f, x, 0.55f, k.red);
        break;
    case Icon::PLoad:
        line(-0.9f, 0.7f, 0.9f, 0.7f, k.fg, th * 1.3f);
        arrow(0.1f, -0.9f, 0.1f, 0.55f, k.red);
        break;
    case Icon::MLoad: {
        line(-0.9f, 0.6f, 0.9f, 0.6f, k.fg, th * 1.3f);
        dl->PathArcTo(P(0, 0.1f), 0.55f * h, -1.1f * kPi, -0.1f * kPi, 16);
        dl->PathStroke(k.red, 0, th);
        arrowHead(dl, P(0.55f * std::cos(-0.1f * kPi), 0.1f + 0.55f * std::sin(-0.1f * kPi)), ImVec2(0.2f, 1), 0.35f * h, k.red);
        break;
    }
    case Icon::Temp:
        dl->AddRect(P(-0.15f, -0.85f), P(0.15f, 0.35f), k.fg, 0.15f * h, 0, th);
        dl->AddCircleFilled(P(0, 0.55f), 0.3f * h, k.red);
        line(0, -0.3f, 0, 0.4f, k.red, th * 1.6f);
        break;
    case Icon::Settle:
        tri(-0.75f, k.green);
        hatch(-0.6f, 0.6f, -0.05f, k.green);
        arrow(0, 0.2f, 0, 0.9f, k.red);
        break;
    case Icon::SelfW:
        dl->AddRectFilled(P(-0.55f, -0.75f), P(0.55f, 0.05f), (k.fg & 0x00FFFFFF) | 0x90000000, 0.1f * h);
        arrow(0, 0.1f, 0, 0.9f, k.red);
        text("g", 0, -0.35f, ImGui::GetColorU32(ImGuiCol_WindowBg), 0.55f);
        break;
    case Icon::Case:
        for (int i = 0; i < 3; ++i)
            dl->AddRectFilled(P(-0.7f + 0.15f * i, -0.6f + 0.35f * i), P(0.5f + 0.15f * i, -0.25f + 0.35f * i),
                              i == 2 ? k.blue : k.dim, 0.08f * h);
        break;
    case Icon::Combo: text("Σ", 0, 0, k.blue, 0.95f); break;
    case Icon::Run: dl->AddTriangleFilled(P(-0.5f, -0.75f), P(0.75f, 0), P(-0.5f, 0.75f), k.green); break;
    case Icon::Buckle: {
        std::vector<ImVec2> pts;
        for (int i = 0; i <= 16; ++i) {
            const float t = (float)i / 16;
            pts.push_back(P(0.45f * std::sin(t * kPi), -0.6f + 1.4f * t));
        }
        dl->AddPolyline(pts.data(), (int)pts.size(), k.blue, 0, th);
        line(0, -0.6f, 0, 0.8f, k.dim, th * 0.6f);
        arrow(0, -1.0f, 0, -0.62f, k.red);
        break;
    }
    case Icon::Modal: {
        std::vector<ImVec2> pts;
        for (int i = 0; i <= 24; ++i) {
            const float t = (float)i / 24;
            pts.push_back(P(-0.9f + 1.8f * t, -0.45f * std::sin(t * 2 * kPi)));
        }
        line(-0.9f, 0, 0.9f, 0, k.dim, th * 0.6f);
        dl->AddPolyline(pts.data(), (int)pts.size(), k.blue, 0, th * 1.2f);
        break;
    }
    case Icon::Influence:
        line(-0.9f, 0.5f, 0.9f, 0.5f, k.fg);
        dl->AddTriangleFilled(P(-0.9f, 0.5f), P(-0.1f, -0.55f), P(0.9f, 0.5f), (k.blue & 0x00FFFFFF) | 0x60000000);
        dl->AddTriangle(P(-0.9f, 0.5f), P(-0.1f, -0.55f), P(0.9f, 0.5f), k.blue, th);
        break;
    case Icon::Schema:
        line(-0.7f, 0.8f, -0.7f, -0.6f, k.fg);
        line(0.7f, 0.8f, 0.7f, -0.6f, k.fg);
        line(-0.7f, -0.6f, 0.7f, -0.6f, k.fg);
        hatch(-0.95f, -0.45f, 0.8f, k.green);
        hatch(0.45f, 0.95f, 0.8f, k.green);
        break;
    case Icon::Deformed: {
        line(-0.7f, 0.8f, -0.7f, -0.6f, k.dim, th * 0.7f);
        line(0.7f, 0.8f, 0.7f, -0.6f, k.dim, th * 0.7f);
        line(-0.7f, -0.6f, 0.7f, -0.6f, k.dim, th * 0.7f);
        dl->AddBezierCubic(P(-0.7f, 0.8f), P(-0.7f, 0.0f), P(-0.4f, -0.5f), P(-0.4f, -0.55f), k.red, th);
        dl->AddBezierCubic(P(-0.4f, -0.55f), P(0, -0.35f), P(0.6f, -0.75f), P(1.0f, -0.55f), k.red, th);
        dl->AddBezierCubic(P(1.0f, -0.55f), P(0.9f, -0.1f), P(0.7f, 0.4f), P(0.7f, 0.8f), k.red, th);
        break;
    }
    case Icon::DiagN:
        diagram(k.blue, 0);
        text("N", 0.6f, 0.6f, k.fg, 0.5f);
        break;
    case Icon::DiagQ:
        diagram(k.red, 1);
        text("Q", 0.6f, -0.65f, k.fg, 0.5f);
        break;
    case Icon::DiagM:
        diagram(k.blue, 2);
        text("M", -0.55f, -0.55f, k.fg, 0.5f);
        break;
    case Icon::Sigma: text("σ", 0, -0.05f, k.blue, 1.0f); break;
    case Icon::Contour:
        dl->AddRectFilledMultiColor(P(-0.8f, -0.7f), P(0.8f, 0.7f), rgb(30, 60, 220), rgb(240, 60, 30), rgb(240, 220, 40),
                                    rgb(40, 200, 120));
        break;
    case Icon::Envelope: {
        line(-0.9f, 0, 0.9f, 0, k.fg);
        std::vector<ImVec2> a, b;
        for (int i = 0; i <= 12; ++i) {
            const float x = -0.9f + 1.8f * i / 12, v = 1 - (x / 0.9f) * (x / 0.9f);
            a.push_back(P(x, 0.7f * v));
            b.push_back(P(x, -0.45f * v + 0.05f));
        }
        dl->AddPolyline(a.data(), (int)a.size(), k.blue, 0, th);
        dl->AddPolyline(b.data(), (int)b.size(), k.red, 0, th);
        break;
    }
    case Icon::Report:
        dl->AddRect(P(-0.6f, -0.85f), P(0.6f, 0.85f), k.fg, 0.08f * h, 0, th);
        for (float y : {-0.45f, -0.15f, 0.15f, 0.45f}) line(-0.35f, y, 0.35f, y, k.dim);
        break;
    case Icon::ZoomExt:
        dl->AddCircle(P(-0.15f, -0.15f), 0.55f * h, k.fg, 0, th);
        line(0.25f, 0.25f, 0.8f, 0.8f, k.fg, th * 1.6f);
        dl->AddRect(P(-0.42f, -0.42f), P(0.12f, 0.12f), k.blue, 0, 0, th);
        break;
    case Icon::Dist:
        arrow(0, 0, -0.85f, 0, k.fg);
        arrow(0, 0, 0.85f, 0, k.fg);
        line(-0.9f, -0.5f, -0.9f, 0.5f, k.blue);
        line(0.9f, -0.5f, 0.9f, 0.5f, k.blue);
        break;
    case Icon::Id:
        line(-0.8f, 0.2f, 0.8f, 0.2f, k.dim);
        line(-0.2f, -0.8f, -0.2f, 0.8f, k.dim);
        dot(-0.2f, 0.2f, k.blue, 0.2f);
        text("xy", 0.45f, -0.45f, k.fg, 0.5f);
        break;
    case Icon::List:
        for (float y : {-0.55f, 0.0f, 0.55f}) {
            dot(-0.65f, y, k.blue, 0.12f);
            line(-0.35f, y, 0.8f, y, k.fg);
        }
        break;
    case Icon::Props:
        dl->AddRect(P(-0.8f, -0.8f), P(0.8f, 0.8f), k.fg, 0, 0, th);
        line(-0.8f, -0.25f, 0.8f, -0.25f, k.fg, th * 0.7f);
        line(-0.8f, 0.3f, 0.8f, 0.3f, k.fg, th * 0.7f);
        line(-0.15f, -0.8f, -0.15f, 0.8f, k.fg, th * 0.7f);
        break;
    case Icon::New:
        dl->AddRect(P(-0.55f, -0.8f), P(0.55f, 0.8f), k.fg, 0, 0, th);
        line(0.15f, -0.8f, 0.55f, -0.4f, k.fg);
        break;
    case Icon::Open:
        dl->AddRect(P(-0.8f, -0.45f), P(0.8f, 0.7f), k.orange, 0.08f * h, 0, th);
        line(-0.8f, -0.45f, -0.35f, -0.7f, k.orange);
        line(-0.35f, -0.7f, 0.1f, -0.45f, k.orange);
        break;
    case Icon::Save:
        dl->AddRect(P(-0.7f, -0.7f), P(0.7f, 0.7f), k.blue, 0.08f * h, 0, th);
        dl->AddRectFilled(P(-0.4f, -0.7f), P(0.4f, -0.2f), k.blue);
        dl->AddRect(P(-0.4f, 0.15f), P(0.4f, 0.7f), k.blue, 0, 0, th);
        break;
    case Icon::Undo:
    case Icon::Redo: {
        const float sg = ic == Icon::Undo ? 1.f : -1.f;
        dl->PathArcTo(P(0, 0.15f), 0.55f * h, sg > 0 ? 1.1f * kPi : -0.1f * kPi, sg > 0 ? 2.05f * kPi : -1.05f * kPi, 16);
        dl->PathStroke(k.fg, 0, th);
        arrowHead(dl, P(-0.55f * sg, 0.12f), ImVec2(0, 1), 0.38f * h, k.fg);
        break;
    }
    case Icon::Sun:
        dl->AddCircleFilled(P(0, 0), 0.35f * h, k.orange);
        for (int i = 0; i < 8; ++i) {
            const float a = i * kPi / 4;
            line(0.5f * std::cos(a), 0.5f * std::sin(a), 0.8f * std::cos(a), 0.8f * std::sin(a), k.orange);
        }
        break;
    case Icon::Moon:
        dl->AddCircleFilled(P(0, 0), 0.6f * h, k.fg);
        dl->AddCircleFilled(P(0.3f, -0.25f), 0.5f * h, ImGui::GetColorU32(ImGuiCol_MenuBarBg));
        break;
    case Icon::Help: text("?", 0, 0, k.blue, 0.9f); break;
    }
}

// ============================================================================ khung ribbon
struct RibbonLayout {
    float top = 0;      // đỉnh vùng nội dung panel
    float content = 0;  // chiều cao vùng nút
    float title = 0;    // chiều cao dải tên panel
    float panelX = 0;
    IconColors ic;
    App* app = nullptr;
} R;

void beginPanel() {
    ImGui::SetCursorScreenPos(ImVec2(ImGui::GetCursorScreenPos().x, R.top));
    R.panelX = ImGui::GetCursorScreenPos().x;
    ImGui::BeginGroup();
}

void endPanel(const char* title) {
    ImGui::EndGroup();
    const float x0 = R.panelX - 4, x1 = std::max(ImGui::GetItemRectMax().x, R.panelX + 40) + 4;
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float ty = R.top + R.content;
    dl->AddRectFilled(ImVec2(x0, ty), ImVec2(x1, ty + R.title), R.app->ui.panelTitle);
    const ImVec2 ts = ImGui::CalcTextSize(title);
    dl->AddText(ImVec2((x0 + x1 - ts.x) * 0.5f, ty + (R.title - ts.y) * 0.5f), R.app->ui.textDim, title);
    dl->AddLine(ImVec2(x1 + 1, R.top + 2), ImVec2(x1 + 1, ty + R.title - 2), R.app->ui.separator);
    ImGui::SetCursorScreenPos(ImVec2(x1 + 10, R.top));
}

bool bigButton(const char* label, Icon icon, const char* tip, bool active = false) {
    ImGui::PushID(label);
    // nhãn có thể có '\n'
    float tw = 0;
    {
        std::string s = label;
        size_t i = 0;
        while (i <= s.size()) {
            size_t j = s.find('\n', i);
            if (j == std::string::npos) j = s.size();
            tw = std::max(tw, ImGui::CalcTextSize(s.substr(i, j - i).c_str()).x);
            i = j + 1;
        }
    }
    const float w = std::max(46.0f, tw + 10), hgt = R.content - 4;
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const bool pressed = ImGui::InvisibleButton("##b", ImVec2(w, hgt));
    const bool hov = ImGui::IsItemHovered(), held = ImGui::IsItemActive();
    ImDrawList* dl = ImGui::GetWindowDrawList();
    if (active || hov || held)
        dl->AddRectFilled(p, ImVec2(p.x + w, p.y + hgt), active ? R.app->ui.active : R.app->ui.hover, 3);
    const float isz = 28;
    drawIcon(dl, icon, ImVec2(p.x + w * 0.5f, p.y + 4 + isz * 0.5f), isz, R.ic);
    // nhãn
    std::string s = label;
    float y = p.y + 8 + isz;
    size_t i = 0;
    while (i <= s.size()) {
        size_t j = s.find('\n', i);
        if (j == std::string::npos) j = s.size();
        const std::string ln = s.substr(i, j - i);
        const ImVec2 ts = ImGui::CalcTextSize(ln.c_str());
        dl->AddText(ImVec2(p.x + (w - ts.x) * 0.5f, y), R.app->ui.text, ln.c_str());
        y += ts.y;
        i = j + 1;
    }
    if (tip && hov) ImGui::SetTooltip("%s", tip);
    ImGui::PopID();
    ImGui::SameLine(0, 2);
    return pressed;
}

bool smallButton(const char* label, Icon icon, const char* tip, bool active = false) {
    ImGui::PushID(label);
    const float isz = 16, hgt = (R.content - 8) / 3;
    const ImVec2 ts = ImGui::CalcTextSize(label);
    const float w = isz + 8 + ts.x + 6;
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const bool pressed = ImGui::InvisibleButton("##s", ImVec2(w, hgt));
    const bool hov = ImGui::IsItemHovered();
    ImDrawList* dl = ImGui::GetWindowDrawList();
    if (active || hov) dl->AddRectFilled(p, ImVec2(p.x + w, p.y + hgt), active ? R.app->ui.active : R.app->ui.hover, 3);
    drawIcon(dl, icon, ImVec2(p.x + 3 + isz * 0.5f, p.y + hgt * 0.5f), isz, R.ic);
    dl->AddText(ImVec2(p.x + isz + 8, p.y + (hgt - ts.y) * 0.5f), R.app->ui.text, label);
    if (tip && hov) ImGui::SetTooltip("%s", tip);
    ImGui::PopID();
    return pressed;
}

// Cột nút nhỏ (tối đa 3 nút xếp dọc)
void beginColumn() {
    ImGui::BeginGroup();
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(2, 2));
}
void endColumn() {
    ImGui::PopStyleVar();
    ImGui::EndGroup();
    ImGui::SameLine(0, 6);
}

bool iconButton(const char* id, Icon icon, const char* tip, float sz = 22) {
    ImGui::PushID(id);
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const bool pressed = ImGui::InvisibleButton("##i", ImVec2(sz, sz));
    if (ImGui::IsItemHovered()) {
        ImGui::GetWindowDrawList()->AddRectFilled(p, ImVec2(p.x + sz, p.y + sz), R.app->ui.hover, 3);
        if (tip) ImGui::SetTooltip("%s", tip);
    }
    drawIcon(ImGui::GetWindowDrawList(), icon, ImVec2(p.x + sz * 0.5f, p.y + sz * 0.5f), sz * 0.72f, R.ic);
    ImGui::PopID();
    return pressed;
}

} // namespace

// ============================================================================ chủ đề
void App::applyTheme() {
    ImGuiStyle& st = ImGui::GetStyle();
    ImVec4* c = st.Colors;
    st.WindowRounding = 0;
    st.ChildRounding = 0;
    st.FrameRounding = 2;
    st.PopupRounding = 2;
    st.TabRounding = 0;
    st.GrabRounding = 2;
    st.WindowBorderSize = 1;
    st.FrameBorderSize = darkTheme ? 0.0f : 1.0f;
    if (darkTheme) {
        ImGui::StyleColorsDark(&st);
        c[ImGuiCol_Text] = v4(232, 235, 240);
        c[ImGuiCol_TextDisabled] = v4(140, 150, 165);
        c[ImGuiCol_WindowBg] = v4(52, 60, 73);
        c[ImGuiCol_ChildBg] = v4(52, 60, 73, 0);
        c[ImGuiCol_PopupBg] = v4(44, 51, 62);
        c[ImGuiCol_Border] = v4(30, 35, 42);
        c[ImGuiCol_FrameBg] = v4(36, 42, 51);
        c[ImGuiCol_FrameBgHovered] = v4(46, 55, 68);
        c[ImGuiCol_FrameBgActive] = v4(54, 66, 82);
        c[ImGuiCol_TitleBg] = v4(40, 46, 56);
        c[ImGuiCol_TitleBgActive] = v4(46, 54, 66);
        c[ImGuiCol_TitleBgCollapsed] = v4(40, 46, 56);
        c[ImGuiCol_MenuBarBg] = v4(40, 46, 56);
        c[ImGuiCol_ScrollbarBg] = v4(44, 51, 62);
        c[ImGuiCol_ScrollbarGrab] = v4(84, 96, 114);
        c[ImGuiCol_ScrollbarGrabHovered] = v4(100, 114, 135);
        c[ImGuiCol_CheckMark] = v4(70, 170, 245);
        c[ImGuiCol_SliderGrab] = v4(70, 170, 245);
        c[ImGuiCol_SliderGrabActive] = v4(100, 190, 255);
        c[ImGuiCol_Button] = v4(66, 77, 93);
        c[ImGuiCol_ButtonHovered] = v4(80, 95, 116);
        c[ImGuiCol_ButtonActive] = v4(0, 120, 200);
        c[ImGuiCol_Header] = v4(0, 120, 200, 0.45f);
        c[ImGuiCol_HeaderHovered] = v4(0, 120, 200, 0.30f);
        c[ImGuiCol_HeaderActive] = v4(0, 120, 200, 0.60f);
        c[ImGuiCol_Separator] = v4(74, 85, 102);
        c[ImGuiCol_Tab] = v4(42, 49, 60);
        c[ImGuiCol_TabHovered] = v4(70, 82, 100);
        c[ImGuiCol_TabSelected] = v4(59, 68, 83);
        c[ImGuiCol_TabSelectedOverline] = v4(0, 140, 230);
        c[ImGuiCol_TabDimmed] = v4(42, 49, 60);
        c[ImGuiCol_TabDimmedSelected] = v4(52, 60, 73);
        c[ImGuiCol_TableHeaderBg] = v4(44, 51, 63);
        c[ImGuiCol_TableBorderStrong] = v4(30, 35, 42);
        c[ImGuiCol_TableBorderLight] = v4(62, 71, 86);
        c[ImGuiCol_TableRowBg] = v4(0, 0, 0, 0);
        c[ImGuiCol_TableRowBgAlt] = v4(255, 255, 255, 0.035f);
        c[ImGuiCol_DockingPreview] = v4(0, 140, 230, 0.6f);
        c[ImGuiCol_DockingEmptyBg] = v4(33, 40, 48);
        c[ImGuiCol_ResizeGrip] = v4(0, 0, 0, 0);
        ui = {rgb(38, 44, 54), rgb(59, 68, 83), rgb(49, 57, 69), rgb(78, 90, 108), rgb(232, 235, 240),
              rgb(165, 174, 188), rgb(78, 92, 112), rgb(0, 112, 192), rgb(0, 140, 230), rgb(36, 42, 51), rgb(44, 51, 62)};
    } else {
        ImGui::StyleColorsLight(&st);
        c[ImGuiCol_Text] = v4(25, 28, 34);
        c[ImGuiCol_TextDisabled] = v4(120, 128, 140);
        c[ImGuiCol_WindowBg] = v4(245, 246, 248);
        c[ImGuiCol_ChildBg] = v4(245, 246, 248, 0);
        c[ImGuiCol_PopupBg] = v4(252, 252, 253);
        c[ImGuiCol_Border] = v4(196, 202, 211);
        c[ImGuiCol_FrameBg] = v4(255, 255, 255);
        c[ImGuiCol_FrameBgHovered] = v4(232, 240, 250);
        c[ImGuiCol_FrameBgActive] = v4(220, 233, 248);
        c[ImGuiCol_TitleBg] = v4(226, 230, 236);
        c[ImGuiCol_TitleBgActive] = v4(214, 220, 229);
        c[ImGuiCol_MenuBarBg] = v4(226, 230, 236);
        c[ImGuiCol_CheckMark] = v4(0, 110, 210);
        c[ImGuiCol_SliderGrab] = v4(0, 110, 210);
        c[ImGuiCol_Button] = v4(230, 233, 238);
        c[ImGuiCol_ButtonHovered] = v4(208, 222, 240);
        c[ImGuiCol_ButtonActive] = v4(0, 110, 210);
        c[ImGuiCol_Header] = v4(0, 110, 210, 0.22f);
        c[ImGuiCol_HeaderHovered] = v4(0, 110, 210, 0.14f);
        c[ImGuiCol_HeaderActive] = v4(0, 110, 210, 0.32f);
        c[ImGuiCol_Separator] = v4(200, 205, 214);
        c[ImGuiCol_Tab] = v4(226, 230, 236);
        c[ImGuiCol_TabHovered] = v4(208, 222, 240);
        c[ImGuiCol_TabSelected] = v4(245, 246, 248);
        c[ImGuiCol_TabSelectedOverline] = v4(0, 110, 210);
        c[ImGuiCol_TabDimmed] = v4(226, 230, 236);
        c[ImGuiCol_TabDimmedSelected] = v4(238, 240, 243);
        c[ImGuiCol_TableHeaderBg] = v4(228, 232, 238);
        c[ImGuiCol_TableBorderStrong] = v4(196, 202, 211);
        c[ImGuiCol_TableBorderLight] = v4(222, 226, 232);
        c[ImGuiCol_TableRowBgAlt] = v4(0, 0, 0, 0.025f);
        c[ImGuiCol_DockingEmptyBg] = v4(255, 255, 255);
        c[ImGuiCol_ResizeGrip] = v4(0, 0, 0, 0);
        ui = {rgb(214, 219, 227), rgb(240, 242, 245), rgb(226, 230, 236), rgb(196, 202, 211), rgb(25, 28, 34),
              rgb(95, 104, 118), rgb(210, 224, 242), rgb(176, 208, 245), rgb(0, 110, 210), rgb(255, 255, 255), rgb(226, 230, 236)};
    }
    pal = makePalette(darkTheme);
}

// ============================================================================ ribbon
void App::appMenu() {
    if (!ImGui::BeginPopup("##appmenu")) return;
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
    if (ImGui::MenuItem("Xuất báo cáo kết quả...")) requestReport();
    ImGui::Separator();
    if (ImGui::MenuItem("Giao diện tối", nullptr, darkTheme)) {
        darkTheme = true;
        applyTheme();
    }
    if (ImGui::MenuItem("Giao diện sáng", nullptr, !darkTheme)) {
        darkTheme = false;
        applyTheme();
    }
    ImGui::Separator();
    if (ImGui::MenuItem("Giới thiệu & quy ước")) showAbout = true;
    if (ImGui::MenuItem("Thoát")) confirmDiscard(Pending::Quit);
    ImGui::EndPopup();
}

void App::ribbon() {
    ImGuiViewport* vp = ImGui::GetMainViewport();
    const float rowH = 30, content = 78, title = 18;
    const float height = rowH + content + title + 6;
    const ImGuiWindowFlags flags = ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings |
                                   ImGuiWindowFlags_NoScrollWithMouse;
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImGui::ColorConvertU32ToFloat4(ui.ribbon));
    const bool open = ImGui::BeginViewportSideBar("##ribbon", vp, ImGuiDir_Up, height, flags);
    ImGui::PopStyleColor();
    ImGui::PopStyleVar(2);
    if (!open) {
        ImGui::End();
        return;
    }
    R.app = this;
    R.ic = iconColors(darkTheme);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 w0 = ImGui::GetWindowPos(), wsz = ImGui::GetWindowSize();
    dl->AddRectFilled(w0, ImVec2(w0.x + wsz.x, w0.y + rowH), ui.tabRow);
    dl->AddLine(ImVec2(w0.x, w0.y + wsz.y - 1), ImVec2(w0.x + wsz.x, w0.y + wsz.y - 1), ui.separator);

    // ---- nút ứng dụng + truy cập nhanh + các thẻ
    ImGui::SetCursorScreenPos(ImVec2(w0.x + 6, w0.y + 3));
    {
        const ImVec2 p = ImGui::GetCursorScreenPos();
        if (ImGui::InvisibleButton("##app", ImVec2(44, 24))) ImGui::OpenPopup("##appmenu");
        const bool hov = ImGui::IsItemHovered();
        dl->AddRectFilled(p, ImVec2(p.x + 44, p.y + 24), hov ? rgb(30, 150, 240) : ui.accent, 3);
        const ImVec2 ts = ImGui::CalcTextSize("FEM");
        dl->AddText(ImVec2(p.x + (44 - ts.x) * 0.5f, p.y + (24 - ts.y) * 0.5f), rgb(255, 255, 255), "FEM");
        if (hov) ImGui::SetTooltip("Menu ứng dụng: tệp, ví dụ, giao diện sáng/tối");
    }
    appMenu();
    ImGui::SameLine(0, 10);
    ImGui::SetCursorScreenPos(ImVec2(ImGui::GetCursorScreenPos().x, w0.y + 4));
    if (iconButton("new", Icon::New, "Mới (Ctrl+N)") && confirmDiscard(Pending::New)) newModel();
    ImGui::SameLine(0, 2);
    if (iconButton("open", Icon::Open, "Mở (Ctrl+O)") && confirmDiscard(Pending::Open)) requestOpen();
    ImGui::SameLine(0, 2);
    if (iconButton("save", Icon::Save, "Lưu (Ctrl+S)")) requestSave(false);
    ImGui::SameLine(0, 2);
    if (iconButton("undo", Icon::Undo, "Hoàn tác (Ctrl+Z)")) undo();
    ImGui::SameLine(0, 2);
    if (iconButton("redo", Icon::Redo, "Làm lại (Ctrl+Y)")) redo();
    ImGui::SameLine(0, 18);

    const char* tabs[] = {"Trang chủ", "Tải trọng", "Phân tích", "Kết quả", "Xem"};
    for (int i = 0; i < 5; ++i) {
        const ImVec2 ts = ImGui::CalcTextSize(tabs[i]);
        const ImVec2 p(ImGui::GetCursorScreenPos().x, w0.y + 2);
        ImGui::SetCursorScreenPos(p);
        ImGui::PushID(i);
        if (ImGui::InvisibleButton("##tab", ImVec2(ts.x + 24, rowH - 2))) ribbonTab = i;
        const bool hov = ImGui::IsItemHovered();
        ImGui::PopID();
        if (i == ribbonTab) {
            dl->AddRectFilled(p, ImVec2(p.x + ts.x + 24, w0.y + rowH), ui.ribbon, 3, ImDrawFlags_RoundCornersTop);
            dl->AddLine(ImVec2(p.x + 2, p.y + 1), ImVec2(p.x + ts.x + 22, p.y + 1), ui.accent, 2);
        } else if (hov) {
            dl->AddRectFilled(p, ImVec2(p.x + ts.x + 24, w0.y + rowH), ui.hover, 3, ImDrawFlags_RoundCornersTop);
        }
        dl->AddText(ImVec2(p.x + 12, p.y + (rowH - 2 - ts.y) * 0.5f), i == ribbonTab ? ui.text : ui.textDim, tabs[i]);
        ImGui::SameLine(0, 2);
    }
    // tiêu đề tệp bên phải hàng thẻ
    {
        const std::string t = windowTitle();
        const ImVec2 ts = ImGui::CalcTextSize(t.c_str());
        const float x = w0.x + wsz.x - ts.x - 44;
        if (x > ImGui::GetCursorScreenPos().x + 20) dl->AddText(ImVec2(x, w0.y + (rowH - ts.y) * 0.5f), ui.textDim, t.c_str());
        ImGui::SetCursorScreenPos(ImVec2(w0.x + wsz.x - 32, w0.y + 4));
        if (iconButton("theme", darkTheme ? Icon::Sun : Icon::Moon, darkTheme ? "Chuyển sang giao diện sáng" : "Chuyển sang giao diện tối")) {
            darkTheme = !darkTheme;
            applyTheme();
        }
    }

    // ---- nội dung thẻ
    R.top = w0.y + rowH + 4;
    R.content = content;
    R.title = title;
    ImGui::SetCursorScreenPos(ImVec2(w0.x + 10, R.top));
    auto run = [&](const char* c) { startCommand(c); };
    const bool cmdIs = (bool)cmd;
    (void)cmdIs;
    switch (ribbonTab) {
    case 0: { // ---------------- Trang chủ
        beginPanel();
        if (bigButton("Thanh\nkhung", Icon::Line, "Vẽ thanh khung (lệnh L)")) run("L");
        if (bigButton("Thanh\ndàn", Icon::Truss, "Vẽ thanh dàn (lệnh DAN)")) run("DAN");
        beginColumn();
        if (smallButton("Nút", Icon::Node, "Thêm nút (lệnh PO)")) run("PO");
        if (smallButton("Lưới phẳng", Icon::Mesh, "Tạo lưới phần tử phẳng (LUOI)")) showMeshDialog = true;
        if (smallButton("Tạo khung", Icon::FrameGen, "Tạo khung nhiều tầng / dầm liên tục (TAOKHUNG)")) showFrameGen = true;
        endColumn();
        beginColumn();
        if (smallButton("Tạo dàn", Icon::TrussGen, "Tạo dàn Pratt/Howe/Warren/tam giác (TAODAN)")) showTrussGen = true;
        endColumn();
        endPanel("Vẽ");

        beginPanel();
        beginColumn();
        if (smallButton("Di chuyển", Icon::Move, "M")) run("M");
        if (smallButton("Sao chép", Icon::Copy, "CO")) run("CO");
        if (smallButton("Xoay", Icon::Rotate, "RO")) run("RO");
        endColumn();
        beginColumn();
        if (smallButton("Đối xứng", Icon::Mirror, "MI")) run("MI");
        if (smallButton("Chia thanh", Icon::Divide, "DIV")) run("DIV");
        if (smallButton("Gộp nút", Icon::Merge, "GOP: gộp các nút trùng nhau")) run("GOP");
        endColumn();
        beginColumn();
        if (smallButton("Xóa", Icon::Erase, "E (Delete)")) run("E");
        endColumn();
        endPanel("Sửa");

        beginPanel();
        if (bigButton("Ngàm", Icon::Fix, "NGAM: khóa ux, uy, φ")) run("NGAM");
        if (bigButton("Gối cố\nđịnh", Icon::Pin, "KHOP: khóa ux, uy")) run("KHOP");
        if (bigButton("Gối di\nđộng", Icon::Roller, "DIDONG: khóa uy")) run("DIDONG");
        beginColumn();
        if (smallButton("Di động ngang", Icon::RollerX, "DIDONGX: khóa ux")) run("DIDONGX");
        if (smallButton("Ngàm trượt", Icon::Fix, "NGAMTRUOT: khóa ux, φ")) run("NGAMTRUOT");
        if (smallButton("Bỏ gối", Icon::Free, "BOGOI")) run("BOGOI");
        endColumn();
        endPanel("Gối tựa");

        beginPanel();
        {
            std::vector<std::string> mats, secs;
            for (auto& m : model.materials) mats.push_back(m.name);
            for (auto& s : model.sections) secs.push_back(s.name);
            ImGui::SetCursorScreenPos(ImVec2(ImGui::GetCursorScreenPos().x, R.top + 10));
            ImGui::BeginGroup();
            ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(6, 8));
            ImGui::AlignTextToFramePadding();
            ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(ui.textDim), "Vật liệu ");
            ImGui::SameLine(70);
            ImGui::SetNextItemWidth(170);
            if (ImGui::BeginCombo("##rmat", curMat < (int)mats.size() ? mats[curMat].c_str() : "—")) {
                for (int i = 0; i < (int)mats.size(); ++i)
                    if (ImGui::Selectable(mats[i].c_str(), i == curMat)) curMat = i;
                ImGui::EndCombo();
            }
            ImGui::AlignTextToFramePadding();
            ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(ui.textDim), "Tiết diện");
            ImGui::SameLine(70);
            ImGui::SetNextItemWidth(170);
            if (ImGui::BeginCombo("##rsec", curSec < (int)secs.size() ? secs[curSec].c_str() : "—")) {
                for (int i = 0; i < (int)secs.size(); ++i)
                    if (ImGui::Selectable(secs[i].c_str(), i == curSec)) curSec = i;
                ImGui::EndCombo();
            }
            ImGui::PopStyleVar();
            ImGui::EndGroup();
            ImGui::SameLine(0, 6);
            ImGui::SetCursorScreenPos(ImVec2(ImGui::GetCursorScreenPos().x, R.top));
            bool anyBar = false;
            for (int e : selElems) anyBar |= model.elements[e].isBar();
            ImGui::BeginDisabled(!anyBar);
            if (bigButton("Gán cho\nthanh chọn", Icon::Props, "Gán vật liệu và tiết diện đang chọn cho các thanh đã chọn")) {
                snapshot();
                for (int e : selElems)
                    if (model.elements[e].isBar()) {
                        model.elements[e].mat = curMat;
                        model.elements[e].sec = curSec;
                    }
                changed();
            }
            ImGui::EndDisabled();
        }
        endPanel("Thuộc tính thanh");

        beginPanel();
        beginColumn();
        if (smallButton("Khoảng cách", Icon::Dist, "DI")) run("DI");
        if (smallButton("Tọa độ", Icon::Id, "ID")) run("ID");
        if (smallButton("Liệt kê", Icon::List, "LI")) run("LI");
        endColumn();
        if (bigButton("Vừa\nkhung", Icon::ZoomExt, "Zoom toàn bộ (Z E, nhấp đúp chuột giữa)")) fitRequested = true;
        endPanel("Tiện ích");
        break;
    }
    case 1: { // ---------------- Tải trọng
        beginPanel();
        {
            ImGui::SetCursorScreenPos(ImVec2(ImGui::GetCursorScreenPos().x, R.top + 6));
            ImGui::BeginGroup();
            ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(ui.textDim), "Trường hợp tải đang soạn");
            ImGui::SetNextItemWidth(200);
            if (ImGui::BeginCombo("##rcase", model.cases[activeCase].name.c_str())) {
                for (int i = 0; i < (int)model.cases.size(); ++i) {
                    ImGui::PushID(i);
                    if (ImGui::Selectable(model.cases[i].name.c_str(), i == activeCase)) activeCase = resultSet = i;
                    ImGui::PopID();
                }
                ImGui::EndCombo();
            }
            if (ImGui::Button("+ Trường hợp tải")) {
                snapshot();
                model.cases.push_back(LoadCase("TH" + std::to_string(model.cases.size() + 1)));
                activeCase = resultSet = (int)model.cases.size() - 1;
                changed();
            }
            ImGui::SameLine();
            if (ImGui::Button("Quản lý...")) showModelPanel = true;
            ImGui::EndGroup();
            ImGui::SameLine(0, 6);
        }
        endPanel("Trường hợp tải");

        beginPanel();
        if (bigButton("Tải nút", Icon::NLoad, "TAINUT: lực Fx, Fy, Mz tại nút chọn")) run("TAINUT");
        if (bigButton("Phân bố", Icon::DLoad, "TAIPB: tải phân bố đều/hình thang trên thanh")) run("TAIPB");
        if (bigButton("Tập trung", Icon::PLoad, "TAITT: lực tập trung, nhấp vị trí trên thanh")) run("TAITT");
        beginColumn();
        if (smallButton("Mômen", Icon::MLoad, "MOMEN: mômen tập trung trên thanh / tại nút")) run("MOMEN");
        if (smallButton("Nhiệt độ, lún...", Icon::Temp, "Nhiệt độ, chuyển vị cưỡng bức, tải cạnh: chọn đối tượng rồi dùng bảng Thuộc tính"))
            showPropsPanel = true;
        if (smallButton("Trọng lượng bản thân", Icon::SelfW, "Bật/tắt trọng lượng bản thân cho trường hợp tải đang soạn")) {
            snapshot();
            double& sw = model.cases[activeCase].selfWeight;
            sw = sw != 0 ? 0 : 1;
            echo(sw != 0 ? "Đã bật trọng lượng bản thân." : "Đã tắt trọng lượng bản thân.");
            changed();
        }
        endColumn();
        endPanel("Tải trọng");

        beginPanel();
        if (bigButton("Sinh tổ\nhợp", Icon::Combo, "Sinh tổ hợp cơ bản 1, 2 theo TCVN 2737:1995")) {
            snapshot();
            std::string note;
            const int n = generateCombos(model, 1.0, 0.9, 300, &note);
            echo("Đã sinh " + std::to_string(n) + " tổ hợp." + (note.empty() ? "" : " " + note));
            changed();
        }
        if (bigButton("Quản lý\ntổ hợp", Icon::Case, "Bảng tổ hợp (thẻ Tổ hợp trong panel Mô hình)")) showModelPanel = true;
        endPanel("Tổ hợp");
        break;
    }
    case 2: { // ---------------- Phân tích
        beginPanel();
        if (bigButton("Tính\n(F5)", Icon::Run, "Phân tích tĩnh tuyến tính")) startCommand("TINH");
        beginColumn();
        ImGui::Checkbox("Tự động tính", &autoAnalyze);
        {
            bool so = model.secondOrder;
            if (ImGui::Checkbox("Sơ đồ biến dạng (P-Δ)", &so)) startCommand("PDELTA");
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Phân tích bậc hai: kể ảnh hưởng của lực dọc lên độ cứng (lệnh PDELTA)");
        }
        ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(ui.textDim), result.ok ? "%d phương trình" : "Chưa tính", result.numEquations);
        endColumn();
        endPanel("Tĩnh học");
        beginPanel();
        if (bigButton("Cấu tạo\nhình học", Icon::Truss, "Bậc tự do W, bậc siêu tĩnh, kết luận bất biến hình (lệnh CTHH)")) startCommand("CTHH");
        endPanel("Cấu tạo");
        beginPanel();
        if (bigButton("Ổn định", Icon::Buckle, "Tải trọng tới hạn, hệ số chiều dài tính toán μ")) startCommand("ONDINH");
        if (bigButton("Dao động\nriêng", Icon::Modal, "Tần số, chu kỳ, dạng dao động")) startCommand("DAODONG");
        endPanel("Trị riêng");
        beginPanel();
        if (bigButton("Đường ảnh\nhưởng", Icon::Influence, "Thiết lập trong panel Kết quả", view == ResultView::Influence)) {
            view = ResultView::Influence;
            showResultPanel = true;
        }
        endPanel("Tải trọng di động");
        break;
    }
    case 3: { // ---------------- Kết quả
        beginPanel();
        {
            std::vector<std::string> names;
            for (auto& c : model.cases) names.push_back("TH: " + c.name);
            for (auto& c : model.combos) names.push_back("THợp: " + c.name);
            ImGui::SetCursorScreenPos(ImVec2(ImGui::GetCursorScreenPos().x, R.top + 6));
            ImGui::BeginGroup();
            ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(ui.textDim), "Bộ kết quả");
            ImGui::SetNextItemWidth(220);
            if (ImGui::BeginCombo("##rset", names.empty() ? "" : names[std::min(resultSet, (int)names.size() - 1)].c_str())) {
                for (int i = 0; i < (int)names.size(); ++i) {
                    ImGui::PushID(i);
                    if (ImGui::Selectable(names[i].c_str(), i == resultSet)) {
                        resultSet = i;
                        if (i < (int)model.cases.size()) activeCase = i;
                    }
                    ImGui::PopID();
                }
                ImGui::EndCombo();
            }
            ImGui::Checkbox("Biểu đồ bao", &showEnvelope);
            ImGui::EndGroup();
            ImGui::SameLine(0, 6);
        }
        endPanel("Bộ kết quả");
        beginPanel();
        struct VB {
            const char* t;
            Icon ic;
            ResultView v;
            const char* tip;
        };
        const VB vbs[] = {{"Sơ đồ", Icon::Schema, ResultView::None, "Chỉ sơ đồ tính"},
                          {"Biến\ndạng", Icon::Deformed, ResultView::Deformed, "Hình biến dạng"},
                          {"N", Icon::DiagN, ResultView::N, "Biểu đồ lực dọc"},
                          {"Q", Icon::DiagQ, ResultView::Q, "Biểu đồ lực cắt"},
                          {"M", Icon::DiagM, ResultView::M, "Biểu đồ mômen (vẽ về phía thớ căng)"},
                          {"σ", Icon::Sigma, ResultView::Sigma, "Ứng suất pháp thớ biên N/A ± M/W"},
                          {"Phổ\nmàu", Icon::Contour, ResultView::Contour, "Phổ màu ứng suất/chuyển vị phần tử phẳng"},
                          {"Đường\nảnh hưởng", Icon::Influence, ResultView::Influence, "Đường ảnh hưởng"},
                          {"Dạng\nriêng", Icon::Modal, ResultView::Mode, "Dạng mất ổn định / dao động"}};
        for (const auto& b : vbs)
            if (bigButton(b.t, b.ic, b.tip, view == b.v)) {
                view = b.v;
                if (view == ResultView::Mode && !(eig.ok && eigVersion == modelVersion)) showEigenDialog = true;
            }
        endPanel("Hiển thị");
        beginPanel();
        {
            ImGui::SetCursorScreenPos(ImVec2(ImGui::GetCursorScreenPos().x, R.top + 4));
            ImGui::BeginGroup();
            ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(6, 3));
            ImGui::SetNextItemWidth(130);
            ImGui::SliderFloat("Biến dạng", &deformFactor, 0.05f, 20.0f, "×%.2f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetNextItemWidth(130);
            ImGui::SliderFloat("Biểu đồ", &diagramFactor, 0.2f, 5.0f, "×%.2f", ImGuiSliderFlags_Logarithmic);
            ImGui::Checkbox("Giá trị", &showValues);
            ImGui::SameLine();
            ImGui::Checkbox("Phản lực", &showReactions);
            ImGui::PopStyleVar();
            ImGui::EndGroup();
            ImGui::SameLine(0, 6);
        }
        endPanel("Tỉ lệ");
        beginPanel();
        if (bigButton("Báo cáo", Icon::Report, "Xuất báo cáo kết quả (.txt)")) requestReport();
        endPanel("Xuất");
        break;
    }
    case 4: { // ---------------- Xem
        beginPanel();
        beginColumn();
        if (smallButton("Thuộc tính (Ctrl+1)", Icon::Props, nullptr, showPropsPanel)) showPropsPanel = !showPropsPanel;
        if (smallButton("Mô hình", Icon::List, nullptr, showModelPanel)) showModelPanel = !showModelPanel;
        if (smallButton("Kết quả", Icon::DiagM, nullptr, showResultPanel)) showResultPanel = !showResultPanel;
        endColumn();
        beginColumn();
        if (smallButton("Nhật ký", Icon::Report, nullptr, showLogPanel)) showLogPanel = !showLogPanel;
        if (smallButton("Khôi phục bố cục", Icon::Props, nullptr)) resetLayout = true;
        endColumn();
        endPanel("Bảng");
        beginPanel();
        beginColumn();
        if (smallButton("Số hiệu nút", Icon::Node, nullptr, showNodeIds)) showNodeIds = !showNodeIds;
        if (smallButton("Số hiệu phần tử", Icon::Line, nullptr, showElemIds)) showElemIds = !showElemIds;
        if (smallButton("Tải trọng", Icon::DLoad, nullptr, showLoads)) showLoads = !showLoads;
        endColumn();
        beginColumn();
        if (smallButton("Gối tựa", Icon::Pin, nullptr, showSupports)) showSupports = !showSupports;
        if (smallButton("Lưới (F7)", Icon::Mesh, nullptr, showGrid)) showGrid = !showGrid;
        if (smallButton("Giá trị", Icon::List, nullptr, showValues)) showValues = !showValues;
        endColumn();
        endPanel("Hiển thị");
        beginPanel();
        if (bigButton("Tối", Icon::Moon, "Giao diện tối", darkTheme)) {
            darkTheme = true;
            applyTheme();
        }
        if (bigButton("Sáng", Icon::Sun, "Giao diện sáng", !darkTheme)) {
            darkTheme = false;
            applyTheme();
        }
        endPanel("Giao diện");
        beginPanel();
        if (bigButton("Vừa\nkhung", Icon::ZoomExt, "Zoom toàn bộ")) fitRequested = true;
        if (bigButton("Trợ giúp", Icon::Help, "Quy ước và danh sách lệnh")) {
            showAbout = true;
            startCommand("?");
        }
        endPanel("Điều hướng");
        break;
    }
    }
    ImGui::End();
}

// ============================================================================ thanh trạng thái
void App::statusBar() {
    ImGuiViewport* vp = ImGui::GetMainViewport();
    const float h = ImGui::GetFrameHeight() + 6;
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(8, 3));
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImGui::ColorConvertU32ToFloat4(ui.statusBg));
    const bool open = ImGui::BeginViewportSideBar("##status", vp, ImGuiDir_Down, h,
                                                  ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings);
    ImGui::PopStyleColor();
    ImGui::PopStyleVar();
    if (open) {
        R.app = this;
        ImGui::AlignTextToFramePadding();
        ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(ui.accent), "MÔ HÌNH");
        ImGui::SameLine(0, 16);
        char b[96];
        if (cursorInCanvas) std::snprintf(b, sizeof b, "%10.4f, %10.4f", cursorWX, cursorWY);
        else std::snprintf(b, sizeof b, "%22s", "");
        ImGui::TextUnformatted(b);
        ImGui::SameLine(0, 16);
        auto toggle = [&](const char* label, bool& v, const char* tip) {
            const ImVec2 ts = ImGui::CalcTextSize(label);
            const ImVec2 p = ImGui::GetCursorScreenPos();
            const float w = ts.x + 14, hh = ImGui::GetFrameHeight();
            ImGui::PushID(label);
            if (ImGui::InvisibleButton("##t", ImVec2(w, hh))) v = !v;
            const bool hov = ImGui::IsItemHovered();
            ImGui::PopID();
            ImDrawList* dl = ImGui::GetWindowDrawList();
            if (v) dl->AddRectFilled(p, ImVec2(p.x + w, p.y + hh), ui.active, 3);
            else if (hov) dl->AddRectFilled(p, ImVec2(p.x + w, p.y + hh), ui.hover, 3);
            dl->AddText(ImVec2(p.x + 7, p.y + (hh - ts.y) * 0.5f), v ? ui.text : ui.textDim, label);
            if (hov && tip) ImGui::SetTooltip("%s", tip);
            ImGui::SameLine(0, 3);
        };
        toggle("GRID", showGrid, "Hiện lưới (F7)");
        toggle("SNAP", snapGrid, "Bắt điểm lưới (F9)");
        toggle("ORTHO", ortho, "Chế độ trực giao: chỉ vẽ ngang/đứng (F8)");
        toggle("OSNAP", osnap, "Bắt điểm đối tượng: nút, trung điểm, giao điểm, vuông góc, gần nhất (F3)");
        toggle("DYN", dynInput, "Nhập liệu động cạnh con trỏ (F12)");
        ImGui::SameLine(0, 14);
        ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(ui.textDim), "Bước lưới");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(60);
        ImGui::InputDouble("##grid", &gridStep, 0, 0, "%.3g");
        gridStep = std::max(gridStep, 1e-3);
        ImGui::SameLine(0, 14);
        ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(ui.textDim), "kN, m");
        // bên phải: trạng thái tính toán
        std::string st;
        ImU32 col = ui.textDim;
        if (!result.ok && !result.error.empty()) {
            st = "Lỗi tính toán";
            col = rgb(230, 70, 60);
        } else if (needAnalyze && !autoAnalyze) {
            st = "Cần tính lại (F5)";
            col = rgb(230, 150, 30);
        } else if (result.ok) {
            std::snprintf(b, sizeof b, "Đã tính: %d phương trình", result.numEquations);
            st = b;
        }
        const ImVec2 ts = ImGui::CalcTextSize(st.c_str());
        const float x = ImGui::GetWindowPos().x + ImGui::GetWindowWidth() - ts.x - 12;
        if (x > ImGui::GetCursorScreenPos().x + 20) {
            ImGui::SameLine();
            ImGui::SetCursorScreenPos(ImVec2(x, ImGui::GetCursorScreenPos().y));
            ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(col), "%s", st.c_str());
        }
    }
    ImGui::End();
}

// ============================================================================ dòng lệnh
std::string App::currentPrompt() {
    if (!cmd) return "Nhập lệnh:";
    return cmd->name + "  " + cmd->prompt(*this);
}

void App::commandLine() {
    ImGuiViewport* vp = ImGui::GetMainViewport();
    const float lh = ImGui::GetTextLineHeightWithSpacing();
    const int nh = cmdExpanded ? 12 : 3;
    const float h = lh * nh + ImGui::GetFrameHeight() + 14;
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(8, 4));
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImGui::ColorConvertU32ToFloat4(ui.cmdBg));
    const bool open = ImGui::BeginViewportSideBar("##cmdline", vp, ImGuiDir_Down, h,
                                                  ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings);
    ImGui::PopStyleColor();
    ImGui::PopStyleVar();
    if (open) {
        ImGui::BeginChild("##hist", ImVec2(0, lh * nh), ImGuiChildFlags_None);
        const int start = std::max(0, (int)log.size() - 200);
        for (int i = start; i < (int)log.size(); ++i) {
            const std::string& l = log[i];
            const bool err = l.rfind("Lỗi", 0) == 0;
            const bool cmdEcho = l.rfind("Lệnh:", 0) == 0;
            if (err) ImGui::PushStyleColor(ImGuiCol_Text, v4(235, 80, 70));
            else if (!cmdEcho) ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(ui.textDim));
            ImGui::TextUnformatted(l.c_str());
            if (err || !cmdEcho) ImGui::PopStyleColor();
        }
        static size_t lastSize = 0;
        if (log.size() != lastSize) {
            ImGui::SetScrollHereY(1.0f);
            lastSize = log.size();
        }
        ImGui::EndChild();
        // dòng nhập lệnh
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const ImVec2 p = ImGui::GetCursorScreenPos();
        const float w = ImGui::GetContentRegionAvail().x, fh = ImGui::GetFrameHeight();
        dl->AddRectFilled(p, ImVec2(p.x + w, p.y + fh), ImGui::GetColorU32(ImGuiCol_FrameBg), 2);
        dl->AddRect(p, ImVec2(p.x + w, p.y + fh), ui.separator, 2);
        const std::string pr = currentPrompt();
        const float ty = p.y + (fh - ImGui::GetTextLineHeight()) * 0.5f;
        dl->AddText(ImVec2(p.x + 6, ty), ui.accent, "▸");
        dl->AddText(ImVec2(p.x + 20, ty), ui.textDim, pr.c_str());
        const float px = p.x + 26 + ImGui::CalcTextSize(pr.c_str()).x;
        dl->AddText(ImVec2(px, ty), ui.text, cmdLine.c_str());
        const float cx = px + ImGui::CalcTextSize(cmdLine.c_str()).x + 1;
        if (std::fmod(ImGui::GetTime(), 1.0) < 0.6) dl->AddLine(ImVec2(cx, ty), ImVec2(cx, ty + ImGui::GetTextLineHeight()), ui.text, 1.5f);
        ImGui::Dummy(ImVec2(w, fh));
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Gõ lệnh ở bất cứ đâu (ví dụ L, E, M, CO, Z, ? ). Enter/Space: thực hiện hoặc lặp lệnh; Esc: hủy; F2: mở rộng lịch sử");
    }
    ImGui::End();
}

// ============================================================================ xử lý lệnh
std::set<int> App::selectionNodes() const {
    std::set<int> s(selNodes.begin(), selNodes.end());
    for (int e : selElems)
        if (e >= 0 && e < (int)model.elements.size())
            for (int k = 0; k < model.elements[e].nodeCount(); ++k) s.insert(model.elements[e].n[k]);
    return s;
}

int App::nodeForSnap(const SnapPoint& p) {
    double x0, y0, x1, y1;
    model.bounds(x0, y0, x1, y1);
    const double tol = 1e-6 * std::max({1.0, x1 - x0, y1 - y0, std::abs(p.x), std::abs(p.y)});
    if (p.node >= 0 && p.node < (int)model.nodes.size() && std::abs(model.nodes[p.node].x - p.x) <= tol &&
        std::abs(model.nodes[p.node].y - p.y) <= tol)
        return p.node;
    const int existing = model.findNode(p.x, p.y, tol);
    if (existing >= 0) return existing;
    // điểm nằm trong thanh → cắt thanh để nối
    std::vector<std::pair<int, double>> hits;
    for (size_t e = 0; e < model.elements.size(); ++e) {
        const Element& el = model.elements[e];
        if (!el.isBar()) continue;
        const Node& A = model.nodes[el.n[0]];
        const Node& B = model.nodes[el.n[1]];
        const double dx = B.x - A.x, dy = B.y - A.y, L2 = dx * dx + dy * dy;
        if (L2 <= 0) continue;
        const double t = ((p.x - A.x) * dx + (p.y - A.y) * dy) / L2;
        if (t <= 1e-9 || t >= 1 - 1e-9) continue;
        if (std::hypot(A.x + t * dx - p.x, A.y + t * dy - p.y) <= tol) hits.push_back({(int)e, t});
    }
    if (hits.empty()) return model.addNode(p.x, p.y, tol);
    int first = -1;
    for (const auto& [e, t] : hits) {
        const int n = splitBarAt(model, e, t);
        if (first < 0) first = n;
    }
    if (hits.size() > 1) mergeNodes(model, tol);
    return first;
}

void App::cancelCommand() {
    if (!cmd) return;
    cmd->onCancel(*this);
    echo("*Hủy*");
    cmd.reset();
}

void App::startCommand(const std::string& raw) {
    const std::string u = upper(trim(raw));
    if (u.empty()) return;
    if (cmd) cancelCommand();
    echo("Lệnh: " + u);
    auto is = [&](std::initializer_list<const char*> names) {
        for (const char* n : names)
            if (u == n) {
                lastCommand = *names.begin();
                return true;
            }
        return false;
    };
    if (is({"U", "UNDO"})) {
        undo();
        return;
    }
    if (is({"REDO"})) {
        redo();
        return;
    }
    if (is({"TINH", "RUN", "CALC"})) {
        analyzeNow(true);
        if (result.ok)
            echo("Tính xong: " + std::to_string(result.numEquations) + " phương trình, " +
                 std::to_string((int)result.timeMs) + " ms.");
        return;
    }
    if (is({"GOP", "OVERKILL", "MERGE"})) {
        double x0, y0, x1, y1;
        model.bounds(x0, y0, x1, y1);
        snapshot();
        const int n = mergeNodes(model, 1e-6 * std::max({1.0, x1 - x0, y1 - y0}));
        echo("Đã gộp " + std::to_string(n) + " nút trùng.");
        if (n) changed();
        else undoStack.pop_back();
        return;
    }
    if (is({"LUOI", "MESH"})) {
        showMeshDialog = true;
        return;
    }
    if (is({"TAOKHUNG", "FRAMEGEN"})) {
        showFrameGen = true;
        return;
    }
    if (is({"TAODAN", "TRUSSGEN"})) {
        showTrussGen = true;
        return;
    }
    if (is({"CTHH", "KINEMATIC"})) {
        kin = kinematicAnalysis(model);
        kinVersion = modelVersion;
        showKinematics = true;
        if (kin.applicable) {
            echo("W = " + std::to_string(kin.W) + ", " + kin.verdict);
            if (!kin.mechanismU.empty()) showMechanism();
        } else echo(kin.note);
        return;
    }
    if (is({"PDELTA"})) {
        snapshot();
        model.secondOrder = !model.secondOrder;
        echo(model.secondOrder ? "Bật tính theo sơ đồ biến dạng (P-Δ)." : "Tắt tính theo sơ đồ biến dạng.");
        changed();
        return;
    }
    if (is({"ONDINH", "BUCKLING"})) {
        eigOpt.kind = EigenKind::Buckling;
        showEigenDialog = true;
        return;
    }
    if (is({"DAODONG", "MODAL"})) {
        eigOpt.kind = EigenKind::Vibration;
        showEigenDialog = true;
        return;
    }
    if (is({"PR", "PROPERTIES", "CH"})) {
        showPropsPanel = !showPropsPanel;
        return;
    }
    if (is({"NEW"})) {
        if (confirmDiscard(Pending::New)) newModel();
        return;
    }
    if (is({"OPEN"})) {
        if (confirmDiscard(Pending::Open)) requestOpen();
        return;
    }
    if (is({"SAVE", "QSAVE"})) {
        requestSave(false);
        return;
    }
    if (is({"SAVEAS"})) {
        requestSave(true);
        return;
    }
    if (is({"?", "HELP", "TROGIUP"})) {
        echo("Các lệnh (không phân biệt hoa thường):");
        for (const auto& [k, d] : commandHelp()) echo("  " + k + " — " + d);
        echo("  Tọa độ: x,y   @dx,dy   @d<góc   hoặc gõ khoảng cách theo hướng con trỏ. F3 OSNAP, F8 ORTHO, F9 SNAP.");
        return;
    }
    std::string canonical;
    cmd = makeCommand(u, canonical);
    if (!cmd) {
        echo("Lỗi: lệnh không xác định \"" + u + "\". Gõ ? để xem danh sách lệnh.");
        return;
    }
    lastCommand = canonical;
    cmd->start(*this);
    if (cmd && cmd->done) cmd.reset();
}

void App::submitCommandLine() {
    const std::string s = trim(cmdLine);
    cmdLine.clear();
    if (!s.empty()) {
        typedHistory.push_back(s);
        typedPos = -1;
    }
    if (!cmd) {
        if (!s.empty()) startCommand(s);
        else if (!lastCommand.empty()) startCommand(lastCommand); // Enter/Space lặp lại lệnh trước
        return;
    }
    if (s.empty()) {
        cmd->onEnter(*this);
    } else {
        echo(cmd->prompt(*this) + " " + s);
        bool ok = cmd->onText(*this, s);
        if (!ok && cmd && cmd->wantsPoint()) {
            double bx = 0, by = 0;
            const bool hb = cmd->base(bx, by);
            double x, y;
            if (parsePoint(s, hb, bx, by, cursorWX - bx, cursorWY - by, x, y)) {
                SnapPoint sp;
                sp.x = x;
                sp.y = y;
                cmd->onPoint(*this, sp);
                ok = true;
            }
        }
        if (!ok) echo("Lỗi: giá trị không hợp lệ.");
    }
    if (cmd && cmd->done) cmd.reset();
}

void App::commandKeys() {
    ImGuiIO& io = ImGui::GetIO();
    // phím chức năng như AutoCAD
    const ImGuiInputFlags g = ImGuiInputFlags_RouteGlobal;
    if (ImGui::Shortcut(ImGuiKey_F2, g)) cmdExpanded = !cmdExpanded;
    if (ImGui::Shortcut(ImGuiKey_F3, g)) osnap = !osnap;
    if (ImGui::Shortcut(ImGuiKey_F7, g)) showGrid = !showGrid;
    if (ImGui::Shortcut(ImGuiKey_F8, g)) ortho = !ortho;
    if (ImGui::Shortcut(ImGuiKey_F9, g)) snapGrid = !snapGrid;
    if (ImGui::Shortcut(ImGuiKey_F12, g)) dynInput = !dynInput;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_1, g)) showPropsPanel = !showPropsPanel;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_A, g)) {
        for (int i = 0; i < (int)model.nodes.size(); ++i) selNodes.insert(i);
        for (int i = 0; i < (int)model.elements.size(); ++i) selElems.insert(i);
    }
    // gõ phím ở bất cứ đâu → dòng lệnh (trừ khi đang nhập vào một ô)
    if (io.WantTextInput || ImGui::GetTopMostPopupModal() != nullptr) return;
    if (io.KeyCtrl || io.KeySuper) return;
    for (ImWchar c : io.InputQueueCharacters) {
        if (c == ' ') submitCommandLine();
        else if (c >= 32 && c != 127) appendUtf8(cmdLine, c);
    }
    if (ImGui::IsKeyPressed(ImGuiKey_Enter) || ImGui::IsKeyPressed(ImGuiKey_KeypadEnter)) submitCommandLine();
    if (ImGui::IsKeyPressed(ImGuiKey_Backspace) && !cmdLine.empty()) popUtf8(cmdLine);
    if (ImGui::IsKeyPressed(ImGuiKey_Escape)) {
        if (!cmdLine.empty()) cmdLine.clear();
        else if (cmd) cancelCommand();
        else {
            selNodes.clear();
            selElems.clear();
        }
    }
    if (ImGui::IsKeyPressed(ImGuiKey_Delete) && !cmd && (!selNodes.empty() || !selElems.empty())) startCommand("E");
    if (ImGui::IsKeyPressed(ImGuiKey_UpArrow) && !typedHistory.empty()) {
        typedPos = typedPos < 0 ? (int)typedHistory.size() - 1 : std::max(0, typedPos - 1);
        cmdLine = typedHistory[typedPos];
    }
    if (ImGui::IsKeyPressed(ImGuiKey_DownArrow) && typedPos >= 0) {
        typedPos = std::min((int)typedHistory.size(), typedPos + 1);
        cmdLine = typedPos < (int)typedHistory.size() ? typedHistory[typedPos] : "";
        if (typedPos >= (int)typedHistory.size()) typedPos = -1;
    }
}
