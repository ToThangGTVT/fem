// Hệ thống lệnh kiểu AutoCAD: dòng lệnh, nhập điểm (chuột hoặc tọa độ), bắt điểm, chọn đối tượng.
#pragma once
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "imgui.h"

struct App;

enum class SnapKind { None, Grid, Node, Mid, Intersection, Perpendicular, Nearest };

struct SnapPoint {
    double x = 0, y = 0;
    SnapKind kind = SnapKind::None;
    int node = -1;          // Node
    int elem = -1;          // điểm nằm trên thanh (Mid, Nearest, Perpendicular, Intersection)
    double t = 0;           // vị trí trên thanh, 0..1
    int elem2 = -1;         // Intersection: thanh thứ hai
    double t2 = 0;
};

const char* snapName(SnapKind k);

// Lệnh đang thực hiện. Mỗi lệnh là một máy trạng thái nhận điểm, chữ, Enter.
struct Command {
    virtual ~Command() = default;
    std::string name;   // tên hiển thị (chữ hoa)
    bool done = false;

    virtual std::string prompt(App& a) = 0;
    virtual void start(App&) {}
    virtual bool wantsPoint() const { return false; }
    virtual bool wantsSelection() const { return false; } // canvas đang ở chế độ chọn đối tượng
    virtual void onPoint(App&, const SnapPoint&) {}
    virtual bool onText(App&, const std::string&) { return false; } // false: không nhận chữ này
    virtual void onEnter(App&) { done = true; }
    virtual void onCancel(App&) { done = true; }
    // điểm gốc cho tọa độ tương đối @, ORTHO, đường cao su
    virtual bool base(double&, double&) const { return false; }
    virtual void preview(App&, ImDrawList*, const SnapPoint&) {}
};

std::unique_ptr<Command> makeCommand(const std::string& upperName, std::string& canonical);
std::vector<std::pair<std::string, std::string>> commandHelp(); // (bí danh, mô tả)

// Phân tích chuỗi tọa độ: "x,y" | "@dx,dy" | "@d<góc" | "d<góc" | "d" (theo hướng con trỏ từ điểm gốc)
bool parsePoint(const std::string& s, bool hasBase, double bx, double by, double dirX, double dirY, double& x, double& y);

std::string upper(const std::string& s);
