#include "core/Model.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fem {

void Model::clear() { *this = Model{}; }

int Model::findNode(double x, double y, double tol) const {
    for (size_t i = 0; i < nodes.size(); ++i)
        if (std::abs(nodes[i].x - x) <= tol && std::abs(nodes[i].y - y) <= tol) return (int)i;
    return -1;
}

int Model::addNode(double x, double y, double tol) {
    int id = findNode(x, y, tol);
    if (id >= 0) return id;
    nodes.push_back({x, y});
    return (int)nodes.size() - 1;
}

int Model::findElement(ElemType type, const std::array<int, 4>& n) const {
    for (size_t i = 0; i < elements.size(); ++i) {
        const Element& e = elements[i];
        if (e.type != type) continue;
        int nc = e.nodeCount();
        std::array<int, 4> a = e.n, b = n;
        std::sort(a.begin(), a.begin() + nc);
        std::sort(b.begin(), b.begin() + nc);
        if (std::equal(a.begin(), a.begin() + nc, b.begin())) return (int)i;
    }
    return -1;
}

Support* Model::supportAt(int node) {
    for (auto& s : supports)
        if (s.node == node) return &s;
    return nullptr;
}

const Support* Model::supportAt(int node) const {
    for (auto& s : supports)
        if (s.node == node) return &s;
    return nullptr;
}

Support& Model::supportFor(int node) {
    if (Support* s = supportAt(node)) return *s;
    Support s;
    s.node = node;
    supports.push_back(s);
    return supports.back();
}

void Model::cleanSupports() {
    supports.erase(std::remove_if(supports.begin(), supports.end(), [](const Support& s) { return s.empty(); }),
                   supports.end());
}

void Model::removeNodes(std::vector<int> ids) {
    if (ids.empty()) return;
    std::vector<int> map(nodes.size(), 0);
    for (int id : ids)
        if (id >= 0 && id < (int)nodes.size()) map[id] = -1;
    // phần tử dùng nút bị xóa cũng bị xóa
    std::vector<int> deadElems;
    for (size_t e = 0; e < elements.size(); ++e)
        for (int k = 0; k < elements[e].nodeCount(); ++k)
            if (map[elements[e].n[k]] < 0) {
                deadElems.push_back((int)e);
                break;
            }
    removeElements(deadElems);

    int next = 0;
    std::vector<Node> kept;
    for (size_t i = 0; i < nodes.size(); ++i)
        if (map[i] >= 0) {
            map[i] = next++;
            kept.push_back(nodes[i]);
        }
    nodes.swap(kept);
    for (auto& e : elements)
        for (int k = 0; k < e.nodeCount(); ++k) e.n[k] = map[e.n[k]];

    auto remapNode = [&](int& n) { n = (n >= 0 && n < (int)map.size()) ? map[n] : -1; };
    for (auto& s : supports) remapNode(s.node);
    supports.erase(std::remove_if(supports.begin(), supports.end(), [](const Support& s) { return s.node < 0; }),
                   supports.end());
    for (auto& lc : cases) {
        for (auto& l : lc.nodal) remapNode(l.node);
        for (auto& d : lc.disp) remapNode(d.node);
        lc.nodal.erase(std::remove_if(lc.nodal.begin(), lc.nodal.end(), [](const NodalLoad& l) { return l.node < 0; }),
                       lc.nodal.end());
        lc.disp.erase(std::remove_if(lc.disp.begin(), lc.disp.end(), [](const PrescribedDisp& d) { return d.node < 0; }),
                      lc.disp.end());
    }
}

void Model::removeElements(std::vector<int> ids) {
    if (ids.empty()) return;
    std::vector<int> map(elements.size(), 0);
    for (int id : ids)
        if (id >= 0 && id < (int)elements.size()) map[id] = -1;
    int next = 0;
    std::vector<Element> kept;
    for (size_t i = 0; i < elements.size(); ++i)
        if (map[i] >= 0) {
            map[i] = next++;
            kept.push_back(elements[i]);
        }
    elements.swap(kept);
    for (auto& lc : cases) {
        for (auto& l : lc.elem) l.elem = (l.elem >= 0 && l.elem < (int)map.size()) ? map[l.elem] : -1;
        lc.elem.erase(std::remove_if(lc.elem.begin(), lc.elem.end(), [](const ElemLoad& l) { return l.elem < 0; }),
                      lc.elem.end());
    }
}

void Model::removeUnusedNodes() {
    std::vector<char> used(nodes.size(), 0);
    for (auto& e : elements)
        for (int k = 0; k < e.nodeCount(); ++k) used[e.n[k]] = 1;
    std::vector<int> dead;
    for (size_t i = 0; i < nodes.size(); ++i)
        if (!used[i]) dead.push_back((int)i);
    removeNodes(dead);
}

double Model::barLength(int e) const {
    const Element& el = elements[e];
    const Node& a = nodes[el.n[0]];
    const Node& b = nodes[el.n[1]];
    return std::hypot(b.x - a.x, b.y - a.y);
}

bool Model::hasPlane() const {
    for (auto& e : elements)
        if (e.isPlane()) return true;
    return false;
}

bool Model::hasBars() const {
    for (auto& e : elements)
        if (e.isBar()) return true;
    return false;
}

void Model::bounds(double& x0, double& y0, double& x1, double& y1) const {
    if (nodes.empty()) {
        x0 = y0 = 0;
        x1 = y1 = 10;
        return;
    }
    x0 = y0 = std::numeric_limits<double>::max();
    x1 = y1 = -std::numeric_limits<double>::max();
    for (auto& n : nodes) {
        x0 = std::min(x0, n.x);
        y0 = std::min(y0, n.y);
        x1 = std::max(x1, n.x);
        y1 = std::max(y1, n.y);
    }
}

const char* elemTypeName(ElemType t) {
    switch (t) {
    case ElemType::Truss: return "Thanh dàn";
    case ElemType::Frame: return "Thanh khung";
    case ElemType::Tri3: return "Tam giác CST";
    case ElemType::Quad4: return "Tứ giác Q4";
    }
    return "?";
}

const char* loadKindName(ElemLoadKind k) {
    switch (k) {
    case ElemLoadKind::Distributed: return "Phân bố";
    case ElemLoadKind::Point: return "Tập trung";
    case ElemLoadKind::Moment: return "Mômen";
    case ElemLoadKind::Temperature: return "Nhiệt độ";
    }
    return "?";
}

const char* loadDirName(LoadDir d) {
    switch (d) {
    case LoadDir::LocalX: return "x cục bộ";
    case LoadDir::LocalY: return "y cục bộ";
    case LoadDir::GlobalX: return "X tổng thể";
    case LoadDir::GlobalY: return "Y tổng thể";
    }
    return "?";
}

const char* dofName(int dof) {
    switch (dof) {
    case UX: return "ux";
    case UY: return "uy";
    case RZ: return "φz";
    }
    return "?";
}

} // namespace fem
