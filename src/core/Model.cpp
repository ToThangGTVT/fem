#include "core/Model.h"
#include "core/I18n.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <bitset>
#include <cstdio>

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

int Model::findElement(ElemType type, const std::array<int, 8>& n) const {
    for (size_t i = 0; i < elements.size(); ++i) {
        const Element& e = elements[i];
        if (e.type != type) continue;
        int nc = e.nodeCount();
        std::array<int, 8> a = e.n, b = n;
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

double FlowDomain::width() const { return std::hypot(bx - ax, by - ay); }

void FlowDomain::direction(double& nx, double& ny) const {
    const double w = width();
    const double tx = w > 0 ? (bx - ax) / w : 0, ty = w > 0 ? (by - ay) / w : 1;
    nx = -ty * side; // pháp tuyến bên trái của A→B nhân với chiều chảy
    ny = tx * side;
}

void FlowDomain::corners(double c[4][2]) const {
    double nx, ny;
    direction(nx, ny);
    c[0][0] = ax, c[0][1] = ay;
    c[1][0] = bx, c[1][1] = by;
    c[2][0] = bx + length * nx, c[2][1] = by + length * ny;
    c[3][0] = ax + length * nx, c[3][1] = ay + length * ny;
}

bool FlowDomain::groundAtB() const {
    double c[4][2];
    corners(c);
    return c[1][1] + c[2][1] < c[0][1] + c[3][1];
}

double FlowDomain::waterLevelY(double depth) const {
    const double w = width();
    const bool gB = groundAtB();
    const double gy = gB ? by : ay, ty = gB ? ay : by;
    return w > 0 ? gy + (ty - gy) / w * depth : gy;
}

void Model::bounds(double& x0, double& y0, double& x1, double& y1) const {
    if (flow.active) {
        double c[4][2];
        flow.corners(c);
        x0 = x1 = c[0][0];
        y0 = y1 = c[0][1];
        for (auto& p : c) {
            x0 = std::min(x0, p[0]);
            y0 = std::min(y0, p[1]);
            x1 = std::max(x1, p[0]);
            y1 = std::max(y1, p[1]);
        }
        for (auto& n : nodes) {
            x0 = std::min(x0, n.x);
            y0 = std::min(y0, n.y);
            x1 = std::max(x1, n.x);
            y1 = std::max(y1, n.y);
        }
        return;
    }
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

int generateCombos(Model& m, double psi1, double psi2, int maxCombos, std::string* note) {
    m.combos.erase(std::remove_if(m.combos.begin(), m.combos.end(), [](const Combo& c) { return c.generated; }),
                   m.combos.end());
    std::vector<int> perm, var;
    for (size_t i = 0; i < m.cases.size(); ++i) (m.cases[i].kind == 0 ? perm : var).push_back((int)i);
    auto fmtF = [](double v) {
        char b[32];
        std::snprintf(b, sizeof b, "%g", v);
        std::string s = b;
        for (char& c : s)
            if (c == '.') c = ',';
        return s;
    };
    std::string permName;
    for (int i : perm) permName += (permName.empty() ? "" : " + ") + m.cases[i].name;
    auto make = [&](const std::vector<int>& vs, double psi, const std::string& prefix) {
        Combo c;
        c.generated = true;
        for (int i : perm) c.terms.push_back({i, m.cases[i].gammaF});
        std::string vn;
        for (int i : vs) {
            c.terms.push_back({i, psi * m.cases[i].gammaF});
            vn += (vn.empty() ? "" : " + ") + m.cases[i].name;
        }
        std::string name = prefix + ": " + permName;
        if (!vs.empty()) {
            if (!permName.empty()) name += " + ";
            name += (psi != 1.0 && vs.size() > 1) ? fmtF(psi) + "(" + vn + ")" : (psi != 1.0 ? fmtF(psi) + "·" : "") + vn;
        }
        c.name = name;
        m.combos.push_back(c);
    };
    int count = 0;
    if (var.empty()) {
        if (!perm.empty()) {
            make({}, 1.0, "CB");
            ++count;
        }
        return count;
    }
    for (int i : var) {
        make({i}, psi1, "CB1");
        ++count;
    }
    const int nv = (int)var.size();
    if (nv > 16) {
        if (note) *note = tr("Quá nhiều tải tạm thời (> 16): chỉ sinh tổ hợp cơ bản 1.");
        return count;
    }
    std::vector<unsigned> masks;
    for (unsigned mask = 0; mask < (1u << nv); ++mask) {
        int pc = 0;
        bool ok = true;
        std::vector<int> groups;
        for (int k = 0; k < nv; ++k)
            if (mask & (1u << k)) {
                ++pc;
                const int g = m.cases[var[k]].group;
                if (g != 0) {
                    if (std::find(groups.begin(), groups.end(), g) != groups.end()) ok = false;
                    groups.push_back(g);
                }
            }
        if (ok && pc >= 2) masks.push_back(mask);
    }
    std::stable_sort(masks.begin(), masks.end(),
                     [](unsigned a, unsigned b) { return std::bitset<32>(a).count() < std::bitset<32>(b).count(); });
    bool truncated = false;
    for (unsigned mask : masks) {
        if (count >= maxCombos) {
            truncated = true;
            break;
        }
        std::vector<int> vs;
        for (int k = 0; k < nv; ++k)
            if (mask & (1u << k)) vs.push_back(var[k]);
        make(vs, psi2, "CB2");
        ++count;
    }
    if (truncated && note) *note = trf("Số tổ hợp vượt giới hạn %d — đã cắt bớt.", (int)maxCombos);
    return count;
}

const char* elemTypeName(ElemType t) {
    switch (t) {
    case ElemType::Truss: return tr("Thanh dàn");
    case ElemType::Frame: return tr("Thanh khung");
    case ElemType::Tri3: return tr("Tam giác CST");
    case ElemType::Quad4: return tr("Tứ giác Q4");
    case ElemType::Tri6: return tr("Tam giác T6");
    case ElemType::Quad8: return tr("Tứ giác Q8");
    }
    return "?";
}

const char* loadKindName(ElemLoadKind k) {
    switch (k) {
    case ElemLoadKind::Distributed: return tr("Phân bố");
    case ElemLoadKind::Point: return tr("Tập trung");
    case ElemLoadKind::Moment: return tr("Mômen");
    case ElemLoadKind::Temperature: return tr("Nhiệt độ");
    case ElemLoadKind::LackOfFit: return tr("Chế tạo ΔL");
    }
    return "?";
}

const char* loadDirName(LoadDir d) {
    switch (d) {
    case LoadDir::LocalX: return tr("x cục bộ");
    case LoadDir::LocalY: return tr("y cục bộ");
    case LoadDir::GlobalX: return tr("X tổng thể");
    case LoadDir::GlobalY: return tr("Y tổng thể");
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
