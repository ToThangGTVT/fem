#include "capi/ModelJson.h"
#include "core/I18n.h"

#include <type_traits>

namespace fem {

namespace {

template <class E, size_t N> const char* enumName(E e, const char* const (&names)[N]) {
    const size_t i = (size_t)e;
    return i < N ? names[i] : names[0];
}

template <class E, size_t N> E enumValue(const std::string& s, const char* const (&names)[N], E def) {
    for (size_t i = 0; i < N; ++i)
        if (s == names[i]) return (E)i;
    return def;
}

const char* const kElemType[] = {"Truss", "Frame", "Tri3", "Quad4", "Tri6", "Quad8"};
const char* const kPlaneMode[] = {"Stress", "Strain", "Axisymmetric"};
const char* const kShapeKind[] = {"None", "Rect", "Tee", "IBeam", "Box", "Pipe", "Circle", "Channel", "Angle", "Polygon"};
const char* const kLoadKind[] = {"Distributed", "Point", "Moment", "Temperature", "LackOfFit"};
const char* const kLoadDir[] = {"LocalX", "LocalY", "GlobalX", "GlobalY"};

template <class A> void readArray(const json::Value& v, const char* k, A& out) {
    const auto& arr = v.getArr(k);
    for (size_t i = 0; i < arr.size() && i < out.size(); ++i) {
        if constexpr (std::is_same_v<typename A::value_type, bool>) out[i] = arr[i].b;
        else out[i] = (typename A::value_type)arr[i].num;
    }
}

void writeShape(json::Writer& w, const SectionShape& s) {
    w.beginObject();
    w.field("kind", enumName(s.kind, kShapeKind));
    w.field("p", std::vector<double>(s.p.begin(), s.p.end()));
    w.key("contours").beginArray();
    for (const auto& c : s.contours) {
        w.beginArray();
        for (const auto& pt : c) w.beginArray().value(pt[0]).value(pt[1]).endArray();
        w.endArray();
    }
    w.endArray();
    w.endObject();
}

void readShape(const json::Value& v, SectionShape& s) {
    s.kind = enumValue(v.getStr("kind"), kShapeKind, ShapeKind::None);
    readArray(v, "p", s.p);
    for (const auto& c : v.getArr("contours")) {
        s.contours.emplace_back();
        for (const auto& pt : c.arr)
            if (pt.arr.size() >= 2) s.contours.back().push_back({pt.arr[0].num, pt.arr[1].num});
    }
}

} // namespace

void writeModel(json::Writer& w, const Model& m) {
    w.beginObject();
    w.field("title", m.title);
    w.field("planeMode", enumName(m.planeMode, kPlaneMode));
    w.field("secondOrder", m.secondOrder);

    w.key("nodes").beginArray();
    for (const auto& n : m.nodes) w.beginObject().field("x", n.x).field("y", n.y).endObject();
    w.endArray();

    w.key("materials").beginArray();
    for (const auto& x : m.materials)
        w.beginObject()
            .field("name", x.name)
            .field("E", x.E)
            .field("nu", x.nu)
            .field("alpha", x.alpha)
            .field("gamma", x.gamma)
            .field("kind", x.kind)
            .field("fd", x.fd)
            .field("ft", x.ft)
            .endObject();
    w.endArray();

    w.key("sections").beginArray();
    for (const auto& x : m.sections) {
        w.beginObject().field("name", x.name).field("A", x.A).field("I", x.I).field("h", x.h);
        w.field("yTop", x.yTop).field("yBot", x.yBot);
        w.key("shape");
        writeShape(w, x.shape);
        w.endObject();
    }
    w.endArray();

    w.key("elements").beginArray();
    for (const auto& e : m.elements) {
        w.beginObject();
        w.field("type", enumName(e.type, kElemType));
        w.key("n").array(e.n.begin(), e.n.begin() + e.nodeCount());
        w.field("mat", e.mat).field("sec", e.sec).field("t", e.t);
        w.field("hingeI", e.hingeI).field("hingeJ", e.hingeJ);
        w.field("kf", e.kf).field("kG", e.kG).field("kfCompOnly", e.kfCompOnly);
        w.field("incompatible", e.incompatible);
        w.endObject();
    }
    w.endArray();

    w.key("supports").beginArray();
    for (const auto& s : m.supports) {
        w.beginObject().field("node", s.node);
        w.key("fix").array(s.fix.begin(), s.fix.end());
        w.key("k").array(s.k.begin(), s.k.end());
        w.endObject();
    }
    w.endArray();

    w.key("cases").beginArray();
    for (const auto& lc : m.cases) {
        w.beginObject();
        w.field("name", lc.name).field("selfWeight", lc.selfWeight).field("kind", lc.kind);
        w.field("group", lc.group).field("gammaF", lc.gammaF).field("category", lc.category);
        w.key("nodal").beginArray();
        for (const auto& l : lc.nodal) {
            w.beginObject().field("node", l.node);
            w.key("F").array(l.F.begin(), l.F.end());
            w.endObject();
        }
        w.endArray();
        w.key("elem").beginArray();
        for (const auto& l : lc.elem)
            w.beginObject()
                .field("elem", l.elem)
                .field("kind", enumName(l.kind, kLoadKind))
                .field("dir", enumName(l.dir, kLoadDir))
                .field("q1", l.q1)
                .field("q2", l.q2)
                .field("a", l.a)
                .field("projected", l.projected)
                .endObject();
        w.endArray();
        w.key("disp").beginArray();
        for (const auto& d : lc.disp)
            w.beginObject().field("node", d.node).field("dof", d.dof).field("value", d.value).endObject();
        w.endArray();
        w.endObject();
    }
    w.endArray();

    w.key("combos").beginArray();
    for (const auto& c : m.combos) {
        w.beginObject().field("name", c.name).field("generated", c.generated);
        w.key("terms").beginArray();
        for (const auto& t : c.terms) w.beginObject().field("loadCase", t.first).field("factor", t.second).endObject();
        w.endArray();
        w.endObject();
    }
    w.endArray();
    w.endObject();
}

bool readModel(const json::Value& v, Model& out, std::string& error) {
    if (v.type != json::Value::Type::Object) {
        error = tr("Mô hình phải là đối tượng JSON");
        return false;
    }
    Model m;
    m.title = v.getStr("title");
    m.planeMode = enumValue(v.getStr("planeMode"), kPlaneMode, PlaneMode::Stress);
    m.secondOrder = v.getBool("secondOrder");

    for (const auto& n : v.getArr("nodes")) m.nodes.push_back({n.getNum("x"), n.getNum("y")});

    for (const auto& j : v.getArr("materials")) {
        Material x;
        x.name = j.getStr("name", x.name);
        x.E = j.getNum("E", x.E);
        x.nu = j.getNum("nu", x.nu);
        x.alpha = j.getNum("alpha", x.alpha);
        x.gamma = j.getNum("gamma", x.gamma);
        x.kind = j.getInt("kind", x.kind);
        x.fd = j.getNum("fd", x.fd);
        x.ft = j.getNum("ft", x.ft);
        m.materials.push_back(x);
    }

    for (const auto& j : v.getArr("sections")) {
        Section x;
        x.name = j.getStr("name", x.name);
        x.A = j.getNum("A", x.A);
        x.I = j.getNum("I", x.I);
        x.h = j.getNum("h", x.h);
        x.yTop = j.getNum("yTop", x.yTop);
        x.yBot = j.getNum("yBot", x.yBot);
        if (const json::Value* s = j.find("shape")) readShape(*s, x.shape);
        m.sections.push_back(x);
    }

    const int nn = (int)m.nodes.size();
    for (const auto& j : v.getArr("elements")) {
        Element e;
        e.type = enumValue(j.getStr("type"), kElemType, ElemType::Frame);
        const auto& n = j.getArr("n");
        if ((int)n.size() < e.nodeCount()) {
            error = trf("Phần tử %zu thiếu nút", m.elements.size());
            return false;
        }
        for (int k = 0; k < e.nodeCount(); ++k) {
            e.n[k] = (int)n[k].num;
            if (e.n[k] < 0 || e.n[k] >= nn) {
                error = trf("Phần tử %zu tham chiếu nút không tồn tại", m.elements.size());
                return false;
            }
        }
        e.mat = j.getInt("mat", e.mat);
        e.sec = j.getInt("sec", e.sec);
        e.t = j.getNum("t", e.t);
        e.hingeI = j.getBool("hingeI", e.hingeI);
        e.hingeJ = j.getBool("hingeJ", e.hingeJ);
        e.kf = j.getNum("kf", e.kf);
        e.kG = j.getNum("kG", e.kG);
        e.kfCompOnly = j.getBool("kfCompOnly", e.kfCompOnly);
        e.incompatible = j.getBool("incompatible", e.incompatible);
        if (e.mat < 0 || e.mat >= (int)m.materials.size() || (e.isBar() && (e.sec < 0 || e.sec >= (int)m.sections.size()))) {
            error = trf("Phần tử %zu tham chiếu vật liệu/tiết diện không tồn tại", m.elements.size());
            return false;
        }
        m.elements.push_back(e);
    }

    for (const auto& j : v.getArr("supports")) {
        Support s;
        s.node = j.getInt("node", -1);
        readArray(j, "fix", s.fix);
        readArray(j, "k", s.k);
        if (s.node >= 0 && s.node < nn) m.supports.push_back(s);
    }

    const int ne = (int)m.elements.size();
    for (const auto& j : v.getArr("cases")) {
        LoadCase lc;
        lc.name = j.getStr("name", lc.name);
        lc.selfWeight = j.getNum("selfWeight", lc.selfWeight);
        lc.kind = j.getInt("kind", lc.kind);
        lc.group = j.getInt("group", lc.group);
        lc.gammaF = j.getNum("gammaF", lc.gammaF);
        lc.category = j.getInt("category", lc.category);
        for (const auto& l : j.getArr("nodal")) {
            NodalLoad x;
            x.node = l.getInt("node", -1);
            readArray(l, "F", x.F);
            if (x.node >= 0 && x.node < nn) lc.nodal.push_back(x);
        }
        for (const auto& l : j.getArr("elem")) {
            ElemLoad x;
            x.elem = l.getInt("elem", -1);
            x.kind = enumValue(l.getStr("kind"), kLoadKind, ElemLoadKind::Distributed);
            x.dir = enumValue(l.getStr("dir"), kLoadDir, LoadDir::GlobalY);
            x.q1 = l.getNum("q1");
            x.q2 = l.getNum("q2");
            x.a = l.getNum("a");
            x.projected = l.getBool("projected");
            if (x.elem >= 0 && x.elem < ne) lc.elem.push_back(x);
        }
        for (const auto& l : j.getArr("disp")) {
            PrescribedDisp x;
            x.node = l.getInt("node", -1);
            x.dof = l.getInt("dof", UY);
            x.value = l.getNum("value");
            if (x.node >= 0 && x.node < nn && x.dof >= 0 && x.dof < 3) lc.disp.push_back(x);
        }
        m.cases.push_back(lc);
    }

    for (const auto& j : v.getArr("combos")) {
        Combo c;
        c.name = j.getStr("name", c.name);
        c.generated = j.getBool("generated");
        for (const auto& t : j.getArr("terms")) {
            const int lc = t.getInt("loadCase", -1);
            if (lc >= 0 && lc < (int)m.cases.size()) c.terms.push_back({lc, t.getNum("factor", 1.0)});
        }
        m.combos.push_back(c);
    }

    out = std::move(m);
    return true;
}

void writeResult(json::Writer& w, const AnalysisResult& r) {
    w.beginObject();
    w.field("ok", r.ok).field("error", r.error);
    w.field("warnings", r.warnings);
    w.field("numEquations", r.numEquations).field("profileSize", r.profileSize).field("timeMs", r.timeMs);
    w.key("dofState").beginArray();
    for (char c : r.dofState) w.value((int)c);
    w.endArray();
    w.key("sets").beginArray();
    for (const auto& s : r.sets) {
        w.beginObject();
        w.field("name", s.name).field("combo", s.combo);
        w.field("U", s.U).field("R", s.R);
        // thanh: lực đầu thanh và các điểm mẫu; phần tử phẳng: null
        w.key("bars").beginArray();
        for (const auto& b : s.bars) {
            if (b.x.empty()) {
                w.null();
                continue;
            }
            w.beginObject();
            w.key("endForces").array(b.endForces.begin(), b.endForces.end());
            w.key("uLocal").array(b.uLocal.begin(), b.uLocal.end());
            w.field("x", b.x).field("N", b.N).field("Q", b.Q).field("M", b.M).field("u", b.u).field("v", b.v);
            w.endObject();
        }
        w.endArray();
        // phần tử phẳng: (σx, σy, τxy, σ_ngoài) tại từng nút, nối liền 8 × 4 giá trị; thanh: mảng rỗng
        w.key("elemStress").beginArray();
        for (size_t e = 0; e < s.elemStress.size(); ++e) {
            const bool plane = e >= s.bars.size() || s.bars[e].x.empty();
            w.beginArray();
            if (plane)
                for (const auto& node : s.elemStress[e])
                    for (double x : node) w.value(x);
            w.endArray();
        }
        w.endArray();
        w.key("nodeStress").beginArray();
        for (const auto& ns : s.nodeStress) w.array(ns.begin(), ns.end());
        w.endArray();
        w.field("nodeStressCount", s.nodeStressCount);
        w.field("sumFx", s.sumFx).field("sumFy", s.sumFy).field("sumRx", s.sumRx).field("sumRy", s.sumRy);
        w.field("foundRx", s.foundRx).field("foundRy", s.foundRy);
        w.endObject();
    }
    w.endArray();
    w.endObject();
}

} // namespace fem
