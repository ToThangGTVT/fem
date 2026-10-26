// Kiểm tra giao diện C (fem_capi) dùng cho ứng dụng di động: JSON khứ hồi, phân tích, xử lý lỗi.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "capi/Json.h"
#include "capi/ModelJson.h"
#include "capi/fem_capi.h"
#include "core/Examples.h"
#include "core/ModelIO.h"

using namespace fem;

static int g_fail = 0, g_pass = 0;

static void expect(bool ok, const std::string& what) {
    if (ok) ++g_pass;
    else {
        ++g_fail;
        std::printf("  FAIL %s\n", what.c_str());
    }
}

// Gọi hàm C, đọc JSON trả về và giải phóng chuỗi
static json::Value call(char* raw) {
    json::Value v;
    std::string err;
    const bool parsed = raw && json::parse(raw, v, err);
    expect(parsed, "JSON trả về hợp lệ: " + err);
    fem_free(raw);
    return v;
}

static void testJson() {
    json::Value v;
    std::string err;
    expect(json::parse(R"({"a":[1,-2.5e3,true,null],"s":"x\"ạ\n","n":NaN})", v, err), "parse: " + err);
    expect(v.getArr("a").size() == 4 && v.getArr("a")[1].num == -2500, "mảng số");
    expect(v.getStr("s") == "x\"ạ\n", "chuỗi thoát + UTF-8");
    expect(!json::parse("{\"a\":1,}", v, err), "từ chối JSON sai");

    json::Writer w;
    w.beginObject().field("t", std::string("Dầm \"A\"")).key("v").beginArray().value(1.5).value(2).endArray().endObject();
    expect(w.str() == R"({"t":"Dầm \"A\"","v":[1.5,2]})", "ghi JSON: " + w.str());
}

// Mọi bài mẫu: mô hình → JSON → mô hình phải giữ nguyên nội dung .fem; phân tích qua C API phải thành công
static void testExamplesRoundTrip() {
    json::Value list = call(fem_examples());
    const auto& ex = list.getArr("examples");
    expect(list.getBool("ok") && ex.size() == examples().size(), "danh sách bài mẫu");

    for (size_t i = 0; i < ex.size(); ++i) {
        const std::string name = ex[i].getStr("name");
        json::Value r = call(fem_example_model((int)i));
        expect(r.getBool("ok"), name + ": lấy mô hình");
        const json::Value* mj = r.find("model");
        Model back;
        std::string err;
        expect(mj && readModel(*mj, back, err), name + ": đọc lại JSON " + err);
        expect(modelToText(back) == modelToText(examples()[i].build()), name + ": khứ hồi JSON giữ nguyên mô hình");

        json::Writer w;
        writeModel(w, back);
        json::Value a = call(fem_analyze(w.str().c_str(), 0));
        const json::Value* res = a.find("result");
        expect(a.getBool("ok") && res && res->getBool("ok"), name + ": phân tích " + (res ? res->getStr("error") : ""));
        if (res) {
            const auto& sets = res->getArr("sets");
            expect(sets.size() == back.cases.size() + back.combos.size(), name + ": số bộ kết quả");
            if (!sets.empty()) {
                expect(sets[0].getArr("U").size() == 3 * back.nodes.size(), name + ": số chuyển vị");
                expect(sets[0].getArr("bars").size() == back.elements.size(), name + ": số phần tử kết quả");
            }
        }

        json::Value f = call(fem_write_fem(w.str().c_str()));
        json::Value p = call(fem_parse_fem(f.getStr("text").c_str()));
        expect(f.getBool("ok") && p.getBool("ok"), name + ": .fem qua C API");
    }
}

static void testSimpleBeamValues() {
    // Dầm đơn giản L = 6 m, q = 10 kN/m: Mmax = qL²/8 = 45 kNm
    const char* model = R"({
      "nodes":[{"x":0,"y":0},{"x":6,"y":0}],
      "materials":[{"name":"Thép","E":2.1e8}], "sections":[{"name":"S","A":0.01,"I":1e-4,"h":0.3}],
      "elements":[{"type":"Frame","n":[0,1]}],
      "supports":[{"node":0,"fix":[true,true,false]},{"node":1,"fix":[false,true,false]}],
      "cases":[{"name":"q","elem":[{"elem":0,"kind":"Distributed","dir":"GlobalY","q1":-10,"q2":-10}]}]
    })";
    json::Value a = call(fem_analyze(model, 40));
    const json::Value* res = a.find("result");
    expect(res && res->getBool("ok"), "dầm đơn giản: phân tích");
    if (!res || res->getArr("sets").empty()) return;
    const json::Value& bar = res->getArr("sets")[0].getArr("bars")[0];
    double mmax = 0;
    for (const auto& v : bar.getArr("M")) mmax = std::max(mmax, v.num);
    expect(std::abs(mmax - 45.0) < 1e-6, "dầm đơn giản: Mmax = 45, nhận " + std::to_string(mmax));
}

static void testErrors() {
    json::Value a = call(fem_analyze("{not json", 0));
    expect(!a.getBool("ok", true) && !a.getStr("error").empty(), "báo lỗi JSON hỏng");
    json::Value b = call(fem_analyze(R"({"nodes":[{"x":0,"y":0}],"elements":[{"type":"Frame","n":[0,5]}]})", 0));
    expect(!b.getBool("ok", true), "báo lỗi nút không tồn tại");
    json::Value c = call(fem_example_model(9999));
    expect(!c.getBool("ok", true), "báo lỗi chỉ số bài mẫu");
    json::Value d = call(fem_analyze(nullptr, 0));
    expect(!d.getBool("ok", true), "báo lỗi con trỏ rỗng");
}

int main() {
    testJson();
    testExamplesRoundTrip();
    testSimpleBeamValues();
    testErrors();
    std::printf("capi: %d đạt, %d lỗi\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
