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

static void testLanguage() {
    fem_set_language("en");
    json::Value ex = call(fem_examples());
    expect(!ex.getArr("examples").empty() && ex.getArr("examples")[0].getStr("name") == "Simply supported beam",
           "tiếng Anh: tên bài mẫu");
    json::Value e = call(fem_example_model(9999));
    expect(e.getStr("error") == "There is no example number 9999", "tiếng Anh: thông báo lỗi có tham số");
    fem_set_language("ja-JP");
    json::Value m = call(fem_example_model(0));
    const json::Value* model = m.find("model");
    expect(model && model->getStr("title").rfind("単純梁", 0) == 0, "tiếng Nhật: tiêu đề mô hình bài mẫu");
    fem_set_language("fr"); // chưa hỗ trợ: về tiếng Việt
    json::Value v = call(fem_example_model(9999));
    expect(v.getStr("error") == "Không có bài mẫu số 9999", "ngôn ngữ chưa hỗ trợ: tiếng Việt");
    fem_set_language("vi");
}

// Hệ số của trường hợp tải [lc] trong tổ hợp [c] (0 nếu không có)
static double comboFactor(const json::Value& c, int lc) {
    for (const auto& t : c.getArr("terms"))
        if ((int)t.getNum("loadCase") == lc) return t.getNum("factor");
    return 0;
}

static void testDesignCodes() {
    json::Value codes = call(fem_design_codes());
    expect(codes.getArr("codes").size() == 5, "5 tiêu chuẩn");
    json::Value lib = call(fem_library("ec"));
    const json::Value* lm = lib.find("model");
    bool ipe = false;
    if (lm)
        for (const auto& s : lm->getArr("sections"))
            if (s.getStr("name") == "IPE 300")
                ipe = std::abs(s.getNum("A") - 53.81e-4) < 0.05e-4 && std::abs(s.getNum("I") - 8356e-8) < 20e-8;
    expect(ipe, "thư viện EC: IPE 300 có A = 53,81 cm², I = 8356 cm⁴");
    expect(lm && lm->getArr("materials").size() == lib.getArr("materialSources").size(), "thư viện: nguồn theo từng vật liệu");
    expect(!call(fem_library("xx")).getBool("ok", true), "báo lỗi tiêu chuẩn không có");

    // G (γf 1,1), hoạt tải (1,2), gió trái, gió phải (1,2): hai tải gió không xuất hiện đồng thời
    const char* model = R"({"cases":[
        {"name":"G","kind":0,"gammaF":1.1,"category":0},
        {"name":"Q","kind":1,"gammaF":1.2,"category":1},
        {"name":"WL","kind":1,"gammaF":1.2,"category":3},
        {"name":"WR","kind":1,"gammaF":1.2,"category":3}]})";
    auto combos = [&](const char* code) {
        json::Value r = call(fem_generate_combos(model, code));
        const json::Value* m = r.find("model");
        return m ? m->getArr("combos") : std::vector<json::Value>{};
    };
    auto find = [](const std::vector<json::Value>& cs, std::initializer_list<std::pair<int, double>> want) {
        for (const auto& c : cs) {
            bool ok = c.getArr("terms").size() == want.size();
            for (const auto& [lc, f] : want) ok = ok && std::abs(comboFactor(c, lc) - f) < 1e-9;
            if (ok) return true;
        }
        return false;
    };
    const auto ec = combos("ec"); // {Q}, {WL}, {WR}, {Q, WL}×2, {Q, WR}×2
    expect(ec.size() == 7, "EN 1990: 7 tổ hợp, nhận " + std::to_string(ec.size()));
    expect(find(ec, {{0, 1.35}, {1, 1.5}, {2, 0.9}}), "EN 1990: 1,35G + 1,5Q + 1,5·0,6·W");
    expect(find(ec, {{0, 1.35}, {1, 1.05}, {2, 1.5}}), "EN 1990: 1,35G + 1,5·0,7·Q + 1,5W");
    const auto tc = combos("tcvn"); // cơ bản 1: 3; cơ bản 2: {Q, WL}, {Q, WR}
    expect(tc.size() == 5, "TCVN: 5 tổ hợp, nhận " + std::to_string(tc.size()));
    expect(find(tc, {{0, 1.1}, {1, 0.9 * 1.2}, {3, 0.9 * 1.2}}), "TCVN: G + 0,9(Q + W) nhân γf");
    const auto sp = combos("sp"); // 3 tổ hợp một tải + mỗi cặp hai thứ tự
    expect(sp.size() == 7, "СП 20: 7 tổ hợp, nhận " + std::to_string(sp.size()));
    expect(find(sp, {{0, 1.1}, {1, 1.2}, {2, 0.9 * 1.2}}), "СП 20: G + Q + 0,9W nhân γf");
    const auto kds = combos("kds"); // 1,4D; 1,2D + 1,6L; 1,2D + 1,3W + 1,0L ×2; 0,9D + 1,3W ×2
    expect(kds.size() == 6, "KDS: 6 tổ hợp, nhận " + std::to_string(kds.size()));
    expect(find(kds, {{0, 1.2}, {1, 1.0}, {3, 1.3}}), "KDS: 1,2D + 1,3W + 1,0L");
    expect(find(kds, {{0, 0.9}, {2, 1.3}}), "KDS: 0,9D + 1,3W");
    const auto jis = combos("jis"); // dài hạn G + P; ngắn hạn G + P + W ×2
    expect(jis.size() == 3, "JIS: 3 tổ hợp, nhận " + std::to_string(jis.size()));
    expect(find(jis, {{0, 1.0}, {1, 1.0}, {2, 1.0}}), "JIS: G + P + W");
}

int main() {
    testJson();
    testExamplesRoundTrip();
    testSimpleBeamValues();
    testErrors();
    testLanguage();
    testDesignCodes();
    std::printf("capi: %d đạt, %d lỗi\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
