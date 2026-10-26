#include "capi/fem_capi.h"

#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>

#include "capi/ModelJson.h"
#include "core/Examples.h"
#include "core/ModelIO.h"
#include "core/Report.h"

using namespace fem;

namespace {

char* dup(const std::string& s) {
    char* p = (char*)std::malloc(s.size() + 1);
    if (p) std::memcpy(p, s.c_str(), s.size() + 1);
    return p;
}

std::string errorStr(const std::string& msg) {
    json::Writer w;
    w.beginObject().field("ok", false).field("error", msg).endObject();
    return w.str();
}

// Bọc phần thân: ngoại lệ C++ không được vượt qua ranh giới C (JNI / Kotlin/Native)
char* guarded(const std::function<std::string()>& body) {
    try {
        return dup(body());
    } catch (const std::exception& e) {
        return dup(errorStr(std::string("Lỗi nội bộ: ") + e.what()));
    } catch (...) {
        return dup(errorStr("Lỗi nội bộ không xác định"));
    }
}

// Đọc mô hình từ JSON; trả về false kèm JSON lỗi trong 'err'
bool modelArg(const char* modelJson, Model& m, std::string& err) {
    if (!modelJson) {
        err = "Thiếu mô hình";
        return false;
    }
    json::Value v;
    return json::parse(modelJson, v, err) && readModel(v, m, err);
}

std::string okText(const char* key, const std::string& text) {
    json::Writer w;
    w.beginObject().field("ok", true).field(key, text).endObject();
    return w.str();
}

std::string okModel(const Model& m) {
    json::Writer w;
    w.beginObject().field("ok", true).key("model");
    writeModel(w, m);
    w.endObject();
    return w.str();
}

} // namespace

extern "C" {

int fem_capi_version(void) { return FEM_CAPI_VERSION; }

char* fem_examples(void) {
    return guarded([] {
        json::Writer w;
        w.beginObject().field("ok", true).key("examples").beginArray();
        for (const auto& e : examples()) w.beginObject().field("name", e.name).field("note", e.note).endObject();
        w.endArray().endObject();
        return w.str();
    });
}

char* fem_example_model(int index) {
    return guarded([index] {
        const auto& ex = examples();
        if (index < 0 || index >= (int)ex.size()) return errorStr("Không có bài mẫu số " + std::to_string(index));
        return okModel(ex[index].build());
    });
}

char* fem_parse_fem(const char* femText) {
    return guarded([femText] {
        Model m;
        std::string err;
        if (!femText || !modelFromText(femText, m, err)) return errorStr(femText ? err : "Thiếu nội dung file");
        return okModel(m);
    });
}

char* fem_write_fem(const char* modelJson) {
    return guarded([modelJson] {
        Model m;
        std::string err;
        if (!modelArg(modelJson, m, err)) return errorStr(err);
        return okText("text", modelToText(m));
    });
}

char* fem_analyze(const char* modelJson, int samplesPerBar) {
    return guarded([modelJson, samplesPerBar] {
        Model m;
        std::string err;
        if (!modelArg(modelJson, m, err)) return errorStr(err);
        const AnalysisResult r = analyze(m, samplesPerBar > 0 ? samplesPerBar : 40);
        json::Writer w;
        w.beginObject().field("ok", true).key("result");
        writeResult(w, r);
        w.endObject();
        return w.str();
    });
}

char* fem_report_text(const char* modelJson) {
    return guarded([modelJson] {
        Model m;
        std::string err;
        if (!modelArg(modelJson, m, err)) return errorStr(err);
        return okText("text", resultReport(m, analyze(m)));
    });
}

char* fem_report_html(const char* modelJson) {
    return guarded([modelJson] {
        Model m;
        std::string err;
        if (!modelArg(modelJson, m, err)) return errorStr(err);
        return okText("html", htmlReport(m, analyze(m)));
    });
}

void fem_free(char* s) { std::free(s); }

} // extern "C"
