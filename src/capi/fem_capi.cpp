#include "capi/fem_capi.h"

#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>

#include "capi/ModelJson.h"
#include "core/Codes.h"
#include "core/Examples.h"
#include "core/I18n.h"
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
        return dup(errorStr(trf("Lỗi nội bộ: %s", e.what())));
    } catch (...) {
        return dup(errorStr(tr("Lỗi nội bộ không xác định")));
    }
}

// Đọc mô hình từ JSON; trả về false kèm JSON lỗi trong 'err'
bool modelArg(const char* modelJson, Model& m, std::string& err) {
    if (!modelJson) {
        err = tr("Thiếu mô hình");
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

void fem_set_language(const char* code) { setLanguage(code); }

char* fem_examples(void) {
    return guarded([] {
        json::Writer w;
        w.beginObject().field("ok", true).key("examples").beginArray();
        for (const auto& e : examples()) w.beginObject().field("name", tr(e.name)).field("note", tr(e.note)).endObject();
        w.endArray().endObject();
        return w.str();
    });
}

char* fem_example_model(int index) {
    return guarded([index] {
        const auto& ex = examples();
        if (index < 0 || index >= (int)ex.size()) return errorStr(trf("Không có bài mẫu số %d", index));
        return okModel(ex[index].build());
    });
}

char* fem_parse_fem(const char* femText) {
    return guarded([femText] {
        Model m;
        std::string err;
        if (!femText || !modelFromText(femText, m, err)) return errorStr(femText ? err : tr("Thiếu nội dung file"));
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

char* fem_design_codes(void) {
    return guarded([] {
        json::Writer w;
        w.beginObject().field("ok", true).key("codes").beginArray();
        for (int i = 0; i < kDesignCodeCount; ++i) {
            const auto c = (DesignCode)i;
            w.beginObject().field("id", codeId(c)).field("name", codeName(c)).key("gammaF").beginArray();
            for (int k = CatDead; k <= CatSeismic; ++k) w.value(defaultGammaF(c, k));
            w.endArray().endObject();
        }
        w.endArray().endObject();
        return w.str();
    });
}

char* fem_library(const char* code) {
    return guarded([code] {
        DesignCode c;
        if (!code || !parseCode(code, c)) return errorStr(trf("Không có tiêu chuẩn \"%s\"", code ? code : ""));
        Model lib;
        for (const auto& p : codeMaterials(c)) lib.materials.push_back(toMaterial(p));
        for (const auto& p : codeProfiles(c)) lib.sections.push_back(toSection(p));
        json::Writer w;
        w.beginObject().field("ok", true).key("model");
        writeModel(w, lib);
        w.key("materialSources").beginArray();
        for (const auto& p : codeMaterials(c)) w.value(p.source);
        w.endArray().key("sectionSources").beginArray();
        for (const auto& p : codeProfiles(c)) w.value(p.source);
        w.endArray().endObject();
        return w.str();
    });
}

char* fem_generate_combos(const char* modelJson, const char* code) {
    return guarded([modelJson, code] {
        Model m;
        std::string err;
        if (!modelArg(modelJson, m, err)) return errorStr(err);
        DesignCode c;
        if (!code || !parseCode(code, c)) return errorStr(trf("Không có tiêu chuẩn \"%s\"", code ? code : ""));
        std::string note;
        const int n = generateCodeCombos(m, c, 300, &note);
        json::Writer w;
        w.beginObject().field("ok", true).key("model");
        writeModel(w, m);
        w.field("count", n).field("note", note).endObject();
        return w.str();
    });
}

void fem_free(char* s) { std::free(s); }

} // extern "C"
