#include "core/I18n.h"

#include <atomic>
#include <cstdarg>
#include <cstdio>
#include <string_view>
#include <unordered_map>

namespace fem {

// Bảng dịch (I18nStrings.cpp): mỗi mục là chữ gốc tiếng Việt và bản dịch en, ja, ko, ru.
struct I18nEntry {
    const char* vi;
    const char* t[4];
};
extern const I18nEntry kI18nEntries[];
extern const int kI18nEntryCount;

namespace {

constexpr const char* kCodes[] = {"vi", "en", "ja", "ko", "ru"};
std::atomic<int> gLang{0};  // chỉ số trong kCodes

const std::unordered_map<std::string_view, const I18nEntry*>& index() {
    static const auto* map = [] {
        auto* m = new std::unordered_map<std::string_view, const I18nEntry*>();
        m->reserve((size_t)kI18nEntryCount);
        for (int i = 0; i < kI18nEntryCount; ++i) m->emplace(kI18nEntries[i].vi, &kI18nEntries[i]);
        return m;
    }();
    return *map;
}

}  // namespace

void setLanguage(const char* code) {
    int lang = 0;
    if (code && code[0] && code[1])
        for (int i = 0; i < 5; ++i)
            if ((code[0] | 0x20) == kCodes[i][0] && (code[1] | 0x20) == kCodes[i][1]) lang = i;
    gLang.store(lang, std::memory_order_relaxed);
}

const char* languageCode() { return kCodes[gLang.load(std::memory_order_relaxed)]; }

const char* tr(const char* vi) {
    const int lang = gLang.load(std::memory_order_relaxed);
    if (lang == 0 || !vi) return vi;
    const auto& map = index();
    const auto it = map.find(vi);
    if (it == map.end()) return vi;
    const char* t = it->second->t[lang - 1];
    return t && t[0] ? t : vi;
}

std::string trf(const char* vi, ...) {
    const char* f = tr(vi);
    va_list ap;
    va_start(ap, vi);
    va_list ap2;
    va_copy(ap2, ap);
    const int n = std::vsnprintf(nullptr, 0, f, ap);
    va_end(ap);
    std::string out;
    if (n > 0) {
        out.resize((size_t)n + 1);
        std::vsnprintf(&out[0], out.size(), f, ap2);
        out.resize((size_t)n);
    }
    va_end(ap2);
    return out;
}

}  // namespace fem
