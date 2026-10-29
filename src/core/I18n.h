// Đa ngôn ngữ cho chuỗi do lõi sinh ra (thông báo lỗi, cảnh báo, báo cáo, bài mẫu).
// Chữ gốc trong code là tiếng Việt, bọc bằng tr("...") hoặc trf("... %d", x); các ngôn ngữ khác tra bảng dịch
// theo đúng chữ gốc (I18nStrings.cpp, sinh từ tools/i18n/core_strings.py — sửa ở đó rồi chạy tools/i18n/gen_core.py).
// Mặc định tiếng Việt: chương trình không gọi setLanguage thì không có gì thay đổi.
#pragma once

#include <string>

namespace fem {

/** "vi", "en", "ja", "ko", "ru" (chỉ xét 2 ký tự đầu; mã khác → tiếng Việt). An toàn khi gọi từ nhiều luồng. */
void setLanguage(const char* code);
const char* languageCode();

/** Bản dịch của chữ gốc [vi] theo ngôn ngữ hiện tại; không có bản dịch thì trả lại chính [vi]. Con trỏ trỏ vào bảng tĩnh. */
const char* tr(const char* vi);

/** Như printf với định dạng [vi] đã dịch. Bản dịch giữ nguyên các %d, %s... theo đúng thứ tự (gen_core.py kiểm tra). */
std::string trf(const char* vi, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

}  // namespace fem
