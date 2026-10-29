/*
 * Giao diện C của lõi tính FEM2D — dùng cho ứng dụng di động (Kotlin Multiplatform):
 *   Android: gọi qua JNI (mobile/native/fem_jni.cpp);  iOS: gọi trực tiếp qua Kotlin/Native cinterop.
 *
 * Mọi hàm nhận và trả chuỗi UTF-8. Kết quả luôn là JSON dạng {"ok": true, ...} hoặc {"ok": false, "error": "..."}.
 * Chuỗi trả về do thư viện cấp phát — người gọi phải giải phóng bằng fem_free().
 * Các hàm không giữ trạng thái toàn cục nên có thể gọi song song từ nhiều luồng.
 *
 * Lược đồ JSON của mô hình/kết quả: xem src/capi/ModelJson.cpp và mobile/shared/.../model.
 */
#ifndef FEM_CAPI_H
#define FEM_CAPI_H

#ifdef __cplusplus
extern "C" {
#endif

/* Phiên bản giao diện C; tăng khi lược đồ JSON thay đổi không tương thích. */
#define FEM_CAPI_VERSION 1

int fem_capi_version(void);

/* Ngôn ngữ của chuỗi do lõi sinh ra (thông báo lỗi, cảnh báo, báo cáo, tên bài mẫu): "vi" (mặc định), "en", "ja", "ko", "ru".
   Áp dụng cho mọi lời gọi sau đó, trên mọi luồng. */
void fem_set_language(const char* code);

/* {"ok":true,"examples":[{"name":..., "note":...}, ...]} */
char* fem_examples(void);

/* Mô hình của bài mẫu thứ index: {"ok":true,"model":{...}} */
char* fem_example_model(int index);

/* Đọc nội dung file .fem: {"ok":true,"model":{...}} */
char* fem_parse_fem(const char* femText);

/* Ghi mô hình (JSON) ra nội dung file .fem: {"ok":true,"text":"..."} */
char* fem_write_fem(const char* modelJson);

/* Phân tích tĩnh: {"ok":true,"result":{...}}. samplesPerBar ≤ 0: mặc định 40. */
char* fem_analyze(const char* modelJson, int samplesPerBar);

/* Báo cáo kết quả dạng văn bản: {"ok":true,"text":"..."} */
char* fem_report_text(const char* modelJson);

/* Báo cáo HTML có hình SVG: {"ok":true,"html":"..."} */
char* fem_report_html(const char* modelJson);

void fem_free(char* s);

#ifdef __cplusplus
}
#endif

#endif /* FEM_CAPI_H */
