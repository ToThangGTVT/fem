// Chuyển đổi mô hình và kết quả phân tích sang/từ JSON (lược đồ khớp với các lớp Kotlin trong mobile/shared).
// Tên trường giữ nguyên như trong core/Model.h; enum ghi bằng tên (ví dụ "Frame", "GlobalY"); chỉ số đánh từ 0.
#pragma once
#include <string>

#include "capi/Json.h"
#include "core/Analysis.h"
#include "core/Model.h"

namespace fem {

void writeModel(json::Writer& w, const Model& m);
bool readModel(const json::Value& v, Model& m, std::string& error);

void writeResult(json::Writer& w, const AnalysisResult& r);

} // namespace fem
