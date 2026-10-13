// Đọc/ghi mô hình dạng văn bản (.fem). Chỉ số nút, phần tử, vật liệu... trong file đánh số từ 1.
#pragma once
#include <string>

#include "core/Analysis.h"
#include "core/Model.h"

namespace fem {

std::string modelToText(const Model& m);
bool modelFromText(const std::string& text, Model& m, std::string& error);
bool saveModel(const std::string& path, const Model& m, std::string& error);
bool loadModel(const std::string& path, Model& m, std::string& error);

// Báo cáo kết quả dạng văn bản (chuyển vị, phản lực, nội lực, ứng suất).
std::string resultReport(const Model& m, const AnalysisResult& r);

} // namespace fem
