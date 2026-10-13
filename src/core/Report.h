// Báo cáo kết quả dạng HTML (có hình SVG), mở bằng trình duyệt và in ra PDF.
#pragma once
#include <string>

#include "core/Analysis.h"

namespace fem {

struct HtmlReportOptions {
    bool diagramM = true, diagramQ = true, diagramN = true, deformed = true;
    bool tables = true;
    int maxRows = 400; // giới hạn số dòng mỗi bảng
};
std::string htmlReport(const Model& m, const AnalysisResult& r, const HtmlReportOptions& o = {});

} // namespace fem
