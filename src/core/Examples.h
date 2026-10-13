// Các bài toán mẫu theo giáo trình Cơ học kết cấu / Lý thuyết đàn hồi.
#pragma once
#include <vector>

#include "core/Model.h"

namespace fem {

struct Example {
    const char* name;
    const char* note;
    Model (*build)();
};

const std::vector<Example>& examples();

} // namespace fem
