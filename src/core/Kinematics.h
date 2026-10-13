// Phân tích cấu tạo hình học hệ thanh phẳng (Cơ học kết cấu, chương 1):
// bậc tự do W theo công thức Чебышёв, bậc siêu tĩnh, kết luận bất biến hình / biến hình / tức thời biến hình.
#pragma once
#include <string>
#include <vector>

#include "core/Model.h"

namespace fem {

struct KinematicResult {
    bool applicable = false;
    std::string note;
    // công thức W = 3D − 3H − 2K − C₀
    int disks = 0;       // D: số miếng cứng (thanh)
    int rigidJoints = 0; // H: số mối hàn đơn
    int hinges = 0;      // K: số khớp đơn
    int supports = 0;    // C₀: số liên kết tựa
    int W = 0;
    // phân tích đại số ma trận cân bằng
    int equations = 0, unknowns = 0, rank = 0;
    int indeterminacy = 0; // bậc siêu tĩnh n = ẩn − hạng
    int mechanisms = 0;    // số bậc tự do biến hình = phương trình − hạng
    std::string verdict;
    std::vector<double> mechanismU; // dạng chuyển động của cơ cấu (3 bậc tự do/nút), rỗng nếu bất biến hình
};

KinematicResult kinematicAnalysis(const Model& m);

} // namespace fem
