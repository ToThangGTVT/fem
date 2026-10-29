// Tiêu chuẩn thiết kế theo nước: thư viện vật liệu, thép hình và quy tắc tự sinh tổ hợp tải.
//   TCVN (Việt Nam) · Eurocode (EN) · JIS / AIJ (Nhật Bản) · KDS / KS (Hàn Quốc) · СП / ГОСТ (Nga)
// Số liệu lấy theo các bảng của tiêu chuẩn tương ứng (ghi ở trường source); thép hình lưu kích thước h, b, tw, tf, r,
// diện tích và mômen quán tính tính từ kích thước (kể bo góc). Cần đối chiếu với bản tiêu chuẩn chính thức khi thiết kế.
#pragma once
#include <string>
#include <vector>

#include "core/Library.h"
#include "core/Model.h"

namespace fem {

enum class DesignCode { TCVN = 0, Eurocode, JIS, KDS, SP };
constexpr int kDesignCodeCount = 5;

const char* codeId(DesignCode c);    // "tcvn", "ec", "jis", "kds", "sp"
const char* codeName(DesignCode c);  // tên hiển thị (đã dịch)
bool parseCode(const std::string& id, DesignCode& out);

// Mác thép, bê tông của tiêu chuẩn. fd: cường độ dùng để thiết kế theo tiêu chuẩn đó (TCVN/СП: R; EC: fyd, fcd;
// JIS: F, Fc; KDS: fy, fck — KDS dùng hệ số giảm cường độ φ cho cấu kiện). ft: cường độ chịu kéo tính toán của bê tông
// nếu tiêu chuẩn định nghĩa (TCVN/СП: Rbt; EC: fctd), không thì 0.
const std::vector<MaterialPreset>& codeMaterials(DesignCode c);
// Thép hình của tiêu chuẩn (A, Ix đã tính từ kích thước với các dãy chữ I/H).
const std::vector<ProfilePreset>& codeProfiles(DesignCode c);

// Tiêu chuẩn tương ứng với ngôn ngữ giao diện ("vi" → TCVN, "en" → Eurocode, "ja" → JIS, "ko" → KDS, "ru" → СП).
DesignCode defaultCodeForLanguage(const std::string& lang);

// Sinh tổ hợp tải theo tiêu chuẩn (thay các tổ hợp tự sinh trước đó), dùng loại tải của từng trường hợp
// (loadCategory). Tải gió cùng loại, động đất cùng loại và tải cùng nhóm loại trừ ≠ 0 không xuất hiện đồng thời.
//   TCVN 2737 / 9386: cơ bản 1: G + Qi; cơ bản 2: G + 0,9·ΣQ; động đất: G + 0,3·Hoạt + E — nhân γf của trường hợp tải.
//   СП 20 / 14: G + Q1 + 0,9·Q2 + 0,7·ΣQ (mỗi tải lần lượt làm Q1, Q2); động đất: 0,9G + 0,8·Hoạt + 0,5·Tuyết + E — nhân γf.
//   EN 1990 (6.10, 6.12b): 1,35G + 1,5Q1 + Σ1,5ψ0·Qi (ψ0: hoạt 0,7, tuyết 0,5, gió 0,6); động đất: G + E + 0,3·Hoạt.
//   KDS 41 10 15: 1,4D; 1,2D + 1,6L + 0,5S; 1,2D + 1,6S + (1,0L | 0,65W); 1,2D + 1,3W + 1,0L + 0,5S;
//                 1,2D + 1,0E + 1,0L + 0,2S; 0,9D + 1,3W; 0,9D + 1,0E.
//   Luật Tiêu chuẩn Xây dựng Nhật (ứng suất cho phép): dài hạn G + P; ngắn hạn G + P + S | W | K.
// Với EN, KDS, JIS hệ số do tiêu chuẩn quy định (bỏ qua γf); chỉ xét trường hợp bất lợi (γG = 1,35 với EN).
// Trả về số tổ hợp đã sinh; vượt maxCombos thì cắt bớt và ghi chú vào note.
int generateCodeCombos(Model& m, DesignCode c, int maxCombos = 300, std::string* note = nullptr);

// γf mặc định của loại tải theo TCVN 2737 / СП 20 (dùng khi người dùng chọn loại tải); EN/KDS/JIS: 1.
double defaultGammaF(DesignCode c, int category);

} // namespace fem
