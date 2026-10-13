// Cửa sổ "Hướng dẫn": cách thao tác, danh sách lệnh, phím tắt, chức năng và quy ước — gom về một chỗ
// để các panel làm việc chỉ còn nội dung chính.
#include <string>

#include "App.h"
#include "Commands.h"
#include "imgui.h"
#include "imgui_internal.h"

namespace {

void bullets(std::initializer_list<const char*> items) {
    for (const char* t : items) {
        ImGui::Bullet();
        ImGui::TextWrapped("%s", t);
    }
}

// Bảng hai cột: phím / lệnh — mô tả
void keyTable(const char* id, std::initializer_list<std::pair<const char*, const char*>> rows, float keyWidth) {
    if (!ImGui::BeginTable(id, 2, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV)) return;
    ImGui::TableSetupColumn("k", ImGuiTableColumnFlags_WidthFixed, keyWidth);
    ImGui::TableSetupColumn("d", ImGuiTableColumnFlags_WidthStretch);
    for (const auto& [k, d] : rows) {
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::TextUnformatted(k);
        ImGui::TableNextColumn();
        ImGui::TextWrapped("%s", d);
    }
    ImGui::EndTable();
}

} // namespace

void App::guidePanel() {
    // Lần đầu xuất hiện (kể cả khi đã có bố cục lưu từ phiên bản cũ): ghép thành thẻ cạnh panel Thuộc tính
    if (ImGuiWindow* props = ImGui::FindWindowByName("Thuộc tính"); props && props->DockId)
        ImGui::SetNextWindowDockID(props->DockId, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(520, 560), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Hướng dẫn", &showGuide)) {
        ImGui::End();
        return;
    }
    if (!ImGui::BeginTabBar("guide", ImGuiTabBarFlags_FittingPolicyScroll)) {
        ImGui::End();
        return;
    }
    const float em = ImGui::GetFontSize();

    if (ImGui::BeginTabItem("Bắt đầu")) {
        ImGui::TextWrapped("FEM2D tính kết cấu phẳng bằng phương pháp phần tử hữu hạn (phương pháp chuyển vị).");
        ImGui::SeparatorText("Các bước");
        bullets({
            "1. Mở một ví dụ (menu Tệp) hoặc dựng sơ đồ: vẽ thanh (L, DAN) hay dùng Tạo khung / Tạo dàn / Tạo vòm / "
            "Lưới phẳng ở thẻ Trang chủ.",
            "2. Khai báo vật liệu, tiết diện trong panel Mô hình; chọn thanh rồi \"Gán cho thanh chọn\".",
            "3. Đặt gối: chọn nút rồi bấm Ngàm / Gối cố định / Gối di động, hoặc gõ NGAM, KHOP, DIDONG.",
            "4. Đặt tải ở thẻ Tải trọng (TAINUT, TAIPB, TAITT, MOMEN). Nhiệt độ, lún gối: chọn đối tượng rồi "
            "sửa trong panel Thuộc tính.",
            "5. Bấm Tính (F5). Xem biểu đồ N, Q, M, chuyển vị, phổ ứng suất ở thẻ Kết quả.",
            "6. Xuất báo cáo HTML (BAOCAO) hoặc bản vẽ DXF (DXFOUT).",
        });
        ImGui::SeparatorText("Các panel");
        keyTable("panels",
                 {{"Thuộc tính", "Sửa nút / phần tử đang chọn: tọa độ, gối, vật liệu, tiết diện, tải. Ctrl+1."},
                  {"Mô hình", "Bảng vật liệu, tiết diện, nút, phần tử, gối, trường hợp tải, tổ hợp."},
                  {"Kết quả", "Chọn bộ kết quả, đại lượng hiển thị, bảng nội lực, đường ảnh hưởng."},
                  {"Nhật ký", "Thông báo của chương trình và kết quả các lệnh đo, liệt kê."}},
                 7 * em);
        ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Vẽ & chọn")) {
        ImGui::SeparatorText("Dòng lệnh");
        bullets({
            "Gõ lệnh ở bất cứ đâu, không cần nhấp vào dòng lệnh; không phân biệt hoa thường. Danh sách ở thẻ Lệnh.",
            "Enter / Space / chuột phải: kết thúc lệnh, hoặc lặp lại lệnh vừa dùng. Esc: hủy.",
            "Mũi tên lên/xuống: gọi lại các dòng đã gõ. F2: mở rộng lịch sử dòng lệnh.",
        });
        ImGui::SeparatorText("Nhập tọa độ");
        keyTable("coords",
                 {{"x,y", "tọa độ tuyệt đối (m)"},
                  {"@dx,dy", "tương đối so với điểm trước"},
                  {"@d<góc", "tương đối theo khoảng cách và góc (độ, ngược chiều kim đồng hồ từ trục X)"},
                  {"d", "chỉ gõ khoảng cách: đặt điểm theo hướng con trỏ chuột"}},
                 6 * em);
        ImGui::SeparatorText("Chọn đối tượng");
        bullets({
            "Nhấp để thêm vào tập chọn, Shift+nhấp để bỏ ra.",
            "Kéo từ trái sang phải: chọn kiểu cửa sổ (nằm trọn trong khung).",
            "Kéo từ phải sang trái: chọn kiểu cắt ngang (chạm khung là chọn).",
            "Nhấp grip (ô vuông xanh) để kéo nút.",
        });
        ImGui::SeparatorText("Xem sơ đồ");
        bullets({
            "Lăn chuột: phóng to / thu nhỏ quanh con trỏ.",
            "Kéo chuột giữa hoặc chuột phải: dịch khung nhìn.",
            "Nhấp đúp chuột giữa (hoặc Z rồi E): vừa khung.",
        });
        ImGui::SeparatorText("Trợ giúp vẽ");
        keyTable("aids",
                 {{"F3 OSNAP", "bắt điểm: nút, trung điểm, giao điểm, chân vuông góc, điểm gần nhất trên thanh"},
                  {"F7", "hiện / ẩn lưới"},
                  {"F8 ORTHO", "chỉ vẽ theo phương ngang, đứng"},
                  {"F9 SNAP", "bắt theo bước lưới (đặt bước lưới ở panel Mô hình, thẻ Chung)"},
                  {"F12 DYN", "ô nhập số cạnh con trỏ"}},
                 6 * em);
        ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Lệnh")) {
        static char filter[64] = "";
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::InputTextWithHint("##filter", "Lọc theo tên lệnh hoặc mô tả...", filter, sizeof filter);
        const std::string f = upper(filter);
        if (ImGui::BeginTable("cmds", 2, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_ScrollY)) {
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableSetupColumn("Lệnh", ImGuiTableColumnFlags_WidthFixed, 11 * em);
            ImGui::TableSetupColumn("Chức năng", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableHeadersRow();
            for (const auto& [k, d] : commandHelp()) {
                if (!f.empty() && upper(k).find(f) == std::string::npos && upper(d).find(f) == std::string::npos) continue;
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TextWrapped("%s", k.c_str());
                ImGui::TableNextColumn();
                ImGui::TextWrapped("%s", d.c_str());
            }
            ImGui::EndTable();
        }
        ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Phím tắt")) {
        keyTable("keys",
                 {{"F1", "mở / đóng Hướng dẫn"},
                  {"F2", "mở rộng lịch sử dòng lệnh"},
                  {"F3 / F8 / F9 / F12", "OSNAP / ORTHO / SNAP / DYN"},
                  {"F5", "tính"},
                  {"F7", "lưới"},
                  {"Ctrl+1", "panel Thuộc tính"},
                  {"Ctrl+Z / Ctrl+Y", "hoàn tác / làm lại"},
                  {"Ctrl+N / O / S", "tệp mới / mở / lưu"},
                  {"Delete", "xóa đối tượng chọn"},
                  {"Esc", "xóa dòng đang gõ, hủy lệnh, rồi bỏ chọn"},
                  {"Enter / Space", "kết thúc hoặc lặp lại lệnh"}},
                 9 * em);
        ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Chức năng")) {
        if (ImGui::CollapsingHeader("Trang chủ — dựng sơ đồ", ImGuiTreeNodeFlags_DefaultOpen))
            bullets({
                "Vẽ: thanh khung, thanh dàn, nút; tạo nhanh khung nhiều tầng / dầm liên tục, dàn Pratt/Howe/Warren, "
                "vòm, cọc chịu tải ngang (phương pháp m), lưới phần tử phẳng.",
                "Sửa: di chuyển, sao chép, xoay, đối xứng, chia thanh, gộp nút trùng, xóa.",
                "Gối tựa: ngàm, gối cố định, gối di động (đứng/ngang), ngàm trượt; lò xo và lún gối ở panel Thuộc tính.",
                "Thuộc tính thanh: gán vật liệu, tiết diện; thiết kế tiết diện (đặc trưng hình học, lõi).",
                "Tiện ích: đo khoảng cách, tọa độ, liệt kê đối tượng.",
            });
        if (ImGui::CollapsingHeader("Tải trọng"))
            bullets({
                "Nhiều trường hợp tải; chọn trường hợp đang soạn ở thanh ribbon.",
                "Tải nút, tải phân bố đều / hình thang, lực và mômen tập trung trên thanh.",
                "Nhiệt độ, sai số chế tạo ΔL, chuyển vị cưỡng bức, trọng lượng bản thân.",
                "Tổ hợp: sinh tổ hợp cơ bản 1, 2 theo TCVN 2737:1995 hoặc tự lập ở panel Mô hình, thẻ Tổ hợp.",
            });
        if (ImGui::CollapsingHeader("Phân tích"))
            bullets({
                "Tính (F5): tĩnh tuyến tính; bật P-Δ (lệnh PDELTA) để tính theo sơ đồ biến dạng.",
                "Cấu tạo hình học (CTHH): W, bậc siêu tĩnh, kết luận bất biến hình, dạng biến hình.",
                "Giáo trình: phương pháp lực từng bước (PPL), phương pháp chuyển vị dạng ma trận (PPCV).",
                "Ổn định: λcr, dạng mất ổn định, hệ số chiều dài tính toán μ. Dao động riêng: tần số, chu kỳ, dạng.",
                "Động lực học (DONGLUC): phổ phản ứng TCVN 9386, cưỡng bức điều hòa, tích phân Newmark.",
                "Đường ảnh hưởng: thiết lập trong panel Kết quả; xếp đoàn xe bất lợi.",
            });
        if (ImGui::CollapsingHeader("Kết quả"))
            bullets({
                "Biểu đồ N, Q, M (M vẽ về phía thớ căng), hình biến dạng, phản lực, phổ màu ứng suất phần tử phẳng.",
                "Biểu đồ bao qua các tổ hợp; chọn bộ kết quả (trường hợp tải / tổ hợp) ở panel Kết quả.",
                "Kiểm tra thiết kế (KTTK): thép TCVN 5575, BTCT TCVN 5574, độ võng.",
                "Ứng suất tiết diện (USTD): σ, τ, vòng Mohr, thuyết bền.",
                "Xuất báo cáo HTML có hình (BAOCAO), bản vẽ DXF (DXFOUT); nhập DXF (DXFIN).",
            });
        if (ImGui::CollapsingHeader("Xem"))
            bullets({
                "Bật / tắt các panel, khôi phục bố cục mặc định.",
                "Hiện số hiệu nút, phần tử, tải trọng, gối tựa, lưới, giá trị trên biểu đồ.",
                "Giao diện tối / sáng.",
            });
        ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Quy ước")) {
        ImGui::SeparatorText("Đơn vị");
        ImGui::BulletText("Lực kN, chiều dài m, mômen kNm, E kN/m², ứng suất hiển thị MPa, chuyển vị mm.");
        ImGui::BulletText("Lò xo gối: kN/m, kNm/rad.");
        ImGui::SeparatorText("Hệ tọa độ & quy ước dấu");
        bullets({
            "Trục X sang phải, Y hướng lên; mômen, góc xoay dương ngược chiều kim đồng hồ.",
            "Trục cục bộ thanh: x từ nút i đến j, y quay x ngược chiều kim đồng hồ 90°.",
            "N > 0: kéo. Q > 0: làm phần thanh quay thuận chiều kim đồng hồ.",
            "M > 0: căng thớ dưới (phía −y cục bộ). Biểu đồ M vẽ về phía thớ căng.",
        });
        ImGui::SeparatorText("Phần tử");
        bullets({
            "Thanh dàn: chỉ chịu kéo nén. Tải trong nhịp quy về nút như dầm đơn giản.",
            "Thanh khung Euler–Bernoulli, khớp đầu thanh, nền Winkler, tải nhiệt độ.",
            "Tam giác CST, T6; tứ giác Q4 (tùy chọn QM6), Q8. Ứng suất phẳng, biến dạng phẳng, đối xứng trục.",
        });
        ImGui::SeparatorText("Tài liệu tham khảo");
        bullets({
            "Lều Thọ Trình — Cơ học kết cấu (tập 1, 2).",
            "Дарков А.В., Шапошников Н.Н. — Строительная механика.",
            "Смирнов А.Ф. и др. — Строительная механика. Стержневые системы.",
            "Варвак П.М. — Метод конечных элементов; Zienkiewicz — The Finite Element Method.",
            "TCVN 5574:2018 (bê tông), TCVN 5575:2012 (thép), ГОСТ 8239-89, ГОСТ 8240-97.",
        });
        ImGui::EndTabItem();
    }

    ImGui::EndTabBar();
    ImGui::End();
}
