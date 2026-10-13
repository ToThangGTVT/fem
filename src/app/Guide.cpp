// Cửa sổ hướng dẫn theo từng chức năng, mở từ thẻ "Hướng dẫn" trên ribbon. Mỗi chủ đề là một cửa sổ nổi riêng:
// chỉ chỗ đặt chức năng trên giao diện, lệnh gõ tương ứng và các bước sử dụng.
#include <string>
#include <vector>

#include "App.h"
#include "Commands.h"
#include "Guide.h"
#include "imgui.h"

namespace {

struct Section {
    const char* title;
    std::vector<const char*> items;
};

struct Topic {
    const char* title;
    const char* tip;
    const char* where; // vị trí trên giao diện
    const char* cmd;   // lệnh mở chức năng (nullptr: không có)
    std::vector<Section> sections;
};

// Thứ tự trùng với enum GuideTopic
const std::vector<Topic>& topics() {
    static const std::vector<Topic> t = {
        // ------------------------------------------------------------------ cơ bản
        {"Bắt đầu",
         "Các bước dựng và tính một kết cấu",
         "Hướng dẫn chung · F1",
         nullptr,
         {{"Các bước",
           {"1. Mở một ví dụ (nút FEM › Ví dụ mẫu) hoặc dựng sơ đồ: vẽ thanh (L, DAN), hoặc dùng Tạo khung / Tạo dàn / "
            "Tạo vòm / Lưới phẳng ở thẻ Trang chủ.",
            "2. Khai báo vật liệu, tiết diện trong panel Mô hình; chọn thanh rồi bấm \"Gán cho thanh chọn\".",
            "3. Đặt gối: chọn nút rồi bấm Ngàm / Gối cố định / Gối di động.",
            "4. Đặt tải ở thẻ Tải trọng. Nhiều trường hợp tải thì lập thêm tổ hợp.",
            "5. Bấm Tính (F5). Xem biểu đồ N, Q, M, chuyển vị ở thẻ Kết quả.",
            "6. Xuất báo cáo HTML hoặc bản vẽ DXF."}},
          {"Các panel",
           {"Thuộc tính: sửa nút / phần tử đang chọn (tọa độ, gối, vật liệu, tiết diện, tải). Ctrl+1.",
            "Mô hình: bảng vật liệu, tiết diện, nút, phần tử, gối, trường hợp tải, tổ hợp.",
            "Kết quả: bộ kết quả, đại lượng hiển thị, bảng chuyển vị / phản lực / nội lực, đường ảnh hưởng.",
            "Nhật ký: thông báo của chương trình, kết quả các lệnh đo và liệt kê.",
            "Panel bị đóng: mở lại ở thẻ Xem › Bảng, hoặc Khôi phục bố cục."}},
          {"Xem thêm",
           {"Mỗi nút trên thẻ Hướng dẫn mở hướng dẫn riêng của một chức năng; có thể mở nhiều cửa sổ cùng lúc.",
            "Trong cửa sổ hướng dẫn, nút \"Chạy lệnh\" mở ngay chức năng đó."}}}},
        {"Vẽ & chọn",
         "Dòng lệnh, nhập tọa độ, chọn đối tượng, xem sơ đồ",
         "Vùng sơ đồ kết cấu và dòng lệnh phía dưới",
         nullptr,
         {{"Dòng lệnh",
           {"Gõ lệnh ở bất cứ đâu, không cần nhấp vào dòng lệnh; không phân biệt hoa thường.",
            "Enter / Space / chuột phải: kết thúc lệnh, hoặc lặp lại lệnh vừa dùng. Esc: hủy.",
            "Mũi tên lên / xuống: gọi lại các dòng đã gõ. F2: mở rộng lịch sử dòng lệnh."}},
          {"Nhập tọa độ (m)",
           {"x,y — tọa độ tuyệt đối.",
            "@dx,dy — tương đối so với điểm trước.",
            "@d<góc — theo khoảng cách và góc (độ, ngược chiều kim đồng hồ từ trục X).",
            "Chỉ gõ một số — đặt điểm cách điểm trước khoảng đó theo hướng con trỏ chuột."}},
          {"Chọn đối tượng",
           {"Nhấp để thêm vào tập chọn, Shift+nhấp để bỏ ra. Esc: bỏ chọn.",
            "Kéo từ trái sang phải: chọn kiểu cửa sổ (nằm trọn trong khung).",
            "Kéo từ phải sang trái: chọn kiểu cắt ngang (chạm khung là chọn).",
            "Nhấp grip (ô vuông xanh) để kéo nút."}},
          {"Xem sơ đồ",
           {"Lăn chuột: phóng to / thu nhỏ quanh con trỏ.",
            "Kéo chuột giữa hoặc chuột phải: dịch khung nhìn.",
            "Nhấp đúp chuột giữa, hoặc nút Vừa khung: xem toàn bộ."}},
          {"Trợ giúp vẽ (thanh trạng thái)",
           {"F3 OSNAP: bắt nút, trung điểm, giao điểm, chân vuông góc, điểm gần nhất trên thanh.",
            "F8 ORTHO: chỉ vẽ theo phương ngang, đứng.",
            "F9 SNAP: bắt theo bước lưới (đặt bước lưới ở panel Mô hình › Chung). F7: hiện / ẩn lưới.",
            "F12 DYN: ô nhập số cạnh con trỏ."}}}},
        {"Danh sách lệnh", "Tất cả lệnh gõ được ở dòng lệnh", "Dòng lệnh (gõ ? để in danh sách vào nhật ký)", nullptr, {}},
        {"Phím tắt", "Phím tắt bàn phím và chuột", "Toàn chương trình", nullptr, {}},
        {"Quy ước & đơn vị",
         "Đơn vị, hệ tọa độ, quy ước dấu nội lực, tài liệu tham khảo",
         "Áp dụng cho toàn bộ số liệu nhập và kết quả",
         nullptr,
         {{"Đơn vị",
           {"Lực kN, chiều dài m, mômen kNm, E kN/m².",
            "Ứng suất hiển thị MPa, chuyển vị mm, góc xoay rad.",
            "Lò xo gối: kN/m, kNm/rad. Kích thước trong Thiết kế tiết diện: cm."}},
          {"Hệ tọa độ & quy ước dấu",
           {"Trục X sang phải, Y hướng lên; mômen, góc xoay dương ngược chiều kim đồng hồ.",
            "Trục cục bộ thanh: x từ nút i đến j, y quay x ngược chiều kim đồng hồ 90°.",
            "N > 0: kéo. Q > 0: làm phần thanh quay thuận chiều kim đồng hồ.",
            "M > 0: căng thớ dưới (phía −y cục bộ). Biểu đồ M vẽ về phía thớ căng."}},
          {"Phần tử",
           {"Thanh dàn: chỉ chịu kéo nén. Tải trong nhịp quy về nút như dầm đơn giản.",
            "Thanh khung Euler–Bernoulli: khớp đầu thanh, nền đàn hồi, tải nhiệt độ.",
            "Phần tử phẳng: tam giác CST, T6; tứ giác Q4 (tùy chọn QM6), Q8.",
            "Bài toán ứng suất phẳng, biến dạng phẳng, đối xứng trục (r = x, z = y)."}},
          {"Tài liệu tham khảo",
           {"Lều Thọ Trình — Cơ học kết cấu (tập 1, 2).",
            "Дарков А.В., Шапошников Н.Н. — Строительная механика.",
            "Смирнов А.Ф. и др. — Строительная механика. Стержневые системы.",
            "Варвак П.М. — Метод конечных элементов; Zienkiewicz — The Finite Element Method.",
            "TCVN 5574:2018 (bê tông), TCVN 5575:2012 (thép), ГОСТ 8239-89, ГОСТ 8240-97."}}}},

        // ------------------------------------------------------------------ mô hình
        {"Tạo khung",
         "Khung nhiều tầng nhiều nhịp, dầm liên tục",
         "Trang chủ › Vẽ › Tạo khung",
         "TAOKHUNG",
         {{"Nhập",
           {"Nhịp (m): các nhịp cách nhau dấu cách; viết 3*6 cho ba nhịp 6 m. Ví dụ: 6 4.5 6.",
            "Chiều cao tầng (m): tương tự, từ dưới lên. Để trống để tạo dầm liên tục.",
            "Liên kết chân cột: ngàm hoặc khớp (gối cố định).",
            "Tiết diện cột, tiết diện dầm: chọn từ danh sách tiết diện của mô hình.",
            "Tải phân bố trên dầm (kN/m) và lực gió tại mỗi sàn đặt ở nút trái (kN)."}},
          {"Tạo",
           {"\"Thay thế mô hình hiện tại\": bỏ chọn để thêm khung vào mô hình đang có.",
            "Bấm Tạo; sau đó có thể sửa từng thanh, nút như mô hình vẽ tay."}}}},
        {"Tạo dàn",
         "Dàn Pratt, Howe, Warren, dàn tam giác",
         "Trang chủ › Vẽ › Tạo dàn",
         "TAODAN",
         {{"Nhập",
           {"Kiểu dàn: Pratt (xiên kéo), Howe (xiên nén), Warren (tam giác), dàn tam giác (mái dốc).",
            "Số khoang (chẵn), chiều dài khoang d (m), chiều cao dàn h (m), tiết diện thanh.",
            "Lực tại mỗi nút (kN), đặt ở cánh dưới hoặc cánh trên (\"Đặt lực ở cánh trên\")."}},
          {"Ghi chú",
           {"Các thanh tạo ra là thanh dàn (chỉ chịu kéo nén).",
            "Xem lực dọc bằng nút N ở thẻ Kết quả: thanh kéo màu xanh, thanh nén màu đỏ."}}}},
        {"Tạo vòm",
         "Vòm parabol hoặc cung tròn, hai khớp / ba khớp / ngàm",
         "Trang chủ › Vẽ › Tạo vòm",
         "VOM",
         {{"Nhập",
           {"Nhịp L (m), mũi tên f (m), số đoạn chia trục vòm.",
            "Trục vòm: parabol y = 4f·x(L − x)/L² hoặc cung tròn.",
            "Liên kết: hai khớp, ba khớp (khớp đỉnh), hai đầu ngàm.",
            "Tải đều theo hình chiếu (kN/m)."}},
          {"Ghi chú",
           {"Vòm parabol ba khớp chịu tải đều theo hình chiếu là vòm có trục hợp lý: M = 0 trên toàn vòm.",
            "Thanh vòm là các đoạn thẳng; tăng số đoạn để kết quả sát vòm cong."}}}},
        {"Lưới phẳng",
         "Lưới phần tử phẳng trên miền chữ nhật: tấm, dầm-tường, tường",
         "Trang chủ › Vẽ › Lưới phẳng",
         "LUOI",
         {{"Miền và lưới",
           {"Góc dưới trái x0, y0; rộng, cao (m); số ô theo X, theo Y.",
            "Loại phần tử: tứ giác Q4, tam giác CST, tứ giác bậc hai Q8, tam giác bậc hai T6.",
            "Q4 có tùy chọn mode không tương thích (QM6): chịu uốn chính xác hơn với lưới thưa.",
            "Chiều dày t (m) và vật liệu."}},
          {"Tùy chọn nhanh",
           {"Ngàm toàn bộ cạnh dưới, hoặc gối cố định + di động tại hai góc dưới.",
            "Tải phân bố mép trên (kN/m)."}},
          {"Ghi chú",
           {"Chọn ứng suất phẳng / biến dạng phẳng / đối xứng trục ở panel Mô hình › Chung.",
            "Nút trùng với nút có sẵn được nối, nên có thể ghép nhiều miền chữ nhật.",
            "Tải trên cạnh biên, áp lực nước: xem hướng dẫn Tải trọng.",
            "Kết quả: nút Phổ màu ở thẻ Kết quả, chọn đại lượng σx, σy, τxy, σ1, σ2, von Mises..."}}}},
        {"Cọc ngang",
         "Cọc chịu tải ngang theo phương pháp m (TCVN 10304:2014)",
         "Trang chủ › Vẽ › Cọc ngang",
         "COC",
         {{"Nhập",
           {"Chiều dài trong đất L, chiều dài tự do l0 (m), đường kính hoặc cạnh d (m), cọc vuông hay tròn.",
            "Hệ số m (kN/m⁴): hệ số nền tăng tuyến tính theo độ sâu C_z = m·z.",
            "Số đoạn chia phần trong đất; lực ngang H (kN) và mômen M (kNm) tại đầu cọc."}},
          {"Ghi chú",
           {"Chương trình tính bề rộng quy ước bp và tạo các thanh trên nền Winkler có hệ số tăng theo độ sâu.",
            "Sau khi tạo, bấm Tính (F5) để xem chuyển vị ngang và mômen dọc thân cọc."}}}},
        {"Vật liệu & tiết diện",
         "Khai báo vật liệu, tiết diện, thiết kế tiết diện",
         "Panel Mô hình › Vật liệu, Tiết diện; Trang chủ › Thuộc tính thanh",
         "TIETDIEN",
         {{"Vật liệu",
           {"Panel Mô hình › Vật liệu: + Thêm, sửa các thông số trong bảng.",
            "Loại vật liệu (thép / bê tông) quyết định cách kiểm tra trong Kiểm tra thiết kế."}},
          {"Tiết diện",
           {"Panel Mô hình › Tiết diện: + Chữ nhật, + Tròn đặc, + Ống, hoặc + Thêm rồi nhập trực tiếp đặc trưng.",
            "Thanh vẽ mới dùng vật liệu, tiết diện chọn ở Trang chủ › Thuộc tính thanh.",
            "Đổi cho thanh có sẵn: chọn thanh, chọn vật liệu / tiết diện rồi bấm \"Gán cho thanh chọn\"."}},
          {"Thiết kế tiết diện (TIETDIEN)",
           {"Chọn dạng tiết diện và kích thước theo cm.",
            "Đa giác: mỗi dòng một chu tuyến \"x,y x,y ...\"; ngược chiều kim đồng hồ là phần đặc, thuận chiều là lỗ.",
            "Kết quả: diện tích, trọng tâm, mômen quán tính, trục chính, lõi tiết diện.",
            "Đưa vào mô hình: \"Thêm tiết diện mới\" hoặc \"Ghi đè\" lên tiết diện có sẵn."}}}},
        {"Gối tựa",
         "Ngàm, khớp, gối di động, lò xo, lún gối",
         "Trang chủ › Gối tựa; panel Thuộc tính khi chọn nút",
         nullptr,
         {{"Đặt gối",
           {"Chọn nút rồi bấm nút gối trên ribbon, hoặc bấm nút gối trước rồi chọn nút và nhấn Enter.",
            "Ngàm (NGAM): khóa ux, uy, φ. Gối cố định (KHOP): khóa ux, uy. Gối di động (DIDONG): khóa uy.",
            "Di động ngang (DIDONGX): khóa ux. Ngàm trượt (NGAMTRUOT): khóa ux, φ. Bỏ gối (BOGOI)."}},
          {"Gối đàn hồi và lún",
           {"Chọn nút, trong panel Thuộc tính bấm \"Gán lò xo\" rồi nhập độ cứng (kN/m, kNm/rad).",
            "\"Chuyển vị cưỡng bức\": lún gối theo ux, uy, φ trong trường hợp tải đang soạn (m, rad).",
            "Bảng tổng hợp ở panel Mô hình › Gối tựa."}}}},
        {"Thuộc tính thanh",
         "Khớp đầu thanh, nền đàn hồi, thanh vút, chia thanh",
         "Panel Thuộc tính khi chọn thanh (Ctrl+1)",
         nullptr,
         {{"Loại và liên kết",
           {"Thanh khung / thanh dàn.",
            "Khớp đầu i, khớp đầu j: giải phóng mômen tại đầu thanh.",
            "Đảo chiều: đổi nút i và j (đổi trục cục bộ, chiều của tải theo trục cục bộ)."}},
          {"Nền đàn hồi",
           {"Nền Winkler k·b (kN/m²): dầm trên nền đàn hồi.",
            "Nền Pasternak kG·b (kN): hệ số thứ hai của nền hai hệ số Pasternak / Власов.",
            "Nền chỉ chịu nén: tính lặp, đoạn dầm vểnh lên được tách khỏi nền."}},
          {"Thanh vút và chia thanh",
           {"Thanh vút (VUT): chọn tiết diện đầu j, số đoạn rồi bấm \"Tạo thanh vút\" để thay bằng chuỗi đoạn "
            "có tiết diện thay đổi dần.",
            "Chia thanh (DIV): chia thanh chọn thành n đoạn bằng nhau."}}}},

        // ------------------------------------------------------------------ tải trọng
        {"Tải trọng",
         "Tải nút, tải trên thanh, nhiệt độ, lún, trọng lượng bản thân, tải cạnh",
         "Thẻ Tải trọng; panel Thuộc tính khi chọn nút / thanh",
         nullptr,
         {{"Đặt tải nhanh trên ribbon",
           {"Chọn trường hợp tải đang soạn ở Tải trọng › Trường hợp tải.",
            "Tải nút (TAINUT): lực Fx, Fy, mômen Mz tại nút chọn.",
            "Phân bố (TAIPB): tải đều / hình thang trên thanh chọn.",
            "Tập trung (TAITT), Mômen (MOMEN): nhấp vị trí trên thanh."}},
          {"Sửa chi tiết trong panel Thuộc tính",
           {"Chọn thanh, bấm \"Thêm tải cho thanh chọn\", chọn Loại: phân bố, tập trung, mômen, nhiệt độ, chế tạo ΔL.",
            "Hướng: x, y cục bộ hoặc X, Y tổng thể; \"Theo hình chiếu\" cho tải tính trên hình chiếu.",
            "Tải phân bố: q đầu i, q đầu j (kN/m). Tải tập trung: vị trí a (m) hoặc a/L.",
            "Nhiệt độ: t (°C) tăng đều và Δt = t dưới − t trên.",
            "Chế tạo ΔL (m): sai số chiều dài chế tạo của thanh."}},
          {"Tải khác",
           {"Trọng lượng bản thân: bật / tắt cho trường hợp đang soạn (hệ số ở panel Mô hình › Tải trọng).",
            "Lún gối: chọn nút, \"Chuyển vị cưỡng bức\" trong panel Thuộc tính.",
            "Phần tử phẳng: chọn các nút trên cạnh biên, \"Đặt tải lên cạnh biên giữa các nút chọn\": lực qx, qy "
            "hoặc áp lực pháp tuyến p = p0 + k·y (áp lực nước); \"Bỏ phần p < 0\" cho phần trên mặt nước."}}}},
        {"Trường hợp & tổ hợp",
         "Nhiều trường hợp tải, tổ hợp, sinh tổ hợp theo TCVN 2737",
         "Panel Mô hình › Tải trọng, Tổ hợp; thẻ Tải trọng › Tổ hợp",
         nullptr,
         {{"Trường hợp tải",
           {"Panel Mô hình › Tải trọng: + TH, Nhân bản, Xóa TH.",
            "Mỗi trường hợp có loại tải: thường xuyên (tĩnh tải) hoặc tạm thời (hoạt tải, gió...).",
            "Tải tạm thời cùng \"nhóm loại trừ\" (≠ 0) không tác dụng đồng thời, ví dụ gió trái và gió phải."}},
          {"Tổ hợp",
           {"Panel Mô hình › Tổ hợp: + Tổ hợp, + thành phần (trường hợp tải và hệ số).",
            "Sinh tổ hợp: tự lập tổ hợp cơ bản 1 và 2 theo TCVN 2737:1995 với hệ số ψ1, ψ2.",
            "Xem kết quả từng tổ hợp: chọn ở Kết quả › Bộ kết quả.",
            "Biểu đồ bao: Kết quả › Bộ kết quả › Biểu đồ bao (max / min qua các tổ hợp)."}}}},

        // ------------------------------------------------------------------ phân tích
        {"Tính toán",
         "Phân tích tĩnh tuyến tính và sơ đồ biến dạng P-Δ",
         "Phân tích › Tĩnh học",
         "TINH",
         {{"Tính",
           {"Bấm Tính hoặc F5. \"Tự động tính\": tính lại mỗi khi mô hình thay đổi.",
            "Panel Kết quả hiện số phương trình, thời gian, tổng lực ΣF và tổng phản lực ΣR để tự kiểm tra cân bằng.",
            "Cảnh báo (⚠) trong panel Kết quả: đọc kỹ; hệ biến hình thì kiểm tra bằng Cấu tạo hình học."}},
          {"Sơ đồ biến dạng (P-Δ)",
           {"Bật \"Sơ đồ biến dạng (P-Δ)\" (lệnh PDELTA) để tính bậc hai: lặp với ma trận độ cứng hình học.",
            "Chương trình báo mất ổn định khi tải vượt tải trọng tới hạn."}}}},
        {"Cấu tạo hình học",
         "Bậc tự do W, bậc siêu tĩnh, bất biến hình",
         "Phân tích › Cấu tạo",
         "CTHH",
         {{"Nội dung",
           {"Công thức Чебышёв: W = 3D − 3H − 2K − C₀ (D miếng cứng, H mối hàn, K khớp đơn, C₀ liên kết tựa).",
            "Phân tích hạng ma trận cân bằng: bậc siêu tĩnh, kết luận bất biến hình / biến hình / tức thời biến hình.",
            "Hệ biến hình: bấm \"Xem dạng biến hình\" để thấy cơ cấu chuyển động trên sơ đồ."}}}},
        {"Phương pháp lực",
         "Giải hệ siêu tĩnh theo phương pháp lực từng bước",
         "Phân tích › Giáo trình › Phương pháp lực",
         "PPL",
         {{"Các bước",
           {"Chọn bộ tải rồi bấm Giải.",
            "1. Bậc siêu tĩnh và hệ cơ bản (bỏ các liên kết thừa, thay bằng ẩn X).",
            "2. Hệ phương trình chính tắc δ·X + Δ_P = 0 với δik, ΔiP tính theo nhân biểu đồ.",
            "3. Nghiệm X.",
            "4. Chọn biểu đồ để vẽ trên sơ đồ: M̄i, M⁰_P, M = M⁰_P + Σ Xi·M̄i.",
            "5. Chuyển vị theo Maxwell–Mohr: nhập nút (hoặc \"Nút đang chọn\") rồi bấm Tính chuyển vị."}}}},
        {"PP chuyển vị ma trận",
         "Phương pháp chuyển vị dạng ma trận K·Z = P",
         "Phân tích › Giáo trình › Chuyển vị ma trận",
         "PPCV",
         {{"Các bước",
           {"Chọn bộ tải rồi bấm Lập.",
            "Ma trận độ cứng phần tử, chuyển sang hệ tổng thể.",
            "Ma trận độ cứng K của hệ, véc tơ tải P và nghiệm Z (chuyển vị nút tự do). Đơn vị kN, m, rad."}}}},
        {"Ổn định",
         "Tải trọng tới hạn, dạng mất ổn định, hệ số chiều dài tính toán μ",
         "Phân tích › Trị riêng & động lực › Ổn định",
         "ONDINH",
         {{"Thiết lập",
           {"Bộ tải gây nén: trường hợp hoặc tổ hợp tạo lực dọc trong thanh.",
            "Số dạng riêng cần tìm; số đoạn chia mỗi thanh (chia nhiều hơn thì dạng mất ổn định chính xác hơn).",
            "Bấm Tính."}},
          {"Kết quả",
           {"λcr: hệ số tải trọng tới hạn (tải tới hạn = λcr × tải đã đặt).",
            "Chọn dạng để xem trên sơ đồ; \"Hoạt họa dạng riêng\" để thấy chuyển động.",
            "Thanh chịu nén ở dạng 1: Ncr = λcr·N, chiều dài tính toán l0 = μ·l."}}}},
        {"Dao động riêng",
         "Tần số, chu kỳ, dạng dao động riêng",
         "Phân tích › Trị riêng & động lực › Dao động riêng",
         "DAODONG",
         {{"Thiết lập",
           {"Hệ số khối lượng bản thân: 1 = kể khối lượng của thanh / phần tử theo trọng lượng riêng.",
            "Quy tải thành khối lượng: chọn bộ tải (thường là tĩnh tải) và hệ số.",
            "Số dạng riêng, số đoạn chia mỗi thanh. Bấm Tính."}},
          {"Kết quả",
           {"Tần số và chu kỳ của từng dạng.",
            "Chọn dạng để xem trên sơ đồ, có hoạt họa."}}}},
        {"Động lực học",
         "Phổ phản ứng TCVN 9386, cưỡng bức điều hòa, tích phân Newmark",
         "Phân tích › Trị riêng & động lực › Động lực học",
         "DONGLUC",
         {{"Khối lượng",
           {"Mục \"Khối lượng và rời rạc hóa\" dùng chung cho cả ba thẻ: khối lượng bản thân, quy tải thành "
            "khối lượng, số dạng, số đoạn chia."}},
          {"Phổ phản ứng TCVN 9386",
           {"Nhập agR (g), hệ số tầm quan trọng γI, hệ số ứng xử q, loại nền; chọn phương X hoặc Y.",
            "Tổ hợp dạng SRSS (mặc định) hoặc CQC.",
            "Bảng dạng dao động: chu kỳ, khối lượng hữu hiệu; chọn một dòng để xem nội lực của dạng đó.",
            "Nội lực tổng hợp hiện trên sơ đồ như kết quả tĩnh."}},
          {"Cưỡng bức điều hòa",
           {"Lực P(t) = P₀·sin θt, P₀ phân bố theo bộ tải chọn; nhập tần số lực kích thích (Hz).",
            "Kết quả: biên độ dao động ổn định (không cản); so tần số kích thích với tần số riêng để tránh cộng hưởng."}},
          {"Theo thời gian (Newmark)",
           {"Dạng tác động: bậc thang, xung chữ nhật, xung nửa sin, điều hòa, gia tốc nền dạng sin, "
            "gia tốc nền tự nhập (mỗi dòng t (s), a (m/s²)).",
            "Bước Δt, thời gian tính, tỉ số cản ζ (cản Rayleigh).",
            "Xem lịch sử chuyển vị tại một nút; kéo thanh thời gian, nút Đỉnh, hoặc Chạy để hoạt họa."}}}},
        {"Đường ảnh hưởng",
         "Đường ảnh hưởng và xếp đoàn tải bất lợi",
         "Phân tích › Tải trọng di động; thiết lập trong panel Kết quả",
         nullptr,
         {{"Đường chạy tải",
           {"Lực P = 1 thẳng đứng, hướng xuống, di chuyển trên đường chạy tải.",
            "Bấm \"Tự động\", hoặc chọn thanh trên sơ đồ rồi \"Dùng thanh đang chọn\"."}},
          {"Đại lượng",
           {"Phản lực gối; nội lực tại tiết diện (thanh và x từ đầu i, có bước nhảy đúng tại tiết diện); "
            "chuyển vị nút.",
            "Kết quả: tung độ max / min, diện tích Ω⁺, Ω⁻. Tải đều q: S = q·Ω; lực P tại tung độ y: S = P·y."}},
          {"Xếp đoàn tải bất lợi",
           {"Chọn đoàn tải mẫu: xe H30, xe tải và xe hai trục HL-93 (kèm tải làn), НК-80; hoặc tự nhập P@khoảng cách.",
            "Chương trình tìm vị trí cho giá trị lớn nhất và nhỏ nhất; \"Vẽ đoàn tải ở vị trí bất lợi\" để xem."}}}},

        // ------------------------------------------------------------------ kết quả
        {"Xem kết quả",
         "Biểu đồ nội lực, biến dạng, phổ màu, bảng kết quả",
         "Thẻ Kết quả; panel Kết quả",
         nullptr,
         {{"Hiển thị",
           {"Bộ kết quả: chọn trường hợp tải hoặc tổ hợp.",
            "Sơ đồ, Biến dạng, N, Q, M, σ (ứng suất thớ biên N/A ± M/W), Phổ màu (phần tử phẳng), Dạng riêng.",
            "Biểu đồ M vẽ về phía thớ căng. Dàn: thanh kéo màu xanh, nén màu đỏ.",
            "Tỉ lệ: phóng đại biến dạng, tỉ lệ biểu đồ; bật / tắt Giá trị, Phản lực."}},
          {"Phổ màu",
           {"Panel Kết quả › Đại lượng: ux, uy, |u|, σx, σy, τxy, σ1, σ2, σz (σθ khi đối xứng trục), von Mises.",
            "\"Trung bình tại nút\" làm trơn phổ màu; \"Trên hình biến dạng\" vẽ phổ trên hình đã biến dạng."}},
          {"Bảng",
           {"Panel Kết quả: Chuyển vị, Phản lực, Nội lực thanh, Ứng suất (phần tử phẳng).",
            "Xuất báo cáo .txt ở panel Kết quả."}}}},
        {"Kiểm tra thiết kế",
         "Kiểm tra cấu kiện thép, BTCT và độ võng",
         "Kết quả › Thiết kế › Kiểm tra thiết kế",
         "KTTK",
         {{"Chuẩn bị",
           {"Cần có kết quả tính (F5). Chọn bộ tải xét.",
            "Loại vật liệu (thép / bê tông) đặt ở panel Mô hình › Vật liệu; BTCT cần tiết diện chữ nhật.",
            "\"Tô màu hệ số sử dụng trên sơ đồ\" để thấy cấu kiện đạt / không đạt."}},
          {"Thép TCVN 5575",
           {"Nhập γc, μ mặc định, hoặc lấy μ từ phân tích ổn định."}},
          {"BTCT TCVN 5574",
           {"Cốt thép dọc dầm, cột (cột đặt thép đối xứng, kể uốn dọc η), cốt đai.",
            "Nhập Rs, Rsc, Rsw (MPa), a (cm), diện tích một lớp đai Asw (cm²), μ cột."}},
          {"Độ võng",
           {"Độ võng so với dây cung nối hai đầu cấu kiện, giới hạn L / n; dùng bộ tải tiêu chuẩn (TTGH II).",
            "\"Chỉ dầm (cấu kiện nằm ngang)\" để bỏ qua cột."}}}},
        {"Ứng suất tiết diện",
         "Ứng suất σ, τ trên tiết diện, vòng Mohr, thuyết bền",
         "Kết quả › Tiết diện › Ứng suất tiết diện",
         "USTD",
         {{"Cách dùng",
           {"Chọn một thanh trên sơ đồ (cần có kết quả tính).",
            "Kéo \"Vị trí tiết diện x/L\" hoặc bấm |M|max, |Q|max để tới tiết diện nguy hiểm.",
            "Biểu đồ σ(y), τ(y) trên chiều cao tiết diện.",
            "Kéo \"Thớ khảo sát\" để xem trạng thái ứng suất tại một thớ: vòng tròn Mohr, ứng suất chính.",
            "Nhập R (MPa) để kiểm tra theo thuyết bền 3, 4."}},
          {"Ghi chú",
           {"Tiết diện không có hình dạng thì τ lấy trung bình Q/A; dùng Thiết kế tiết diện để có phân bố τ."}}}},
        {"Báo cáo & DXF",
         "Báo cáo HTML có hình, xuất và nhập bản vẽ DXF",
         "Kết quả › Xuất; nút FEM (Nhập bản vẽ DXF)",
         "BAOCAO",
         {{"Báo cáo",
           {"Báo cáo HTML (BAOCAO): số liệu, kết quả và hình biểu đồ; mở bằng trình duyệt, in ra PDF.",
            "Báo cáo văn bản .txt: nút FEM, hoặc panel Kết quả."}},
          {"Xuất DXF (DXFOUT)",
           {"Xuất sơ đồ và biểu đồ đang xem ra bản vẽ DXF để mở bằng AutoCAD."}},
          {"Nhập DXF (DXFIN)",
           {"Đoạn thẳng (LINE, LWPOLYLINE, POLYLINE) thành thanh; nút trùng được gộp, thanh tự cắt tại điểm nối chữ T.",
            "Chọn lớp (layer), đơn vị bản vẽ, nhập thành thanh khung hay thanh dàn, rồi bấm Nhập."}}}},
    };
    return t;
}

void bullets(const std::vector<const char*>& items) {
    for (const char* t : items) {
        ImGui::Bullet();
        ImGui::TextWrapped("%s", t);
    }
}

void keyTable() {
    static const std::pair<const char*, const char*> rows[] = {
        {"F1", "hướng dẫn Bắt đầu"},
        {"F2", "mở rộng lịch sử dòng lệnh"},
        {"F3 / F8 / F9 / F12", "OSNAP / ORTHO / SNAP / DYN"},
        {"F5", "tính"},
        {"F7", "lưới"},
        {"Ctrl+1", "panel Thuộc tính"},
        {"Ctrl+Z / Ctrl+Y", "hoàn tác / làm lại"},
        {"Ctrl+N / O / S", "tệp mới / mở / lưu"},
        {"Delete", "xóa đối tượng chọn"},
        {"Esc", "xóa dòng đang gõ, hủy lệnh, rồi bỏ chọn"},
        {"Enter / Space", "kết thúc hoặc lặp lại lệnh"},
        {"Lăn chuột", "phóng to / thu nhỏ quanh con trỏ"},
        {"Kéo chuột giữa / phải", "dịch khung nhìn"},
        {"Nhấp đúp chuột giữa", "vừa khung"},
        {"Shift+nhấp", "bỏ đối tượng khỏi tập chọn"},
    };
    if (!ImGui::BeginTable("keys", 2, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV)) return;
    ImGui::TableSetupColumn("k", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 10);
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

void commandTable() {
    static char filter[64] = "";
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputTextWithHint("##filter", "Lọc theo tên lệnh hoặc mô tả...", filter, sizeof filter);
    const std::string f = upper(filter);
    if (!ImGui::BeginTable("cmds", 2, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_ScrollY)) return;
    ImGui::TableSetupScrollFreeze(0, 1);
    ImGui::TableSetupColumn("Lệnh", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFontSize() * 11);
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

} // namespace

const char* guideTip(GuideTopic t) {
    IM_ASSERT(topics().size() == (size_t)GuideTopic::Count);
    return topics()[(int)t].tip;
}

void App::openGuide(GuideTopic t) {
    openGuides.insert((int)t);
    focusGuide = (int)t;
}

void App::guideWindows() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    for (auto it = openGuides.begin(); it != openGuides.end();) {
        const int id = *it;
        const Topic& tp = topics()[id];
        // xếp so le để các cửa sổ mở liên tiếp không che khít nhau
        const float off = 28.0f * (id % 8);
        ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x * 0.5f - 260 + off, vp->WorkPos.y + 90 + off),
                                ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(520, 560), ImGuiCond_FirstUseEver);
        if (focusGuide == id) {
            ImGui::SetNextWindowFocus();
            ImGui::SetNextWindowCollapsed(false);
            focusGuide = -1;
        }
        bool open = true;
        const std::string title = std::string("Hướng dẫn: ") + tp.title + "###guide" + std::to_string(id);
        if (ImGui::Begin(title.c_str(), &open)) {
            ImGui::TextDisabled("Vị trí:");
            ImGui::SameLine();
            ImGui::TextWrapped("%s", tp.where);
            if (tp.cmd) {
                if (ImGui::Button((std::string("Chạy lệnh ") + tp.cmd).c_str())) startCommand(tp.cmd);
            } else if (id == (int)GuideTopic::Influence) {
                if (ImGui::Button("Mở đường ảnh hưởng")) {
                    view = ResultView::Influence;
                    showResultPanel = true;
                }
            }
            ImGui::Spacing();
            if (id == (int)GuideTopic::Commands) commandTable();
            if (id == (int)GuideTopic::Keys) keyTable();
            for (const Section& s : tp.sections) {
                ImGui::SeparatorText(s.title);
                bullets(s.items);
            }
        }
        ImGui::End();
        it = open ? std::next(it) : openGuides.erase(it);
    }
}
