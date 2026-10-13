# FEM2D — tính toán kết cấu phẳng bằng phần tử hữu hạn

Chương trình C++17 tính kết cấu 2D theo phương pháp chuyển vị, giao diện **SDL3 + Dear ImGui** (docking).
Các phần tử, quy ước dấu và cách trình bày biểu đồ theo giáo trình Cơ học kết cấu Việt Nam / Liên Xô
(biểu đồ M vẽ về phía thớ căng, có gạch tung độ). Phần 3D sẽ bổ sung sau.

## Biên dịch

Yêu cầu: CMake ≥ 3.20, trình biên dịch C++17, SDL3 (nếu hệ thống không có, CMake sẽ tự tải mã nguồn).
Dear ImGui (v1.92.9b-docking) được CMake tải tự động.

```sh
# macOS
brew install cmake ninja sdl3
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/fem_tests          # các phép kiểm tra so với lời giải giải tích
./build/fem2d              # mở giao diện
./build/fem2d examples/khung_2_tang_2_nhip.fem
```

Tùy chọn CMake: `-DFEM_BUILD_GUI=OFF` (chỉ lõi tính + test), `-DFEM_IMGUI_DIR=/đường/dẫn/imgui` (dùng ImGui có sẵn, không cần mạng).

## Tính năng

**Phần tử**
- Thanh dàn (chỉ kéo nén; tải trong nhịp quy về nút như dầm đơn giản).
- Thanh khung Euler–Bernoulli: khớp đầu thanh (giải phóng mômen bằng khử tĩnh), dầm trên nền đàn hồi Winkler
  (ma trận nền tương thích), tải nhiệt độ (t đều và chênh lệch Δt).
- Tam giác biến dạng hằng (CST), tứ giác Q4 đẳng tham số có tùy chọn mode không tương thích QM6 (Wilson–Taylor).
- Bài toán ứng suất phẳng / biến dạng phẳng.

**Tải trọng**: tải nút (Fx, Fy, Mz); trên thanh: phân bố hình thang (theo trục cục bộ, tổng thể, hoặc trên hình chiếu),
tập trung, mômen tập trung, nhiệt độ, sai số chế tạo ΔL; trọng lượng bản thân; chuyển vị cưỡng bức (lún gối); gối đàn hồi.

**Nền**: Winkler, hai hệ số Pasternak / Власов, nền chỉ chịu nén (phân tích lặp, đoạn dầm vểnh lên tách khỏi nền).

**Phân tích bậc hai** (tính theo sơ đồ biến dạng, P-Δ và P-δ): lặp với ma trận độ cứng hình học đến khi hội tụ,
báo mất ổn định khi tải vượt tải trọng tới hạn.

**Cấu tạo hình học** (lệnh `CTHH`): W = 3D − 3H − 2K − C₀, hạng ma trận cân bằng, bậc siêu tĩnh, kết luận
bất biến hình / biến hình / tức thời biến hình và hiển thị dạng chuyển động của cơ cấu.
Nhiều trường hợp tải và tổ hợp tải tuyến tính.

**Kết quả**: chuyển vị, phản lực, biểu đồ N, Q, M, ứng suất pháp thớ biên N/A ± M/W, hình biến dạng (đường đàn hồi
chính xác trong nhịp), phổ màu ux, uy, |u|, σx, σy, τxy, σ1, σ2, σ von Mises; xuất báo cáo `.txt`.

**Biểu đồ bao** N, Q, M qua các tổ hợp (ví dụ hoạt tải cách nhịp của dầm liên tục).

**Đường ảnh hưởng** (đah) do lực P = 1 thẳng đứng di chuyển trên đường chạy tải (tự động hoặc chọn thanh):
phản lực gối, nội lực N/Q/M tại tiết diện bất kỳ (có bước nhảy đúng tại tiết diện), chuyển vị nút; kèm tung độ
max/min và diện tích Ω⁺, Ω⁻ để tính S = q·Ω, S = P·y.

**Ổn định** (phương pháp chuyển vị với ma trận độ cứng hình học): hệ số tải trọng tới hạn λcr, các dạng mất
ổn định, lực tới hạn Ncr và hệ số chiều dài tính toán μ (l0 = μ·l) của từng thanh chịu nén.

**Dao động riêng**: tần số ω, f, chu kỳ T và dạng dao động (có hoạt họa); khối lượng từ trọng lượng bản thân
và/hoặc quy từ tải trọng thẳng đứng của một bộ tải (m = |Fy|/g).

**Tổ hợp tải trọng tự động theo TCVN 2737:1995**: phân loại tải thường xuyên/tạm thời, nhóm loại trừ
(gió trái/phải, các sơ đồ chất hoạt tải), hệ số độ tin cậy γf; sinh tổ hợp cơ bản 1 (ψ = 1) và cơ bản 2 (ψ = 0,9).

**Tải trên cạnh phần tử phẳng**: lực phân bố qx, qy hoặc áp lực pháp tuyến biến đổi tuyến tính p = p0 + k·y
(ví dụ áp lực thủy tĩnh), quy về lực nút tương đương.

**Giao diện kiểu AutoCAD** (sáng/tối): ribbon có biểu tượng, dòng lệnh, thanh trạng thái GRID/SNAP/ORTHO/OSNAP/DYN,
con trỏ chữ thập, bắt điểm (nút, trung điểm, giao điểm, vuông góc, gần nhất — vẽ vào giữa thanh sẽ tự cắt thanh để nối),
chọn cửa sổ (trái→phải) / cắt ngang (phải→trái), grip để kéo nút, nhập liệu động cạnh con trỏ, bảng thuộc tính,
hoàn tác/làm lại, tự động tính lại khi sửa, bộ tạo nhanh khung nhiều tầng, dàn, lưới phần tử phẳng.
Thư viện vật liệu TCVN 5574:2018, TCVN 5575:2012 và thép hình ГОСТ 8239-89 (chữ I), ГОСТ 8240-97 (chữ [).

## Quy ước

- Đơn vị: kN, m, °C; E nhập bằng kN/m²; chuyển vị hiển thị mm, ứng suất MPa.
- Trục X sang phải, Y hướng lên; mômen và góc xoay dương ngược chiều kim đồng hồ.
- Trục cục bộ thanh: x từ nút i → j, y quay x ngược chiều kim đồng hồ 90°.
- N > 0 kéo; Q > 0 làm phần thanh quay thuận chiều kim đồng hồ; M > 0 căng thớ dưới (phía −y cục bộ).

## Thao tác (giống AutoCAD)

Gõ lệnh ở bất cứ đâu, không cần nhấp vào dòng lệnh. **Enter / Space / chuột phải** để kết thúc hoặc lặp lại lệnh,
**Esc** để hủy. Lệnh không phân biệt hoa thường; gõ `?` để xem danh sách.

| Lệnh | Chức năng |
|---|---|
| `L` / `DAN` / `PO` | vẽ thanh khung / thanh dàn / nút |
| `E`, `M`, `CO`, `RO`, `MI` | xóa, di chuyển, sao chép, xoay, đối xứng |
| `DIV`, `GOP` | chia thanh thành n đoạn, gộp nút trùng |
| `NGAM`, `KHOP`, `DIDONG`, `DIDONGX`, `NGAMTRUOT`, `BOGOI` | gối tựa |
| `TAINUT`, `TAIPB`, `TAITT`, `MOMEN` | tải nút, tải phân bố, lực tập trung, mômen |
| `TINH` (F5), `PDELTA`, `CTHH` | tính tĩnh, bật/tắt P-Δ, cấu tạo hình học |
| `ONDINH`, `DAODONG` | ổn định, dao động riêng |
| `DI`, `ID`, `LI`, `Z` | đo khoảng cách, tọa độ, liệt kê, zoom (`E` toàn bộ, `W` cửa sổ) |
| `U`, `REDO` | hoàn tác, làm lại |

Nhập tọa độ: `x,y` (tuyệt đối), `@dx,dy` (tương đối), `@d<góc` (cực), hoặc gõ một số để đi theo hướng con trỏ.

| Phím / chuột | Chức năng |
|---|---|
| F3 / F7 / F8 / F9 / F12 | OSNAP / lưới / ORTHO / SNAP / nhập liệu động |
| F2 | mở rộng lịch sử dòng lệnh |
| Ctrl+1 | bảng thuộc tính |
| Ctrl+Z / Ctrl+Y, Ctrl+S, Delete | hoàn tác / làm lại, lưu, xóa |
| Lăn chuột | zoom quanh con trỏ; kéo chuột giữa/phải: pan; nhấp đúp chuột giữa: vừa khung |
| Nhấp / Shift+nhấp | thêm / bỏ đối tượng khỏi tập chọn; nhấp vào grip để kéo nút |

## Định dạng file `.fem`

File văn bản, mỗi dòng một khai báo, chỉ số đánh từ 1, `#` là chú thích. Xem các ví dụ trong `examples/`.

```
MATERIAL "Bê tông B20" 27500000 0.2 1e-05 25      # tên E nu alpha gamma
SECTION "CN 30x60" 0.18 0.0054 0.6                 # tên A I h
NODE 0 0
FRAME 1 2 1 1 0 0 0                                # n1 n2 vật_liệu tiết_diện khớp_i khớp_j kf
SUPPORT 1 1 1 1 0 0 0                              # nút khóa_x khóa_y khóa_φ kx ky kφ
LOADCASE "Tĩnh tải" 1                              # tên hệ_số_trọng_lượng_bản_thân
  ELOAD 2 DIST GY -15 -15 0                        # phần_tử DIST hướng q1 q2 theo_hình_chiếu
  NLOAD 2 10 0 0                                   # nút Fx Fy Mz
COMBO "TT + 0,9Gió" 1 1.0 2 0.9
```

## Cấu trúc mã nguồn

```
src/core/   lõi tính toán, không phụ thuộc giao diện
  Model         dữ liệu mô hình
  Elements      ma trận độ cứng, tải tương đương, nội lực phần tử
  Skyline       ma trận skyline + phân tích LDLᵀ (COLSOL)
  Analysis      lắp ghép, đánh số lại RCM, điều kiện biên, giải, hậu xử lý, phát hiện hệ biến hình
  PostProcess   biểu đồ bao, đường ảnh hưởng
  Eigen         ổn định, dao động riêng (Householder + QL cho ma trận đối xứng)
  Kinematics    phân tích cấu tạo hình học (ma trận cân bằng, hạng, không gian rỗng)
  Mesher        tạo lưới, chia thanh
  Library       vật liệu, thép hình
  ModelIO       đọc/ghi .fem, báo cáo
  Examples      bài toán mẫu
src/app/    giao diện SDL3 + ImGui (Ribbon: ribbon, dòng lệnh, chủ đề; Commands: các lệnh; View: vùng vẽ; Panels: bảng)
tests/      kiểm chứng với lời giải giải tích
```

Chế độ dành cho phát triển: `fem2d --shot ảnh.bmp <số ví dụ> <chế độ xem> [đại lượng màu] [e | m<thanh>:<x>]` chụp màn hình rồi thoát;
`fem2d --smoke` giả lập gõ lệnh và thao tác chuột như người dùng để kiểm tra giao diện.

## Tài liệu tham khảo

- Lều Thọ Trình — Cơ học kết cấu, tập 1–2.
- Дарков А.В., Шапошников Н.Н. — Строительная механика.
- Смирнов А.Ф. и др. — Строительная механика. Стержневые системы.
- Варвак П.М. и др. — Метод конечных элементов; Zienkiewicz, Taylor — The Finite Element Method.
- Жемочкин Б.Н. — Расчёт балок на упругом основании.
- Смирнов А.Ф. — Устойчивость и колебания сооружений; СП 16.13330 (chiều dài tính toán cột khung).
- TCVN 2737:1995 — Tải trọng và tác động.
