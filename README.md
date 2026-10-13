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
tập trung, mômen tập trung, nhiệt độ; trọng lượng bản thân; chuyển vị cưỡng bức (lún gối); gối đàn hồi.
Nhiều trường hợp tải và tổ hợp tải tuyến tính.

**Kết quả**: chuyển vị, phản lực, biểu đồ N, Q, M, ứng suất pháp thớ biên N/A ± M/W, hình biến dạng (đường đàn hồi
chính xác trong nhịp), phổ màu ux, uy, |u|, σx, σy, τxy, σ1, σ2, σ von Mises; xuất báo cáo `.txt`.

**Biểu đồ bao** N, Q, M qua các tổ hợp (ví dụ hoạt tải cách nhịp của dầm liên tục).

**Đường ảnh hưởng** (đah) do lực P = 1 thẳng đứng di chuyển trên đường chạy tải (tự động hoặc chọn thanh):
phản lực gối, nội lực N/Q/M tại tiết diện bất kỳ (có bước nhảy đúng tại tiết diện), chuyển vị nút; kèm tung độ
max/min và diện tích Ω⁺, Ω⁻ để tính S = q·Ω, S = P·y.

**Tải trên cạnh phần tử phẳng**: lực phân bố qx, qy hoặc áp lực pháp tuyến biến đổi tuyến tính p = p0 + k·y
(ví dụ áp lực thủy tĩnh), quy về lực nút tương đương.

**Giao diện**: vẽ nút/thanh bằng chuột (bắt lưới, bắt nút), chọn bằng khung, bảng nhập liệu, hoàn tác/làm lại,
tự động tính lại khi sửa, bộ tạo nhanh khung nhiều tầng, dàn (Pratt, Howe, Warren, tam giác), lưới phần tử phẳng.
Thư viện vật liệu TCVN 5574:2018, TCVN 5575:2012 và thép hình ГОСТ 8239-89 (chữ I), ГОСТ 8240-97 (chữ [).

## Quy ước

- Đơn vị: kN, m, °C; E nhập bằng kN/m²; chuyển vị hiển thị mm, ứng suất MPa.
- Trục X sang phải, Y hướng lên; mômen và góc xoay dương ngược chiều kim đồng hồ.
- Trục cục bộ thanh: x từ nút i → j, y quay x ngược chiều kim đồng hồ 90°.
- N > 0 kéo; Q > 0 làm phần thanh quay thuận chiều kim đồng hồ; M > 0 căng thớ dưới (phía −y cục bộ).

## Thao tác

| Thao tác | Cách làm |
|---|---|
| Phóng to/thu nhỏ | lăn chuột (quanh con trỏ) |
| Di chuyển khung nhìn | kéo chuột phải hoặc chuột giữa |
| Vẽ thanh | chọn "+ Thanh khung"/"+ Thanh dàn", nhấp liên tiếp; chuột phải hoặc Esc để dừng |
| Chọn | nhấp, hoặc kéo khung; Shift để chọn thêm |
| Gán gối, tải | chọn nút/thanh → panel **Thuộc tính**; hoặc chuột phải → menu nhanh |
| Phím tắt | Ctrl+Z/Ctrl+Y hoàn tác/làm lại, Delete xóa, F vừa khung, F5 tính, Ctrl+S lưu |

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
  Mesher        tạo lưới, chia thanh
  Library       vật liệu, thép hình
  ModelIO       đọc/ghi .fem, báo cáo
  Examples      bài toán mẫu
src/app/    giao diện SDL3 + ImGui (View: vẽ sơ đồ; Panels: bảng nhập liệu, kết quả)
tests/      kiểm chứng với lời giải giải tích
```

Chế độ dành cho phát triển: `fem2d --shot ảnh.bmp <số ví dụ> <chế độ xem> [đại lượng màu] [e | m<thanh>:<x>]` chụp màn hình rồi thoát;
`fem2d --smoke` giả lập thao tác chuột để kiểm tra công cụ vẽ.

## Tài liệu tham khảo

- Lều Thọ Trình — Cơ học kết cấu, tập 1–2.
- Дарков А.В., Шапошников Н.Н. — Строительная механика.
- Смирнов А.Ф. и др. — Строительная механика. Стержневые системы.
- Варвак П.М. и др. — Метод конечных элементов; Zienkiewicz, Taylor — The Finite Element Method.
- Жемочкин Б.Н. — Расчёт балок на упругом основании.
