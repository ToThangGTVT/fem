// Mô hình kết cấu phẳng (2D) cho phương pháp phần tử hữu hạn.
// Đơn vị thống nhất: kN, m, °C. Hệ tọa độ tổng thể: X sang phải, Y hướng lên, góc xoay ngược chiều kim đồng hồ dương.
// Mỗi nút có 3 bậc tự do: ux, uy, φz (các bậc tự do không có độ cứng sẽ được tự động khóa khi tính).
#pragma once
#include <array>
#include <string>
#include <utility>
#include <vector>

namespace fem {

enum class ElemType { Truss, Frame, Tri3, Quad4, Tri6, Quad8 };
enum class PlaneMode { Stress, Strain, Axisymmetric }; // ứng suất phẳng / biến dạng phẳng / đối xứng trục (x = r, y = z)

enum Dof { UX = 0, UY = 1, RZ = 2 };

struct Node {
    double x = 0, y = 0;
};

struct Material {
    std::string name = "Vật liệu";
    double E = 2.1e8;      // mô đun đàn hồi, kN/m²
    double nu = 0.3;       // hệ số Poisson
    double alpha = 1.2e-5; // hệ số giãn nở nhiệt, 1/°C
    double gamma = 78.5;   // trọng lượng riêng, kN/m³
    int kind = 0;          // 0: khác, 1: thép, 2: bê tông (dùng khi kiểm tra thiết kế)
    double fd = 0;         // cường độ tính toán, MPa: thép f (TCVN 5575); bê tông Rb (TCVN 5574)
    double ft = 0;         // bê tông: cường độ chịu kéo tính toán Rbt, MPa
};

// Hình dạng tiết diện (để tính đặc trưng hình học, ứng suất tiếp, vẽ). Kích thước theo m.
enum class ShapeKind { None = 0, Rect, Tee, IBeam, Box, Pipe, Circle, Channel, Angle, Polygon };
struct SectionShape {
    ShapeKind kind = ShapeKind::None;
    // Rect: b, h | Tee: b_cánh, t_cánh, t_sườn, h | IBeam: b_trên, t_trên, t_bụng, h, b_dưới, t_dưới
    // Box: b, h, t | Pipe: D, t | Circle: D | Channel: b, h, t_bụng, t_cánh | Angle: b, h, t
    std::array<double, 6> p{0, 0, 0, 0, 0, 0};
    // Polygon: các chu tuyến (x, y), ngược chiều kim đồng hồ là phần đặc, thuận chiều là lỗ
    std::vector<std::vector<std::array<double, 2>>> contours;
};

struct Section {
    std::string name = "Tiết diện";
    double A = 0.01;  // diện tích, m²
    double I = 1e-5;  // mômen quán tính với trục trung hòa nằm ngang, m⁴
    double h = 0.2;   // chiều cao (tải nhiệt), m
    double yTop = 0;  // khoảng cách trọng tâm → thớ trên (0: lấy h/2)
    double yBot = 0;  // khoảng cách trọng tâm → thớ dưới (0: lấy h/2)
    SectionShape shape;
    double top() const { return yTop > 0 ? yTop : h / 2; }
    double bot() const { return yBot > 0 ? yBot : h / 2; }
};

struct Element {
    ElemType type = ElemType::Frame;
    // nút: các góc theo chiều ngược kim đồng hồ, sau đó các nút giữa cạnh (T6: 1-2, 2-3, 3-1; Q8: 1-2, 2-3, 3-4, 4-1)
    std::array<int, 8> n{-1, -1, -1, -1, -1, -1, -1, -1};
    int mat = 0;
    int sec = 0;              // thanh: chỉ số tiết diện
    double t = 0.1;           // phần tử phẳng: chiều dày, m
    bool hingeI = false;      // khớp tại đầu i (giải phóng mômen)
    bool hingeJ = false;      // khớp tại đầu j
    double kf = 0;            // hệ số nền Winkler × bề rộng, kN/m² (dầm trên nền đàn hồi)
    double kG = 0;            // hệ số trượt của nền Pasternak × bề rộng, kN (mô hình nền hai hệ số)
    bool kfCompOnly = false;  // nền chỉ chịu nén (tách khỏi nền khi thanh bị nâng lên, phía −y cục bộ là nền)
    bool incompatible = true; // Quad4: dùng mode không tương thích (Wilson/Taylor QM6)

    int nodeCount() const {
        switch (type) {
        case ElemType::Truss:
        case ElemType::Frame: return 2;
        case ElemType::Tri3: return 3;
        case ElemType::Quad4: return 4;
        case ElemType::Tri6: return 6;
        case ElemType::Quad8: return 8;
        }
        return 0;
    }
    int cornerCount() const { return type == ElemType::Tri3 || type == ElemType::Tri6 ? 3 : type == ElemType::Quad4 || type == ElemType::Quad8 ? 4 : 2; }
    // Thứ tự nút đi quanh biên (góc – giữa cạnh – góc ...); trả về số nút
    int outline(int out[8]) const {
        static const int t3[] = {0, 1, 2}, q4[] = {0, 1, 2, 3}, t6[] = {0, 3, 1, 4, 2, 5}, q8[] = {0, 4, 1, 5, 2, 6, 3, 7};
        const int* src = t3;
        int n = 0;
        switch (type) {
        case ElemType::Tri3: src = t3, n = 3; break;
        case ElemType::Quad4: src = q4, n = 4; break;
        case ElemType::Tri6: src = t6, n = 6; break;
        case ElemType::Quad8: src = q8, n = 8; break;
        default: out[0] = 0, out[1] = 1; return 2;
        }
        for (int i = 0; i < n; ++i) out[i] = src[i];
        return n;
    }
    bool isBar() const { return type == ElemType::Truss || type == ElemType::Frame; }
    bool isPlane() const { return !isBar(); }
};

struct Support {
    int node = -1;
    std::array<bool, 3> fix{false, false, false}; // khóa ux, uy, φz
    std::array<double, 3> k{0, 0, 0};             // lò xo (gối đàn hồi): kN/m, kN/m, kNm/rad
    bool empty() const { return !fix[0] && !fix[1] && !fix[2] && k[0] == 0 && k[1] == 0 && k[2] == 0; }
};

struct NodalLoad {
    int node = -1;
    std::array<double, 3> F{0, 0, 0}; // Fx, Fy (kN), Mz (kNm)
};

enum class ElemLoadKind { Distributed, Point, Moment, Temperature, LackOfFit };
enum class LoadDir { LocalX, LocalY, GlobalX, GlobalY };

struct ElemLoad {
    int elem = -1;
    ElemLoadKind kind = ElemLoadKind::Distributed;
    LoadDir dir = LoadDir::GlobalY;
    // Distributed: q1, q2 = cường độ tại đầu i, j (kN/m), hình thang tuyến tính trên toàn thanh
    // Point:       q1 = P (kN), a = khoảng cách từ đầu i (m)
    // Moment:      q1 = M (kNm, ngược chiều kim đồng hồ dương), a = khoảng cách từ đầu i
    // Temperature: q1 = độ tăng nhiệt trung bình t (°C), q2 = t_dưới − t_trên (°C)
    // LackOfFit:   q1 = sai số chế tạo ΔL (m), dương: thanh dài hơn thiết kế
    double q1 = 0, q2 = 0;
    double a = 0;
    bool projected = false; // tải phân bố theo phương tổng thể tính trên chiều dài hình chiếu
};

struct PrescribedDisp {
    int node = -1;
    int dof = UY;
    double value = 0; // m hoặc rad
};

struct LoadCase {
    LoadCase() = default;
    LoadCase(std::string n, double sw = 0) : name(std::move(n)), selfWeight(sw) {}
    std::string name = "TH1";
    double selfWeight = 0; // hệ số trọng lượng bản thân (0 = không xét)
    int kind = 0;          // 0: tải thường xuyên (tĩnh tải), 1: tải tạm thời (hoạt tải, gió...)
    int group = 0;         // nhóm loại trừ: các tải tạm thời cùng nhóm ≠ 0 không xuất hiện đồng thời (gió trái/phải)
    double gammaF = 1.0;   // hệ số độ tin cậy dùng khi tự sinh tổ hợp
    std::vector<NodalLoad> nodal;
    std::vector<ElemLoad> elem;
    std::vector<PrescribedDisp> disp; // chuyển vị cưỡng bức của gối
};

struct Combo {
    Combo() = default;
    Combo(std::string n, std::vector<std::pair<int, double>> t) : name(std::move(n)), terms(std::move(t)) {}
    std::string name = "THop1";
    std::vector<std::pair<int, double>> terms; // (chỉ số trường hợp tải, hệ số)
    bool generated = false;                    // được sinh tự động (sẽ bị thay khi sinh lại)
};

// Dòng chảy chất lỏng: miền chữ nhật dựng trên cạnh cửa vào AB, chảy vuông góc với AB, dài L. Vật cản là các
// đường khép kín (thanh) vẽ bên trong miền. Đơn vị SI: m/s, kg/m³, Pa·s.
struct FlowDomain {
    bool active = false;
    // 0: mặt bằng (nhìn từ trên xuống, hai biên bên là bờ); 1: mặt cắt đứng (biên dưới là mặt đất / đáy, biên trên là
    // mặt thoáng hoặc trần, có trọng lực)
    int section = 0;
    double ax = 0, ay = 0, bx = 0, by = 1; // cạnh cửa vào A → B
    double length = 5;      // chiều dài miền theo hướng chảy (m)
    int side = -1;          // +1: chảy về phía bên trái của A→B, −1: bên phải
    double U = 1;           // vận tốc trung bình ở cửa vào (m/s)
    bool parabolic = false; // phân bố parabol (vận tốc lớn nhất 1,5·U) thay vì đều
    int medium = 0;         // 0 nước, 1 không khí, 2 dầu, 3 glyxerin, 4 tùy chọn
    double rho = 998.2, mu = 1.002e-3;
    bool slipWalls = false; // mặt bằng: hai bờ trượt (dòng tự do) thay vì thành cứng
    bool topWall = false;   // mặt cắt đứng: biên trên là trần cứng (cống, ống) thay vì mặt thoáng / không khí tự do
    bool hydrostatic = true; // mặt cắt đứng: cộng áp suất thủy tĩnh ρ·g·(y_mặt thoáng − y) vào áp suất và lực
    double h = 0, hFar = 0; // kích thước ô sát vật cản / xa nhất (m); 0: tự chọn
    bool transient = false;
    double dt = 0.01, endTime = 5;
    int turbulence = 0;          // 0: tự động (k-ε khi Re > 2300), 1: chảy tầng, 2: chảy rối k-ε
    double turbIntensity = 0.05; // cường độ rối ở cửa vào
    // mặt cắt đứng có mặt thoáng tự do (hai pha nước – không khí, VOF): AB gồm cả lớp không khí phía trên mặt nước
    bool freeSurface = false;
    double waterDepth = 0;  // độ sâu nước ở cửa vào, đo từ đáy dọc AB (m)
    double initDepth = -1;  // độ sâu nước ban đầu trong miền (m); < 0: bằng waterDepth
    double gate = 0;        // nước ban đầu chỉ trong đoạn dài này từ cửa vào (vỡ đập, xả cửa cống); 0: toàn miền
    // mặt cắt đứng: mặt đất / đáy thấm nước — nước thoát ra qua đáy với tốc độ thấm drainRate (m/s, ≈ hệ số thấm K khi
    // dòng thấm thẳng đứng, gradient thủy lực 1); có mặt thoáng tự do: chỉ phần ngập nước thoát
    bool drainGround = false;
    double drainRate = 0.01;
    int drainMode = 0; // 0: thấm với tốc độ drainRate; 1: thoát tự do, không giới hạn (cần mặt thoáng tự do)

    double width() const;
    void direction(double& nx, double& ny) const; // hướng chảy (véc tơ đơn vị)
    void corners(double c[4][2]) const;           // A, B, B + L·n, A + L·n
    // mặt cắt đứng: đầu A hay B của cửa vào nằm ở đáy (cạnh bên thấp hơn); cao độ mặt nước ở độ sâu depth đo dọc AB
    bool groundAtB() const;
    double waterLevelY(double depth) const;
};

struct Model {
    std::string title;
    FlowDomain flow;
    PlaneMode planeMode = PlaneMode::Stress;
    bool secondOrder = false; // tính theo sơ đồ biến dạng (P-Δ, phân tích bậc hai)
    std::vector<Node> nodes;
    std::vector<Material> materials;
    std::vector<Section> sections;
    std::vector<Element> elements;
    std::vector<Support> supports;
    std::vector<LoadCase> cases;
    std::vector<Combo> combos;

    void clear();
    int findNode(double x, double y, double tol) const;
    int addNode(double x, double y, double tol = 1e-9); // gộp nút trùng
    int findElement(ElemType type, const std::array<int, 8>& n) const;
    Support* supportAt(int node);
    const Support* supportAt(int node) const;
    Support& supportFor(int node);
    void cleanSupports();
    // Xóa nút (kèm các phần tử, gối, tải liên quan) và đánh số lại.
    void removeNodes(std::vector<int> ids);
    // Xóa phần tử (kèm tải trên phần tử) và đánh số lại.
    void removeElements(std::vector<int> ids);
    void removeUnusedNodes();
    double barLength(int e) const;
    bool hasPlane() const;
    bool hasBars() const;
    void bounds(double& x0, double& y0, double& x1, double& y1) const;
};

// Sinh tổ hợp tải trọng theo TCVN 2737:1995 (thay các tổ hợp tự sinh trước đó):
//   Tổ hợp cơ bản 1: Σ tải thường xuyên + 1 tải tạm thời (ψ1 = 1,0);
//   Tổ hợp cơ bản 2: Σ tải thường xuyên + ψ2·Σ (≥ 2 tải tạm thời), ψ2 = 0,9; mỗi nhóm loại trừ tối đa một tải.
// Hệ số của mỗi trường hợp tải nhân thêm γf. Trả về số tổ hợp đã sinh.
int generateCombos(Model& m, double psi1 = 1.0, double psi2 = 0.9, int maxCombos = 300, std::string* note = nullptr);

const char* elemTypeName(ElemType t);
const char* loadKindName(ElemLoadKind k);
const char* loadDirName(LoadDir d);
const char* dofName(int dof);

} // namespace fem
