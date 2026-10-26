// Các chủ đề hướng dẫn: mỗi chủ đề mở thành một cửa sổ riêng từ thẻ "Hướng dẫn" trên ribbon.
#pragma once

enum class GuideTopic {
    Start, Draw, Commands, Keys, Conventions,
    FrameGen, TrussGen, ArchGen, MeshGen, PileGen, MatSec, Supports, Bars,
    Loads, Combos,
    Static, Kinematics, ForceMethod, MatrixMethod, Buckling, Vibration, Dynamics, Influence, Fluid,
    Results, Design, SectionStress, Export,
    Count
};

const char* guideTip(GuideTopic t); // mô tả một dòng, dùng làm tooltip trên ribbon
