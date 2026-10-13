// Vẽ biểu tượng nét (stroke) từ dữ liệu SVG 24×24 kiểu Lucide / Feather bằng ImDrawList.
#pragma once

#include "imgui.h"

enum class SvgIconId { FilePlus, FolderOpen, Save, Undo, Redo, Sun, Moon, Count };

// Vẽ biểu tượng vào ô vuông cạnh `size` tâm `c`; nét dày 2/24 cạnh như bản gốc.
void drawSvgIcon(ImDrawList* dl, SvgIconId id, ImVec2 c, float size, ImU32 col);
