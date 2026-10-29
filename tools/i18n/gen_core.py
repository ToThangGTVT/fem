#!/usr/bin/env python3
"""Sinh src/core/I18nStrings.cpp từ core_strings.py và kiểm tra bảng dịch của lõi.

Kiểm tra: mọi khóa tr()/trf() trong src có bản dịch; mỗi bản dịch giữ đúng các %d/%s... theo thứ tự (trf dùng printf)
và đúng các thẻ HTML của chữ gốc.
Dùng: python3 tools/i18n/gen_core.py            (sinh file + kiểm tra)
      python3 tools/i18n/gen_core.py --check    (chỉ kiểm tra, mã lỗi 1 nếu có lỗi)
"""
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from core_strings import ROWS  # noqa: E402
from extract import used_keys  # noqa: E402

LANGS = ["en", "ja", "ko", "ru"]
PRINTF = re.compile(r"%(?:%|[-+ #0]*\d*(?:\.\d+)?(?:hh|h|ll|l|z|j|t|L)?[diouxXeEfgGcsp])")
TAG = re.compile(r"<[^>]*>")


def c_literal(s):
    out = []
    for ch in s:
        if ch == "\\":
            out.append("\\\\")
        elif ch == '"':
            out.append('\\"')
        elif ch == "\n":
            out.append("\\n")
        elif ch == "\t":
            out.append("\\t")
        else:
            out.append(ch)
    return '"' + "".join(out) + '"'


def check():
    errors = []
    vi = [r[0] for r in ROWS]
    for k in sorted({k for k in vi if vi.count(k) > 1}):
        errors.append(f"trùng khóa: {k!r}")
    used = used_keys()
    for k in used:
        if k not in vi:
            errors.append(f"thiếu bản dịch: {k!r}")
    for k in vi:
        if k not in used:
            print(f"khóa không còn dùng: {k!r}")
    for row in ROWS:
        if len(row) != 5:
            errors.append(f"dòng thiếu cột: {row[0]!r}")
            continue
        want_fmt = [m.group(0) for m in PRINTF.finditer(row[0])]
        want_tag = TAG.findall(row[0])
        for lang, t in zip(LANGS, row[1:]):
            got_fmt = [m.group(0) for m in PRINTF.finditer(t)]
            if got_fmt != want_fmt:
                errors.append(f"{lang}: định dạng {got_fmt} ≠ {want_fmt} trong {row[0]!r}")
            if TAG.findall(t) != want_tag:
                errors.append(f"{lang}: thẻ HTML khác chữ gốc trong {row[0]!r}")
    return errors


def generate():
    lines = [
        "// Sinh tự động bởi tools/i18n/gen_core.py từ tools/i18n/core_strings.py — sửa ở đó, không sửa file này.",
        '#include "core/I18n.h"',
        "",
        "namespace fem {",
        "",
        "struct I18nEntry {",
        "    const char* vi;",
        "    const char* t[4];  // en, ja, ko, ru",
        "};",
        "",
        "extern const I18nEntry kI18nEntries[] = {",
    ]
    for r in ROWS:
        lines.append("    {" + c_literal(r[0]) + ", {" + ", ".join(c_literal(t) for t in r[1:]) + "}},")
    lines += [
        "};",
        "extern const int kI18nEntryCount = (int)(sizeof(kI18nEntries) / sizeof(kI18nEntries[0]));",
        "",
        "}  // namespace fem",
        "",
    ]
    (ROOT / "src/core/I18nStrings.cpp").write_text("\n".join(lines))


def main():
    errors = check()
    for e in errors:
        print(e)
    if "--check" not in sys.argv:
        generate()
        print(f"Đã sinh src/core/I18nStrings.cpp: {len(ROWS)} chuỗi × {len(LANGS)} ngôn ngữ.")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
