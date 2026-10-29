"""Khóa dịch dùng trong lõi: tr("…")/trf("…"), tên và mô tả bài mẫu (Examples.cpp), tên vật liệu có sẵn (Library.cpp)."""
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
LIT = r'"((?:[^"\\]|\\.)*)"'


def unescape(s):
    return re.sub(r'\\(.)', lambda m: {'n': '\n', 't': '\t'}.get(m.group(1), m.group(1)), s)


def used_keys():
    keys = []
    for p in sorted((ROOT / "src").rglob("*.cpp")):
        if p.name == "I18nStrings.cpp":
            continue
        t = p.read_text()
        keys += [unescape(m.group(1)) for m in re.finditer(r"\btrf?\(" + LIT, t)]
        if p.name == "Examples.cpp":
            keys += [unescape(x) for m in re.finditer(r'^\s*\{' + LIT + r', ' + LIT + r', \w+\},', t, re.M) for x in m.groups()]
        if p.name == "Library.cpp":
            body = t[t.index("materialPresets()"):t.index("profilePresets()")]
            keys += [unescape(m.group(1)) for m in re.finditer(r'^\s*\{' + LIT + r', ' + LIT + r', [\d.e]+', body, re.M)]
    seen, out = set(), []
    for k in keys:
        if k not in seen:
            seen.add(k)
            out.append(k)
    return out


if __name__ == "__main__":
    for k in used_keys():
        print(repr(k))
