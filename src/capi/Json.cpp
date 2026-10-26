#include "capi/Json.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace fem::json {

// ---------------------------------------------------------------- ghi

void Writer::sep() {
    if (afterKey_) {
        afterKey_ = false;
        return;
    }
    if (!first_.empty()) {
        if (!first_.back()) out_ += ',';
        first_.back() = false;
    }
}

Writer& Writer::beginObject() {
    sep();
    out_ += '{';
    first_.push_back(true);
    return *this;
}

Writer& Writer::endObject() {
    out_ += '}';
    first_.pop_back();
    return *this;
}

Writer& Writer::beginArray() {
    sep();
    out_ += '[';
    first_.push_back(true);
    return *this;
}

Writer& Writer::endArray() {
    out_ += ']';
    first_.pop_back();
    return *this;
}

Writer& Writer::key(const char* k) {
    sep();
    quoted(k);
    out_ += ':';
    afterKey_ = true;
    return *this;
}

Writer& Writer::value(double v) {
    sep();
    // NaN/Infinity ghi dạng token để phía Kotlin đọc với allowSpecialFloatingPointValues
    if (std::isnan(v)) out_ += "NaN";
    else if (std::isinf(v)) out_ += v > 0 ? "Infinity" : "-Infinity";
    else {
        char buf[32];
        std::snprintf(buf, sizeof buf, "%.17g", v);
        out_ += buf;
    }
    return *this;
}

Writer& Writer::value(int v) { return value((long long)v); }

Writer& Writer::value(long long v) {
    sep();
    out_ += std::to_string(v);
    return *this;
}

Writer& Writer::value(bool v) {
    sep();
    out_ += v ? "true" : "false";
    return *this;
}

Writer& Writer::value(const std::string& v) { return value(v.c_str()); }

Writer& Writer::value(const char* v) {
    sep();
    quoted(v);
    return *this;
}

void Writer::quoted(const char* v) {
    out_ += '"';
    for (const unsigned char* p = (const unsigned char*)v; *p; ++p) {
        switch (*p) {
        case '"': out_ += "\\\""; break;
        case '\\': out_ += "\\\\"; break;
        case '\n': out_ += "\\n"; break;
        case '\r': out_ += "\\r"; break;
        case '\t': out_ += "\\t"; break;
        default:
            if (*p < 0x20) {
                char buf[8];
                std::snprintf(buf, sizeof buf, "\\u%04x", *p);
                out_ += buf;
            } else {
                out_ += (char)*p; // UTF-8 giữ nguyên
            }
        }
    }
    out_ += '"';
}

Writer& Writer::null() {
    sep();
    out_ += "null";
    return *this;
}

// ---------------------------------------------------------------- đọc

const Value* Value::find(const char* k) const {
    if (type != Type::Object) return nullptr;
    auto it = obj.find(k);
    return it == obj.end() ? nullptr : &it->second;
}

double Value::getNum(const char* k, double def) const {
    const Value* v = find(k);
    return v && v->type == Type::Number ? v->num : def;
}

int Value::getInt(const char* k, int def) const {
    const Value* v = find(k);
    return v && v->type == Type::Number ? (int)std::lround(v->num) : def;
}

bool Value::getBool(const char* k, bool def) const {
    const Value* v = find(k);
    return v && v->type == Type::Bool ? v->b : def;
}

std::string Value::getStr(const char* k, const std::string& def) const {
    const Value* v = find(k);
    return v && v->type == Type::String ? v->str : def;
}

const std::vector<Value>& Value::getArr(const char* k) const {
    static const std::vector<Value> empty;
    const Value* v = find(k);
    return v && v->type == Type::Array ? v->arr : empty;
}

namespace {

struct Parser {
    explicit Parser(const std::string& text) : s(text) {}
    const std::string& s;
    size_t i = 0;
    std::string err;

    void ws() {
        while (i < s.size() && (s[i] == ' ' || s[i] == '\n' || s[i] == '\r' || s[i] == '\t')) ++i;
    }
    bool fail(const char* what) {
        if (err.empty()) err = std::string(what) + " tại vị trí " + std::to_string(i);
        return false;
    }
    bool lit(const char* w) {
        size_t n = 0;
        while (w[n]) ++n;
        if (s.compare(i, n, w) != 0) return false;
        i += n;
        return true;
    }
    static void utf8(std::string& o, unsigned cp) {
        if (cp < 0x80) o += (char)cp;
        else if (cp < 0x800) o += (char)(0xC0 | (cp >> 6)), o += (char)(0x80 | (cp & 0x3F));
        else if (cp < 0x10000)
            o += (char)(0xE0 | (cp >> 12)), o += (char)(0x80 | ((cp >> 6) & 0x3F)), o += (char)(0x80 | (cp & 0x3F));
        else
            o += (char)(0xF0 | (cp >> 18)), o += (char)(0x80 | ((cp >> 12) & 0x3F)),
                o += (char)(0x80 | ((cp >> 6) & 0x3F)), o += (char)(0x80 | (cp & 0x3F));
    }
    bool hex4(unsigned& cp) {
        if (i + 4 > s.size()) return fail("thiếu mã \\u");
        cp = (unsigned)std::strtoul(s.substr(i, 4).c_str(), nullptr, 16);
        i += 4;
        return true;
    }
    bool string(std::string& o) {
        ++i; // "
        while (i < s.size() && s[i] != '"') {
            char c = s[i++];
            if (c != '\\') {
                o += c;
                continue;
            }
            if (i >= s.size()) return fail("chuỗi chưa đóng");
            char e = s[i++];
            switch (e) {
            case 'n': o += '\n'; break;
            case 'r': o += '\r'; break;
            case 't': o += '\t'; break;
            case 'b': o += '\b'; break;
            case 'f': o += '\f'; break;
            case 'u': {
                unsigned cp;
                if (!hex4(cp)) return false;
                if (cp >= 0xD800 && cp < 0xDC00 && s.compare(i, 2, "\\u") == 0) { // cặp thay thế
                    i += 2;
                    unsigned lo;
                    if (!hex4(lo)) return false;
                    cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                }
                utf8(o, cp);
                break;
            }
            default: o += e;
            }
        }
        if (i >= s.size()) return fail("chuỗi chưa đóng");
        ++i;
        return true;
    }
    bool value(Value& v, int depth) {
        if (depth > 64) return fail("lồng nhau quá sâu");
        ws();
        if (i >= s.size()) return fail("thiếu giá trị");
        char c = s[i];
        if (c == '{') {
            v.type = Value::Type::Object;
            ++i;
            ws();
            if (i < s.size() && s[i] == '}') return ++i, true;
            for (;;) {
                ws();
                if (i >= s.size() || s[i] != '"') return fail("cần tên khóa");
                std::string k;
                if (!string(k)) return false;
                ws();
                if (i >= s.size() || s[i] != ':') return fail("cần ':'");
                ++i;
                if (!value(v.obj[k], depth + 1)) return false;
                ws();
                if (i < s.size() && s[i] == ',') {
                    ++i;
                    continue;
                }
                if (i < s.size() && s[i] == '}') return ++i, true;
                return fail("cần ',' hoặc '}'");
            }
        }
        if (c == '[') {
            v.type = Value::Type::Array;
            ++i;
            ws();
            if (i < s.size() && s[i] == ']') return ++i, true;
            for (;;) {
                v.arr.emplace_back();
                if (!value(v.arr.back(), depth + 1)) return false;
                ws();
                if (i < s.size() && s[i] == ',') {
                    ++i;
                    continue;
                }
                if (i < s.size() && s[i] == ']') return ++i, true;
                return fail("cần ',' hoặc ']'");
            }
        }
        if (c == '"') {
            v.type = Value::Type::String;
            return string(v.str);
        }
        if (lit("true")) return v.type = Value::Type::Bool, v.b = true, true;
        if (lit("false")) return v.type = Value::Type::Bool, v.b = false, true;
        if (lit("null")) return v.type = Value::Type::Null, true;
        v.type = Value::Type::Number;
        if (lit("NaN")) return v.num = NAN, true;
        if (lit("Infinity")) return v.num = INFINITY, true;
        if (lit("-Infinity")) return v.num = -INFINITY, true;
        const char* start = s.c_str() + i;
        char* end = nullptr;
        v.num = std::strtod(start, &end);
        if (end == start) return fail("giá trị không hợp lệ");
        i += (size_t)(end - start);
        return true;
    }
};

} // namespace

bool parse(const std::string& text, Value& out, std::string& error) {
    Parser p(text);
    out = Value{};
    if (!p.value(out, 0)) {
        error = "JSON lỗi: " + p.err;
        return false;
    }
    p.ws();
    if (p.i != text.size()) {
        error = "JSON lỗi: dư ký tự tại vị trí " + std::to_string(p.i);
        return false;
    }
    return true;
}

} // namespace fem::json
