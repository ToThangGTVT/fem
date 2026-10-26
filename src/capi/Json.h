// JSON tối giản cho lớp giao tiếp C (fem_capi): bộ ghi dạng luồng và bộ đọc DOM.
// Chỉ dùng nội bộ để trao đổi dữ liệu với ứng dụng di động (Kotlin Multiplatform).
#pragma once
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace fem::json {

// ---------------------------------------------------------------- ghi
class Writer {
public:
    Writer& beginObject();
    Writer& endObject();
    Writer& beginArray();
    Writer& endArray();
    Writer& key(const char* k);
    Writer& value(double v);
    Writer& value(int v);
    Writer& value(long long v);
    Writer& value(bool v);
    Writer& value(const std::string& v);
    Writer& value(const char* v);
    Writer& null();
    template <class T> Writer& field(const char* k, const T& v) { return key(k).value(v); }
    template <class It> Writer& array(It first, It last) {
        beginArray();
        for (; first != last; ++first) value(*first);
        return endArray();
    }
    template <class C> Writer& field(const char* k, const std::vector<C>& v) { return key(k).array(v.begin(), v.end()); }
    const std::string& str() const { return out_; }

private:
    void sep();
    void quoted(const char* v);
    std::string out_;
    std::vector<bool> first_; // phần tử đầu tiên của từng mức lồng nhau
    bool afterKey_ = false;
};

// ---------------------------------------------------------------- đọc
struct Value {
    enum class Type { Null, Bool, Number, String, Array, Object } type = Type::Null;
    bool b = false;
    double num = 0;
    std::string str;
    std::vector<Value> arr;
    std::map<std::string, Value> obj;

    bool isNull() const { return type == Type::Null; }
    const Value* find(const char* k) const;
    // Truy cập có giá trị mặc định khi thiếu khóa hoặc sai kiểu
    double getNum(const char* k, double def = 0) const;
    int getInt(const char* k, int def = 0) const;
    bool getBool(const char* k, bool def = false) const;
    std::string getStr(const char* k, const std::string& def = {}) const;
    const std::vector<Value>& getArr(const char* k) const;
};

bool parse(const std::string& text, Value& out, std::string& error);

} // namespace fem::json
