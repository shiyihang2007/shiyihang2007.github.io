---
title: C++ 代码模板
---

ToDo:

1. SublimeText 适配
2. 编译命令行模板

## 目录

1. [VSCode 一键导入](#vscode-一键导入)
2. [常用模板](#常用模板)

## VSCode 一键导入

```json
{
	"base": {
		"prefix": ["main", "base"],
		"body": [
			"#include <algorithm>",
			"#include <cstdio>",
			"#include <iostream>",
			"using namespace std;",
			"int main()",
			"{",
			"    return 0;",
			"}"
		],
		"description": "基本代码模板",
		"isFileTemplate": true
	},
	"quickIO": {
		"prefix": ["quickIO", "qio", "fio", "qread", "qprint", "fread", "fprint"],
		"body": [
			"const size_t SIZE_IOBUF = 1e5 + 5;",
			"char bufIn[SIZE_IOBUF], *nowIn = bufIn, *endIn = bufIn;",
			"char getc()",
			"{",
			"    if (nowIn == endIn) {",
			"        endIn = (nowIn = bufIn) + fread(bufIn, 1, SIZE_IOBUF, stdin);",
			"        if (nowIn == endIn) {",
			"            return EOF;",
			"        }",
			"    }",
			"    return *nowIn++;",
			"}",
			"// char getc() { return (nowIn == endIn && (endIn = (nowIn = bufIn) + fread(bufIn, 1, SIZE_IOBUF, stdin), endIn == nowIn) ? EOF : *(nowIn++)); }",
			"char bufOut[SIZE_IOBUF], *nowOut = bufOut;",
			"void flush()",
			"{",
			"    fwrite(bufOut, 1, nowOut - bufOut, stdout);",
			"    nowOut = bufOut;",
			"}",
			"void putc(const char &ch)",
			"{",
			"    if (nowOut == bufOut + SIZE_IOBUF) {",
			"        flush();",
			"    }",
			"    *(nowOut++) = ch;",
			"}",
			"template <typename _Tp>",
			"void read(_Tp &x)",
			"{",
			"    x = 0;",
			"    short f = 1;",
			"    char ch = getc();",
			"    while (!isdigit(ch)) {",
			"        f = ch == '-' ? -1 : f;",
			"        ch = getc();",
			"    }",
			"    while (isdigit(ch)) {",
			"        x = x * 10 + ch - '0';",
			"        ch = getc();",
			"    }",
			"    x *= f;",
			"}",
			"template <typename _Tp>",
			"_Tp read(void)",
			"{",
			"    _Tp x = 0;",
			"    short f = 1;",
			"    char ch = getc();",
			"    while (!isdigit(ch)) {",
			"        f = ch == '-' ? -1 : f;",
			"        ch = getc();",
			"    }",
			"    while (isdigit(ch)) {",
			"        x = x * 10 + ch - '0';",
			"        ch = getc();",
			"    }",
			"    x *= f;",
			"    return x;",
			"}",
			"template <typename _Tp>",
			"void print(_Tp x)",
			"{",
			"    if (x < 0) {",
			"        putc('-');",
			"        x = -x;",
			"    }",
			"    if (x > 9) {",
			"        print(x / 10);",
			"        x %= 10;",
			"    }",
			"    putc(x + '0');",
			"}",
			"template <typename First, typename... Rest>",
			"void read(First &first, Rest &...rest)",
			"{",
			"    return read(first), read(rest...);",
			"}",
			"template <typename First, typename... Rest>",
			"void print(First first, Rest... rest)",
			"{",
			"    return print(first), print(rest...);",
			"}"
		],
		"description": "快读/快写模板"
	},
	"Sparse Table": {
		"prefix": ["st", "sparseTable", "rmq"],
		"body": [
			"namespace FstLog2 {",
			"    const int MAXN = 1e5 + 5;",
			"    int lg[MAXN];",
			"    void flogInit()",
			"    {",
			"        for (int i = 2; i < MAXN; ++i) {",
			"            lg[i] = lg[i / 2] + 1;",
			"        }",
			"    }",
			"    int flog(int x)",
			"    {",
			"        return x < MAXN ? lg[x] : std::ceil(std::log2(x));",
			"    }",
			"}",
			"template <typename T>",
			"class SparseTable {",
			"    using VT = vector<T>;",
			"    using VVT = vector<VT>;",
			"    using func_type = function<T(const T &, const T &)>;",
			"    VVT ST;",
			"    static T default_func(const T &t1, const T &t2) { return max(t1, t2); }",
			"    func_type op;",
			"",
			"public:",
			"    SparseTable(const vector<T> &v, func_type _func = default_func)",
			"    {",
			"        op = _func;",
			"        int len = v.size(), l1 = FstLog2::flog(len) + 1;",
			"        ST.assign(len, VT(l1, 0));",
			"        for (int i = 0; i < len; ++i) {",
			"            ST[i][0] = v[i];",
			"        }",
			"        for (int j = 1; j < l1; ++j) {",
			"            int pj = (1 << (j - 1));",
			"            for (int i = 0; i + pj < len; ++i) {",
			"                ST[i][j] = op(ST[i][j - 1], ST[i + (1 << (j - 1))][j - 1]);",
			"            }",
			"        }",
			"    }",
			"    T query(int l, int r)",
			"    {",
			"        int lt = r - l + 1;",
			"        int q = FstLog2::flog(lt) - 1;",
			"        return op(ST[l][q], ST[r - (1 << q) + 1][q]);",
			"    }",
			"};"
		],
		"description": "ST表模板"
	},
	"Dsu": {
		"prefix": ["dsu", "merge"],
		"body": [
			"struct dsu {",
			"    std::vector<size_t> pa, size;",
			"    explicit dsu(size_t size_) : pa(size_), size(size_, 1)",
			"    {",
			"        for (int i = 0; i < size_; ++i) {",
			"            pa[i] = i;",
			"        }",
			"    }",
			"    int find(int x)",
			"    {",
			"        return pa[x] == x ? x : pa[x] = find(pa[x]);",
			"    }",
			"    void unite(size_t x, size_t y)",
			"    {",
			"        x = find(x), y = find(y);",
			"        if (x == y) {",
			"            return;",
			"        }",
			"        if (size[x] < size[y]) {",
			"            std::swap(x, y);",
			"        }",
			"        pa[y] = x;",
			"        size[x] += size[y];",
			"    }",
			"};"
		],
		"description": "并查集"
	},
	"Binary Indexed Tree": {
		"prefix": ["bit", "lowbitTree", "BinaryIndexedTree"],
		"body": [],
		"description": "树状数组"
	},
	"qpow": {
		"prefix": ["qpow", "fpow"],
		"body": [
			"long long fpow(long long x, long long y, long long mod){x %= mod;long long res = 1;while (y) {if (y & 1) {res = res * x % mod;}x = x * x % mod;y >>= 1;}return res;}"
		],
		"description": "快速幂模板"
	},
	"geometry": {
		"prefix": ["geo", "math"],
		"body": [
			"const double eps = 1e-7;class Point {public:};class Vector {public:virtual int size() const = 0;virtual double length() const = 0;virtual double operator[](const int &index) const = 0;virtual double &operator[](const int &index) = 0;};class Point2 : public Point {public:double x, y;Point2(double x = 0.0, double y = 0.0) : x(x), y(y) {}Point2(const Point2 &x) : x(x.x), y(x.y) {}~Point2() {}};class Point3 : public Point {public:double x, y, z;Point3(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}Point3(const Point3 &x) : x(x.x), y(x.y), z(x.z) {}~Point3() {}};class Vector2 : public Vector {public:double x, y;Vector2(double a = 0.0, double b = 0.0) : x(a), y(b) {}explicit Vector2(double a) : x(a), y(a) {}Vector2(const Vector2 &other) : x(other.x), y(other.y) {}Vector2(const Point2 &other) : x(other.x), y(other.y) {}~Vector2() = default;friend Vector2 operator+(const Vector2 &a, const Vector2 &b) { return Vector2(a.x + b.x, a.y + b.y); }friend Vector2 operator-(const Vector2 &a, const Vector2 &b) { return Vector2(a.x - b.x, a.y - b.y); }friend Vector2 operator*(const Vector2 &a, const double &b) { return Vector2(a.x * b, a.y * b); }friend Vector2 operator/(const Vector2 &a, const double &b) { return b == 0 ? Vector2() : Vector2(a.x / b, a.y / b); }friend Vector2 operator-(const Vector2 &a) { return Vector2(-a.x, -a.y); }friend Vector2 operator*(const double &a, const Vector2 &b) { return b * a; }Vector2 &operator=(const Vector2 &b) { return x = b.x, y = b.y, *this; }Vector2 &operator+=(const Vector2 &b) { return *this = *this + b; }Vector2 &operator-=(const Vector2 &b) { return *this = *this - b; }Vector2 &operator*=(const double &b) { return *this = *this * b; }Vector2 &operator/=(const double &b) { return *this = *this / b; }bool operator==(const Vector2 &b) { return x == b.x && y == b.y; }bool operator!=(const Vector2 &b) { return x != b.x || y != b.y; }double operator[](const int &index) const { return (&x)[index]; }double &operator[](const int &index) { return (&x)[index]; }double length() const { return ::sqrtl(x * x + y * y); }Vector2 &normalize(){double scale = 1.0 / length();x *= scale, y *= scale;return *this;}Vector2 &clear(){x = 0, y = 0;return *this;}int size() const { return 2; }friend double operator*(const Vector2 &a, const Vector2 &b) { return a.x * b.x + a.y * b.y; }friend double cross(const Vector2 &a, const Vector2 &b) { return a.x * b.y - a.y * b.x; }double cross(const Vector2 &b) const { return x * b.y - y * b.x; }};class Vector3 : public Vector {public:double x, y, z;Vector3(double a = 0.0, double b = 0.0, double c = 0.0) : x(a), y(b), z(c) {}explicit Vector3(double a) : x(a), y(a), z(a) {}Vector3(const Vector3 &other) : x(other.x), y(other.y), z(other.z) {}~Vector3() = default;friend Vector3 operator+(const Vector3 &a, const Vector3 &b) { return Vector3(a.x + b.x, a.y + b.y, a.z + b.z); }friend Vector3 operator-(const Vector3 &a, const Vector3 &b) { return Vector3(a.x - b.x, a.y - b.y, a.z - b.z); }friend Vector3 operator*(const Vector3 &a, const double &b) { return Vector3(a.x * b, a.y * b, a.z * b); }friend Vector3 operator/(const Vector3 &a, const double &b) { return b == 0 ? Vector3() : Vector3(a.x / b, a.y / b, a.z / b); }friend Vector3 operator-(const Vector3 &a) { return Vector3(-a.x, -a.y, -a.z); }friend Vector3 operator*(const double &a, const Vector3 &b) { return b * a; }Vector3 &operator=(const Vector3 &b) { return x = b.x, y = b.y, z = b.z, *this; }Vector3 &operator+=(const Vector3 &b) { return *this = *this + b; }Vector3 &operator-=(const Vector3 &b) { return *this = *this - b; }Vector3 &operator*=(const double &b) { return *this = *this * b; }Vector3 &operator/=(const double &b) { return *this = *this / b; }bool operator==(const Vector3 &b) { return (fabs(x - b.x) < 0.1 && fabs(y - b.y) < 0.1 && fabs(z - b.z) < 0.1); }bool operator!=(const Vector3 &b) { return (fabs(x - b.x) > 0.5 || fabs(y - b.y) > 0.5 || fabs(z - b.z) > 0.5); }double operator[](const int &index) const { return (&x)[index]; }double &operator[](const int &index) { return (&x)[index]; }inline double length() const { return ::sqrtl(x * x + y * y + z * z); }inline Vector3 &normalize(){double scale = 1.0 / length();x *= scale, y *= scale, z *= scale;return *this;}inline Vector3 &clear(){x = 0, y = 0, z = 0;return *this;}inline int size() const { return 3; }friend double operator*(const Vector3 &a, const Vector3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }friend Vector3 cross(const Vector3 &a, const Vector3 &b) { return Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }Vector3 cross(const Vector3 &b) const { return Vector3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }};struct Edge2 {Vector2 start, end;Edge2(const Vector2 &a = {}, const Vector2 &b = {}) : start(a), end(b) {}Edge2(const Edge2 &a) : start(a.start), end(a.end) {}Edge2 &operator=(const Edge2 &a) { return this->start = a.start, this->end = a.end, *this; }friend bool operator==(const Edge2 &a, const Edge2 &b){double t1 = atan2((a.end - a.start).x, (a.end - a.start).y);double t2 = atan2((b.end - b.start).x, (b.end - b.start).y);return fabs(t1 - t2) < eps;}friend bool operator!=(const Edge2 &a, const Edge2 &b){return !(a == b);}};"
		],
		"description": "计算几何模板"
	},
	"andrew": {
		"prefix": ["geoAndrew", "mathAndrew"],
		"body": [
			"std::vector<Vector2> andrew(std::vector<Vector2> mp){static std::vector<int> stk;static std::vector<Vector2> res;static std::vector<bool> used;sort(mp.begin(), mp.end(), [](const Vector2 &x, const Vector2 &y) { return ((x.x != y.x) ? (x.x < y.x) : (x.y < y.y)); });stk.clear();stk.push_back(0);used.resize(mp.size());used.clear();for (int i = 1; i < mp.size(); ++i) {while (stk.size() > 1 && cross((mp[stk.back()] - mp[*(stk.end() - 2)]), (mp[i] - mp[stk.back()])) <= 0) {used[stk.back()] = 0;stk.pop_back();}used[i] = 1;stk.push_back(i);}int tmp = stk.size();for (int i = mp.size() - 2; i >= 0; --i) {if (used[i]) {continue;}while ((int)stk.size() > tmp && cross((mp[stk.back()] - mp[*(stk.end() - 2)]), (mp[i] - mp[stk.back()])) <= 0) {used[stk.back()] = 0;stk.pop_back();}used[i] = 1;stk.push_back(i);}res.clear();for (int i = 1; i < stk.size(); ++i) {res.push_back(mp[stk[i]]);}return res;}"
		],
		"description": "Andrew 求凸包"
	},
	"rotatingCalipers": {
		"prefix": ["geoGetLongest", "mathGetLongest"],
		"body": [
			"double getLongest(std::vector<Vector2> mp){double mx = 0;if (mp.size() < 3) {mx = (mp[0] - mp[1]).length();return mx;}mp.push_back(mp[0]);for (int i = 0, j = 2; i < mp.size() - 1; ++i) {while (abs(cross((mp[i + 1] - mp[i]), (mp[j] - mp[i + 1]))) <= abs(cross((mp[i + 1] - mp[i]), (mp[j % mp.size()] - mp[i + 1])))) {j = j % mp.size();}mx = std::max(mx, std::max((mp[i + 1] - mp[j]).length(), (mp[i] - mp[j]).length()));}return mx;}"
		],
		"description": "旋转卡壳 求直径"
	},
	"getPoint": {
		"prefix": ["geoEdgePoint", "mathEdgePoint"],
		"body": [
			"Vector2 getPoint(const Edge2 &a, const Edge2 &b){Vector2 p1 = a.start, p2 = b.start, v1 = a.end, v2 = b.end;v1 = v1 - p1;v2 = v2 - p2;Vector2 u = p2 - p1;Vector2 p = p2 + (cross(u, v1) / cross(v1, v2)) * v2;return p;}"
		],
		"description": "求向量交点"
	},
	"segtree": {
		"prefix": ["segmentTree"],
		"body": [
			"template <typename Tp> struct operAdd {auto operator()(const Tp &a, const Tp &b) -> Tp { return a + b; }};template <typename Tp, typename SizeT> struct operRAdd {auto operator()(const Tp &a, const SizeT &l, const SizeT &r) -> Tp {return a * (r - l + 1);}};",
			"/* 汎用高性能セグメントツリーテンプレート!",
			" *",
			" * 模板参数 (template arguments):",
			" * - `Tp = int                            // 节点值类型`",
			" * - `SizeT = int                         // 区间值类型`",
			" * - `OperPushup = operAdd<Tp>            // 节点合并运算`",
			" * - `OperModify = operAdd<Tp>            // 节点修改运算`",
			" * - `OperRange = operRAdd<Tp, SizeT>     // 区间修改运算`",
			" * - `defaultVal = 0                      // 节点默认值`",
			" * - `defaultTag = 0                      // 标记默认值`",
			" *",
			" * 静态方法 (static method):",
			" * - `NewNode() -> segtree<Tp, SizeT, OperPushup, OperModify, OperRange,",
			" * defaultVal, defaultTag> *              // 返回一个指向新节点的指针`",
			" *",
			" * 实例方法 (instance method):",
			" * - `build(SizeT, SizeT, [optional]Tp*)  // 传入区间大小，建立线段树`",
			" * - `modify(SizeT, Tp)                   // 传入位置和值，修改节点`",
			" * - `query(SizeT) -> Tp                  // 传入位置，获取节点的值`",
			" */template <typename Tp = int, typename SizeT = int,  class OperPushup = operAdd<Tp>, class OperModify = operAdd<Tp>,  class OperRange = operRAdd<Tp, SizeT>, int defaultVal = 0,  int defaultTag = 0>struct segtree {  private:SizeT l, r;Tp val;Tp tag;segtree<Tp, SizeT, OperPushup, OperModify, OperRange, defaultVal,defaultTag> *son[2];",
			"#define ls son[0]",
			"#define rs son[1]",
			"#define mid (((l) + (r)) >> 1)",
			"/* 注意：这个函数不会进行边界检查 */void pushdown() {if (this->tag == defaultTag) {return;}ls->val = OperModify()(ls->val, OperRange()(this->tag, ls->l, ls->r));rs->val = OperModify()(rs->val, OperRange()(this->tag, rs->l, rs->r));ls->tag = ls->tag == defaultTag ? this->tag: OperModify()(ls->tag, this->tag);rs->tag = rs->tag == defaultTag ? this->tag: OperModify()(rs->tag, this->tag);this->tag = defaultTag;}/* 注意：这个函数不会进行边界检查 */void pushup() { this->val = OperPushup()(ls->val, rs->val); }",
			" public:segtree(): l(0), r(0), val(defaultVal), tag(defaultTag), son{nullptr, nullptr} {}static auto newNode() -> segtree<Tp, SizeT, OperPushup, OperModify, OperRange, defaultVal, defaultTag> * {constexpr SizeT bufSize = 1e5;static segtree<Tp, SizeT, OperPushup, OperModify, OperRange,   defaultVal, defaultTag> *now =new segtree<Tp, SizeT, OperPushup, OperModify, OperRange,defaultVal, defaultTag>[bufSize];static segtree<Tp, SizeT, OperPushup, OperModify, OperRange,   defaultVal, defaultTag> *end = now + bufSize;return now == end   ? (now = new segtree<Tp, SizeT, OperPushup, OperModify,OperRange, defaultVal,defaultTag>[bufSize],  end = now + bufSize, now)   : ++now;}/* 递归建树，传入的 `l` 应当小于等于 `r` */void build(const SizeT &l, const SizeT &r, const Tp *a = nullptr) {this->l = l;this->r = r;this->tag = defaultTag;if (l >= r) {this->val = a != nullptr ? a[l] : defaultVal;return;}ls = newNode();rs = newNode();ls->build(l, mid, a);rs->build(mid + 1, r, a);this->pushup();}/* 单点修改，传入的 `pos` 应在该节点的区间范围内 */void modify(const SizeT &pos, const Tp &val) {if (this->l == this->r) {this->val = OperModify()(this->val, val);return;}this->pushdown();if (ls->r >= pos) {ls->modify(pos, val);}else {rs->modify(pos, val);}this->pushup();}/* 单点查询，传入的 `pos` 应在该节点的区间范围内 */[[nodiscard]] auto query(const SizeT &pos) -> const Tp & {if (this->l == this->r) {return this->val;}this->pushdown();if (ls->r >= pos) {return ls->query(pos);}return rs->query(pos);}/* 区间修改，传入的 `l`, `r` 应在该节点的区间范围内 */void modify(const SizeT &l, const SizeT &r, const Tp &val) {if (l <= this->l && this->r <= r) {this->tag =this->tag == defaultTag ? val : OperModify()(this->tag, val);this->val =OperModify()(this->val, OperRange()(val, this->l, this->r));return;}this->pushdown();if (ls->r >= l) {ls->modify(l, r, val);}if (rs->l <= r) {rs->modify(l, r, val);}this->pushup();}/* 区间查询，传入的 `l`, `r` 应在该节点的区间范围内 */[[nodiscard]] auto query(const SizeT &l, const SizeT &r) -> Tp {if (l <= this->l && this->r <= r) {return this->val;}this->pushdown();Tp res = defaultVal;if (ls->r >= l) {res = OperPushup()(res, ls->query(l, r));}if (rs->l <= r) {res = OperPushup()(res, rs->query(l, r));}return res;}",
			"#undef mid",
			"#undef ls",
			"#undef rs",
			"};"
		],
		"description": "泛用型高性能线段树模板"
	},
	"modint": {
		"prefix": ["modint"],
		"body": [
			"constexpr int mod = 1e9 + 7;struct getInv {auto operator()(int x) -> int {auto fpow = [](long long x, long long y, long long mod = ::mod) {x %= mod;long long res = 1;while (y != 0) {if ((y & 1) != 0) {res = res * x % mod;}x = x * x % mod;y >>= 1;}return res;};return static_cast<int>(fpow(x, mod - 2));}};template <typename Tp, int mod = ::mod, class getInv = getInv> class modInt {private:Tp val;public:",
			"constexpr modInt() { static_cast<int>(this->val = 0); } /* NOLINT */constexpr modInt(const Tp &val) { /* NOLINT */static_cast<int>(this->val = (val % mod + mod) % mod); /* NOLINT */}constexpr modInt(const modInt<Tp> &x) { /* NOLINT */static_cast<int>(this->val = (x.val % mod + mod) % mod); /* NOLINT */}constexpr modInt(const modInt<Tp> &&x) { /* NOLINT */static_cast<int>(this->val = (x.val % mod + mod) % mod); /* NOLINT */}",
			"auto operator=(const Tp &x) -> modInt & {this->val = (x % mod + mod) % mod;return *this;}auto operator=(const modInt<Tp> &x) -> modInt & {this->val = (x.val % mod + mod) % mod;return *this;}template <class rTp> explicit operator rTp() {return static_cast<rTp>(this->val);}explicit operator bool() { return this->val != 0; }constexpr auto operator-() const -> modInt { return (-this->val) % mod; }constexpr auto operator+(const modInt<Tp> &x) const -> modInt {return ((this->val + x.val) % mod + mod) % mod;}constexpr auto operator-(const modInt<Tp> &x) const -> modInt {return ((this->val - x.val) % mod + mod) % mod;}constexpr auto operator*(const modInt<Tp> &x) const -> modInt {return ((this->val * x.val) % mod + mod) % mod;}constexpr auto operator/(const modInt<Tp> &x) const -> modInt {return ((this->val / x.val) % mod + mod) % mod;}constexpr auto operator%(const modInt<Tp> &x) const -> modInt {return ((this->val % x.val) % mod + mod) % mod;}[[nodiscard]] auto safeDiv(const modInt<Tp> &x) const -> modInt {return (this->val * getInv()(x.val)) % mod;}auto operator+=(const modInt<Tp> &x) -> modInt & {return *this = *this + x;}auto operator-=(const modInt<Tp> &x) -> modInt & {return *this = *this - x;}auto operator*=(const modInt<Tp> &x) -> modInt & {return *this = *this * x;}auto operator/=(const modInt<Tp> &x) -> modInt & {return *this = *this / x;}auto operator%=(const modInt<Tp> &x) -> modInt & {return *this = *this % x;}[[maybe_unused]] auto safeDivEqual(const modInt<Tp> &x) -> modInt & {return *this = this->safeDiv(x);}constexpr auto operator<(const modInt<Tp> &x) const -> bool {return this->val < x.val;}constexpr auto operator==(const modInt<Tp> &x) const -> bool {return this->val == x.val;}constexpr auto operator>(const modInt<Tp> &x) const -> bool {return this->val > x.val;}friend auto operator>>(std::istream &is, modInt<Tp> &x) -> std::istream & {return is >> x.val;}friend auto operator<<(std::ostream &os, const modInt<Tp> &x)-> std::ostream & {return os << x.val;}}; using modint = modInt<int>;"
		]
	},
	"debug": {
		"prefix": ["debug"],
		"body": [
			"namespace __debug { /* NOLINT */",
			"\ttemplate <typename T>",
			"\tinline void _debug(const char *format, T t)",
			"\t{",
			"\tstd::cerr << format << '=' << t << \"\\n\";",
			"\t}",
			"\ttemplate <class First, class... Rest>",
			"\tinline void _debug(const char *format, First first, Rest... rest) {",
			"\twhile (*format != ',') {",
			"\tstd::cerr << *format++;",
			"\t}",
			"\tstd::cerr << '=' << first << \",\";",
			"\t_debug(format + 1, rest...);",
			"\t}",
			"\ttemplate <typename T>",
			"\tauto operator<<(std::ostream &os, const std::vector<T> &V)",
			"\t-> std::ostream & {",
			"\tos << \"[\";",
			"\tfor (const auto &vv : V) {",
			"\tos << vv << \",\";",
			"\t}",
			"\tos << \"]\";",
			"\treturn os;",
			"\t}",
			"",
			"#define debug(...) _debug(#__VA_ARGS__, __VA_ARGS__)",
			"} // namespace __debug"
		]
	},
	"base2": {
		"prefix": ["main2", "base2"],
		"body": [
			"// You know I'm not Alone...",
			"#include <algorithm>",
			"#include <cstdio>",
			"#include <iostream>",
			"#include <string>",
			"#include <vector>",
			"#define STANDING_IO      // 当需要文件输入输出时，注释这一行",
			"#define NO_DEBUG_INFO // 当需要调试信息时，注释这一行",
			"// clang-format off",
			"class fastIO{private:char ibuf[50007],*p1=ibuf,*p2=ibuf,obuf[50007],*p3=obuf,sta[50];bool file_end=false;char get(){return p1==p2&&(p2=(p1=ibuf)+fread(ibuf,1,50007,stdin),p1==p2)?(file_end=true),char(EOF):*p1++;}void put(const char x){p3-obuf<50007?*p3++=x:(fwrite(obuf,p3-obuf,1,stdout),p3=obuf,*p3++=x);}public:explicit operator bool(){return!file_end;}size_t flush(){size_t f=fwrite(obuf,p3-obuf,1,stdout);p3=obuf;*p3=0;return f;}fastIO&operator>>(char&t){for(t=get();!isgraph(t);t=get()){;}return*this;}template<typename any>typename std::enable_if<std::is_same<any,char>::value,any>::type tpval(){char t;for(t=get();!isgraph(t);t=get()){;}return t;}fastIO&operator>>(char*t){char c;for(c=get();!isgraph(c);c=get()){;}for(;isgraph(c);c=get())*t=c,t++;*t=0;return*this;}fastIO&operator>>(std::string&t){t.clear();char c;for(c=get();!isgraph(c);c=get()){;}for(;isgraph(c);c=get())t+=c;return*this;}template<typename any>typename std::enable_if<std::is_same<any,std::string>::value,any>::type tpval(){std::string t;char c;for(c=get();!isgraph(c);c=get()){;}for(;isgraph(c);c=get())t+=c;return t;}template<typename any>typename std::enable_if<(std::is_signed<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__int128_t>::value,fastIO>::type&operator>>(any&t){t=0;bool y=0;char c=get();for(;!isdigit(c);c=get())if(c==45)y=true;for(;isdigit(c);c=get())t=t*10+c-48;if(y==1)t=-t;return*this;}template<typename any>typename std::enable_if<(std::is_signed<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__int128_t>::value,any>::type tpval(){any t=0;bool y=0;char c=get();for(;!isdigit(c);c=get())if(c==45)y=true;for(;isdigit(c);c=get())t=t*10+c-48;if(y==1)t=-t;return t;}template<typename any>typename std::enable_if<(std::is_unsigned<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__uint128_t>::value,fastIO>::type&operator>>(any&t){t=0;char c=get();for(;!isdigit(c);c=get()){;}for(;isdigit(c);c=get())t=t*10+c-48;return*this;}template<typename any>typename std::enable_if<(std::is_unsigned<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__uint128_t>::value,any>::type tpval(){any t=0;char c=get();for(;!isdigit(c);c=get()){;}for(;isdigit(c);c=get())t=t*10+c-48;return t;}template<typename any1,typename any2>fastIO&operator>>(std::pair<any1,any2>&t){return*this>>t.first>>t.second;}template<typename any1,typename any2>std::pair<any1,any2>tpval(){return std::pair<any1,any2>(tpval<any1>(),tpval<any2>());}template<typename any>fastIO&read(any&t){return*this>>t;}fastIO&read(char*t){char c;for(c=get();!isgraph(c);c=get()){;}for(;isgraph(c);c=get())*t=c,t++;*t=0;return*this;}template<typename any,typename...args>fastIO&read(any&t1,args&...t2){return(*this>>t1).read(t2...);}fastIO&operator<<(const char t){put(t);return*this;}fastIO&operator<<(const char*t){for(;*t;t++)put(*t);return*this;}fastIO&operator<<(const std::string&t){for(const char it:t)put(it);return*this;}template<typename any>typename std::enable_if<(std::is_signed<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__int128_t>::value,fastIO>::type&operator<<(any t){if(!t){put(48);return*this;}int len=0;if(t<0)t=-t,put(45);while(t)sta[len++]=char(t%10+48),t/=10;while(len--)put(sta[len]);return*this;}template<typename any>typename std::enable_if<(std::is_unsigned<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__uint128_t>::value,fastIO>::type&operator<<(any t){if(!t){put(48);return*this;}int len=0;while(t)sta[len++]=char(t%10+48),t/=10;while(len--)put(sta[len]);return*this;}template<typename any1,typename any2>fastIO&operator<<(const std::pair<any1,any2>&t){return*this<<t.first<<' '<<t.second;}template<typename any>fastIO&write(const any&t){return*this<<t;}template<typename any,typename...args>fastIO&write(const any&t1,const args&...t2){return(*this<<t1).write(t2...);}~fastIO(){fwrite(obuf,p3-obuf,1,stdout);}}fio; // NOLINT",
			"#ifndef STANDING_IO",
			"#ifdef _WIN32",
			"#define FILENAME []() -> string { std::string tmp(__FILE__); std::size_t pos = tmp.find_last_of('\\\\\\\\') + 1; return tmp.substr(pos, tmp.find_first_of('.') - pos);}()",
			"#else",
			"#define FILENAME []() -> string { std::string tmp(__FILE__); std::size_t pos = tmp.find_last_of('/') + 1; return tmp.substr(pos, tmp.find_first_of('.') - pos);}()",
			"#endif",
			"#define INPUT_FILENAME (std::string(FILENAME) + \".in\").c_str()",
			"#define OUTPUT_FILENAME (std::string(FILENAME) + \".out\").c_str()",
			"#endif",
			"#ifndef NO_DEBUG_INFO",
			"template <typename T> inline void debugValueOutput(const char *format, T t) {std::cerr << \"\\e[96m\" << format << \"\\e[0m\" << '=' << \"\\e[97;1m\" << t << \"\\e[0m\\n\";}template <class First, class... Rest>inline void debugValueOutput(const char *format, First first, Rest... rest) { std::cerr << \"\\e[96m\"; while (*format != ',') {std::cerr << *format++;}std::cerr << \"\\e[0m\" << '=' << \"\\e[97;1m\" << first << \"\\e[0m,\";debugValueOutput(format + 1, rest...);}template <typename T>auto operator<<(std::ostream &os, const std::vector<T> &V) -> std::ostream & {os << \"[\";for (const auto &vv : V) {os << vv << \",\";}os << \"]\";return os;} // NOLINT",
			"#define PREFIX_DEBUG cerr << \"\\e[1m\" << (__FILE__) << \":\" << (__LINE__) << \"-(\" << __FUNCTION__ << \")-[\\e[33mdebug\\e[0;1m]: \\e[0m\" // NOLINT",
			"#define debug(...) PREFIX_DEBUG, debugValueOutput(#__VA_ARGS__, __VA_ARGS__) // NOLINT",
			"#else",
			"#define debug(...)",
			"#endif",
			"#define runOnce(arg) [&]{ static bool flg = false; if (!flg) { arg(); flg = true; }}()",
			"// clang-format on",
			"using namespace std;",
			"auto main() -> int {",
			"#ifndef STANDING_IO",
			"    freopen(INPUT_FILENAME, \"r\", stdin);",
			"    freopen(OUTPUT_FILENAME, \"w\", stdout);",
			"#endif",
			"    // Your Code Hare...",
			"#ifndef STANDING_IO",
			"    fio.flush();",
			"    fclose(stdin);",
			"    fclose(stdout);",
			"#endif",
			"    return 0;",
			"}",
			""
		]
	},
	"msort": {
		"prefix": ["msort"],
		"body": [
			"template <typename Tp>",
			"auto msort(",
			"    Tp *begin, Tp *end,",
			"    bool (*Compare)(const Tp &, const Tp &) =",
			"        [](const Tp &x, const Tp &y) -> bool { return x < y; }) -> size_t {",
			"    if (begin >= end - 1) {",
			"        return 0;",
			"    }",
			"    Tp *mid = (begin + (end - begin) / 2);",
			"    size_t res = msort(begin, mid) + msort(mid, end);",
			"    Tp *buffer = new Tp[end - begin];",
			"    auto left = begin;",
			"    auto right = mid;",
			"    auto pos = buffer;",
			"    while (left < mid && right < end) {",
			"        if (*left == *right || Compare(*left, *right)) {",
			"            *(pos++) = *(left++);",
			"        }",
			"        else {",
			"            *(pos++) = *(right++);",
			"            res += mid - left;",
			"        }",
			"    }",
			"    while (left < mid) {",
			"        *(pos++) = *(left++);",
			"    }",
			"    while (right < end) {",
			"        *(pos++) = *(right++);",
			"    }",
			"    memmove(begin, buffer, sizeof(Tp) * (end - begin));",
			"    delete[] buffer;",
			"    return res;",
			"}"
		]
	}
}
```

## 常用模板

### main

```cpp
#include <algorithm>
#include <cstdio>
#include <iostream>
using namespace std;
int main()
{
	return 0;
}
```

### main2

```cpp
// You know I'm not Alone...
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#define STANDING_IO   // 当需要文件输入输出时，注释这一行
#define NO_DEBUG_INFO // 当需要调试信息时，注释这一行
// clang-format off
class fastIO{private:char ibuf[50007],*p1=ibuf,*p2=ibuf,obuf[50007],*p3=obuf,sta[50];bool file_end=false;char get(){return p1==p2&&(p2=(p1=ibuf)+fread(ibuf,1,50007,stdin),p1==p2)?(file_end=true),char(EOF):*p1++;}void put(const char x){p3-obuf<50007?*p3++=x:(fwrite(obuf,p3-obuf,1,stdout),p3=obuf,*p3++=x);}public:explicit operator bool(){return!file_end;}size_t flush(){size_t f=fwrite(obuf,p3-obuf,1,stdout);p3=obuf;*p3=0;return f;}fastIO&operator>>(char&t){for(t=get();!isgraph(t);t=get()){;}return*this;}template<typename any>typename std::enable_if<std::is_same<any,char>::value,any>::type tpval(){char t;for(t=get();!isgraph(t);t=get()){;}return t;}fastIO&operator>>(char*t){char c;for(c=get();!isgraph(c);c=get()){;}for(;isgraph(c);c=get())*t=c,t++;*t=0;return*this;}fastIO&operator>>(std::string&t){t.clear();char c;for(c=get();!isgraph(c);c=get()){;}for(;isgraph(c);c=get())t+=c;return*this;}template<typename any>typename std::enable_if<std::is_same<any,std::string>::value,any>::type tpval(){std::string t;char c;for(c=get();!isgraph(c);c=get()){;}for(;isgraph(c);c=get())t+=c;return t;}template<typename any>typename std::enable_if<(std::is_signed<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__int128_t>::value,fastIO>::type&operator>>(any&t){t=0;bool y=0;char c=get();for(;!isdigit(c);c=get())if(c==45)y=true;for(;isdigit(c);c=get())t=t*10+c-48;if(y==1)t=-t;return*this;}template<typename any>typename std::enable_if<(std::is_signed<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__int128_t>::value,any>::type tpval(){any t=0;bool y=0;char c=get();for(;!isdigit(c);c=get())if(c==45)y=true;for(;isdigit(c);c=get())t=t*10+c-48;if(y==1)t=-t;return t;}template<typename any>typename std::enable_if<(std::is_unsigned<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__uint128_t>::value,fastIO>::type&operator>>(any&t){t=0;char c=get();for(;!isdigit(c);c=get()){;}for(;isdigit(c);c=get())t=t*10+c-48;return*this;}template<typename any>typename std::enable_if<(std::is_unsigned<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__uint128_t>::value,any>::type tpval(){any t=0;char c=get();for(;!isdigit(c);c=get()){;}for(;isdigit(c);c=get())t=t*10+c-48;return t;}template<typename any1,typename any2>fastIO&operator>>(std::pair<any1,any2>&t){return*this>>t.first>>t.second;}template<typename any1,typename any2>std::pair<any1,any2>tpval(){return std::pair<any1,any2>(tpval<any1>(),tpval<any2>());}template<typename any>fastIO&read(any&t){return*this>>t;}fastIO&read(char*t){char c;for(c=get();!isgraph(c);c=get()){;}for(;isgraph(c);c=get())*t=c,t++;*t=0;return*this;}template<typename any,typename...args>fastIO&read(any&t1,args&...t2){return(*this>>t1).read(t2...);}fastIO&operator<<(const char t){put(t);return*this;}fastIO&operator<<(const char*t){for(;*t;t++)put(*t);return*this;}fastIO&operator<<(const std::string&t){for(const char it:t)put(it);return*this;}template<typename any>typename std::enable_if<(std::is_signed<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__int128_t>::value,fastIO>::type&operator<<(any t){if(!t){put(48);return*this;}int len=0;if(t<0)t=-t,put(45);while(t)sta[len++]=char(t%10+48),t/=10;while(len--)put(sta[len]);return*this;}template<typename any>typename std::enable_if<(std::is_unsigned<any>::value&&std::is_integral<any>::value&&!std::is_same<any,char>::value)||std::is_same<any,__uint128_t>::value,fastIO>::type&operator<<(any t){if(!t){put(48);return*this;}int len=0;while(t)sta[len++]=char(t%10+48),t/=10;while(len--)put(sta[len]);return*this;}template<typename any1,typename any2>fastIO&operator<<(const std::pair<any1,any2>&t){return*this<<t.first<<' '<<t.second;}template<typename any>fastIO&write(const any&t){return*this<<t;}template<typename any,typename...args>fastIO&write(const any&t1,const args&...t2){return(*this<<t1).write(t2...);}~fastIO(){fwrite(obuf,p3-obuf,1,stdout);}}fio; // NOLINT
#ifndef STANDING_IO
#ifdef _WIN32
#define FILENAME []() -> string { std::string tmp(__FILE__); std::size_t pos = tmp.find_last_of('\\') + 1; return tmp.substr(pos, tmp.find_first_of('.') - pos);}()
#else
#define FILENAME []() -> string { std::string tmp(__FILE__); std::size_t pos = tmp.find_last_of('/') + 1; return tmp.substr(pos, tmp.find_first_of('.') - pos);}()
#endif
#define INPUT_FILENAME (std::string(FILENAME) + ".in").c_str()
#define OUTPUT_FILENAME (std::string(FILENAME) + ".out").c_str()
#endif
#ifndef NO_DEBUG_INFO
template <typename T> inline void debugValueOutput(const char *format, T t) {std::cerr << "\e[96m" << format << "\e[0m" << '=' << "\e[97;1m" << t << "\e[0m\n";}template <class First, class... Rest>inline void debugValueOutput(const char *format, First first, Rest... rest) { std::cerr << "\e[96m"; while (*format != ',') {std::cerr << *format++;}std::cerr << "\e[0m" << '=' << "\e[97;1m" << first << "\e[0m,";debugValueOutput(format + 1, rest...);}template <typename T>auto operator<<(std::ostream &os, const std::vector<T> &V) -> std::ostream & {os << "[";for (const auto &vv : V) {os << vv << ",";}os << "]";return os;} // NOLINT
#define PREFIX_DEBUG cerr << "\e[1m" << (__FILE__) << ":" << (__LINE__) << "-(" << __FUNCTION__ << ")-[\e[33mdebug\e[0;1m]: \e[0m" // NOLINT
#define debug(...) PREFIX_DEBUG, debugValueOutput(#__VA_ARGS__, __VA_ARGS__) // NOLINT
#else
#define debug(...)
#endif
#define runOnce(arg) [&]{ static bool flg = false; if (!flg) { arg(); flg = true; }}()
// clang-format on
using namespace std;
auto main() -> int {
#ifndef STANDING_IO
	freopen(INPUT_FILENAME, "r", stdin);
	freopen(OUTPUT_FILENAME, "w", stdout);
#endif
	// Your Code Hare...
#ifndef STANDING_IO
	fio.flush();
	fclose(stdin);
	fclose(stdout);
#endif
	return 0;
}
```

### fpow

```cpp
long long fpow(long long x, long long y, long long mod){x %= mod;long long res = 1;while (y) {if (y & 1) {res = res * x % mod;}x = x * x % mod;y >>= 1;}return res;}
```

### math

```cpp
const double eps = 1e-7;class Point {public:};class Vector {public:virtual int size() const = 0;virtual double length() const = 0;virtual double operator[](const int &index) const = 0;virtual double &operator[](const int &index) = 0;};class Point2 : public Point {public:double x, y;Point2(double x = 0.0, double y = 0.0) : x(x), y(y) {}Point2(const Point2 &x) : x(x.x), y(x.y) {}~Point2() {}};class Point3 : public Point {public:double x, y, z;Point3(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}Point3(const Point3 &x) : x(x.x), y(x.y), z(x.z) {}~Point3() {}};class Vector2 : public Vector {public:double x, y;Vector2(double a = 0.0, double b = 0.0) : x(a), y(b) {}explicit Vector2(double a) : x(a), y(a) {}Vector2(const Vector2 &other) : x(other.x), y(other.y) {}Vector2(const Point2 &other) : x(other.x), y(other.y) {}~Vector2() = default;friend Vector2 operator+(const Vector2 &a, const Vector2 &b) { return Vector2(a.x + b.x, a.y + b.y); }friend Vector2 operator-(const Vector2 &a, const Vector2 &b) { return Vector2(a.x - b.x, a.y - b.y); }friend Vector2 operator*(const Vector2 &a, const double &b) { return Vector2(a.x * b, a.y * b); }friend Vector2 operator/(const Vector2 &a, const double &b) { return b == 0 ? Vector2() : Vector2(a.x / b, a.y / b); }friend Vector2 operator-(const Vector2 &a) { return Vector2(-a.x, -a.y); }friend Vector2 operator*(const double &a, const Vector2 &b) { return b * a; }Vector2 &operator=(const Vector2 &b) { return x = b.x, y = b.y, *this; }Vector2 &operator+=(const Vector2 &b) { return *this = *this + b; }Vector2 &operator-=(const Vector2 &b) { return *this = *this - b; }Vector2 &operator*=(const double &b) { return *this = *this * b; }Vector2 &operator/=(const double &b) { return *this = *this / b; }bool operator==(const Vector2 &b) { return x == b.x && y == b.y; }bool operator!=(const Vector2 &b) { return x != b.x || y != b.y; }double operator[](const int &index) const { return (&x)[index]; }double &operator[](const int &index) { return (&x)[index]; }double length() const { return ::sqrtl(x * x + y * y); }Vector2 &normalize(){double scale = 1.0 / length();x *= scale, y *= scale;return *this;}Vector2 &clear(){x = 0, y = 0;return *this;}int size() const { return 2; }friend double operator*(const Vector2 &a, const Vector2 &b) { return a.x * b.x + a.y * b.y; }friend double cross(const Vector2 &a, const Vector2 &b) { return a.x * b.y - a.y * b.x; }double cross(const Vector2 &b) const { return x * b.y - y * b.x; }};class Vector3 : public Vector {public:double x, y, z;Vector3(double a = 0.0, double b = 0.0, double c = 0.0) : x(a), y(b), z(c) {}explicit Vector3(double a) : x(a), y(a), z(a) {}Vector3(const Vector3 &other) : x(other.x), y(other.y), z(other.z) {}~Vector3() = default;friend Vector3 operator+(const Vector3 &a, const Vector3 &b) { return Vector3(a.x + b.x, a.y + b.y, a.z + b.z); }friend Vector3 operator-(const Vector3 &a, const Vector3 &b) { return Vector3(a.x - b.x, a.y - b.y, a.z - b.z); }friend Vector3 operator*(const Vector3 &a, const double &b) { return Vector3(a.x * b, a.y * b, a.z * b); }friend Vector3 operator/(const Vector3 &a, const double &b) { return b == 0 ? Vector3() : Vector3(a.x / b, a.y / b, a.z / b); }friend Vector3 operator-(const Vector3 &a) { return Vector3(-a.x, -a.y, -a.z); }friend Vector3 operator*(const double &a, const Vector3 &b) { return b * a; }Vector3 &operator=(const Vector3 &b) { return x = b.x, y = b.y, z = b.z, *this; }Vector3 &operator+=(const Vector3 &b) { return *this = *this + b; }Vector3 &operator-=(const Vector3 &b) { return *this = *this - b; }Vector3 &operator*=(const double &b) { return *this = *this * b; }Vector3 &operator/=(const double &b) { return *this = *this / b; }bool operator==(const Vector3 &b) { return (fabs(x - b.x) < 0.1 && fabs(y - b.y) < 0.1 && fabs(z - b.z) < 0.1); }bool operator!=(const Vector3 &b) { return (fabs(x - b.x) > 0.5 || fabs(y - b.y) > 0.5 || fabs(z - b.z) > 0.5); }double operator[](const int &index) const { return (&x)[index]; }double &operator[](const int &index) { return (&x)[index]; }inline double length() const { return ::sqrtl(x * x + y * y + z * z); }inline Vector3 &normalize(){double scale = 1.0 / length();x *= scale, y *= scale, z *= scale;return *this;}inline Vector3 &clear(){x = 0, y = 0, z = 0;return *this;}inline int size() const { return 3; }friend double operator*(const Vector3 &a, const Vector3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }friend Vector3 cross(const Vector3 &a, const Vector3 &b) { return Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }Vector3 cross(const Vector3 &b) const { return Vector3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }};struct Edge2 {Vector2 start, end;Edge2(const Vector2 &a = {}, const Vector2 &b = {}) : start(a), end(b) {}Edge2(const Edge2 &a) : start(a.start), end(a.end) {}Edge2 &operator=(const Edge2 &a) { return this->start = a.start, this->end = a.end, *this; }friend bool operator==(const Edge2 &a, const Edge2 &b){double t1 = atan2((a.end - a.start).x, (a.end - a.start).y);double t2 = atan2((b.end - b.start).x, (b.end - b.start).y);return fabs(t1 - t2) < eps;}friend bool operator!=(const Edge2 &a, const Edge2 &b){return !(a == b);}};
```

### mathAndrew

```cpp
std::vector<Vector2> andrew(std::vector<Vector2> mp){static std::vector<int> stk;static std::vector<Vector2> res;static std::vector<bool> used;sort(mp.begin(), mp.end(), [](const Vector2 &x, const Vector2 &y) { return ((x.x != y.x) ? (x.x < y.x) : (x.y < y.y)); });stk.clear();stk.push_back(0);used.resize(mp.size());used.clear();for (int i = 1; i < mp.size(); ++i) {while (stk.size() > 1 && cross((mp[stk.back()] - mp[*(stk.end() - 2)]), (mp[i] - mp[stk.back()])) <= 0) {used[stk.back()] = 0;stk.pop_back();}used[i] = 1;stk.push_back(i);}int tmp = stk.size();for (int i = mp.size() - 2; i >= 0; --i) {if (used[i]) {continue;}while ((int)stk.size() > tmp && cross((mp[stk.back()] - mp[*(stk.end() - 2)]), (mp[i] - mp[stk.back()])) <= 0) {used[stk.back()] = 0;stk.pop_back();}used[i] = 1;stk.push_back(i);}res.clear();for (int i = 1; i < stk.size(); ++i) {res.push_back(mp[stk[i]]);}return res;}
```

### mathEdgePoint

```cpp
Vector2 getPoint(const Edge2 &a, const Edge2 &b){Vector2 p1 = a.start, p2 = b.start, v1 = a.end, v2 = b.end;v1 = v1 - p1;v2 = v2 - p2;Vector2 u = p2 - p1;Vector2 p = p2 + (cross(u, v1) / cross(v1, v2)) * v2;return p;}
```

### mathGetLongest

```cpp
double getLongest(std::vector<Vector2> mp){double mx = 0;if (mp.size() < 3) {mx = (mp[0] - mp[1]).length();return mx;}mp.push_back(mp[0]);for (int i = 0, j = 2; i < mp.size() - 1; ++i) {while (abs(cross((mp[i + 1] - mp[i]), (mp[j] - mp[i + 1]))) <= abs(cross((mp[i + 1] - mp[i]), (mp[j % mp.size()] - mp[i + 1])))) {j = j % mp.size();}mx = std::max(mx, std::max((mp[i + 1] - mp[j]).length(), (mp[i] - mp[j]).length()));}return mx;}
```

### segmentTree

```cpp
template <typename Tp> struct operAdd {auto operator()(const Tp &a, const Tp &b) -> Tp { return a + b; }};template <typename Tp, typename SizeT> struct operRAdd {auto operator()(const Tp &a, const SizeT &l, const SizeT &r) -> Tp {return a * (r - l + 1);}};
template <typename Tp = int, typename SizeT = int,  class OperPushup = operAdd<Tp>, class OperModify = operAdd<Tp>,  class OperRange = operRAdd<Tp, SizeT>, int defaultVal = 0,  int defaultTag = 0>struct segtree {  private:SizeT l, r;Tp val;Tp tag;segtree<Tp, SizeT, OperPushup, OperModify, OperRange, defaultVal,defaultTag> *son[2];
#define ls son[0]
#define rs son[1]
#define mid (((l) + (r)) >> 1)
/* 注意：这个函数不会进行边界检查 */void pushdown() {if (this->tag == defaultTag) {return;}ls->val = OperModify()(ls->val, OperRange()(this->tag, ls->l, ls->r));rs->val = OperModify()(rs->val, OperRange()(this->tag, rs->l, rs->r));ls->tag = ls->tag == defaultTag ? this->tag: OperModify()(ls->tag, this->tag);rs->tag = rs->tag == defaultTag ? this->tag: OperModify()(rs->tag, this->tag);this->tag = defaultTag;}/* 注意：这个函数不会进行边界检查 */void pushup() { this->val = OperPushup()(ls->val, rs->val); }
 public:segtree(): l(0), r(0), val(defaultVal), tag(defaultTag), son{nullptr, nullptr} {}static auto newNode() -> segtree<Tp, SizeT, OperPushup, OperModify, OperRange, defaultVal, defaultTag> * {constexpr SizeT bufSize = 1e5;static segtree<Tp, SizeT, OperPushup, OperModify, OperRange,   defaultVal, defaultTag> *now =new segtree<Tp, SizeT, OperPushup, OperModify, OperRange,defaultVal, defaultTag>[bufSize];static segtree<Tp, SizeT, OperPushup, OperModify, OperRange,   defaultVal, defaultTag> *end = now + bufSize;return now == end   ? (now = new segtree<Tp, SizeT, OperPushup, OperModify,OperRange, defaultVal,defaultTag>[bufSize],  end = now + bufSize, now)   : ++now;}/* 递归建树，传入的 `l` 应当小于等于 `r` */void build(const SizeT &l, const SizeT &r, const Tp *a = nullptr) {this->l = l;this->r = r;this->tag = defaultTag;if (l >= r) {this->val = a != nullptr ? a[l] : defaultVal;return;}ls = newNode();rs = newNode();ls->build(l, mid, a);rs->build(mid + 1, r, a);this->pushup();}/* 单点修改，传入的 `pos` 应在该节点的区间范围内 */void modify(const SizeT &pos, const Tp &val) {if (this->l == this->r) {this->val = OperModify()(this->val, val);return;}this->pushdown();if (ls->r >= pos) {ls->modify(pos, val);}else {rs->modify(pos, val);}this->pushup();}/* 单点查询，传入的 `pos` 应在该节点的区间范围内 */[[nodiscard]] auto query(const SizeT &pos) -> const Tp & {if (this->l == this->r) {return this->val;}this->pushdown();if (ls->r >= pos) {return ls->query(pos);}return rs->query(pos);}/* 区间修改，传入的 `l`, `r` 应在该节点的区间范围内 */void modify(const SizeT &l, const SizeT &r, const Tp &val) {if (l <= this->l && this->r <= r) {this->tag =this->tag == defaultTag ? val : OperModify()(this->tag, val);this->val =OperModify()(this->val, OperRange()(val, this->l, this->r));return;}this->pushdown();if (ls->r >= l) {ls->modify(l, r, val);}if (rs->l <= r) {rs->modify(l, r, val);}this->pushup();}/* 区间查询，传入的 `l`, `r` 应在该节点的区间范围内 */[[nodiscard]] auto query(const SizeT &l, const SizeT &r) -> Tp {if (l <= this->l && this->r <= r) {return this->val;}this->pushdown();Tp res = defaultVal;if (ls->r >= l) {res = OperPushup()(res, ls->query(l, r));}if (rs->l <= r) {res = OperPushup()(res, rs->query(l, r));}return res;}
#undef mid
#undef ls
#undef rs
};
```

