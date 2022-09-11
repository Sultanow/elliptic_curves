#include <cmath>
#include <cstdint>
#include <thread>
#include <future>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <map>
#include <set>
#include <chrono>
#include <mutex>

#define ASSERT_MSG(cond, msg) { if (!(cond)) throw std::runtime_error("Assertion (" #cond ") failed at line " + std::to_string(__LINE__) + "! Msg '" + std::string(msg) + "'"); }
#define ASSERT(cond) ASSERT_MSG(cond, "")
#define COUT(code) { if constexpr(1) { auto lock = glog.Locker(); lock.Stream() code; } }

using i8 = int8_t;
using u8 = uint8_t;
using i16 = int16_t;
using u16 = uint16_t;
using i32 = int32_t;
using u32 = uint32_t;
using i64 = int64_t;
using u64 = uint64_t;
using i128 = signed __int128;
using u128 = unsigned __int128;

static i64 constexpr limit_A = 1 << 16, limit_B = 1 << 11;
//static i64 constexpr limit_A = 1 << 16, limit_B = 1 << 11;

double Time() {
    static auto const gtb = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - gtb)
        .count();
}

double constexpr flush_time = 15.0;

class Log {
public:
    Log(std::string const & fname = "ec_pq_find_ab_cases.log")
        : file_(fname, std::ios::binary | std::ios::app), file_flush_time_(Time()) {}
    ~Log() { Flush(); }
    void Flush(bool force = true) {
        std::lock_guard lock(mux_);
        
        std::cout << ss_.view().substr(done_cout_.size()) << std::flush;
        done_cout_ += ss_.view().substr(done_cout_.size());
        
        if (force || Time() - file_flush_time_ >= flush_time) {
            file_ << ss_.view().substr(done_file_.size()) << std::flush;
            done_file_ += ss_.view().substr(done_file_.size());
            file_flush_time_ = Time();
        }
    }
    class LockerT {
    public:
        LockerT(Log & _log)
            : log_(_log), lock_(std::lock_guard(log_.mux_)) {}
        ~LockerT() { log_.Flush(false); }
        auto & Stream() { return log_.ss_; }
    private:
        Log & log_;
        std::lock_guard<std::recursive_mutex> lock_;
    };
    auto Locker() { return LockerT(*this); }
private:
    std::string done_cout_, done_file_;
    std::ostringstream ss_;
    std::ofstream file_;
    double file_flush_time_ = 0;
    std::recursive_mutex mux_;
    friend class LockerT;
};

Log glog;

template <typename T>
std::vector<T> const & Factor_TrialDiv(T n) {
    thread_local std::vector<T> fs0;
    auto & fs = fs0;
    fs.clear();
    if (n <= 1)
        return fs;
    while ((n & 1) == 0) {
        fs.push_back(2);
        n >>= 1;
    }
    for (T d = 3; d * d <= n; d += 2)
        while (n % d == 0) {
            fs.push_back(d);
            n /= d;
        }
    if (n > 1)
        fs.push_back(n);
    return fs;
}

#define Factor Factor_TrialDiv

template <typename T>
T GCD(T a, T b) {
    using UT = std::make_unsigned_t<T>;
    ASSERT(a >= 0 && b >= 0);
    while (b != 0)
        std::tie(a, b) = std::make_tuple(b, UT(a) % UT(b));
    return a;
}

template <typename T>
bool Coprime(T a, T b) {
    return GCD(a, b) == 1;
}

template <typename T>
auto Task(T const p, T const q) {
    using MT = i128;
    MT const pq = MT(p) * q;
    static auto const ab = []{
        std::vector<std::pair<i32, i16/*, i64, i64*/>> r;
        for (i32 a = /*-limit_A + 1*/ 1; a < limit_A; ++a) {
            if (a % 15000 == 0)
                COUT(<< std::fixed << "a:" << std::llround(double(a /*- (-limit_A + 1)*/ - 1) * 100 / (limit_A)) << "% ");
            if (a == 0)
                continue;
            for (i16 b = 1; b < limit_B; ++b) {
                if (!Coprime(u32(std::abs(a)), u32(b)))
                    continue;
                r.push_back({a, b/*, i64(a) * a, i64(b) * b * b * b*/});
            }
        }
        COUT(<< std::endl);
        return r;
    }();
    static auto const bs = [&]{
        std::vector<std::vector<std::tuple<i16, i16, i64, i64>>> r(limit_B);
        for (i16 n = 1; n < r.size(); ++n)
            for (i16 b1 = 1; b1 <= n; ++b1) {
                if (n % b1 != 0)
                    continue;
                T const b2 = n / b1;
                r[n].push_back({b1, b2, i64(b1) * b1 * b1 * b1, i64(b2) * b2 * b2 * b2});
            }
        COUT(<< "Threads started..." << std::endl);
        return r;
    }();
    auto Pow = [](auto a, size_t b){
        ASSERT(b >= 2);
        MT r = MT(a) * a;
        for (size_t i = 2; i < b; ++i)
            r *= a;
        return r;
    };
    std::map<size_t, std::vector<std::pair<T, T>>> cases;
    for (auto const [_a, _b/*, _ap2, _bp4*/]: ab) {
        T a = 0, b = 0;
        MT ap2 = 0, bp4 = 0;
        //std::tie(a, b, ap2, bp4) = std::make_tuple(_a, _b, _ap2, _bp4);
        std::tie(a, b) = std::make_tuple(_a, _b);
        std::tie(ap2, bp4) = std::make_tuple(MT(_a) * _a, MT(_b) * _b * _b * _b);

        for (i32 a_sign: {1/*, -1*/}) {
            a *= a_sign;

            auto Case = [&](size_t icase, bool cond){
                if (!cond)
                    return;
                cases[icase].push_back({a, b});
            };
            
            Case(1, q == ap2 + p * bp4);
            Case(2, pq == ap2 + bp4);
            Case(3, pq * 4 * bp4 == ap2 + 1);
            Case(4, p == q * bp4 - ap2);
            Case(5, q == p * bp4 - ap2);
            Case(6, pq == bp4 - ap2);
            
            for (auto const [b1, b2, _b1p4, _b2p4]: bs.at(b)) {
                MT b1p4 = 0, b2p4 = 0;
                std::tie(b1p4, b2p4) = std::make_tuple(_b1p4, _b2p4);
                //if (b1 == 1) continue;
                Case(8,  b1p4 * q == ap2 + p * b2p4);
                Case(9,  b1p4 * pq == ap2 + b2p4);
                
                if (b1 < b2) {
                    Case(7,  b1p4 * pq == b2p4 - ap2);
                    Case(10, b1p4 * q == b2p4 * p - ap2);
                }
            }
        }
    };
    return cases;
}

template <typename T>
auto const & PQs() {
    static auto const pqs = []{
        std::vector<std::pair<T, T>> pqs;
        for (T pq = 1; pq < (1 << 21); ++pq) {
            auto fs = Factor(pq);
            if (fs.size() != 2)
                continue;
            if (fs[0] == fs[1])
                continue;
            if (fs[0] == 2 || fs[1] == 2)
                continue;
            if (fs[0] > fs[1])
                std::swap(fs[0], fs[1]);
            pqs.push_back({fs[0], fs[1]});
        }
        return pqs;
    }();
    return pqs;
}

template <typename T>
void Main() {
    COUT(<< "limit_A " << limit_A << ", limit_B " << limit_B << std::endl);
    size_t const nthr = std::thread::hardware_concurrency();
    auto const & pqs = PQs<T>();
    std::vector<std::future<void>> tasks;
    std::set<std::pair<T, T>> const limit_pqs = {
         // /*{11, 29},*/ {11, 43}, /*{5, 163},*/ {11, 79}, /*! {5, 193},*/ /*! {29, 41},*/ {7, 223}, {7, 227}, /*{37, 43},*/
    };
    for (size_t ipq = 0; ipq < pqs.size(); ++ipq) {
        T p = 0, q = 0;
        std::tie(p, q) = pqs.at(ipq);
        if (!limit_pqs.empty() && !limit_pqs.count({p, q}))
            continue;
        while (tasks.size() >= nthr) {
            bool erased = false;
            for (ptrdiff_t i = ptrdiff_t(tasks.size()) - 1; i >= 0; --i)
                if (tasks[i].wait_for(std::chrono::milliseconds(1)) == std::future_status::ready) {
                    tasks[i].get();
                    tasks.erase(tasks.begin() + i);
                    erased = true;
                }
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
            std::this_thread::yield();
        }
        tasks.push_back(std::async(std::launch::async, [&, ipq, p, q]{
            auto res = Task<T>(p, q);
            std::ostringstream ss;
            ss << std::setw(5) << ipq << ": " << std::setw(5) << p << ", " << std::setw(5) << q << "    ";
            if (!res.empty()) {
                ss << "{";
                size_t ires = 0;
                for (auto const & [icase, ab]: res) {
                    ss << icase << ": [";
                    size_t iab = 0;
                    for (auto const [a, b]: ab) {
                        ss << "(" << a << ", " << b << ")";
                        if (iab + 1 < ab.size())
                            ss << ", ";
                        ++iab;
                    }
                    ss << "]";
                    if (ires + 1 < res.size())
                        ss << ", ";
                    ++ires;
                }
                ss << "}";
            }
            ss << std::endl;
            COUT(<< ss.view());
        }));
    }
    glog.Flush();
}

int main() {
    try {
        Main<i64>();
        return 0;
    } catch (std::exception const & ex) {
        COUT(<< "Exception: " << ex.what() << std::endl);
        return -1;
    }
}