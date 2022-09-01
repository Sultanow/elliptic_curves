#include <cstring>
#include <cstdint>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <tuple>
#include <bit>
#include <memory>
#include <set>
#include <map>
#include <stdexcept>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <functional>
#include <future>
#include <thread>
#include <sstream>
#include <fstream>

#define ASSERT_MSG(cond, msg) { if (!(cond)) { throw std::runtime_error("Assertion (" #cond ") failed at line " + std::to_string(__LINE__) + "! Msg '" + std::string(msg) + "'"); } }
#define ASSERT(cond) ASSERT_MSG(cond, "")
#define TRY(code) [&]{ try { return (code); } catch (std::exception const & ex) { throw std::runtime_error("Nested exception for code '" #code "' at line " + std::to_string(__LINE__) + ", inner exception:\n" + std::string(ex.what())); } }()
#define LN //{ COUT(<< "LN " << __LINE__ << std::endl << std::flush); }
#define COUT(code) if (1) { auto lock = glog.Locker(); lock.Stream() code; }
#define bit_sizeof(x) (sizeof(x) * 8)

using u8  = uint8_t;
using u32 = uint32_t;
using i64 = int64_t;
using u64 = uint64_t;
using i128 = signed __int128;
using u128 = unsigned __int128;

class AtomicMutex {
public:
    auto Locker() { return std::lock_guard<AtomicMutex>(*this); }
    auto & Flag() { return f_; }
    
    void lock() { while (f_.test_and_set(std::memory_order_acquire)) {} }
    void unlock() { f_.clear(std::memory_order_release); }
private:
    std::atomic_flag f_ = ATOMIC_FLAG_INIT;
};

double Time() {
    static auto const gtb = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - gtb)
        .count();
}

static double constexpr flush_time = 15.0;

class Log {
public:
    Log(std::string const & fname = "ec_pq_find_rational_xy.log")
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

class Timer {
public:
    Timer(std::string const & name = "", bool enable = true)
        : enable_(enable), name_(name), start_(Time()) {}
    double Elapsed() const { return Time() - start_; }
    ~Timer() {
        if (enable_)
            COUT(<< "'" << name_ << "' " << std::fixed << std::setprecision(3)
                << std::setw(7) << Elapsed() << " sec" << std::endl << std::flush);
    }
private:
    bool enable_ = true;
    std::string name_;
    double start_ = 0;
};

class RunInDestr {
public:
    RunInDestr(std::function<void()> const & f) : f_(f) {}
    ~RunInDestr() { f_(); }
private:
    std::function<void()> f_;
};

template <typename T = u32>
T GetPrime(size_t pos) {
    static std::vector<T> primes = {2, 3};
    static AtomicMutex mux;
    if (!mux.Flag().test() && pos < primes.size())
        return primes[pos];
    auto lock = mux.Locker();
    while (pos >= primes.size())
        for (T p = primes.back() + 2;; p += 2) {
            bool is_prime = true;
            for (auto d: primes) {
                if (d * d > p)
                    break;
                if (p % d == 0) {
                    is_prime = false;
                    break;
                }
            }
            if (is_prime) {
                primes.push_back(p);
                break;
            }
        }
    return primes[pos];
}

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

template <typename T>
std::vector<T> const & Factor_TrialDiv_Primes(T n) {
    thread_local std::vector<T> fs0;
    auto & fs = fs0;
    fs.clear();
    if (n <= 1)
        return fs;
    while ((n & 1) == 0) {
        fs.push_back(2);
        n >>= 1;
    }
    for (size_t i = 1;; ++i) {
        auto d = GetPrime(i);
        if (d * d > n)
            break;
        while (n % d == 0) {
            fs.push_back(d);
            n /= d;
        }
    }
    if (n > 1)
        fs.push_back(n);
    return fs;
}

#define Factor Factor_TrialDiv_Primes

class BitVector {
    using Word = u64;
    static size_t constexpr bits = bit_sizeof(Word);
public:
    size_t Size() const { return size_; }
    BitVector & Resize(size_t size, bool fillv = false) {
        v.resize((size + bits - 1) / bits, fillv ? Word(-1) : Word(0));
        size_ = size;
        return *this;
    }
    bool Get(size_t i) const {
        return u8(v[i / bits] >> (i % bits)) & u8(1);
    }
    bool GetC(size_t i) const {
        TRY(CheckInRange(i));
        return Get(i);
    }
    BitVector & Set(size_t i, bool val = true) {
        if (val)
            v[i / bits] |=   Word(1) << (i % bits);
        else
            v[i / bits] &= ~(Word(1) << (i % bits));
        return *this;
    }
    BitVector & SetC(size_t i, bool val = true) {
        TRY(CheckInRange(i));
        return Set(i);
    }
    bool InRange(size_t i) const {
        return i < Size();
    }
private:
    void CheckInRange(size_t i) const {
        ASSERT_MSG(InRange(i), "i " + std::to_string(i) + " Size() " + std::to_string(Size()));
    }
    
    std::vector<Word> v;
    size_t size_ = 0;
};

template <typename To, typename From>
inline To BitCast(From const & f) {
    return std::bit_cast<To>(f);
    static_assert(sizeof(To) == sizeof(From));
    To dst;
    std::memcpy(&dst, &f, sizeof(f));
    return dst;
}

size_t Log2I(u64 n) {
    return ((BitCast<u64>(double(n)) >> 52) & 0x7FF) - 1023;
}

template <typename T>
std::vector<T> Uniquize(std::vector<T> v) {
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return std::move(v);
}

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

BitVector const & CoprimeNums(u64 n, size_t end) {
    static std::unordered_map<u64, std::shared_ptr<BitVector>> m;
    static AtomicMutex mux;
    size_t const num_bits = ((Log2I(std::max<size_t>(end, 1 << 15)) + 1) + 2) / 3 * 3;
    ASSERT(num_bits < 256 && n < (u64(1) << 56));
    u64 const key = n | (u64(num_bits) << 56);
    {
        auto lock = mux.Locker();
        auto it = m.find(key);
        if (it != m.end())
            return *it->second;
    }
    auto const ufs = Uniquize(Factor(n));
    u64 un = 1;
    for (auto f: ufs)
        un *= f;
    u64 const ukey = un | (u64(num_bits) << 56);
    {
        auto lock = mux.Locker();
        auto it = m.find(ukey);
        if (it != m.end())
            return *m.insert({key, it->second}).first->second;
    }
    auto bv = std::make_shared<BitVector>();
    bv->Resize(1 << num_bits, true);
    ASSERT(end <= bv->Size());
    for (auto p: ufs)
        for (size_t i = 0; i < bv->Size(); i += p)
            bv->Set(i, false);
    if constexpr(0)
        for (size_t i = 1; i < bv->Size(); ++i)
            ASSERT_MSG(bv->Get(i) == Coprime(un, i), "i " + std::to_string(i) + " un " + std::to_string(un));
    auto lock = mux.Locker();
    return *m.insert({key, std::move(bv)}).first->second;
}

template <typename T, typename F>
T BinarySearch(F const & f, T begin, T end) {
    auto const end0 = end;
    while (begin + 2 < end) {
        auto const mid = (begin + end - 1) / 2;
        if (f(mid))
            end = mid + 1;
        else
            begin = mid + 1;
    }
    for (auto i = begin; i < end; ++i)
        if (f(i))
            return i;
    return end0;
}

template <typename T, typename F>
T ExpSearch(F const & f, T begin, T end, bool back = false) {
    auto const begin0 = begin, end0 = end;
    T l = 0;
    if (!back) {
        for (l = 1;; l <<= 1) {
            if (begin0 + l >= end0)
                break;
            if (f(begin0 + l - 1))
                break;
            begin = begin0 + l;
        }
        end = std::min(begin0 + l, end0);
    } else {
        for (l = 1;; l <<= 1) {
            if (begin0 + l >= end0)
                break;
            if (!f(end0 - l))
                break;
            end = end0 - l + 1;
        }
        begin = std::max(begin0 + l, end0) - l;
    }
    if (begin >= end || !f(end - 1))
        return end0;
    return BinarySearch(f, begin, end);
}

template <typename T>
struct Result {
    std::set<std::tuple<T, T, T, T>> xy;
    double rtime = 0, xA_coprime_ratio = 0;
    size_t complexity = 0;
};

template <typename T>
Result<T> Task(T pq, T A_begin, T A_end, T B_begin, T B_end, T X_begin, T X_end) {
    // A^3 * y^2 = B^2 * x^3 - A^2 * B^2 * p * q * x
    using UT = std::make_unsigned_t<T>;
    using MT = i128;
    Result<T> result;
    size_t complexity = 0, xA_coprime = 0, xA_noncoprime = 0;
    auto tim = Time();
    RunInDestr finalize([&]{
        tim = Time() - tim;
        result.rtime = tim;
        result.complexity = complexity;
        result.xA_coprime_ratio = xA_coprime / std::max<double>(1, xA_coprime + xA_noncoprime);
    });
    ASSERT(A_begin >= 0 && A_begin < A_end && B_begin >= 0 && B_begin < B_end && X_begin < X_end);
    UT const X_coprime_max = std::max<UT>(1 << 11, std::max<UT>(std::abs(X_begin), std::abs(X_end)));
    //Y_coprime_max = (1ULL << std::llround((std::log2(X_coprime_max) + 0.5) * 2)) - 2;
    
    for (T A = A_begin; A < A_end; ++A) {
        if (A == 0)
            continue;
        MT const A2 = MT(A) * A, A3 = A2 * A;
        auto const & A_coprime = CoprimeNums(UT(A), X_coprime_max);
        
        for (T B = B_begin; B < B_end; ++B) {
            if (B == 0)
                continue;
            MT const B2 = MT(B) * B, A2B2 = A2 * B2, A2B2pq = A2B2 * pq;
            
            for (auto [X_begin0, X_end0]: {
                std::pair<T, T>{X_begin < 0 ? X_begin : 0, X_end < 0 ? X_end : 0},
                std::pair<T, T>{X_begin < 0 ? 0 : X_begin, X_end < 0 ? 0 : X_end},
            }) {
                if (X_begin0 >= X_end0)
                    continue;
                bool const xsign = X_begin0 < 0;
                // https://www.wolframalpha.com/input?i=plot+100+-10+*+x%5E3+%2B+100+*+x
                X_begin0 = BinarySearch<T>([&](auto x){ return - B2 * x * x * x + A2B2pq * x <= 0; }, X_begin0, X_end0);
                if (X_begin0 >= X_end0)
                    continue;
                T constexpr Y_max = 1ULL << 62;
                T ylast = 1, flast = 0;
                bool yfirst = true;
                for (T x = X_begin0; x < X_end0; ++x) {
                    if (x == 0)
                        continue;
                    if (!TRY(A_coprime.GetC(std::abs(x)))) {
                        ++xA_noncoprime;
                        continue;
                    }
                    ++xA_coprime;
                    ++complexity;
                    MT const X3 = MT(x) * x * x, right = - B2 * X3 + A2B2pq * x;
                    auto f = [&](auto y){ return A3 * y * y + right; };
                    T y = 0;
                    if (xsign) {
                        if (yfirst || f(ylast) <= flast) {
                            y = ExpSearch<T>([&](auto y){ return f(y) >= 0; }, ylast, Y_max);
                            ASSERT(y != Y_max);
                        } else {
                            y = ExpSearch<T>([&](auto y){ return f(y) >= 0; }, 1, ylast + 1, true);
                            ASSERT(y != ylast + 1);
                        }
                    } else {
                        if (!yfirst)
                            ASSERT(f(ylast) <= flast);
                        y = ExpSearch<T>([&](auto y){ return f(y) >= 0; }, ylast, Y_max);
                        ASSERT(y != Y_max);
                    }
                    auto const fy = f(y);
                    ASSERT(fy >= 0);
                    ASSERT(y >= 1);
                    if (fy == 0 && Coprime(y, B)) {
                        ASSERT(MT(A) * A * A * y * y == MT(B) * B * x * x * x - MT(A) * A * B * B * pq * x);
                        ASSERT_MSG(Coprime(std::abs(x), A), "x " + std::to_string(x) + " A " + std::to_string(A));
                        ASSERT_MSG(Coprime(y, B), "y " + std::to_string(y) + " B " + std::to_string(B));
                        result.xy.insert({x, A, y, B});
                    }
                    yfirst = false;
                    ylast = y;
                    flast = fy;
                }
            }
        }
    }
    return result;
}

auto const & PQs() {
    static std::vector<std::pair<u64, u64>> pqs;// = {{3, 11}}; return pqs;
    static std::mutex mux;
    std::lock_guard lock(mux);
    if (!pqs.empty())
        return pqs;
    for (size_t pq = 1; pq < (1 << 21); ++pq) {
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
}

size_t NThreads() {
    return std::thread::hardware_concurrency();
}

void Main() {
    COUT(<< "============  Program Start  ============" << std::endl);
    
    size_t constexpr block = 1 << 7;
    double const blocks_per_pq = 0.02;
    size_t const nthr = NThreads();
    auto const gstart_time = Time();
    auto const & pqs = PQs();
    std::vector<std::future<void>> tasks;
    std::set<std::tuple<u64, u64, u64, u64, i64>> added;
    std::set<std::pair<u64, u64>> done_pq;
    std::recursive_mutex gmux;
    size_t done_iPQ = 0;

    std::map<std::pair<u64, u64>, size_t> pqs_set;
    for (size_t i = 0; i < pqs.size(); ++i)
        pqs_set[pqs[i]] = i;
    
    {
        std::ifstream fsols("ec_pq_find_rational_xy_solutions.txt");
        u64 p = 0, q = 0, iA = 0, iB = 0, xy_cnt = 0;
        i64 iX = 0, tmp = 0;
        while (fsols >> p >> q >> iA >> iB >> iX >> xy_cnt) {
            for (size_t i = 0; i < xy_cnt; ++i)
                fsols >> tmp >> tmp >> tmp >> tmp;
            ASSERT(fsols);
            added.insert({p, q, iA, iB, iX});
            if (xy_cnt > 0)
                done_pq.insert({p, q});
            done_iPQ = std::max(done_iPQ, pqs_set.at({p, q}) + 1);
        }
    }
    
    std::ofstream fsols("ec_pq_find_rational_xy_solutions.txt", std::ios::app);
    std::string fsols_buf;
    double fsols_flush_time = Time();
    
    auto MAdd = [&](u64 iPQ, u64 iA, u64 iB, i64 iX) -> bool {
        u64 p = 0, q = 0;
        std::tie(p, q) = pqs.at(iPQ);
        auto const key = std::make_tuple(p, q, iA, iB, iX);
        if (added.count(key))
            return false;
        {
            std::lock_guard lock(gmux);
            if (done_pq.count({p, q}))
                return false;
        }
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
        tasks.push_back(std::async(std::launch::async, [&gmux, &done_pq, &fsols, &fsols_buf, &fsols_flush_time, gstart_time, p, q, iPQ, iA, iB, iX]{
            auto const res = Task<i64>(p * q, block * iA, block * (iA + 1), block * iB, block * (iB + 1),
                block * block * iX, block * block * (iX + 1));
            std::ostringstream ss;
            ss  << std::fixed << std::setw(6) << std::llround(Time() - gstart_time) << " sec <" << std::setw(3) << iPQ << ": p = "
                << std::setw(3) << p << ", q = " << std::setw(3) << q
                << ">[" << std::setw(4) << iA << ", " << std::setw(4) << iB << ", " << std::setw(5) << iX << "]: time " << std::fixed
                << std::setprecision(3) << std::setw(7) << res.rtime << " sec, complexity " << std::setw(7)
                << (std::ostringstream{} << "2^" << std::fixed << std::setprecision(2) << std::setw(5) << std::log2(std::max<size_t>(1, res.complexity))).str()
                << ", xA_coprime_r " << std::setprecision(3) << std::setw(5) << res.xA_coprime_ratio << ", [";
            for (auto const [x, A, y, B]: res.xy)
                ss << "(x = " << x << "/" << A << ", y = " << y << "/" << B << "), ";
            ss << "]";
            if (!res.xy.empty()) {
                std::lock_guard lock(gmux);
                done_pq.insert({p, q});
            }
            COUT(<< ss.str() << std::endl << std::flush);
            
            {
                std::ostringstream ss;
                ss << p << " " << q << " " << iA << " " << iB << " " << iX << " " << res.xy.size();
                for (auto const [x, A, y, B]: res.xy)
                    ss << " " << x << " " << A << " " << y << " " << B;
                ss << std::endl;
                std::lock_guard lock(gmux);
                fsols_buf += ss.str();
                if (Time() - fsols_flush_time >= flush_time) {
                    fsols << fsols_buf << std::flush;
                    fsols_buf.clear();
                    fsols_flush_time = Time();
                }
            }
        }));
        added.insert(key);
        return true;
    };
    
    auto NeedBlocks = [&](auto iPQ){ return std::max<size_t>(1, std::llround((done_iPQ - iPQ + 1) * blocks_per_pq + 0.5)); };
    auto AllSolved = [&]{
        std::lock_guard lock(gmux);
        bool all_solved = true;
        for (size_t i = 0; i < pqs.size(); ++i)
            if (!done_pq.count(pqs[i])) {
                all_solved = false;
                break;
            }
        return all_solved;
    };
    while (!AllSolved()) {
        bool modified = false;
        for (size_t iPQ = 0;; ++iPQ) {
            size_t const need_blocks = NeedBlocks(iPQ);
            for (size_t iA = 0;; ++iA) {
                if (iA >= need_blocks)
                    break;
                for (size_t iB = 0;; ++iB) {
                    if (iB >= need_blocks)
                        break;
                    for (size_t iX = 0;; ++iX) {
                        if (iX >= need_blocks * need_blocks)
                            break;
                        for (bool xsign: {false, true})
                            if (MAdd(iPQ, iA, iB, xsign ? (-i64(iX) - 1) : i64(iX)))
                                modified = true;
                    }
                }
            }
            if (iPQ >= done_iPQ || iPQ + 1 >= pqs.size()) {
                ASSERT(iPQ <= done_iPQ);
                ++done_iPQ;
                if (!modified)
                    break;
                std::lock_guard lock(gmux);
                size_t const done_pqs = std::min(done_iPQ, pqs.size());
                size_t iworst = size_t(-1);
                for (size_t i = 0; i < done_pqs; ++i)
                    if (!done_pq.count(pqs.at(i))) {
                        iworst = i;
                        break;
                    }
                std::ostringstream pq_map0, pq_map;
                for (size_t i = 0; i < done_pqs; ++i) {
                    bool const is_done = done_pq.count(pqs.at(i));
                    pq_map0 << (is_done ? "+" : ".");
                    pq_map << (is_done ? "+" : ".")
                        << pqs.at(i).first << "," << pqs.at(i).second;
                    if (!is_done)
                        pq_map << "(" << NeedBlocks(i) << ")";
                    pq_map << " ";
                }
                COUT(<< std::endl << "---------------  Total PQs " << done_pqs << " (last " << pqs.at(done_pqs - 1).first
                    << "," << pqs.at(done_pqs - 1).second << "), Solved " << done_pq.size() << " ("
                    << std::fixed << std::setprecision(1) << done_pq.size() * 100.0 / done_pqs
                    << "%), Worst " << (iworst == size_t(-1) ? 0 : pqs.at(iworst).first)
                    << "," << (iworst == size_t(-1) ? 0 : pqs.at(iworst).second)
                    << " (" << (iworst == size_t(-1) ? 0 : NeedBlocks(iworst)) << " blocks) ---------------" << std::endl
                    << pq_map0.str() << std::endl << pq_map.str() << std::endl << std::endl << std::flush);
                break;
            }
        }
    }
    
    for (auto & e: tasks)
        e.get();
    tasks.clear();
    fsols << fsols_buf << std::flush;
    fsols_buf.clear();
    glog.Flush();
}

int main() {
    try {
        Main();
        return 0;
    } catch (std::exception const & ex) {
        COUT(<< "Exception:\n" << ex.what() << std::endl);
        return -1;
    }
}