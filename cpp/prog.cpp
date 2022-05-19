#include <iostream>
#include <iomanip>
#include <string>
#include <future>
#include <thread>
#include <vector>
#include <mutex>
#include <chrono>

#include <gmpxx.h>

inline double Time() {
    static auto const gtb = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::high_resolution_clock::now() - gtb).count();
}

inline std::string mpz_to_str(mpz_class n, size_t base = 10) {
    size_t const bit_size_hi = mpz_size(n.get_mpz_t()) * GMP_LIMB_BITS;
    std::string s(bit_size_hi + 4, '\x00');
    mpz_get_str(s.data(), base, n.get_mpz_t());
    return s.c_str();
}

int main(int argc, char ** argv) {
    uint64_t const gbegin = std::stoll(argc >= 2 ? argv[1] : "0"),
        gend = std::stoll(argc >= 3 ? argv[2] : "0");
    size_t const nthreads = std::thread::hardware_concurrency();
    std::cout << "Num threads: " << nthreads << std::endl
        << "Begin: " << gbegin << ", End: " << gend << std::endl;

    auto stime = Time();
    if (gbegin < gend) {
        std::mutex cout_mux;
        std::vector<std::future<void>> threads;
        for (size_t ithread = 0; ithread < nthreads; ++ithread)
            threads.emplace_back(std::async(std::launch::async, [&, ithread]{
                uint64_t const block = (gend - gbegin + nthreads - 1) / nthreads,
                    begin = gbegin + ithread * block,
                    end = std::min<size_t>(gbegin + (ithread + 1) * block, gend);
                mpz_class x, t0, t1, y, root, rem;
                for (uint64_t ix = begin; ix < end; ++ix) {
                    x = uint32_t(ix >> 32); x <<= 32; x |= uint32_t(ix);
                    mpz_pow_ui(t0.get_mpz_t(), x.get_mpz_t(), 6);
                    mpz_pow_ui(t1.get_mpz_t(), x.get_mpz_t(), 2);
                    y = t0 - 4 * t1 + 4;
                    mpz_sqrtrem(root.get_mpz_t(), rem.get_mpz_t(), y.get_mpz_t());
                    if (rem == 0) {
                        std::lock_guard<std::mutex> lock(cout_mux);
                        std::cout << "[" << mpz_to_str(x) << ", " << mpz_to_str(root)
                            << ", " << mpz_to_str(y) << "]" << std::endl;
                    }
                }
            }));
    }
    stime = Time() - stime;
    std::cout << "Time " << std::fixed << std::setprecision(3)
        << stime << " sec" << std::endl;
}