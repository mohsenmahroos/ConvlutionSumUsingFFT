// Accepted solution of FOUR ARRAYS problem in May 2022 Long-time Competitive Programming Contest 1 at CodeChef.com
#include "FFT.h"
#include <iostream>
#include <numeric>

int main() {
    using i_limits = std::numeric_limits<size_t>;
    const size_t M = 4;
    ivector array_size(M), f_min(M,i_limits::max()), f_max(M,i_limits::min());
    std::cin.tie(nullptr)->sync_with_stdio(false);
    for (size_t i = 0; i < M; ++i)
        std::cin >> array_size[i];
    using ull = unsigned long long;
    ull K, ans = std::numeric_limits<ull>::max(); std::cin >> K;
    const FFT_engine<18> FFT;
    std::vector<ivector> f(M,ivector(FFT.N));
    for (size_t i = 0; i < M; ++i)
        for (size_t elem; array_size[i] > 0; ++f[i][elem])
            std::cin >> elem, --array_size[i],
            f_min[i] = std::min(f_min[i],elem),
            f_max[i] = std::max(f_max[i],elem);
    auto g = FFT.convolution_sum(f[0],f[1]);
    auto h = FFT.convolution_sum(f[2],f[3]);
    for (size_t y = 1; y < FFT.N; ++y)
        h[y] += h[y-1];
    const auto y_index = [&](ull m, ull x) {
        const ull y = FFT.N-1;
        return x == 0 ? y: std::min(y,m/x); };
    const auto p_count = [&](ull m) {
        ull n = 0;
        for (size_t x = 0; x < FFT.N and n < K; ++x)
            if (g[x] > 0)
                n += 1ll*g[x]*h[y_index(m,x)];
        return n; };
    const auto product_of_sums = [&](const ivector &x) {
        return 1ull*(x[0]+x[1])*(x[2]+x[3]); };
    for (ull l = product_of_sums(f_min), r = product_of_sums(f_max); l <= r; ) {
        const ull m = (l+r)/2, n = p_count(m);
        if (n < K)
            l = m+1;
        else
            r = m-1, ans = std::min(ans,m); }
    std::cout << ans;
	return 0; }
