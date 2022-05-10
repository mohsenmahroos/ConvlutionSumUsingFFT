// Accepted solution of FOUR ARRAYS problem 
// May 2022 Long-time Competitive Programming Contest 1
// at CodeChef.com

#include "FFT.h"
#include <iostream>
#include <numeric>

using namespace std;
using ilimits = numeric_limits<size_t>;
using ivector = vector<size_t>;
using ull = unsigned long long;

const size_t M = 4;
const FFT_engine<18,ivector> FFT;

size_t array_size[M], f_min[M], f_max[M];
ull K, ans = numeric_limits<ull>::max();
vector<ivector> f(M,ivector(FFT.N));

int main () {
    cin.tie(nullptr)->sync_with_stdio(false);
    for (size_t min = ilimits::min(), max = ilimits::max(), i = 0; i < M; ++i)
        cin >> array_size[i], f_min[i] = max, f_max[i] = min;
    cin >> K;
    for (size_t i = 0; i < M; ++i)
        for (size_t elem; array_size[i] > 0; ++f[i][elem])
            cin >> elem, --array_size[i],
            f_min[i] = min(f_min[i],elem),
            f_max[i] = max(f_max[i],elem);
    auto g = FFT.convolution_sum(f[0],f[1]);
    auto h = FFT.convolution_sum(f[2],f[3]);
    for (size_t y = 1; y < FFT.N; ++y)
        h[y] += h[y-1];
    const auto y_index = [&](ull m, ull x) {
        const ull y = FFT.N-1;
        return x == 0 ? y: min(y,m/x); };
    const auto p_count = [&](ull m) {
        ull n = 0;
        for (size_t x = 0; x < FFT.N and n < K; ++x)
            if (g[x] > 0)
                n += 1ll*g[x]*h[y_index(m,x)];
        return n; };
    const auto prod_of_sums = [&](size_t x[]) {
        return 1ull*(x[0]+x[1])*(x[2]+x[3]); };
    for (ull m, l = prod_of_sums(f_min), r = prod_of_sums(f_max); l <= r; )
        if (m = (l+r)/2, p_count(m) < K)
            l = m+1;
        else
            r = m-1, ans = min(ans,m);
    cout << ans;
    return 0; }
