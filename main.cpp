#include "FFT.h"
#include <iostream>
#include <array>
#include <numeric>
using namespace std;

int main() {
    using ull = unsigned long long;
    using imatrix = vector<ivector>;
    FFT_engine<18> FFT;
    const size_t M = 4, N = FFT.N;
    array<size_t,M> array_size;
    cin.tie(nullptr)->sync_with_stdio(false);
    for (size_t i = 0; i < M; ++i)
        cin >> array_size[i];
    ull K; cin >> K;
    imatrix F(M,ivector(N));
    ivector f_min(M,numeric_limits<size_t>::max());
    ivector f_max(M,numeric_limits<size_t>::min());
    for (size_t i = 0; i < M; ++i)
        for (size_t elem; array_size[i] > 0; ++F[i][elem])
            cin >> elem, --array_size[i],
            f_min[i] = min(f_min[i],elem),
            f_max[i] = max(f_max[i],elem);
    auto G = FFT.convolution_sum(F[0],F[1]);
    auto H = FFT.convolution_sum(F[2],F[3]);
    for (size_t y = 1; y < N; ++y)
        H[y] += H[y-1];
    const auto y_index = [&](ull m, ull x) {
        const ull y = N-1;
        return x == 0 ? y: min(y,m/x); };
    const auto p_count = [&](ull m) {
        ull n = 0;
        for (size_t x = 0; x < N and n < K; ++x)
            if (G[x] > 0)
                n += 1ll*G[x]*H[y_index(m,x)];
        return n; };
    const auto product_of_sums = [&](const ivector& x) {
        return 1ull*(x[0]+x[1])*(x[2]+x[3]); };
    ull ans = numeric_limits<ull>::max();
    for (ull l = product_of_sums(f_min), r = product_of_sums(f_max); l <= r; ) {
        const ull m = (l+r)/2, n = p_count(m);
        if (n < K)
            l = m+1;
        else
            r = m-1, ans = min(ans,m); }
    cout << ans;
	return 0; }
