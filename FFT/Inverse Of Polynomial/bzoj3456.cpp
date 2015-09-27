#include <cstdio>
#include <cassert>
#include <cstring>
#include <algorithm>

const int mod = 1004535809;
const int prt = 3;
const int MAXN = 800001;

int n, a[MAXN], b[MAXN], aR[MAXN], c[MAXN], fac[MAXN], invfac[MAXN], inv[MAXN];

int fpm(int a, int b, int p) {
	assert(b >= 0);
	int ret = 1;
	for (; b; b >>= 1) {
		if (b & 1) ret = (long long)ret * a % p;
		a = (long long)a * a % p;
	}
	return ret;
}

int prepare(int n) {
	int len = 1;
	for (; len <= 2 * n; len <<= 1);
	return len;
}

void DFT(int *a, int n, int f) {
	for (int i = 0, j = 0; i < n; i++) {
		if (i > j) std::swap(a[i], a[j]);
		for (int t = n >> 1; (j ^= t) < t; t >>= 1);
	}
	for (int i = 2; i <= n; i <<= 1) {
		int z = (f == 0) ? ((mod - 1)) / i : ((mod - 1) / i * (i - 1));
		int wm = fpm(prt, z, mod);
		for (int j = 0; j < n; j += i) {
			for (int k = 0, w = 1; k < (i >> 1); k++) {
				int A = a[j + k];
				int B = (long long)a[j + k + (i >> 1)] * w % mod;
				a[j + k] = (A + B) % mod;
				a[j + k + (i >> 1)] = (A - B + mod) % mod;
				w = (long long)w * wm % mod;
			}
		}
	}
	if (f == 1) {
		long long rev = fpm(n, mod - 2, mod);
		for (int i = 0; i < n; i++) {
			a[i] = (long long)a[i] * rev % mod;
		}
	}
}

void getInv(int *a, int *b, int n) {
	static int tmp[MAXN];
	std::fill(b, b + n, 0);
	b[0] = fpm(a[0], mod - 2, mod);
	for (int c = 1; c <= n; c <<= 1) {
		for (int i = 0; i < c; i++) tmp[i] = a[i];
		std::fill(b + c, b + (c << 1), 0);
		std::fill(tmp + c, tmp + (c << 1), 0);
		DFT(tmp, c << 1, 0);
		DFT(b, c << 1, 0);
		for (int i = 0; i < (c << 1); i++) {
			b[i] = (long long)(2 - (long long)tmp[i] * b[i] % mod + mod) * b[i] % mod;
		}
		DFT(b, c << 1, 1);
		std::fill(b + c, b + (c << 1), 0);
	}
}

int main() {
	freopen("input.txt", "r", stdin);
	scanf("%d", &n);
	fac[0] = invfac[0] = a[0] = inv[1] = 1;
	for (int i = 2; i <= n; i++) {
		inv[i] = (long long)(mod - mod / i) * inv[mod % i] % mod;
	}
	for (int i = 1; i <= n; i++) {
		fac[i] = (long long)fac[i - 1] * i % mod;
		invfac[i] = (long long)invfac[i - 1] * inv[i] % mod;
		int q = fpm(2, ((long long)i * (i - 1) >> 1) % (mod - 1), mod);
		a[i] = (long long)q * invfac[i] % mod;
		b[i] = (long long)q * invfac[i - 1] % mod;
	}
	int len = prepare(n);
	getInv(a, aR, len);
	DFT(aR, len, 0);
	DFT(b, len, 0);
	for (int i = 0; i < len; i++) c[i] = (long long)aR[i] * b[i] % mod;
	DFT(c, len, 1);
	int ans = (long long)c[n] * fac[n - 1] % mod;
	printf("%d\n", ans);
	return 0;
}

