#include <cstdio>
#include <cstring>
#include <algorithm>

const int mod = 998244353;
const int prt = 3;
const int MAXN = 400001;

int n, a[MAXN], b[MAXN], c[MAXN];

int fpm(int a, int b, int p) {
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
	for (int i = 2; i <= n; i <<= 1)
		for (int j = 0; j < n; j += i) {
			int w = 1, z = (f == 0) ? ((mod - 1)) / i : ((mod - 1) / i * (i - 1));
			int wm = fpm(prt, z, mod);
			for (int k = 0; k < (i >> 1); k++) {
				int A = a[j + k];
				int B = (long long)a[j + k + (i >> 1)] * w % mod;
				a[j + k] = (A + B) % mod;
				a[j + k + (i >> 1)] = (A - B + mod) % mod;
				w = (long long)w * wm % mod;
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
	scanf("%d", &n);
	for (int i = 0; i < n; i++) scanf("%d", a + i);
	n = prepare(n);
	getInv(a, b, n);
	for (int i = 0; i < 2 * n; i++) printf("%d%c", a[i], " \n"[i == 2 * n - 1]);
	for (int i = 0; i < 2 * n; i++) printf("%d%c", b[i], " \n"[i == 2 * n - 1]);
	n *= 2;
	DFT(a, n, 0);
	DFT(b, n, 0);
	for (int i = 0; i < n; i++) c[i] = (long long)a[i] * b[i] % mod;
	DFT(c, n, 1);
	for (int i = 0; i < n; i++) printf("%d%c", c[i], " \n"[i == n - 1]);
	return 0;
}
