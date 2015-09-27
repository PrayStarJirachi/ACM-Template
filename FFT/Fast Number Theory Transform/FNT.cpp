#include <cstdio>
#include <algorithm>

const int mod = 998244353;
const int prt = 3;
const int MAXN = 400001;

int n, e[2][MAXN], a[MAXN], b[MAXN], c[MAXN];

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
	for (int i = 0; i <= len; i++) {
		e[0][i] = fpm(prt, (mod - 1) / len * i, mod);
		e[1][i] = fpm(prt, (mod - 1) / len * (len - i), mod);
	}
	return len;
}

void DFT(int *a, int n, int f) {
	for (int i = 0, j = 0; i < n; i++) {
		if (i > j) std::swap(a[i], a[j]);
		for (int t = n >> 1; (j ^= t) < t; t >>= 1);
	}
	for (int i = 2; i <= n; i <<= 1)
		for (int j = 0; j < n; j += i)
			for (int k = 0; k < (i >> 1); k++) {
				int A = a[j + k];
				int B = (long long)a[j + k + (i >> 1)] * e[f][n / i * k] % mod;
				a[j + k] = (A + B) % mod;
				a[j + k + (i >> 1)] = (A - B + mod) % mod;
			}
	if (f == 1) {
		long long rev = fpm(n, mod - 2, mod);
		for (int i = 0; i < n; i++) {
			a[i] = (long long)a[i] * rev % mod;
		}
	}
}

int main() {
	freopen("input.txt", "r", stdin);
	scanf("%d", &n);
	for (int i = 0; i < n; i++) scanf("%d%d", a + i, b + i);
	std::reverse(a, a + n);
	int len = prepare(n);
	DFT(a, len, 0);
	DFT(b, len, 0);
	for (int i = 0; i < len; i++) c[i] = (long long)a[i] * b[i] % mod;
	DFT(c, len, 1);
	for (int i = n - 1; i >= 0; i--) printf("%d\n", c[i]);
	return 0;
}
