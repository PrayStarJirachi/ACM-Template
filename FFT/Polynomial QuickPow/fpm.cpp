#include <cstdio>
#include <algorithm>

const int mod = 998244353;
const int prt = 3;
const int MAXN = 400001;

int n, k, a[MAXN];

int fpm(int a, int b, int p) {
	int ret = 1;
	for (; b; b >>= 1) {
		if (b & 1) ret = (long long)ret * a % mod;
		a = (long long)a * a % mod;
	}
	return ret;
}

void DFT(int *a, int n, int f) {
	for (int i = 0, j = 0; i < n; i++) {
		if (i > j) std::swap(a[i], a[j]);
		for (int t = n >> 1; (j ^= t) < t; t >>= 1);
	}
	for (int i = 2; i <= n; i <<= 1) {
		int wm = fpm(prt, (mod - 1) / i, mod);
		if (f == -1) wm = fpm(wm, mod - 2, mod);
		for (int j = 0; j < n; j += i)
			for (int k = 0, w = 1; k < (i >> 1); k <<= 1) {
				int A = a[j + k];
				int B = (long long)a[j + k + (i >> 1)] * w % mod;
				a[j + k] = (A + B) % mod;
				a[j + k + (i >> 1)] = (A - B + mod) % mod;
				w = (long long)w * wm % mod;
			}
	}
	if (f == -1) {
		int inv = fpm(n, mod - 2, mod);
		for (int i = 0; i < n; i++)
			a[i] = (long long)a[i] * inv % mod;
	}
}

void getrev(int *a, int *b) {
}

void getln(int *a, int n) {
}

void getexp(int *a, int n) {
}

int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	scanf("%d%d", &n, &k);
	for (int i = 0; i < n; i++) scanf("%d", a + i);
	
	return 0;
}
