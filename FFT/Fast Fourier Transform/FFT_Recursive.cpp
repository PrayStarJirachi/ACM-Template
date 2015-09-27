#include <cmath>
#include <cstdio>
#include <algorithm>

const int MAXN = 200001;
const double pi = (double)3.1415926535897932384626;

struct Complex{
	double a, b;
	Complex () {}
	Complex (double a, double b) : a(a), b(b) {}
	Complex operator +(const Complex &p)const {return Complex(a + p.a, b + p.b);}
	Complex operator -(const Complex &p)const {return Complex(a - p.a, b - p.b);}
    Complex operator *(const Complex &p)const {return Complex(a * p.a - b * p.b, a * p.b + b * p.a);}
	bool read() {return scanf("%lf%lf", &a, &b) == 2;}
	void print() {printf("A = %7.3f B = %7.3f\n", a, b);}
}e[MAXN], a[MAXN], b[MAXN], aR[MAXN], bR[MAXN], c[MAXN], cR[MAXN];

int n;

int prepare(int n) {
	int len = 1;
	for (; len <= 2 * (n + 1); len <<= 1);
	return len;
}

void getUnitRoot(int n, int sgn) {
	for (int i = 0; i < n; i++) {
		e[i] = Complex(cos(sgn * 2 * pi * i / n), sin(sgn * 2 * pi * i / n));
	}
}

void DFT(Complex *a, Complex *b, int s, int l, int c) {
	if (l == 1) {b[s] = a[s]; return;}
	int m = l >> 1, h = s + m;
	for (int i = 0; i < m; i++) {
		b[s + i] = a[s + (i << 1)];
		b[h + i] = a[s + (i << 1 ^ 1)];
	}
	DFT(b, a, s, m, c << 1);
	DFT(b, a, h, m, c << 1);
	for (int i = 0; i < m; i++) {
		b[s + i] = a[s + i] + e[i * c] * a[h + i];
		b[h + i] = a[s + i] - e[i * c] * a[h + i];
	}
}

int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	scanf("%d", &n);
	for (int i = 0; i < n; i++) {
		int x; scanf("%d", &x);
		a[i] = Complex(x, 0);
	}
	for (int i = 0; i < n; i++) {
		int x; scanf("%d", &x);
		b[i] = Complex(x, 0);
	}
	n = prepare(n);
	getUnitRoot(n, 1);
	DFT(a, aR, 0, n, 1);
	DFT(b, bR, 0, n, 1);
	for (int i = 0; i < n; i++) cR[i] = aR[i] * bR[i];
	getUnitRoot(n, -1);
	DFT(cR, c, 0, n, 1);
	for (int i = 0; i < n; i++) {
		printf("%.0f%c", c[i].a / n, " \n"[i == n - 1]);
	}
	return 0;
}
