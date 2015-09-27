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
}e[2][MAXN], a[MAXN], b[MAXN], c[MAXN];

int n;

int prepare(int n) {
	int len = 1;
	for (; len <= 2 * n; len <<= 1);
	for (int i = 0; i < len; i++) {
		e[0][i] = Complex(cos(2 * pi * i / len), sin(2 * pi * i / len));
		e[1][i] = Complex(cos(2 * pi * i / len), -sin(2 * pi * i / len));
	}
	return len;
}

void DFT(Complex *a, int n, int f) {
	for (int i = 0, j = 0; i < n; i++) {
		if (i > j) std::swap(a[i], a[j]);
		for (int t = n >> 1; (j ^= t) < t; t >>= 1);
	}
	for (int i = 2; i <= n; i <<= 1)
		for (int j = 0; j < n; j += i)
			for (int k = 0; k < (i >> 1); k++) {
				Complex A = a[j + k];
				Complex B = e[f][n / i * k] * a[j + k + (i >> 1)];
				a[j + k] = A + B;
				a[j + k + (i >> 1)] = A - B;
			}
	if (f == 1) {
		for (int i = 0; i < n; i++)
			a[i].a /= n;
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
	DFT(a, n, 0);
	DFT(b, n, 0);
	for (int i = 0; i < n; i++) c[i] = a[i] * b[i];
	DFT(c, n, 1);
	for (int i = 0; i < n; i++) {
		printf("%.0f%c", c[i].a, " \n"[i == n - 1]);
	}
	return 0;
}

