#include <cmath>
#include <cstdio>
#include <algorithm>

const int G = 10;
const double eps = 1e-9;
const double pi = (double)3.1415926535897932384626;

const int MAXN = 200001;

int dcmp(const double &x) {return fabs(x) < eps ? 0 : (x > 0 ? 1 : -1);}

struct Point{
	double x, y;
	Point() {}
	Point(double x, double y) : x(x), y(y) {}
	Point operator +(const Point &p)const {return Point(x + p.x, y + p.y);}
	Point operator -(const Point &p)const {return Point(x - p.x, y - p.y);}
	Point operator *(const double &p)const {return Point(x * p, y * p);}
	Point operator /(const double &p)const {return Point(x / p, y / p);}
	bool read() {return scanf("%lf%lf", &x, &y) == 2;}
}c, v, p[MAXN];

double dot(const Point &a, const Point &b) {return a.x * b.x + a.y * b.y;}
double det(const Point &a, const Point &b) {return a.x * b.y - b.x * a.y;}
double sqrdist(const Point &a) {return a.x * a.x + a.y * a.y;}
double sqrdist(const Point &a, const Point &b) {return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);}
double dist(const Point &a, const Point &b) {return sqrt(sqrdist(a, b));}
double dist(const Point &a) {return sqrt(sqrdist(a));}

int n;
double r, h;

double getSectorArea(const Point &a, const Point &b, const double &r) {
	double c = (2.0 * r * r - sqrdist(a, b)) / (2.0 * r * r);
	double alpha = acos(c);
	return r * r * alpha / 2.0;
}

std::pair<double, double> getSolution(const double &a, const double &b, const double &c) {
	double delta = b * b - 4.0 * a * c;
	if (dcmp(delta) < 0) return std::make_pair(0, 0);
	else return std::make_pair((-b - sqrt(delta)) / (2.0 * a), (-b + sqrt(delta)) / (2.0 * a));
}

std::pair<Point, Point> getIntersection(const Point &a, const Point &b, const double &r) {
	Point d = b - a;
	double A = dot(d, d);
	double B = 2.0 * dot(d, a);
	double C = dot(a, a) - r * r;
	std::pair<double, double> s = getSolution(A, B, C);
	return std::make_pair(a + d * s.first, a + d * s.second);
}

double getPointDist(const Point &a, const Point &b) {
	Point d = b - a;
	int sA = dcmp(dot(a, d)), sB = dcmp(dot(b, d));
	if (sA * sB <= 0) return det(a, b) / dist(a, b);
	else return std::min(dist(a), dist(b));
}

double getArea(const Point &a, const Point &b, const double &r) {
	double dA = dot(a, a), dB = dot(b, b), dC = getPointDist(a, b), ans = 0.0;
	if (dcmp(dA - r * r) <= 0 && dcmp(dB - r * r) <= 0) return det(a, b) / 2.0;
	Point tA = a / dist(a) * r;
	Point tB = b / dist(b) * r;
	if (dcmp(dC - r) > 0) return getSectorArea(tA, tB, r);
	std::pair<Point, Point> ret = getIntersection(a, b, r);
	if (dcmp(dA - r * r) > 0 && dcmp(dB - r * r) > 0) {
		ans += getSectorArea(tA, ret.first, r);
		ans += det(ret.first, ret.second) / 2.0;
		ans += getSectorArea(ret.second, tB, r);
		return ans;
	}
	if (dcmp(dA - r * r) > 0) return det(ret.first, b) / 2.0 + getSectorArea(tA, ret.first, r);
	else return det(a, ret.second) / 2.0 + getSectorArea(ret.second, tB, r);
}

double getArea(int n, Point *p, const Point &c, const double r)  {
	double ret = 0.0;
	for (int i = 0; i < n; i++) {
		int sgn = dcmp(det(p[i] - c, p[(i + 1) % n] - c));
		if (sgn > 0) ret += getArea(p[i] - c, p[(i + 1) % n] - c, r);
		else ret -= getArea(p[(i + 1) % n] - c, p[i] - c, r);
	}
	return fabs(ret);
}

int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	while (c.read()) {
		scanf("%lf", &h);
		v.read();
		c = c + v * sqrt(2 * h / G);
		scanf("%lf", &r);
		scanf("%d", &n);
		for (int i = 0; i < n; i++) p[i].read();
		printf("%.2f\n", getArea(n, p, c, r));
	}
	return 0;
}
