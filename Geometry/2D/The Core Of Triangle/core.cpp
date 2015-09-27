#include <cmath>
#include <cstdio>

const double eps = 1e-9;

double dcmp(const double &x) {return fabs(x) < eps ? 0 : (x > 0 ? 1 : -1);}

struct Point{
	double x, y;
	Point() {}
	Point(double x, double y) : x(x), y(y) {}
	Point operator +(const Point &p)const {return Point(x + p.x, y + p.y);}
	Point operator -(const Point &p)const {return Point(x - p.x, y - p.y);}
	Point operator *(const double &p)const {return Point(x * p, y * p);}
	Point operator /(const double &p)const {return Point(x / p, y / p);}
	double length() {return sqrt(x * x + y * y);}
	bool read() {return scanf("%lf%lf", &x, &y) == 2;}
};

double dot(const Point &a, const Point &b) {return a.x * b.x + a.y * b.y;}
double det(const Point &a, const Point &b) {return a.x * b.y - b.x * a.y;}
double sqrdist(const Point &a, const Point &b) {return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);}
double dist(const Point &a, const Point &b) {return sqrt(sqrdist(a, b));}

Point getIncenter(const Point &a, const Point &b, const Point &c) {
	double p = (a - b).length() + (b - c).length() + (c - a).length();
	return (a * (b - c).length() + b * (c - a).length() + c * (a - b).length()) / p;
}

Point getCircumcenter(const Point &a, const Point &b, const Point &c) {
	Point p = b - a, q = c - a, s(dot(p, p) / 2, dot(q, q) / 2);
	double d = det(p, q);
	return a + Point(det(s, Point(p.y, q.y)), det(Point(p.x, q.x), s)) / d;
}

Point getOrthocenter(const Point &a, const Point &b, const Point &c) {
	return a + b + c - getCircumcenter(a, b, c) * 2.0;
}

int main() {
}
