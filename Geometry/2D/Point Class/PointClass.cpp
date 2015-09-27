#include <cmath>
#include <cstdio>

double dcmp(const double &x) {return fabs(x) < eps ? 0 : (x > 0 ? 1 : -1);}

struct Point{
	double x, y;
	Point() {}
	Point(double x, double y) : x(x), y(y) {}
	Point operator +(const Point &p)const {return Point(x + p.x, y + p.y);}
	Point operator -(const Point &p)const {return Point(x - p.x, y - p.y);}
	Point operator *(const double &p)const {return Point(x * p, y * p);}
	Point operator /(const double &p)const {return Point(x / p, y / p);}
	bool read() {return scanf("%lf%lf", &x, &y) == 2;}
};

double dot(const Point &a, const Point &b) {return a.x * b.x + a.y * b.y;}
double det(const Point &a, const Point &b) {return a.x * b.y - b.x * a.y;}
double sqrdist(const Point &a, const Point &b) {return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);}
double dist(const Point &a, const Point &b) {return sqrt(sqrdist(a, b));}

int main() {
}
