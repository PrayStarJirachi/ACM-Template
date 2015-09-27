#include <ctime>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;

const int MAXN = 100001;
const int INF = ~0u >> 2;
const double eps = 1e-12;
const double pi = 3.1415926535897932384626;
const int mod = 0;

struct Point{
	double x, y;
	Point() {}
	Point(double x, double y):x(x), y(y) {}
	Point operator +(const Point &p)const {return Point(x + p.x, y + p.y);}
	Point operator -(const Point &p)const {return Point(x - p.x, y - p.y);}
	Point operator *(const double &p)const {return Point(x * p, y * p);}
	Point operator /(const double &p)const {return Point(x / p, y / p);}
	int read() {return scanf("%lf%lf", &x, &y);}
};

struct Line{
	Point a, b;
	Line() {}
	Line(Point a, Point b):a(a), b(b) {}
};

template <typename numtype>
numtype gcd(numtype a, numtype b) {
	while (b != 0) {
		numtype r = a % b;
		a = b; b = r;
	}
	return a;
}

int dcmp(double x) {return fabs(x) < eps ? 0 : (x > 0 ? 1 : -1);}
double det(Point a, Point b) {return a.x * b.y - b.x * a.y;}
double dot(Point a, Point b) {return a.x * b.x + a.y * b.y;}
Point getLeftNormalVector(const Point &a) {return Point(-a.y, a.x);}
Point getRightNormalVector(const Point &a) {return Point(a.y, -a.x);}
bool operator ==(const Point &a, const Point &b) {return dcmp(a.x - b.x) && dcmp(a.y - b.y);}
double sqrdist(Point a, Point b) {return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);}
double dist(Point a, Point b) {return sqrt(sqrdist(a, b));}
double sqrdist(Point a) {return a.x * a.x + a.y * a.y;}
double dist(Point a) {return sqrt(sqrdist(a));}
double dist(Point x, Line l) {return fabs(det(x - l.a, l.b - l.a)) / dist(l.b - l.a);}

bool is_OnSegment(const Point &x, const Line &l) {
	Point vA = l.a - x, vB = l.b - x;
	if (!dcmp(det(vA, vB))) return false;
	return dcmp(dot(vA, vB)) <= 0;
}

bool is_OnLine(const Point &x, const Line &l) {
	Point vA = l.a - x, vB = l.b - x;
	return dcmp(det(vA, vB));
}

bool is_Colinear(const Line &l, const Line &r) {
	double d = det(r.b - r.a, l.b - l.a);
	return dcmp(d);
}

bool is_Line_Same(const Line &l, const Line &r) {
	Point lv = l.b - l.a;
	Point rv = r.b - r.a;
	return dcmp(lv.x - rv.x) && dcmp(lv.y - rv.y);
}

bool is_Segment_Intersect(const Line &l, const Line &r, bool STANDARD_FLAG = false) {
	double dRA_L = det(r.a - l.a, l.b - l.a);
	double dRB_L = det(r.b - l.a, l.b - l.a);
	double dLA_R = det(l.a - r.a, r.b - r.a);
	double dLB_R = det(l.b - r.a, r.b - r.a);
	if (dcmp(dRA_L * dRB_L) < 0 && dcmp(dLA_R * dLB_R) < 0) return true;
	if (STANDARD_FLAG) return false;
	if (dcmp(dRA_L) == 0 && is_OnSegment(r.a, l)) return true;
	if (dcmp(dRB_L) == 0 && is_OnSegment(r.b, l)) return true;
	if (dcmp(dLA_R) == 0 && is_OnSegment(l.a, r)) return true;
	if (dcmp(dLB_R) == 0 && is_OnSegment(l.b, r)) return true;
	return false;
}

Point getIntersection(const Line &l, const Line &r) { // 两直线交点
	double sA = det(r.a - l.a, l.b - l.a);
	double sB = det(l.b - l.a, r.b - l.a);
	return r.a + (r.b - r.a) * (sA / (sA + sB));
}

int getQuad(const Point &a) {
	if (dcmp(a.x) > 0 && dcmp(a.y) >= 0) return 1;
	if (dcmp(a.x) <= 0 && dcmp(a.y) > 0) return 2;
	if (dcmp(a.x) < 0 && dcmp(a.y) <= 0) return 3;
	if (dcmp(a.x) >= 0 && dcmp(a.y) < 0) return 4;
	return 0;
}

bool Polar_Angle_Comp_Point(const Point &a, const Point &b) {  // 极角排序
	int aQ = getQuad(a), bQ = getQuad(b);
	if (aQ < bQ) return true;
	if (aQ > bQ) return false;
	double d = det(a, b);
	if (dcmp(d) > 0) return true;
	if (dcmp(d) < 0) return false;
	return dist(a) < dist(b);
}

bool Polar_Angle_Comp_Line(const Line &l, const Line &r) {  // 将按直线的向量对直线极角排序
	return Polar_Angle_Comp_Point(l.b - l.a, r.b - r.a);
}

bool Pair_Comp(const Point &a, const Point &b) {
	if (dcmp(a.x - b.x) < 0) return true;
	if (dcmp(a.x - b.x) > 0) return false;
	return dcmp(a.y - b.y) < 0;
}

int Convex_Hull(int n, Point *P, Point *C, bool COLINEAR_LIMIT_FLAG = true) {
	sort(P, P + n, Pair_Comp);
	int top = 0;
	for (int i = 0; i < n; i++) {
		if (COLINEAR_LIMIT_FLAG) while (top >= 2 && dcmp(det(C[top - 1] - C[top - 2], P[i] - C[top - 2])) <= 0) top--;
		else while (top >= 2 && dcmp(det(C[top - 1] - C[top - 2], P[i] - C[top - 2])) < 0) top--;
		C[top++] = P[i];
	}
	int lasttop = top;
	for (int i = n - 1; i >= 0; i--) {
		if (COLINEAR_LIMIT_FLAG) while (top > lasttop && dcmp(det(C[top - 1] - C[top - 2], P[i] - C[top - 2])) <= 0) top--;
		else while (top > lasttop && dcmp(det(C[top - 1] - C[top - 2], P[i] - C[top - 2])) < 0) top--;
		C[top++] = P[i];
	}
	return top;
}

// Half Plane Intersection Begin

bool isOnLeft(const Point &x, const Line &l) {
	double d = det(x - l.a, l.b - l.a);
	return dcmp(d) <= 0;
}

int getIntersectionOfHalfPlane(int n, Line *L, Line *A) {
	Line *q = new Line[n + 1];
	Point *p = new Point[n + 1];
	sort(L, L + n, Polar_Angle_Comp_Line);
	int l = 1, r = 0;
	for (int i = 0; i < n; i++) {
		while (l < r && !isOnLeft(p[r - 1], L[i])) r--;
		while (l < r && !isOnLeft(p[l], L[i])) l++;
		q[++r] = L[i];
		if (l < r && is_Colinear(q[r], q[r - 1])) {
			r--;
			if (isOnLeft(L[i].a, q[r])) q[r] = L[i];
		}
		if (l < r) p[r - 1] = getIntersection(q[r - 1], q[r]);
	}
	while (l < r && !isOnLeft(p[r - 1], q[l])) r--;
	if (r - l + 1 <= 2) return 0;
	int tot = 0;
	for (int i = l; i <= r; i++) A[tot++] = q[i];
	return tot;
}

// Half Plane Intersection End

double getArea(int n, Point *P) {  // 求面积
	double ret = 0.0;
	if (n == 0) return ret;
	for (int i = 0; i < n - 1; i++) ret += det(P[i], P[i + 1]);
	ret += det(P[n - 1], P[0]);
	return fabs(ret / 2.0);
}

bool isClockwise(int n, Point *p) {  // 判断点集P所对应的简单多边形是否是顺时针的
	double ret = 0.0;
	for (int i = 1; i < n - 1; i++)
		ret += det(p[i] - p[0], p[i + 1] - p[0]);
	return dcmp(ret) < 0;
}

bool isInside(const Point &a, int n, Point *P, int MAX_RANGE = 10000) { // 判断点在多边形内
    bool ret = false;
    Line ray = Line(Point(-(MAX_RANGE + 1), a.y), a);
    for (int i = 0; i < n; i++) {
    	Line now_line = Line(P[i], P[(i + 1) % n]);
		ret ^= (dcmp(min(P[i].y, P[(i + 1) % n].y) - a.y) && is_Segment_Intersect(now_line, ray));
	}
	return ret;
}

bool isConvex(int n, Point *p) { // 判断是否是凸包，但要求点集p以逆时针传入
	for (int i = 2; i < n; i++)
		if (dcmp(det(p[i] - p[i - 2], p[i - 1] - p[i - 2])) > 0) return false;
	return true;
}

int getInside(int n, Point *P) {  // 求多边形P内有多少个整数点
	int OnEdge = n;
	double area = getArea(n, P);
	for (int i = 0; i < n - 1; i++) {
		Point now = P[i + 1] - P[i];
		int y = (int)now.y, x = (int)now.x;
		OnEdge += abs(gcd(x, y)) - 1;
	}
	Point now = P[0] - P[n - 1];
	int y = (int)now.y, x = (int)now.x;
	OnEdge += abs(gcd(x, y)) - 1;
	double ret = area - (double)OnEdge / 2 + 1;
	return (int)ret;
}

int main() {
}
