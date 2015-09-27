#include <cmath>
#include <queue>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>

const int MAXN = 5001;
const int MAXF = 100001;
const double eps = 1e-9;

int dcmp(const double &x) {
	return fabs(x) < eps ? 0 : (x > 0 ? 1 : -1);
}

struct TPoint{
	double x, y, z;
	TPoint() {}
	TPoint(double x, double y, double z) : x(x), y(y), z(z) {}
	TPoint operator +(const TPoint &p)const {
		return TPoint(x + p.x, y + p.y, z + p.z);
	}
	TPoint operator -(const TPoint &p)const {
		return TPoint(x - p.x, y - p.y, z - p.z);
	}
	TPoint operator *(const double &p)const {
		return TPoint(x * p, y * p, z * p);
	}
	TPoint operator /(const double &p)const {
		return TPoint(x / p, y / p, z / p);
	}
	bool operator <(const TPoint &p)const {
		int dX = dcmp(x - p.x), dY = dcmp(y - p.y), dZ = dcmp(z - p.z);
		return dX < 0 || (dX == 0 && (dY < 0 || (dY == 0 && dZ < 0)));
	}
	bool read() {
		return scanf("%lf%lf%lf", &x, &y, &z) == 3;
	}
};

double sqrdist(const TPoint &a) {
	double ret = 0;
	ret += a.x * a.x;
	ret += a.y * a.y;
	ret += a.z * a.z;
	return ret;
}
double sqrdist(const TPoint &a, const TPoint &b) {
	double ret = 0;
	ret += (a.x - b.x) * (a.x - b.x);
	ret += (a.y - b.y) * (a.y - b.y);
	ret += (a.z - b.z) * (a.z - b.z);
	return ret;
}
double dist(const TPoint &a) {
	return sqrt(sqrdist(a));
}
double dist(const TPoint &a, const TPoint &b) {
	return sqrt(sqrdist(a, b));
}
TPoint det(const TPoint &a, const TPoint &b) {
	TPoint ret;
	ret.x = a.y * b.z - b.y * a.z;
	ret.y = a.z * b.x - b.z * a.x;
	ret.z = a.x * b.y - b.x * a.y;
	return ret;
}
double dot(const TPoint &a, const TPoint &b) {
	double ret = 0;
	ret += a.x * b.x;
	ret += a.y * b.y;
	ret += a.z * b.z;
	return ret;
}
double detdot(const TPoint &a, const TPoint &b, const TPoint &c, const TPoint &d) {
	return dot(det(b - a, c - a), d - a);
}

struct Triangle{
	TPoint a, b, c;
	Triangle() {}
	Triangle(TPoint a, TPoint b, TPoint c) : a(a), b(b), c(c) {}
	double getArea() {
		TPoint ret = det(b - a, c - a);
		return dist(ret) / 2.0;
	}
};

namespace Convex_Hull {
	struct Face{
		int a, b, c;
		bool isOnConvex;
		Face() {}
		Face(int a, int b, int c) : a(a), b(b), c(c) {}
	};

	int nFace, left, right, whe[MAXN][MAXN];
	Face queue[MAXF], tmp[MAXF];

	bool isVisible(const std::vector<TPoint> &p, const Face &f, const TPoint &a) {
		return dcmp(detdot(p[f.a], p[f.b], p[f.c], a)) > 0;
	}

	bool init(std::vector<TPoint> &p) {
		bool check = false;
		for (int i = 1; i < (int)p.size(); i++) {
			if (dcmp(sqrdist(p[0], p[i]))) {
				std::swap(p[1], p[i]);
				check = true;
				break;
			}
		}
		if (!check) return false;
		check = false;
		for (int i = 2; i < (int)p.size(); i++) {
			if (dcmp(sqrdist(det(p[i] - p[0], p[1] - p[0])))) {
				std::swap(p[2], p[i]);
				check = true;
				break;
			}
		}
		if (!check) return false;
		check = false;
		for (int i = 3; i < (int)p.size(); i++) {
			if (dcmp(detdot(p[0], p[1], p[2], p[i]))) {
				std::swap(p[3], p[i]);
				check = true;
				break;
			}
		}
		if (!check) return false;
		for (int i = 0; i < (int)p.size(); i++)
			for (int j = 0; j < (int)p.size(); j++) {
				whe[i][j] = -1;
			}
		return true;
	}

	void pushface(const int &a, const int &b, const int &c) {
		nFace++;
		tmp[nFace] = Face(a, b, c);
		tmp[nFace].isOnConvex = true;
		whe[a][b] = nFace;
		whe[b][c] = nFace;
		whe[c][a] = nFace;
	}

	bool deal(const std::vector<TPoint> &p, const std::pair<int, int> &now, const TPoint &base) {
		int id = whe[now.second][now.first];
		if (!tmp[id].isOnConvex) return true;
		if (isVisible(p, tmp[id], base)) {
			queue[++right] = tmp[id];
			tmp[id].isOnConvex = false;
			return true;
		}
		return false;
	}

	std::vector<Triangle > getConvex(std::vector<TPoint> &p) {
		static std::vector<Triangle> ret;
		ret.clear();
		if (!init(p)) return ret;
		if (!isVisible(p, Face(0, 1, 2), p[3])) pushface(0, 1, 2); else pushface(0, 2, 1);
		if (!isVisible(p, Face(0, 1, 3), p[2])) pushface(0, 1, 3); else pushface(0, 3, 1);
		if (!isVisible(p, Face(0, 2, 3), p[1])) pushface(0, 2, 3); else pushface(0, 3, 2);
		if (!isVisible(p, Face(1, 2, 3), p[0])) pushface(1, 2, 3); else pushface(1, 3, 2);
		for (int a = 4; a < (int)p.size(); a++) {
			TPoint base = p[a];
			for (int i = 1; i <= nFace; i++) {
				if (tmp[i].isOnConvex && isVisible(p, tmp[i], base)) {
					left = 0, right = 0;
					queue[++right] = tmp[i];
					tmp[i].isOnConvex = false;
					while (left < right) {
						Face now = queue[++left];
						if (!deal(p, std::make_pair(now.a, now.b), base)) pushface(now.a, now.b, a);
						if (!deal(p, std::make_pair(now.b, now.c), base)) pushface(now.b, now.c, a);
						if (!deal(p, std::make_pair(now.c, now.a), base)) pushface(now.c, now.a, a);
					}
					break;
				}
			}
		}
		for (int i = 1; i <= nFace; i++) {
			Face now = tmp[i];
			if (now.isOnConvex) {
				ret.push_back(Triangle(p[now.a], p[now.b], p[now.c]));
			}
		}
		return ret;
	}
};

int n;
std::vector<TPoint> p;
std::vector<Triangle> answer;

int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	scanf("%d", &n);
	for (int i = 1; i <= n; i++) {
		TPoint a;
		a.read();
		p.push_back(a);
	}
	answer = Convex_Hull::getConvex(p);
	double areaCounter = 0.0;
	for (int i = 0; i < (int)answer.size(); i++) {
		areaCounter += answer[i].getArea();
	}
	printf("%.3f\n", areaCounter);
	return 0;
}
