#include <map>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

const int MAXS = 10001;
const int MAXN = 100001;
const int MAXM = 200001;
const double eps = 1e-9;

__inline int dcmp(const int &x) {return fabs(x) < eps ? 0 : (x > 0 ? 1 : -1);}

struct Point{
	int x, y;
	Point() {}
	Point(int x, int y) : x(x), y(y) {}
	Point operator +(const Point &p)const {return Point(x + p.x, y + p.y);}
	Point operator -(const Point &p)const {return Point(x - p.x, y - p.y);}
	Point operator *(const int &p)const {return Point(x * p, y * p);}
	Point operator /(const int &p)const {return Point(x / p, y / p);}
	bool operator ==(const Point &p)const {return dcmp(x - p.x) == 0 && dcmp(y - p.y) == 0;}
	bool operator <(const Point &p)const {
		if (dcmp(x - p.x) < 0) return true;
		if (dcmp(x - p.x) > 0) return false;
		return dcmp(y - p.y) < 0;
	}
	bool read() {return scanf("%d%d", &x, &y) == 2;}
}z[MAXN], p[MAXN], base;

int dot(const Point &a, const Point &b) {return a.x * b.x + a.y * b.y;}
int det(const Point &a, const Point &b) {return a.x * b.y - b.x * a.y;}
int sqrdist(const Point &a) {return a.x * a.x + a.y * a.y;}
int sqrdist(const Point &a, const Point &b) {return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);}

struct Line{
	Point a, b;
	Line() {}
	Line(Point a, Point b) : a(a), b(b) {}
	bool operator <(const Line &p)const {return a < p.a || (a == p.a && b < p.b);}
	bool operator ==(const Line &p)const {return a == p.a && b == p.b;}
}l[MAXN];

bool iscross(const Line &l, const Line &r) {
	int dRA_L = dcmp(det(r.a - l.a, l.b - l.a));
	int dRB_L = dcmp(det(r.b - l.a, l.b - l.a));
	int dLA_R = dcmp(det(l.a - r.a, r.b - r.a));
	int dLB_R = dcmp(det(l.b - r.a, r.b - r.a));
	return dRA_L * dRB_L <= 0 && dLA_R * dLB_R <= 0;
}

struct Polygon{
	std::vector<Point> p;
	int area;
	bool operator <(const Polygon &z)const {return area < z.area;}
	void init() {
		p.clear();
		area = 0;
	}
	void getArea() {
		for (int i = 0, s = (int)p.size(); i < s; i++) {
			area += det(p[(i + 1) % s], p[i]);
		}
	}
    bool isInside(const Point &x) {
        bool ret = false;
        Line ray = Line(Point(-10001, x.y), x);
        for (int i = 0, s = (int)p.size(); i < s; i++) {
          Line now_line = Line(p[i], p[(i + 1) % s]);
          ret ^= (dcmp(std::min(p[i].y, p[(i + 1) % s].y) - x.y) && iscross(now_line,ray));
        }
        return ret;
    }
}tmp, s[MAXS];

struct Edge{
	int node, next;
}e[MAXM];

int n, m, idP, idL, nState, t, o[MAXM], sL[MAXM], h[MAXS], w[MAXS], c[MAXN], b[MAXN];
std::map<Point, int> mP;
std::map<Line, int> mL, pL;
std::vector<int> g[MAXN];
bool flag[MAXS];

void addedge(int x, int y) {
	//printf("Addedge(%d, %d)\n", x, y);
	t++; e[t] = (Edge){y, h[x]}; h[x] = t;
}

int getQuad(const Point &a) {
	if (dcmp(a.x) > 0 && dcmp(a.y) >= 0) return 1;
	if (dcmp(a.x) <= 0 && dcmp(a.y) > 0) return 2;
	if (dcmp(a.x) < 0 && dcmp(a.y) <= 0) return 3;
	if (dcmp(a.x) >= 0 && dcmp(a.y) < 0) return 4;
	return 0;
}

bool cmp(const int &iA, const int &iB) {
	Point a = p[iA] - base, b = p[iB] - base;
	int aQ = getQuad(a), bQ = getQuad(b);
	if (aQ < bQ) return true;
	if (aQ > bQ) return false;
	int d = det(a, b);
	if (dcmp(d) > 0) return true;
	if (dcmp(d) < 0) return false;
	return sqrdist(a) < sqrdist(b);
}

void getPolygon(int nEdge) {
	static bool v[MAXM];
	std::fill(v + 1, v + nEdge + 1, false);
	for (int i = 1; i <= nEdge; i++) {
		if (v[i]) continue;
		tmp.init();
		Line now = l[i];
		v[mL[now]] = true;
		tmp.p.push_back(now.a);
		while (true) {
			int tot = g[mP[now.b]].size();
			Line rev = Line(now.b, now.a);
			Point nxt = p[g[mP[now.b]][(pL[rev] + 1) % tot]];
			now = Line(now.b, nxt);
			v[mL[now]] = true;
			if (now == l[i]) break;
			tmp.p.push_back(now.a);
		}
		tmp.getArea();
		if (dcmp(tmp.area) > 0) s[++nState] = tmp;
	}
}

int main() {
	freopen("risk.in", "r", stdin);
	freopen("risk.out", "w", stdout);
	scanf("%d%d", &n, &m);
	for (int i = 1; i <= n; i++) z[i].read();
	for (int i = 1; i <= m; i++) {
		Point a, b;
		a.read();
		b.read();
		int &iA = mP[a], &iB = mP[b];
		if (!iA) p[iA = ++idP] = a;
		if (!iB) p[iB = ++idP] = b;
		l[++idL] = Line(a, b); mL[l[idL]] = idL;
		l[++idL] = Line(b, a); mL[l[idL]] = idL;
		g[iA].push_back(iB);
		g[iB].push_back(iA);
	}
	for (int i = 1; i <= idP; i++) {
		base = p[i];
		std::sort(g[i].begin(), g[i].end(), cmp);
		for (int j = 0; j < (int)g[i].size(); j++) {
			pL[Line(p[i], p[g[i][j]])] = j;
		}
	}
	getPolygon(idL);
	std::sort(s + 1, s + nState + 1);
	for (int i = 1; i <= nState; i++)
		for (int j = 0, sz = (int)s[i].p.size(); j < sz; j++) {
			Line now = Line(s[i].p[j], s[i].p[(j + 1) % sz]);
			sL[mL[now]] = i;
		}
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= nState; j++) {
			if (s[j].isInside(z[i])) {
				w[j] = i;
				b[i] = j;
				break;
			}
		}
	}
	for (int i = 1; i <= idL; i++) {
		if (!sL[i]) continue;
		Line rev = Line(l[i].b, l[i].a);
		if (sL[mL[rev]]) {
			addedge(sL[i], sL[mL[rev]]);
			addedge(sL[mL[rev]], sL[i]);
		}
	}
	for (int i = 1; i <= nState; i++) {
		bool isEmpty = false;
		for (int j = 0, sz = (int)s[i].p.size(); j < sz; j++) {
			Line rev = Line(s[i].p[(j + 1) % sz], s[i].p[j]);
			if (!sL[mL[rev]]) isEmpty = true;
			else addedge(i, sL[mL[rev]]);
		}
		for (int j = 1; j <= nState; j++) {
			if (i == j) continue;
			bool check = false;
			if (!check && isEmpty && !flag[i]) {
				if (s[j].isInside(z[w[i]])) {
					flag[i] = true;
					addedge(i, j);
					addedge(j, i);
				}
			}
		}
	}
	for (int i = 1; i <= n; i++) {
		int tot = 0;
		for (int k = h[b[i]]; k; k = e[k].next) c[++tot] = w[e[k].node];
		std::sort(c + 1, c + tot + 1);
		tot = std::unique(c + 1, c + tot + 1) - c - 1;
		//printf("%d\n", b[i]);
		printf("%d", tot);
		for (int j = 1; j <= tot; j++) printf(" %d", c[j]);
		printf("\n");
	}
	return 0;
}
