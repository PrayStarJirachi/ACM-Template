#include <cstdio>
#include <cstring>
#include <algorithm>

const int MAXN = 100011;
const int MAXT = 100011;

struct Palindromic_Tree{
	int nTree, nStr, last, c[MAXT][26], fail[MAXT], r[MAXN], l[MAXN], s[MAXN];
	int allocate(int len) {
		l[nTree] = len;
		r[nTree] = 0;
		fail[nTree] = 0;
		memset(c[nTree], 0, sizeof(c[nTree]));
		return nTree++;
	}
	void init() {
		nTree = nStr = 0;
		int newEven = allocate(0);
		int newOdd = allocate(-1);
		last = newEven;
		fail[newEven] = newOdd;
		fail[newOdd] = newEven;
		s[0] = -1;
	}
	void add(int x) {
		s[++nStr] = x;
		int nownode = last;
		while (s[nStr - l[nownode] - 1] != s[nStr]) nownode = fail[nownode];
		if (!c[nownode][x]) {
			int newnode = allocate(l[nownode] + 2), &newfail = fail[newnode];
			newfail = fail[nownode];
			while (s[nStr - l[newfail] - 1] != s[nStr]) newfail = fail[newfail];
			newfail = c[newfail][x];
			c[nownode][x] = newnode;
		}
		last = c[nownode][x];
		r[last]++;
	}
	void count() {
		for (int i = nTree - 1; i >= 0; i--) {
			r[fail[i]] += r[i];
		}
	}
}ord, rev;

int cs, pA[MAXN], pB[MAXN];
char s[MAXN];

int main() {
	freopen("input.txt", "r", stdin);
	while (scanf("%s", s + 1) == 1) {
		ord.init();
		rev.init();
		int n = strlen(s + 1);
		//printf("%s\n", s + 1);
		for (int i = 1; i <= n; i++) {
			ord.add(s[i] - 'a');
			pA[i] = ord.l[ord.last];
		}
		for (int i = n; i >= 1; i--) {
			rev.add(s[i] - 'a');
			pB[i] = rev.l[rev.last];
		}
		int ans = 0;
		for (int i = 1; i <= n; i++) {
			ans = std::max(ans, pA[i] + pB[i + 1]);
		}
		printf("%d\n", ans);
	}
	return 0;
}
