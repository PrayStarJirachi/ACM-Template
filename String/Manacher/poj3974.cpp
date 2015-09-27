#include <cstdio>
#include <cstring>
#include <algorithm>

const int MAXN = 2000010;

int cs, r[MAXN * 2];
char s[MAXN];

int manacher(char *s, int *r) {
	static char t[MAXN * 2];
	int len = strlen(s + 1);
	for (int i = 1; i <= len; i++) {
		t[2 * i - 1] =  s[i];
		t[2 * i] = '#';
	}
	t[2 * len] = '\0';
	for (int i = 1, p = 0; i <= 2 * len - 1; i++) {
		if (p + r[p] - 1 < i) r[i] = 1;
		else r[i] = std::min(r[2 * p - i],  p + r[p] - i);
		while (i + r[i] <= 2 * len - 1 && i - r[i] >= 1 && t[i + r[i]] == t[i - r[i]]) r[i]++;
		if (p + r[p] - 1 < i + r[i] - 1) p = i;
	}
	int ret = 0;
	for (int i = 1; i <= 2 * len - 1; i++) {
		int tmp = r[i];
		if ((i & 1) && (r[i] & 1 ^ 1)) tmp--;
		if ((i & 1 ^ 1) && (r[i] & 1)) tmp--;
		ret = std::max(ret, tmp);
	}
	printf("Case %d: %d\n", ++cs, ret);
}

int main() {
	freopen("input.txt", "r", stdin);
	while (scanf("%s", s + 1) == 1) {
		if (!strcmp(s + 1, "END")) break;
		manacher(s, r);
	}
	return 0;
}
