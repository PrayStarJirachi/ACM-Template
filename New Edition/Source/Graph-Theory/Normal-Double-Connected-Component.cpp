void tarjan(int x){
	dfn[x] = low[x] = ++ind2;
	v[x] = 1;
	for (int i = nt[x]; pt[i]; i = nt[i])
		if (!dfn[pt[i]]){
			tarjan(pt[i]);
			low[x] = min(low[x], low[pt[i]]);
			if (dfn[x] <= low[pt[i]])
				++v[x];
		}
		else
			low[x] = min(low[x], dfn[pt[i]]);
}
int main(){
	for (; ; ){
		scanf("%d%d", &n, &m);
		if (n == 0 && m == 0)
			return 0;
		for (int i = 1; i <= ind; ++i)
			nt[i] = pt[i] = 0;
		ind = n;
		for (int i = 1; i <= ind; ++i)
			last[i] = i;
		for (int i = 1; i <= m; ++i){
			scanf("%d%d", &x, &y);
			++x, ++y;
			edge(x, y), edge(y, x);
		}
		memset(dfn, 0, sizeof(dfn));
		memset(v, 0, sizeof(v));
		ans = num = ind2 = 0;
		for (int i = 1; i <= n; ++i)
			if (!dfn[i]){
				root = i;
				size = 0;
				++num;
				tarjan(i);
				--v[root];
			}
		for (int i = 1; i <= n; ++i)
			if (v[i] + num - 1 > ans)
				ans = v[i] + num - 1;
		printf("%d\n",ans);
	}
	return 0;
}
