#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
using namespace std;
#define MAXN 100100
#define MAXK 26
#define palin(x, id, c) (s[id - 1 - l[x]] - 'a' == c)
#define MOD 1000000007
#define LL long long

char str1[MAXN], str2[MAXN];
int p1[MAXN], p2[MAXN];

struct Trie
{
    int next[MAXN][MAXK], fail[MAXN], l[MAXN];
    int sz, ro, re, last;
    void Init(char *s)
    {
        sz = 0;
        ro = ++sz, l[ro] = -1, fail[ro] = ro;
        memset(next[sz], 0, sizeof(next[sz]));
        re = ++sz, l[re] =  0, fail[re] = ro;
        memset(next[sz], 0, sizeof(next[sz]));
        last = ro;
        s[0] = '$';
    }
    void Add(char *s, int c, int id, int p[])
    {
        while(!palin(last, id, c)) last = fail[last];
        if (next[last][c]) last = next[last][c];
        else
        {
            int x = last;
            ++sz;
            memset(next[sz], 0, sizeof(next[sz]));
            next[x][c] = sz, l[sz] = l[x] + 2;
            if (x == ro) fail[sz] = re;
            else
            {
                x = fail[x];
                while(!palin(x, id, c)) x = fail[x];
                fail[sz] = next[x][c];
            }
            last = sz;
        }
        p[id] = id - l[last] + 1;
    }
}pt;

int main()
{
	freopen("input.txt", "r", stdin);
    while(gets(str1 + 1))
    {
        int len = strlen(str1 + 1);
        for(int i = 1; i <= len; i++) str2[len - i + 1] = str1[i];

        pt.Init(str1);
        for(int i = 1; i <= len; i++) pt.Add(str1, str1[i] - 'a', i, p1);
        pt.Init(str2);
        for(int i = 1; i <= len; i++) pt.Add(str2, str2[i] - 'a', i, p2);

        int ans = 2;
        for(int i = 1; i < len; i++)
        {
            int len1 = i - p1[i] + 1;
            int len2 = (len - i) - p2[len - i] + 1;
            ans = max(ans, len1 + len2);
        }
        printf("%d\n", ans);
    }
    return 0;
}
