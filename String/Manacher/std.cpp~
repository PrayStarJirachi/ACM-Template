#include <iostream>  
#include <cstring>  
#include <cstdio>  
#include <cstdlib>  
using namespace std;  
const int MAX = 1000010;  
char s[MAX];  
char ss[MAX<<1];  
int dp[MAX<<1];  
  
int solve(int len)  
{  
    int ans = 0;  
    int right = -1;  
    int id = -1;  
    for(int i=0; i<len; i++)  
    {  
        int r = 1;  
        if(right >= i)  
            r = max(r, min(right-i+1, dp[2*id-i]));  
        while((i-r+1>=0&&i+r-1<len)&&(ss[i-r+1]==ss[i+r-1]))  
            r++;  
        r--;  
        if(i+r-1 > right)  
        {  
            right = i+r-1;  
            id = i;  
        }  
        dp[i] = r;  
        if(ans < r)  
            ans = r;  
    }  
    return ans - 1;  
}  
  
  
int main()  
{ 
	freopen("input.txt", "r", stdin);
    int cas = 1;  
    while(scanf("%s", s) != EOF)  
    {  
        if(strcmp(s, "END") == 0)  
            break;  
        int len = strlen(s);  
        int cnt = 0;  
        for(int i=0; i<len; i++)  
        {  
            ss[cnt++] = '#';  
            ss[cnt++] = s[i];  
        }  
        ss[cnt++] = '#';  
  
        printf("Case %d: %d\n", cas++, solve(cnt));  
    }  
    return 0;  
}  
