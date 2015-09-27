#include <cstdio>
#include <cstdlib>
#include <ctime>

int main() {
	srand(time(NULL));
	freopen("input.txt", "w", stdout);
	int n = 5000;
	printf("%d\n", n);
	for (int i = 1; i <= n; i++) {
		printf("%d %d %d\n", rand() % 10000, rand() % 10000, rand() % 10000);
	}
	return 0;
}
