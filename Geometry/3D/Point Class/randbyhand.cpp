#include <cstdio>
#include <ctime>

int ran() {
	static int x = 18273817;
	return (x += (x << 10) + 187281) & (~0u >> 1);
}

int main() {
	for (int i = 1; i <= 100; i++) {
		printf("%d\n", ran() & 1);
	}
	return 0;
}
