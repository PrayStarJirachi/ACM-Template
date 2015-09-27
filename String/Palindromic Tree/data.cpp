#include <cstdio>
#include <cstdlib>
#include <ctime>

int main() {
	srand(time(NULL));
	freopen("input.txt", "w", stdout);
	int len = 10;
	for (int i = 1; i <= len; i++) {
		putchar(rand() % 3 + 'a');
	}
	puts("");
	return 0;
}
