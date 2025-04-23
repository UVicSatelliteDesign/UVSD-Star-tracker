#include "trees.h"

void print_binary_int64(unsigned long long int N) {
	for (int i = 63; i >= 0; i--) {
		std::cout << ((N >> i) & 1);
	}
}
