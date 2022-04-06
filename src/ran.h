#include <limits.h>

#define URAN_MAX UINT_MAX
#define IRAN_MAX INT_MAX

void init_ran(int iseed);	// Initialize

unsigned int uran();		// Integer 0...URAN_MAX
int iran();			// Integer 0...IRAN_MAX
int iran_sign();		// Integer INT_MIN...IRAN_MAX
double dran();			// Double [0, 1)
double dran_sign();		// Double [-1, 1)
