#include <stdlib.h>
#include <stdio.h>
#include "../complex.h"

int main(int agrc, char** argv)
{
	char* op = argv[1];
	double x = atof(argv[2]);
	double y = atof(argv[3]);

	printf("op %s, x %lf, y %lf\n", op, x, y);

	struct Complex z = init(x, y);

	if (*op == 'd')
	{
		printf("dist_sq = %lf\n", dist_sq(z));
	} else if (*op == 's')  
	{
		print_complex(sq(z));
	} else if (*op == 'a')
	{
		print_complex(add(z, z));
	} else
	{
		printf("no matching condition\n");
	}
}
