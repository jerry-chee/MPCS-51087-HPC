#include <stdio.h>
#include <stdlib.h>


int main(int argc, char** argv)
{
	FILE* fp = fopen("test.txt", "w");
	fprintf(fp, "testing\n");	

}
