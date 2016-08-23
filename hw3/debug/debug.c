#include <stdlib.h>
#include <stdio.h>
#include "debug_transport.c"
#include "debug_utils.c"
#include "debug_buffer.c"
#include <math.h>

/*
#include "../src/prng.c"
#include "../src/tally.c"
#include "../src/eigenvalue.c"
#include "../src/io.c"
*/

	/*
		 int mype_loc = 20;
		 int n_workers = 64;
		 // transform local coordinates to global coordinates
		 int x_lvl = x_lvl_from_mype(mype_loc, n_workers); 
		 int y_lvl = y_lvl_from_mype(mype_loc, n_workers); 
		 int z_lvl = z_lvl_from_mype(mype_loc, n_workers); 

		 printf("mype_loc:%d, n_workers:%d \n x_lvl:%d, y_lvl:%d, z_lvl:%d\n\n",
		 mype_loc, n_workers, x_lvl, y_lvl, z_lvl);

		 int x_lvl = 1;
		 int y_lvl = 2;
		 int z_lvl = 3;
		 int n_workers = 64;
		 int mype = mype_from_lvl(x_lvl, y_lvl, z_lvl, n_workers);
		 printf("mype:%d, n_workers:%d \n\nx_lvl:%d, y_lvl:%d, z_lvl:%d\n",
		 mype, n_workers, x_lvl, y_lvl, z_lvl);

		 int mype = 61;
		 int n_workers = 64;
		 int surface = Y0;
		 printf("mype:%d, n_workers:%d\n", mype, n_workers);
		 printf("surface:%d, exterior?:%d\n", 
		 surface, exterior_surface(mype, n_workers, surface));	
		 Buffer *buff = init_buffer(-1);
		 print_buffer(buff);
		 putchar('\n');

		 int i;
		 for (i = 1; i <= 1000; i++)
		 {
		 add_particle(buff, &i, (i%10));
		 }
		 print_buffer(buff);	

		 Buffer *search = search_buffer(buff, 11);
		 if (search == NULL)
		 printf("NULL search for mype 11\n");
		 */
int main(int argc, char** argv)
{
}
