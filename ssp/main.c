
/*
 * Subset-Sum Problem (SSP)
 *
 * The main C file
 *
 * last update: September 22nd, 2024
 *
 * AM
 */

#include "ssp.h"

// main
int main()
{
   size_t n = 30;  // total number of integers
   ssp_set_seed(n + 123456);
   double density = 1.0;
   //ssp_load("filename.txt");
   //ssp_t *instance = ssp_gen_free(n,100);
   ssp_t *instance = ssp_gen_hard(n,density);
   ssp_print(instance);
   //ssp_print_matlab(instance);
   printf("Density is %lf\n",ssp_density(instance));
   unsigned long total = ssp_total_sum(instance);
   printf("Total sum is %lu\n",total);
   printf("The instance is feasible: ");
   if (ssp_isfeasible(instance))
      printf("YES!\n");
   else
      printf("NO :-(\n");
   printf("Sum of current subset : %lu\n",ssp_current_sum(instance));
   printf("Selected integers in solution : %lu / %lu\n",ssp_count_selected(instance),n);
   printf("Number of integers larger than target : %lu\n",ssp_count_larger_target(instance));
   ssp_erasesol(instance);
   printf("\n");

   // solving by branch-and-prune
   printf("Branch-and-prune\n");
   time_t start = clock();
   bool found = bp(instance);
   time_t end = clock();
   double time = compute_time(start,end);
   printf("Time is %lf\n",time);
   printf("\n");

   // ending
   ssp_free(instance);
   return 0;
};

