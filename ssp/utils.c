
/*
 * Subset-Sum Problem (SSP)
 *
 * Utilities
 *
 * last update: September 22nd, 2024
 *
 * AM
 */

#include "ssp.h"

// drand (double-precision random number in [0,1])
double drand()
{
   return ((double)rand() / RAND_MAX);
};

// computing CPU clock time
float compute_time(time_t start,time_t end)
{
   return ((float)((int)end - (int)start))/CLOCKS_PER_SEC;
};

// ln with bit shifts
unsigned long my_ln(unsigned long number)
{
   unsigned long log = 0UL;
   while (number > 0UL)
   {
      number = number >> 1;
      if (number > 0UL)  log++;
   };
   return log;
};

// exp with bit shifts
unsigned long my_exp(unsigned long exponent)
{
   unsigned long exp = 1UL;
   if (exponent > 64UL)
   {
      fprintf(stderr,"This wont work with unsigned long precision (exponent is %lu)\n",exponent);
      return 0UL;
   };
   return (exp << exponent);
};

