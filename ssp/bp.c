
/*
 * Subset-Sum Problem (SSP)
 *
 * Branch-and-Prune style tree exploration
 *
 * Pruning is performed by verifying whether the current sum is larger than the target,
 * and whether the current sum + the remaining integers can actually reach the target or not.
 *
 * last update: September 22nd, 2024
 *
 * AM
 */

#include "ssp.h"

// how many solutions?
bool ONLYONE = false;

// counter of found solutions
size_t nsols = 0UL;

// sum of remaining integers
unsigned long rest = 0UL;

// recursive calls to branch-and-prune
void bp_rec(ssp_t* instance,size_t current,unsigned long psum)
{
   unsigned long sum = psum;
   rest = rest - instance->set[current];

   // suppose the current integer is NOT part of the subset
   instance->solution[current] = false;
   if (sum + rest >= instance->target)
   {
      if (current + 1 < instance->n)
      {
         bp_rec(instance,current + 1,sum);
      }
      else if (sum == instance->target)
      {
         nsols++;
         printf("#%lu: ",nsols);
         unsigned long tmp = 0UL;
         for (size_t i = 0; i < instance->n; i++)
         {
            if (instance->solution[i])
            {
               printf("%lu ",instance->set[i]);
               tmp = tmp + instance->set[i];
            };
         };
         printf("[%lu]\n",tmp);
      };
   };
   if (ONLYONE)  if (nsols > 0UL)  goto BACK;

   // suppose the current integer is part of the subset
   instance->solution[current] = true;
   sum = sum + instance->set[current];
   if (sum <= instance->target)
   {
      if (current + 1 < instance->n)
      {
         bp_rec(instance,current + 1,sum);
      }
      else if (sum == instance->target)
      {
         nsols++;
         printf("#%lu: ",nsols);
         unsigned long tmp = 0UL;
         for (size_t i = 0; i < instance->n; i++)
         {  
            if (instance->solution[i])
            {  
               printf("%lu ",instance->set[i]);
               tmp = tmp + instance->set[i];
            };
         };
         printf("[%lu]\n",tmp);
      };
   };

BACK:
   // updating variables for backtracking
   sum = sum - instance->set[current];
   rest = rest + instance->set[current];
};

// SSP by branch-and-prune
bool bp(ssp_t* instance)
{
   // erasing any previous solutions
   ssp_erasesol(instance);
   nsols = 0UL;

   // computing the total sum
   rest = ssp_total_sum(instance);

   // two possible choices for the very first element
   unsigned long sum = 0UL;
   rest = rest - instance->set[0];
   instance->solution[0] = false;
   if (rest >= instance->target)  bp_rec(instance,1,sum);
   if (ONLYONE)  if (nsols > 0UL)  return true;
   instance->solution[0] = true;
   sum = sum + instance->set[0];
   if (sum < instance->target)  bp_rec(instance,1,sum);

   // ending
   return nsols > 0UL;
};

