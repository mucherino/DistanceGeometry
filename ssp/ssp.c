
/*
 * Subset-Sum Problem (SSP)
 *
 * SSP (ssp_t type) main functions
 *
 * last update: September 22nd, 2024
 *
 * AM
 */

#include "ssp.h"

// global variable (seed for random numbers)
unsigned int ssp_seed = 0;

// SSP seed definition (should not be 0)
void ssp_set_seed(unsigned int seed)
{
   assert(seed > 0U);
   srand(seed);
   ssp_seed = seed;
};

// SSP empty initialization
ssp_t* ssp_init(size_t n)
{
   assert(n > 2);
   ssp_t *instance = NULL;
   instance = (ssp_t*)malloc(sizeof(ssp_t));
   if (instance == NULL)  return NULL;
   instance->n = n;
   instance->solution = NULL;
   instance->solution = (bool*)calloc(n,sizeof(bool));
   if (instance->solution == NULL)
   {
      free(instance);
      return NULL;
   };
   for (size_t i = 0; i < n; i++)  instance->solution[i] = false;
   instance->target = 0UL;
   instance->set = NULL;
   instance->set = (unsigned long*)calloc(n,sizeof(unsigned long));
   if (instance->set == NULL)
   {
      free(instance->solution);
      free(instance);
      return NULL;
   };
   for (size_t i = 0; i < n; i++)  instance->set[i] = 0UL;
   return instance;
};

// SSP instance loader from text file
ssp_t* ssp_load(char *filename)
{
   // operning text file
   FILE *input = NULL;
   input = fopen(filename,"r");
   if (input == NULL)
   {
      fprintf(stderr,"Impossible to read input text file\n");
      return NULL;
   };

   // reading the size and ssp init
   size_t n = 0;
   if (fscanf(input,"%lu",&n) == EOF)
   {
      fprintf(stderr,"Error while reading instance size from input file\n");
      return NULL;
   };
   ssp_t *instance = ssp_init(n);
   if (instance == NULL)
   {
      fprintf(stderr,"Error while allocating memory for SSP instance\n");
      return NULL;
   };

   // reading the target
   if (fscanf(input,"%lu",&instance->target) == EOF)
   {
      fprintf(stderr,"Error while reading target from input file\n");
      return NULL;
   };

   // reading the integers forming the SSP instance
   unsigned long element = 0UL;
   for (size_t i = 0; i < n; i++)
   {
      if (fscanf(input,"%lu",&element) == EOF)
      {
         fprintf(stderr,"Error while reading integers forming the SSP instance\n");
         return NULL;
      };
      instance->set[i] = element;
   };
 
   // we are done
   fclose(input);
   return instance;
};

// SSP "free" instance generator (simply random integers)
ssp_t* ssp_gen_free(size_t n,unsigned long max)
{
   assert(n > 2);

   // is seed available?
   if (ssp_seed == 0)  srand(time(NULL));  // time-based seed

   // initializing the instance
   ssp_t *instance = ssp_init(n);
   if (instance == NULL)  return NULL;

   // random choice on the number of integers forming a solution
   double p = 0.2 + 0.6*drand();

   // creating the instance and one of its solutions
   for (size_t i = 0; i < n; i++)
   {
      instance->set[i] = 1UL + rand()%(max - 1);
      instance->solution[i] = false;
      if (drand() < p)
      {
         instance->target = instance->target + instance->set[i];
         instance->solution[i] = true;
      }
   };

   return instance;
};

// SSP instance generator for hard instances 
// (controlled density, target between 2^n-1 (n/2 - sqrt(n)) and 2^n-1 (n/2 + sqrt(n)),
//  where 2^n is the max, and 50% of selected integers in the given solution)
ssp_t* ssp_gen_hard(size_t n,double density)
{
   assert(n > 2);
   assert(density > 0.0);

   // is seed available?
   if (ssp_seed == 0)  srand(time(NULL));  // time-based seed

   // initializing the instance
   ssp_t *instance = NULL;
   instance = ssp_init(n);
   if (instance == NULL)  return NULL;

   // computing the max value to satisfy the given density
   unsigned long ln = n;
   if (density != 1.0)
   {
      ln = (unsigned long) n/density;
      if (ln == 0UL)  ln++;
   }
   unsigned long max = my_exp(ln);
 
   // computing the bounds for the target
   unsigned long sq = sqrt(n);
   unsigned long lb = (max >> 1) * (n/2 - sq);
   unsigned long ub = (max >> 1) * (n/2 + sq);

   // randomly selecting the target
   instance->target = lb + rand()%(ub - lb);
   unsigned long current = instance->target;

   // the max is currently placed in the first position
   instance->set[0] = max;
   instance->solution[0] = false;
   size_t count = 0;
   size_t i = 1;

   /* creating the instance and one of its solutions */

   // first attempt with purely random generation
   for (; i < n; i++)
   {
      instance->set[i] = 1UL + rand()%(max - 1UL);
      instance->solution[i] = false;
      if (current > 0UL)  if (rand()%2)
      {
         if (instance->set[i] <= current)
         {
            current = current - instance->set[i];
         }
         else
         {
            instance->set[i] = current;
            current = 0UL;
         };
         instance->solution[i] = true;
         count++;
      };
   };

   // is the current solution equal to the target?
   if (current > 0UL)
   {
      // adding parts of 'current' to integers belonging to the solution
      while (current > 0UL)
      {
         size_t k = 0UL;
         while (!instance->solution[k])  k = 1UL + rand()%(n - 1UL);
         unsigned long adjust = current;
         if (instance->set[k] + adjust > max)
         {
            unsigned long remainder = max - instance->set[k];
            if (remainder > 1UL)
            {
               adjust = rand()%remainder;
            }
            else  // low probability case
            {     // where we could get stuck
               adjust = 0UL;
               instance->target--; 
               current--;
            };
         };
         instance->set[k] = instance->set[k] + adjust;
         current = current - adjust;
      };
   };

   // is the number of integers in the solution close to 50%?
   double tp = (n/2)*(1.0/n);  // target percentage
   double p = (double)(count) / n;
   if (p < tp)
   {   
      while (p < tp)  // case "too few"
      {
         // look for a selected integer, and a non-selected one
         size_t a = n - 1;
         while (!instance->solution[a])  a = 1UL + rand()%(n - 1UL);
         size_t b = a;
         while (a == b || instance->solution[b])  b = 1UL + rand()%(n - 1UL);

         // share the value of a with b
         instance->set[b] = 1UL + rand()%(instance->set[a] - 1UL);
         instance->solution[b] = true;
         instance->set[a] = instance->set[a] - instance->set[b];
         count++;

         // updating percentage
         p = (double)(count) / n;
      };
   }
   else if (p > 1.0 - tp)
   {
      while (p > 1.0 - tp)  // case "too many"
      {
         // look for two selected integers
         size_t a = n - 1;
         while (!instance->solution[a])  a = 1UL + rand()%(n - 1UL);

         // transferring parts of set[a] to other integers
         while (instance->set[a] != 0UL)
         {
            size_t b = a;
            while (a == b || !instance->solution[b] || instance->set[b] == max)  b = 1UL + rand()%(n - 1UL);

            // transfer part of set[a] to set[b]
            unsigned long part = instance->set[a];
            while (part + instance->set[b] > max)  part = 1UL + rand()%(instance->set[a] - 1UL);
            instance->set[a] = instance->set[a] - part;
            instance->set[b] = instance->set[b] + part;
            if (instance->set[a] == 0UL)
            {
               instance->solution[a] = false;
               count--;
            };
         };

         // updating percentage
         p = (double)(count) / n;
      };
   };

   // the max is placed at a random position
   size_t imax = n - 1;
   while (instance->solution[imax])  imax = 1UL + rand()%(n - 1UL);
   instance->set[0] = instance->set[imax];
   instance->set[imax] = max;

   // the hard instance is finally ready
   return instance;
};

// SSP clone
ssp_t* ssp_clone(ssp_t *instance)
{
   size_t n = instance->n;
   ssp_t *clone = (ssp_t*)malloc(sizeof(ssp_t));
   clone->n = n;
   clone->solution = (bool*)calloc(n,sizeof(bool));
   for (size_t i = 0; i < n; i++)  clone->solution[i] = instance->solution[i];
   clone->target = instance->target;
   clone->set = (unsigned long*)calloc(n,sizeof(unsigned long));
   for (size_t i = 0; i < n; i++)  clone->set[i] = instance->solution[i];
   return clone;
};

// SSP solution loader from char string of 0 and 1's
bool ssp_load_solution(ssp_t *instance,char *solution)
{
   for (size_t i = 0; i < instance->n; i++)
   {
      if (solution[i] != '0' || solution[i] != '1')  return false;
      instance->solution[i] = false;
      if (solution[i] == '1')  instance->solution[i] = true;
   };
   return true;
};

// SSP solution loader from text file
bool ssp_load_solution_from_file(ssp_t *instance,char *filename)
{
   // operning text file
   FILE *input = NULL;
   input = fopen(filename,"r");
   if (input == NULL)  return false;

   // loading the solution
   unsigned short sol = false;
   for (size_t i = 0; i < instance->n; i++)
   {
      if (fscanf(input,"%hu",&sol) == EOF)
      {
         fclose(input);
         return false;
      };
      instance->solution[i] = sol;
   };

   // ending
   fclose(input);
   return true;
};

// SSP max integer
unsigned long ssp_max_integer(ssp_t *instance)
{
   unsigned long max = 0UL;
   for (size_t i = 0; i < instance->n; i++)
   {
      if (max < instance->set[i])  max = instance->set[i];
   };
   return max;
};

// SSP total sum
unsigned long ssp_total_sum(ssp_t *instance)
{
   unsigned long sum = 0UL;
   for (size_t i = 0; i < instance->n; i++)  sum = sum + instance->set[i];
   return sum;
};

// SSP current subset sum up to a given element
unsigned long ssp_sum_upto(ssp_t *instance,size_t upto)
{
   assert(upto <= instance->n);
   unsigned long sum = 0UL;
   for (size_t i = 0; i < upto; i++)  sum = sum + instance->set[i];
   return sum;
}; 

// SSP current subset sum
unsigned long ssp_current_sum(ssp_t *instance)
{
   unsigned long sum = 0UL;
   for (size_t i = 0; i < instance->n; i++)  if (instance->solution[i])  sum = sum + instance->set[i];
   return sum;
};

// SSP count selected integers in current solution
size_t ssp_count_selected(ssp_t *instance)
{
   size_t count = 0UL;
   for (size_t i = 0; i < instance->n; i++)  if (instance->solution[i])  count++;
   return count;
};

// SSP count number of integers larger than target
size_t ssp_count_larger_target(ssp_t *instance)
{
   size_t count = 0UL;
   for (size_t i = 0; i < instance->n; i++)  if (instance->set[i] > instance->target)  count++;
   return count;
};

// SSP density
double ssp_density(ssp_t *instance)
{
   unsigned long max = ssp_max_integer(instance);
   unsigned long nbits = my_ln(max);
   return (double) instance->n / (double) nbits;
};

// SSP error (difference between current sum and target)
unsigned long ssp_error(ssp_t *instance)
{
   unsigned long diff = 0UL;
   unsigned long sum = ssp_current_sum(instance);
   if (sum > instance->target)
      diff = sum - instance->target;
   else
      diff = instance->target - sum;
   return diff;
};

// SSP feasibility
bool ssp_isfeasible(ssp_t *instance)
{
   return ssp_total_sum(instance) >= instance->target;
};

// SSP solution (checks whether the internal solution is OK or not)
bool ssp_issolution(ssp_t *instance)
{
   return ssp_current_sum(instance) == instance->target;
};

// SSP solution eraser
void ssp_erasesol(ssp_t *instance)
{
   for (size_t i = 0; i < instance->n; i++)  instance->solution[i] = false;
};

// SSP printer
void ssp_print(ssp_t *instance)
{
   printf("Subset-Sum instance (n = %lu, target = %lu)\n",instance->n,instance->target);
};

// SSP printer for Matlab (and other similar software tools)
void ssp_print_matlab(ssp_t *instance)
{
   printf("n = %lu;\n",instance->n);
   printf("target = %lu;\n",instance->target);
   printf("set = [%lu",instance->set[0]);
   for (size_t i = 1; i < instance->n; i++)  printf(" %lu",instance->set[i]);
   printf("];\n");
};

// SSP printer for current solution in binary format
void ssp_print_solution(ssp_t *instance)
{
   for (size_t i = 0; i < instance->n; i++)  printf("%d",instance->solution[i]);
   printf("\n");
};

// freeing SSP internal memory
void ssp_free(ssp_t *instance)
{
   free(instance->solution);
   free(instance->set);
   free(instance);
};

