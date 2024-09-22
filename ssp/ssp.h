
/*
 * Subset-Sum Problem (SSP)
 *
 * Header file
 *
 * last update: September 22nd, 2024
 *
 * AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <time.h>

// SSP data structure
// - "n" is the size of the instance (# of elements in the set, size_t type)
// - "solution" is an array of unsigned long that may contain one of the known solutions (bool type)
// - "target" is the (unsigned long) integer target
// - "set" is the pointer to the array of (unsigned long) integer values
struct ssp
{
   size_t n;
   unsigned long target;
   unsigned long *set;
   bool *solution;
};
typedef struct ssp ssp_t;

// ssp_t
void ssp_set_seed(unsigned int seed);
ssp_t* ssp_init(size_t n);
ssp_t* ssp_load(char *filename);
ssp_t* ssp_gen_free(size_t n,unsigned long max);
ssp_t* ssp_gen_hard(size_t n,double density);
ssp_t* ssp_clone(ssp_t *instance);
bool ssp_load_solution(ssp_t *instance,char *solution);
bool ssp_load_solution_from_file(ssp_t *instance,char *filename);
unsigned long ssp_max_integer(ssp_t *instance);
unsigned long ssp_total_sum(ssp_t *instance);
unsigned long ssp_sum_upto(ssp_t *instance,size_t upto);
unsigned long ssp_current_sum(ssp_t *instance);
size_t ssp_count_selected(ssp_t *instance);
size_t ssp_count_larger_target(ssp_t *instance);
double ssp_density(ssp_t *instance);
unsigned long ssp_error(ssp_t *instance);
bool ssp_isfeasible(ssp_t *instance);
bool ssp_issolution(ssp_t *instance);
void ssp_erasesol(ssp_t *instance);
void ssp_print(ssp_t *instance);
void ssp_print_matlab(ssp_t *instance);
void ssp_print_solution(ssp_t *instance);
void ssp_free(ssp_t *instance);

// Branch-and-Prune
void bp_rec(ssp_t* instance,size_t current,unsigned long psum);
bool bp(ssp_t* instance);

// utils
double drand();
float compute_time(time_t start,time_t end);
unsigned long my_ln(unsigned long number);
unsigned long my_exp(unsigned long exponent);

