//
//  differentialEvolution.c
//  Inversr
//
//  Created by Rafael Verduzco on 8/17/11.
//  Copyright 2011 Universidad Autónoma del Estado de Morelos. All rights reserved.
//

#if __STDC_VERSION__ < 199901L
#error "This code requires a C99 compiler. Try cc -std=c99"
#else

#include "differentialEvolution.h"

void 
differentialEvolution(const size_t maxGenerations,
                      const size_t NP,
                      const size_t D,
                      const double CR,
                      const double F,
                      size_t seed,
                      double (* evaluate) (double [], void *),
                      double * minBounds,
                      double * maxBounds,
                      double * bestVector,
                      double * bestCost,
                      void * userData)
{
    
    double currentGen[NP][D];
    //double newGen[NP][D];
    double u[NP][D];
    double uCosts[NP];
    double currentCosts[NP];
    
    // Initialize random numbers
    gsl_rng * myRNG;
    myRNG = gsl_rng_alloc(gsl_rng_ranlxs2);
    gsl_rng_set(myRNG, seed);
    
    // Initialize the population with random values
    // chosen from within the given boudaries
    // Uniform distribution is a safe first guess
    // one could use Gaussian distribution, etc
    for (size_t agent_i = 0; agent_i < NP; agent_i++) {
        for (size_t param_j = 0; param_j < D; param_j++) {
            currentGen[agent_i][param_j] = 
            gsl_rng_uniform(myRNG) * (maxBounds[param_j] - minBounds[param_j])
            + minBounds[param_j];
        }
        currentCosts[agent_i] = 5000.0;
    }

    size_t generationNumber = 0;
    // Main loop
    while (generationNumber < maxGenerations) {
        for (size_t agent_i = 0; agent_i < NP; agent_i++) {
            // First, chose three different agents
            size_t r1, r2, r3;
            do {r1 = gsl_rng_uniform_int(myRNG, NP);}while (r1 == agent_i);
            do {r2 = gsl_rng_uniform_int(myRNG, NP);}while (r2 == agent_i || r2 == r1);
            do {r3 = gsl_rng_uniform_int(myRNG, NP);}while (r3 == agent_i || r3 == r2 || r3 == r1);
            size_t k = gsl_rng_uniform_int(myRNG, D);
            // Generate the trial population
            for (size_t param_j = 0; param_j < D; param_j++) {
				if (gsl_rng_uniform(myRNG) <= CR || param_j == k) {
                    u[agent_i][param_j] = currentGen[r3][param_j] + 
                    F * (currentGen[r1][param_j] - currentGen[r2][param_j]);
				}
				else u[agent_i][param_j] = currentGen[agent_i][param_j];
			}
            // Handle boundary constraints according to Lampinen J. 2002 IEEE
            // Essentially, replace an out-of-boundaries parameter value
            // with a random value lying within boundaries
            for (size_t param_j = 0; param_j < D; param_j++) {
                if (u[agent_i][param_j] < minBounds[param_j] || 
                    u[agent_i][param_j] > maxBounds[param_j]) {
                    u[agent_i][param_j] = gsl_rng_uniform(myRNG) * (maxBounds[param_j] - minBounds[param_j])
                    + minBounds[param_j];
                }
            }
            
            // Evaluate
            for (size_t cost_i = 0; cost_i < NP; cost_i++) {
                uCosts[cost_i] = evaluate(u[agent_i], userData); 
            }
            // Decide whether to replace the population or not
            if (uCosts[agent_i] < currentCosts[agent_i]) {
                // Substitute the currentGen's agent_i with u[aent_i] as well as its cost
                for (size_t param_j = 0; param_j < D; param_j++) {
                    currentGen[agent_i][param_j] = u[agent_i][param_j];
                }
                currentCosts[agent_i] = uCosts[agent_i];
            } // Otherwise leave it like that
        }//for
        generationNumber++;
    }//while
    
    // Get the smallest cost
    size_t bestCostIndex;
    gsl_sort_smallest_index(&bestCostIndex, 1, currentCosts, 1, NP); 
    (*bestCost) = currentCosts[bestCostIndex];
    // Copy the parameters corresponding to that cost into the output vector
    for (size_t param_j = 0; param_j < D; param_j++) {
        bestVector[param_j] = currentGen[bestCostIndex][param_j];
    }
    gsl_rng_free(myRNG);
}

#endif
