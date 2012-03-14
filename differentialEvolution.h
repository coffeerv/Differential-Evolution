//
//  differentialEvolution.h
//  Inversr
//  
//  Created by Rafael Verduzco on 8/5/11.
//  Copyright 2011 Universidad Autónoma del Estado de Morelos. All rights reserved.
//

#ifndef Inversr_differentialEvolution_h
#define Inversr_differentialEvolution_h

#if __STDC_VERSION__ < 199901L
#error "This code requires a C99 compiler. Try cc -std=c99"
#else

#include <stdio.h>
#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>
//#include <gsl/gsl_matrix.h>

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
                      void * userData);

#endif // C language version guards

#endif // Include guards
