/**@file   conic_knapsack_utility.h
 * @brief  Conic Knapsack utility
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#pragma once
#include "probdata_cbp.h"
#include "scip/scip_numerics.h"
#include "estimator.h"

#ifndef tol
#define tol 1e-6
#endif

/* type of solution of conic knapsack problem */
enum SolType
{ 
	Unknown = 0,
	Infeasible  = 1,
	Feasible_Heur  = 2, 
	Feasible_Exact  = 3,  
	Optimal = 4, 
};
//typedef enum SolType SOLTYPE_CKNAP;
typedef SolType SOLTYPE_CKNAP;
