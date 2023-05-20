/**@file   conic_knapsack_heur.h
 * @brief  Conic Knapsack heuristics
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#pragma once
#include "utilities.h"
using namespace std;

/** find a feasible solution for conic knapsack problem by heuristics
*/
void solve_conic_knap_heur(
	const vector<SCIP_Real> & objs, /** objective coefficients of items */
	const vector<SCIP_Real> & mus, /** mu of items */
	const vector<SCIP_Real> & bs, /** b of items */
	SCIP_Real Dalpha, /** Dalpha */
	SCIP_Real capacity, /** capacity */
	int numitems, /** number of items */
	const vector<pair<int,int>>& items_diff, /** items in same constraints */
	const conflict_graph & conflict, /** the conflict graph*/
	vector<int> & items_bin, /* items in the maximal bin, assumed to be sorted */
	SCIP_Real & sol_val, /** solution value */
	SCIP_Real time_limit ,  /** solving time left */
	SCIP_Real target_lb = 1.0 /* the tagret lower bound bound  */
);

