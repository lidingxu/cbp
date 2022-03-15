/**@file   conic_knapsack.h
 * @brief  Conic Knapsack solver
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <vector>
#include <utility>
#include "conic_knapsack_utility.h"
using namespace std;


/** solve conic knapsack problem before timie_limit,
*   find a solution with value greater than the tagret lower bound, and store the solution
*   possible return values for:
*  - Sol_Ge_Bd : a solution with value greater than the tagret lower bound is found
*  - No_Sol_Ge_Bd :  no solution with value greater than the tagret lower bound
*  - Unsolved: the problem greater than the tagret lower bound is unsolved.
*/
void solve_conic_knap(
	SCIP * scip, /** scip pointer */
	stat & stat_prc, /** statistic of pricing */
	const vector<SCIP_Real> & objs, /** objective coefficients of items */
	const vector<SCIP_Real> & mus, /** mu of items */
	const vector<SCIP_Real> & bs, /** b of items */
	SCIP_Real Dalpha, /** Dalpha */
	SCIP_Real capacity, /** capacity */
	int numitems, /** number of items */
	estimator & init_estimator, /** initial estimator*/
	Pr_param & pr_param, 
	const vector<pair<int,int>>& items_same, /** items in same constraints */
	const vector<pair<int,int>>& items_diff, /** items in same constraints */
	const conflict_graph & conflict, /** the conflict graph*/
	vector<SCIP_Real> & stable_center, /** the stablization center */
	const conf & algo_conf, /** algorithm configuration*/
    list<list<int>> & sol_pool, /* solutions pools, solution items are assumed to be sorted */
	SCIP_Real & sol_val, /** solution value */
	SCIP_Real & sol_ub, /** solution value */
	SOLTYPE_CKNAP & sol_type, /** solution type */
	SCIP_Real stop_pricing_obj, // if the optimal  value of pricing solution is <= stop_pricing_obj, early stops
    SCIP_Real time_limit = 3600,  /** solving time left */
    SCIP_Real target_lb = 1.0 /* the tagret lower bound bound  */
);
