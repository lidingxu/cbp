/**@file   conic_knapsack_pwlbc.h
 * @brief  Conic Knapsack pwl bc alogirthm
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <vector>
#include <list>
#include <tuple>
#include <utility>
#include "utilities.h"
#include "estimator.h"
using namespace std;

/** solve conic knapsack problem before timie_limit,
*   either find a solution with the value greater than the tagret lower bound,  or conclude that there is no such solution
*/
void solve_conic_knap_pwlbc(
	SCIP * scip,
	stat & stat_prc, /** statistic of pricing */
	const vector<SCIP_Real> & objs, /** objective coefficients of items */
	const vector<SCIP_Real> & mus, /** mu of items */
	const vector<SCIP_Real> & bs, /** b of items */
	SCIP_Real Dalpha, /** Dalpha */
	SCIP_Real capacity, /** capacity */
	int numitems, /** number of items */
    Estimator & init_estimator, /** initial estimator*/
	const conf & algo_conf, /** algorithm configuration*/
	SCIP_Real heur_mu, /* heur mu*/
	const vector<pair<int,int>>& items_diff, /** items in same constraints */
	const conflict_graph & conflict, /** the conflict graph*/
    list<list<int>> & sol_pool, /* solutions pools, solution items are assumed to be sorted */
	SCIP_Real & sol_val, /** solution value */
	SCIP_Real & sol_ub, /** solution value */
	SOLTYPE_CKNAP & sol_type, /** solution type */
    SCIP_Real time_limit = 3600,  /** solving time left */
    SCIP_Real target_lb = 1.0 /* the tagret lower bound bound  */
);

/**
 * @brief tightenning the bound of milp piece-wise linear relaxation
 */
void rel_milp_bd(
	SCIP * scip,
	const vector<SCIP_Real> & mus, /** mu of items */
	const vector<SCIP_Real> & bs, /** b of items */
	SCIP_Real Dalpha, /** Dalpha */
	SCIP_Real capacity, /** capacity */
	int numitems, /** number of items */
	const vector<pair<int,int>>& items_diff, /** items in same constraints */
	SCIP_Real & lb, /*valid lower bound*/
	SCIP_Real & ub, /*valid upper bound*/
	int max_iter, /** the maximum iteration of bound tightenning*/
	SCIP_Real rel_eps_bd, /** relative tolerance of bound tightenning */
	SCIP_Real time_limit /** the time limit of bound tightenning*/
);