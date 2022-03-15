/**@file   conic_knapsack_heur.cpp
 * @brief  Conic Knapsack heuristics
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <vector>
#include <list>
#include <limits>
#include <ilcplex/ilocplex.h>
#include "conic_knapsack_heur.h"


/** find a feasible solution for conic knapsack problem by the best fit heuristic,
 *   each time, the method forces one of the items to be in the solution and storing the best solution found.
*/
void solve_conic_knap_best_fit(
	const vector<SCIP_Real> & objs, /** objective coefficients of items */
	const vector<SCIP_Real> & mus, /** mu of items */
	const vector<SCIP_Real> & bs, /** b of items */
	SCIP_Real Dalpha, /** Dalpha */
	SCIP_Real capacity, /** capacity */
	int numitems, /** number of items */
	const conflict_graph & conflict, /** the conflict graph*/
    vector<int> & items_bin, /* items in the maximal bin, assumed to be sorted */
	SCIP_Real & sol_val /** solution value */
){
	vector<short> status(numitems); // -1 unknown, 0 inpacked, 1 packed
	SCIP_Real best_sol_val  = -1;
	list<int> best_sol;
	vector<bool>never_try(numitems, false);
	int num_never_try = 0;
	for(int item1 = 0; item1 < numitems; item1++){
			SCIP_Real mu1 = mus[item1];
			SCIP_Real b1 = bs[item1];
			SCIP_Real obj1 = objs[item1];
			for(int item2 = item1 + 1; item2 < numitems; item2++){
				if(mu1 >= mus[item2] && b1 >= bs[item2] && obj1<= objs[item2]){
					never_try[item1] = true;
					num_never_try++;
					break;
				}
		}	
	}
	for(int fix = 0; fix < numitems; fix++){
		fill(status.begin(), status.end(), -1);
		SCIP_Real sumb = bs[fix];
		SCIP_Real summu = mus[fix];
		SCIP_Real obj = objs[fix];
		if(objs[fix] < tol || never_try[fix]){
			continue;
		}
		status[fix] = 1; // add to the solution
		SCIP_Real lhs = summu +  Dalpha * sqrt(sumb); 
		if(lhs > capacity){ // if not packable cont.
			continue;
		}
		for(int item: conflict.get_diffs(fix)){ // fixed by conflict
			status[item] = 0;
		}
		while(true){
			SCIP_Real best_ratio  = -1;
			int candidate = -1;
			for(int item = 0; item < numitems; item++){
				if( status[item] != -1){
					continue;
				}
				SCIP_Real lhs_ = summu + mus[item] + Dalpha * sqrt(sumb + bs[item]);
				if(lhs_ > capacity){
					status[item] = 0;
					continue;
				}
				SCIP_Real ratio =  objs[item] / (lhs_ - lhs);
				if(ratio > best_ratio){
					best_ratio = ratio;
					candidate = item;
				}
			}
			if(candidate != -1){
				summu += mus[candidate];
				sumb += bs[candidate];
				obj += objs[candidate];
				lhs = summu + Dalpha * sqrt(sumb);
				for(int item_: conflict.get_diffs(candidate)){ // fixed by conflict
					status[item_] = 0;
				}
				status[candidate] = 1;
			}
			else{
				break;
			}
		}
		// an improving solution found, update
		if(obj > best_sol_val){
			best_sol_val = obj;
			best_sol.clear();
			for(int item = 0; item < numitems; item++){
				if(status[item] == 1){
					best_sol.push_back(item);
				}
			}
		}
	}
	sol_val = best_sol_val;
	items_bin = vector<int>(best_sol.begin(), best_sol.end());
};


/** find a feasible solution for conic knapsack problem by heuristics, 
 * Notice: the heuristics only applies for the merged items.
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
	SCIP_Real target_lb /* the tagret lower bound bound  */
){
   vector<int> items_bin_best_fit; /* items in the maximal bin, assumed to be sorted */
   SCIP_Real  sol_val_best_fit = 0; /** solution value */
    //SCIPdebugMessage("enter heur,%lf\n", time_limit);
    solve_conic_knap_best_fit(objs, mus, bs, Dalpha, capacity,
    numitems, conflict, items_bin_best_fit, sol_val_best_fit);
	//SCIPdebugMessage("quit heur,%lf\n", sol_val_best_fit);
	if(sol_val_best_fit > target_lb +tol ){
		sol_val = sol_val_best_fit;
		items_bin = items_bin_best_fit;
		//printf("best_f\n");
	}
	
};