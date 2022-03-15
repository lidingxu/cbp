/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pricer_cknap.cpp
 * @brief  Conic Knapsack pricer
 * @author Liding Xu
 * 
 * This file implements the variable pricer which check if variables exist with negative reduced cost. See
 * @ref BINPACKING_PRICER for more details.
 *
 * @page BINPACKING_PRICER Pricing new variables
 *
 * The task of the pricer is to search for new variables with negative reduced costs. For this, the following integer
 * program is solved:
 *
 *  \f[
 *  \begin{array}[t]{rll}
 *       \max & \displaystyle \sum_{i=1}^n (\lambda_S)_i y^\star_i\\
 *        & \\
 *        subject \ to & \displaystyle \sum_{i=0}^n (\lambda_S)_i s_i \leq \kappa \\
 *        & \\
 *        & (\lambda_S)_i \in \{0,1\} & \quad \forall i \in \{ 1, \dots , n \} \\
 *  \end{array}
 * \f]
 *
 * where \f$ (\lambda_S)_i \f$ for \f$i\in\{1,\dots,n\}\f$ are binary variables and \f$y^\star_i\f$ given by the dual
 * solution of the restricted master problem. See the \ref BINPACKING_PROBLEM "problem description" for more details.
 *
 * To solve the above integer program, we create a new SCIP instance within SCIP and use the usual functions to create
 * variables and constraints. Besides, we need the current dual solutions to all set covering constraints (each stands
 * for one item) which are the objective coefficients of the binary variables. Therefore, we use the function
 * SCIPgetDualsolSetppc() which returns the dual solutions for the given set covering constraint.
 *
 * Since we also want to generate new variables during search, we have to care that we do not generate variables over
 * and over again. For example, if we branched or fixed a certain packing to zero, we have to make sure that we do not
 * generate the corresponding variables at that node again. For this, we have to add constraints forbidding to generate
 * variables which are locally fixed to zero. See the function addFixedVarsConss() for more details. While using the
 * \ref BINPACKING_BRANCHING "Ryan/Foster branching", we also have to ensure that these branching decisions are respected. This is
 * realized within the function addBranchingDecisionConss().
 *
 * @note In case of this binpacking example, the master LP should not get infeasible after branching, because of the way
 *       branching is performed. Therefore, the Farkas pricing is not implemented.
 *       1. In case of Ryan/Foster branching, the two items are selected in a way such that the sum of the LP values of
 *          all columns/packings containing both items is fractional. Hence, it exists at least one column/packing which
 *          contains both items and also at least one column/packing for each item containing this but not the other
 *          item. That means, branching in the "same" direction stays LP feasible since there exists at least one
 *          column/packing with both items and branching in the "differ" direction stays LP feasible since there exists
 *          at least one column/packing containing one item, but not the other.
 *       2. In case of variable branching, we only branch on fractional variables. If a variable is fixed to one, there
 *          is no issue.  If a variable is fixed to zero, then we know that for each item which is part of that
 *          column/packing, there exists at least one other column/packing containing this particular item due to the
 *          covering constraints.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <iostream>
#include <utility>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <set>


#include "objscip/objscip.h"
#include "cons_samediff.h"
#include "pricer_cknap.h"
#include "probdata_cbp.h"
#include "scip/cons_setppc.h"
#include "conic_knapsack.h"
#include "conic_knapsack_approx.h"

class UnionFindSet{
public:
    unordered_map<int, int> father;   // value表示父节点
    unordered_map<int, int> size;
    explicit UnionFindSet(const vector<int> & data){
        father.clear();
        size.clear();
        for(auto var : data){
            father.insert({var, var});
            size.insert({var, 1});
        }
    }

    int findRep(int node){
        int f = father[node];
        if (f != node){
            f = findRep(f);
            father[node] = f;
        }
        return f;
    }

    bool isSameSet(int a, int b){
        return findRep(a) == findRep(b);
    }

    void Union(int a, int b){
        int aRep = findRep(a);
        int bRep = findRep(b);
        if(aRep != bRep){
            int aSize = size[aRep];
            int bSize = size[bRep];
            if(aSize > bSize){
                father[bRep] = aRep;
                size[aRep] += size[bRep];
            }
            else{
                father[aRep] = bRep;
                size[bRep] += size[aRep];
            }
        }
    }
};

/**
 * @brief merge the items
 * 
 */
int merge(
   const vector<pair<int,int>> & items_same, // items in the same constraints
   const vector<pair<int,int>> & items_differ, // items in the differ constraints
   const vector<SCIP_Real> & mus, // items' mus
   const vector<SCIP_Real> & bs, // items' bs
   const int numitems, // number of the items
   vector<int> & belongs, // items corresponding to the to new items
   vector<vector<int>> & new_items, // new items that contains the original items 
   vector<SCIP_Real> & new_mus, // new items' mus
   vector<SCIP_Real> & new_bs, // new items' bs
   vector<pair<int, int>> & new_differ // new different constraints 
){
   assert(mus.size() == belongs.size());
   //iota(belongs.begin(), belongs,end(), 0); // fill the belongs to itself
   for(int i = 0; i < belongs.size(); i++){
      belongs[i] = i;
   }
   UnionFindSet findset = UnionFindSet(belongs);

   for(auto p: items_same){ // assign the belongs to the lexicograph  less item
      if(!findset.isSameSet(p.first, p.second)){
         findset.Union(p.first, p.second);
      }
   } 

   int num_new_items = 0;
   unordered_map<int, int> indexmap;
   new_items.clear();
   for(int item = 0; item < numitems; item++){
      int belong = findset.findRep(item);
      if(indexmap.find(belong) != indexmap.end()){ // add to new item
         int new_item = indexmap[belong];
         belongs[item] = new_item;
         new_items[new_item].push_back(item);
      }
      else{
         // create and add to new item
         indexmap[belong] = num_new_items;
         belongs[item] = num_new_items;
         new_items.push_back(vector<int>());
         new_items[num_new_items].push_back(item);
         num_new_items++;
      }
   }


   
   new_mus = vector<SCIP_Real>(num_new_items, 0);
   new_bs = vector<SCIP_Real>(num_new_items, 0);
   for(int item = 0; item < numitems; item++){
      int belong = belongs[item];
      new_mus[belong] += mus[item];
      new_bs[belong] += bs[item];
   }
   // new diff constraints;
   set<pair<int, int>> differ_dict;
   new_differ.clear();
   for(auto pair_item: items_differ){
      int new_item1 = belongs[pair_item.first];
      int new_item2 = belongs[pair_item.second];
      assert(new_item1 != new_item2);
      pair<int, int> s = new_item1 < new_item2 ? make_pair(new_item1, new_item2): make_pair(new_item2, new_item1);
      if(differ_dict.find(s) == differ_dict.end()){
         new_differ.push_back(s);
         differ_dict.insert(s);
      }
   }
   return num_new_items;
}


// apply greedy heuristics to find a maximal number of items in one bin
int greedy_heuristic(
   vector<SCIP_Real> & mus, // mus 
   vector<SCIP_Real> & bs, // bs
   SCIP_Real Dalpha, // dalpha
   SCIP_Real capacity, // capacity
    int numitems, // the number of items
    conflict_graph & conflict // conflict
){
   vector<short>  status(numitems, -1);
	SCIP_Real sumb = 0;
	SCIP_Real summu = 0;
	SCIP_Real lhs = 0;
   int bin_size = 0;
   while(true){
		SCIP_Real best_use = capacity + 1;
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
         SCIP_Real cap_use = lhs_ - lhs;
         if(cap_use < best_use){
            best_use = cap_use;
				candidate = item;
         }
      }
      if(candidate != -1){
         summu += mus[candidate];
         sumb += bs[candidate];
         lhs = summu + Dalpha * sqrt(sumb);
         for(int item_: conflict.get_diffs(candidate)){ // fixed by conflict
            status[item_] = 0;
         }
         status[candidate] = 1;
         bin_size++;
      }
      else{
         break;
      }
   }
   return bin_size;
}
// solve the pricing problem
void solve_pricing(
	SCIP * scip, /** scip pointer */
   int numitems, // the number of items
   vector<SCIP_Real> & objs, // objective values
   ProbDataCBP * probdata, // problem data
   SCIP_Real stop_pricing_obj, // if the optimal  value of pricing solution is <= stop_pricing_obj, early stops
   SCIP_Real timelimit, // timie limit 
   list<vector<int>> & sol_pool, /* solutions pools, solution items are assumed to be sorted */
   SCIP_Real & sol_val, // solution value
   SCIP_Real & sol_ub, // solution upper bound
   SOLTYPE_CKNAP & sol_type, // solution_type
   SCIP_Real target_bd // target bound 
){
   
   SCIP_Real timelimit_bd = (numitems) * 0.015;
   SCIP_Real rel_eps_bd  = 5e-5;
   SCIP_Real max_iter_bd = 200;
   long long nodenumber= SCIPnodeGetNumber(SCIPgetFocusNode(scip));
   //SCIPdebugMessage("%d\n", int(nodenumber));
   if(SCIPinDive(scip)){
      vector<pair<int, int>> items_same = getItemsInSame(scip); // items in same constraints
	   vector<pair<int, int>> items_differ = getItemsInDiffer(scip); //items in differ constraints
      vector<int> belongs(numitems); // items corresponding to the to new items
      vector<vector<int>> new_items; // new items that contains the original items 
      vector<SCIP_Real> new_mus; // new items' mus
      vector<SCIP_Real> new_bs; // new items' bs
      vector<pair<int, int>> new_differ; // new different constraints 
      int num_new_items = merge(items_same, items_differ, probdata->mus, probdata->bs, numitems, belongs, new_items, new_mus, new_bs, new_differ);
      vector<SCIP_Real> new_objs(num_new_items, 0); // new items' objs
      for(int item = 0; item < numitems; item++){
         int belong = belongs[item];
         new_objs[belong] += objs[item];
      }
      vector<int> new_items_bin;
      vector<pair<int,int>> new_same(0);
      conflict_graph conflict = conflict_graph(num_new_items, new_differ);

      // construct initial quadratic estimator
      int full_piece_size = greedy_heuristic(new_mus,  new_bs, probdata->Dalpha, probdata->capacity,  num_new_items,  conflict); // compute the piece size
      //normalization(new_mus,   probdata->capacity,  num_new_items);  // normalize constraints
      SCIP_Real lb_milp = 0, ub_milp = probdata->capacity;
      int piece_size = ceil((ub_milp - lb_milp) /  probdata->capacity *  full_piece_size);
      estimator  init_estimator(probdata->capacity , lb_milp, ub_milp, piece_size);
      SCIP_Real cbp_time = numitems * 0.015;

      list<list<int>> new_sol_pool;
      vector<SCIP_Real> stable_center(num_new_items, 0);
      solve_conic_knap(scip, probdata->stat_pr, new_objs, new_mus, new_bs,   probdata->Dalpha, probdata->capacity, num_new_items,  init_estimator,
        new_same,  new_differ, conflict, stable_center, probdata->algo_conf, new_sol_pool,  sol_val,  sol_ub, sol_type,   stop_pricing_obj, 
         cbp_time < timelimit?  cbp_time : timelimit, target_bd);
      for(auto it =  new_sol_pool.begin(); it != new_sol_pool.end(); it++){
         sol_pool.push_back(vector<int>());
         vector<int> & items_bin = sol_pool.back();
         for(int new_item: *it){
            items_bin.insert(items_bin.end(), new_items[new_item].begin(), new_items[new_item].end());
         }
         it->clear();
      }
      new_sol_pool.clear();
   }
   else if(nodenumber != probdata->currentnode){
      vector<pair<int, int>> items_same = getItemsInSame(scip); // items in same constraints
	   vector<pair<int, int>> items_differ = getItemsInDiffer(scip); //items in differ constraints
      probdata->currentnode = nodenumber;
      //SCIPdebugMessage("1.1\n");
      probdata->num_new_items = merge(items_same, items_differ, probdata->mus, probdata->bs, numitems, probdata->belongs, probdata->new_items, probdata->new_mus, 
      probdata->new_bs, probdata->new_differ);
      vector<SCIP_Real> new_objs(probdata->num_new_items, 0); // new items' objs
      for(int item = 0; item < numitems; item++){
         int belong = probdata->belongs[item];

         new_objs[belong]+= objs[item];
      }

      vector<int> new_items_bin;
      vector<pair<int,int>> new_same(0);
      probdata->conflict = conflict_graph(probdata->num_new_items, probdata->new_differ);

      // construct initial quadratic estimator
      probdata->cbp_time = (probdata->num_new_items) * 0.022;
      if(nodenumber == (long long)1 && !probdata->algo_conf.is_misocp && probdata->algo_conf.is_bd_tight){
         int full_piece_size = greedy_heuristic(probdata->new_mus,  probdata->new_bs, probdata->Dalpha, probdata->capacity, 
         probdata->num_new_items,  probdata->conflict); // compute the piece size
         SCIP_Real lb_milp = 0, ub_milp = probdata->capacity;
         rel_milp_bd( scip, probdata->new_mus, probdata->new_bs,  probdata->Dalpha, probdata->capacity , probdata->num_new_items, 
         probdata->new_differ,  lb_milp,  ub_milp,  max_iter_bd,  rel_eps_bd, timelimit_bd);
         SCIPdebugMessage("%lf %lf\n", lb_milp, ub_milp);
         int piece_size = ceil((ub_milp - lb_milp) /  probdata->capacity *  full_piece_size);
         probdata->init_estimator = estimator(probdata->capacity , lb_milp, ub_milp, piece_size);
         probdata->cbp_time += (log(piece_size) + 2) * 0.022;
      }

      list<list<int>> new_sol_pool;
      probdata->stable_center = vector<SCIP_Real> (probdata->num_new_items, 0);
      solve_conic_knap(scip, probdata->stat_pr,  new_objs, probdata->new_mus, probdata->new_bs, probdata->Dalpha, probdata->capacity, probdata->num_new_items, 
       probdata->init_estimator, new_same,  probdata->new_differ, probdata->conflict, probdata->stable_center, probdata->algo_conf, new_sol_pool,  sol_val, sol_ub, sol_type,   stop_pricing_obj, 
       probdata->cbp_time < timelimit?  probdata->cbp_time : timelimit, target_bd);
      for(auto it =  new_sol_pool.begin(); it != new_sol_pool.end(); it++){
         sol_pool.push_back(vector<int>());
         vector<int> & items_bin = sol_pool.back();
         for(int new_item: *it){
            items_bin.insert(items_bin.end(), probdata->new_items[new_item].begin(), probdata->new_items[new_item].end());
         }
         it->clear();
      }
      new_sol_pool.clear();
   }
   else{
      vector<SCIP_Real> new_objs(probdata->num_new_items, 0); // new items' objs
      for(int item = 0; item < numitems; item++){
         int belong = probdata->belongs[item];
         //assert(belong <= probdata->num_new_items);
         new_objs[belong]+= objs[item];
      }
      vector<int> new_items_bin;
      vector<pair<int,int>> new_same(0);
      list<list<int>> new_sol_pool;
      solve_conic_knap(scip, probdata->stat_pr, new_objs, probdata->new_mus, probdata->new_bs,  probdata->Dalpha, probdata->capacity,  probdata->num_new_items, 
       probdata->init_estimator,  new_same,  probdata->new_differ, probdata->conflict, probdata->stable_center, probdata->algo_conf, new_sol_pool,  sol_val, sol_ub, sol_type,   stop_pricing_obj, 
        probdata->cbp_time < timelimit?  probdata->cbp_time : timelimit, target_bd);
      //SCIPdebugMessage("2.3\n");
      for(auto it =  new_sol_pool.begin(); it != new_sol_pool.end(); it++){
         sol_pool.push_back(vector<int>());
         vector<int> & items_bin = sol_pool.back();
         for(int new_item: *it){
            items_bin.insert(items_bin.end(), probdata->new_items[new_item].begin(), probdata->new_items[new_item].end());
         }
         it->clear();
      }
      new_sol_pool.clear();
      //SCIPdebugMessage("2.4\n");
   }
   //SCIPdebugMessage("%d %d %f %d %f\n", int(sol_pool.size()) , probdata->stat_pr.col_exact, probdata->stat_pr.time_exact, probdata->stat_pr.col_heur, probdata->stat_pr.time_heur );
}
/** Pricing of additional variables if LP is feasible.
 *
 *  - get the values of the dual variables you need
 *  - if this tour has negative reduced cost, add it to the LP
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : at least one improving variable was found, or it is ensured that no such variable exists
 *  - SCIP_DIDNOTRUN  : the pricing process was aborted by the pricer, there is no guarantee that the current LP solution is optimal
 */
SCIP_DECL_PRICERREDCOST(PricerConicKnap::scip_redcost)
{  /*lint --e{715}*/

   assert(scip != NULL);

   (*result) = SCIP_DIDNOTRUN;

   ProbDataCBP * probdata = NULL;
	probdata = dynamic_cast<ProbDataCBP *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);


   SCIP_Real timelimit;

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if(timelimit - SCIPgetSolvingTime(scip) < 100 && SCIPnodeGetNumber(SCIPgetFocusNode(scip)) == SCIPnodeGetNumber(SCIPgetRootNode(scip)) ){
      *lowerbound = probdata->global_lb;
      *stopearly = TRUE;
      return SCIP_OKAY;
   }
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   
   
   
   int numitems = probdata->numitems;
   vector<SCIP_Real> objs(numitems); // objective values of items given by dual variables
   for(int item = 0; item < numitems; item++){
      SCIP_Real dual = SCIPgetDualsolSetppc(scip, probdata->sc_conss[item].sc_cons); 
      if(!SCIPisGE(scip, dual, 0)){
         SCIPdebugMessage("%lf\n" , dual);

      }
      objs[item] = dual;
   }


   list<vector<int>> sol_pool;
   SCIP_Real sol_val = 0, sol_ub = MAXFLOAT;
   SOLTYPE_CKNAP sol_type = Unknown;
   //assert(!SCIPinDive(scip));
   SCIP_Real  lp_obj = SCIPgetLPObjval(scip); // LP objective
   SCIP_Real  primal_bd = SCIPgetPrimalbound(scip); // primal bound
   SCIP_Real  lp_ub = SCIPceil(scip,  lp_obj); // LP objective ceil
   SCIP_Real stop_bd = primal_bd < lp_ub ? primal_bd : lp_ub; // stopping bound
   assert(stop_bd-  1 > 0.5);
   SCIP_Real stop_pricing_obj = lp_obj / (stop_bd - 1); 
   SCIP_Real stop_pricing_obj2 = lp_obj / probdata->global_lb;
   stop_pricing_obj = stop_pricing_obj > stop_pricing_obj2 ? stop_pricing_obj : stop_pricing_obj2;

   solve_pricing(scip,  numitems,  objs,  probdata, stop_pricing_obj, timelimit, sol_pool, sol_val, sol_ub, sol_type, 1 );

   
      SCIP_Real Farley_Bd_Ori =  lp_obj / sol_ub; 
      //SCIPdebugMessage("%lf %lf %lf %lf %lf %lf\n", (Farley_Bd_Ori - lp_obj) / lp_obj , SCIPgetLPObjval(scip), sol_val, sol_ub, stop_pricing_obj, fabs(stop_pricing_obj - sol_ub) / sol_ub );

   if(sol_type == Optimal || sol_type == Feasible_Exact ){
      //SCIPdebugMessage("%lf\n", sol_val);
      (*result) = SCIP_SUCCESS;
      if(SCIPisPositive(scip, sol_val - 1)){ // add the optimal pack
         SCIP_Real Farley_Bd_Ori =  lp_obj / sol_ub; 
         SCIP_Real Farley_Bd = SCIPceil(scip, Farley_Bd_Ori); 
         SCIP_Real disp_Bd = Farley_Bd_Ori > probdata->global_lb? Farley_Bd_Ori: probdata->global_lb;
         SCIP_Bool prune_bound = SCIPisGE(scip, Farley_Bd - primal_bd, 0);
         SCIP_Bool prune_improve = SCIPisEQ(scip, Farley_Bd, lp_ub);
         if(prune_bound || prune_improve){ // early termination
            *lowerbound = disp_Bd;
            probdata->global_lb = disp_Bd;
            (*result) = SCIP_DIDNOTRUN;
            (*stopearly) = TRUE;
         }
         else{
            *lowerbound = disp_Bd;
            probdata->global_lb = disp_Bd;
            for(auto it = sol_pool.begin(); it != sol_pool.end(); it++){
               sort(it->begin(), it->end());
               SCIP_CALL(probdata->addPackVar(scip, (*it), TRUE));
            }
         }
      }
   }
   else if(sol_type == Feasible_Heur && SCIPisPositive(scip, sol_val - 1)){
         for(auto it = sol_pool.begin(); it != sol_pool.end(); it++){
            sort(it->begin(), it->end());
            SCIP_CALL(probdata->addPackVar(scip, (*it), TRUE));
         }
      (*result) = SCIP_SUCCESS;
   }
   else if(sol_type == Infeasible){ // stop the pricing iteration
      //*lowerbound = SCIPceil(scip,  SCIPgetLPObjval(scip));
      probdata->global_lb  = lp_obj > probdata->global_lb? lp_obj  : probdata->global_lb;
      (*result) = SCIP_SUCCESS;
      //SCIPdebugMessage("termination by success");
   }
   else{
      assert(!SCIPisPositive(scip, sol_val - 1));
      SCIPdebugMessage("Unkown %lf\n", sol_val);
   }
   return SCIP_OKAY;
}





/** Pricing of additional variables if LP is infeasible.
 *
 *  - get the values of the dual Farks multipliers you need
 *  - construct the reduced-cost arc lengths from these values
 *  - find the shortest admissible path with respect to these lengths
 *  - if this tour has negative reduced cost, add it to the LP
 */
SCIP_DECL_PRICERFARKAS(PricerConicKnap::scip_farkas)
{

   assert(scip != NULL);

   (*result) = SCIP_DIDNOTRUN;

   ProbDataCBP * probdata = NULL;
	probdata = dynamic_cast<ProbDataCBP *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);


   SCIP_Real timelimit;

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   

   int numitems = probdata->numitems;
   vector<SCIP_Real> objs(numitems); // objective values of items given by dual variables
   for(int item = 0; item < numitems; item++){
      SCIP_Real dual = SCIPgetDualfarkasSetppc(scip, probdata->sc_conss[item].sc_cons); 
      assert(SCIPisGE(scip, dual, 0));
      objs[item] = dual;
   }


   vector<int> items_bin;
   SCIP_Real sol_val = 0, sol_ub = MAXFLOAT;
   SOLTYPE_CKNAP sol_type = Unknown;
   //assert(!SCIPinDive(scip));
   list<vector<int>> sol_pool;
   solve_pricing(scip,  numitems,  objs,  probdata,  SCIP_DEFAULT_INFINITY, timelimit, sol_pool, sol_val, sol_ub,  sol_type,  0 );


   if(sol_type == Optimal){
      (*result) = SCIP_SUCCESS;
      if(SCIPisPositive(scip, sol_val)){ // add the optimal pack
            for(auto it = sol_pool.begin(); it != sol_pool.end(); it++){
               sort(it->begin(), it->end());
               SCIP_CALL(probdata->addPackVar(scip, (*it), TRUE));
            }
      }
   }
   else if((sol_type == Feasible_Exact || sol_type == Feasible_Heur) && SCIPisPositive(scip, sol_val)){
            for(auto it = sol_pool.begin(); it != sol_pool.end(); it++){
               sort(it->begin(), it->end());
               SCIP_CALL(probdata->addPackVar(scip, (*it), TRUE));
            }
      (*result) = SCIP_SUCCESS;
   }
   else if(sol_type == Infeasible){ // stop the pricing iteration
      (*result) = SCIP_SUCCESS;
      //SCIPdebugMessage("termination by success");
   }
   else{
      SCIPdebugMessage("Unkown %lf\n", sol_val);
   }
   return SCIP_OKAY;
} 

/**@} */
