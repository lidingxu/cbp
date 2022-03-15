/**@file   conic_knapsack.cpp
 * @brief  Conic Knapsack solver
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <ilcplex/ilocplex.h>
#include "conic_knapsack.h"
#include "conic_knapsack_heur.h"
#include "conic_knapsack_approx.h"


/** solve conic knapsack problem by cplex before timie_limit,
*   find a solution with value greater than the tagret lower bound, and store the solution
*   possible return values for:
*  - Sol_Ge_Bd : a solution with value greater than the tagret lower bound is found
*  - No_Sol_Ge_Bd :  no solution with value greater than the tagret lower bound
*  - Unsolved: the problem greater than the tagret lower bound is unsolved.
*/
void solve_micp(
	const vector<SCIP_Real> & objs, /** objective coefficients of items */
	const vector<SCIP_Real> & mus, /** mu of items */
	const vector<SCIP_Real> & bs, /** b of items */
	SCIP_Real Dalpha, /** Dalpha */
	SCIP_Real capacity, /** capacity */
	int numitems, /** number of items */
   conf algo_conf, /** algorithm configuration */
	const vector<pair<int,int>>& items_diff, /** items in same constraints */
	list<list<int>> & sol_pool, /* items in the maximal bin, assumed to be sorted */
   SCIP_Real & sol_val, /** solution value */
   SCIP_Real & sol_ub, /** solution value upper bound*/
	SOLTYPE_CKNAP & sol_type, /** solution type */
   SCIP_Real time_limit,  /** solving time left */
   SCIP_Real target_lb /* the tagret lower bound bound  */
){
   IloEnv env;


   try {
         // Create the model, populate by row
         IloModel model(env);
         IloCplex cplex(env);
         IloNumVarArray x_vars(env, numitems);
         IloNumVarArray xc_vars(env, numitems);
         IloNumVar z(env, 0.0, IloInfinity, ILOFLOAT); 
         

         IloExpr expr_sum_mu_x(env);
         IloExpr expr_sum_xc2(env);
         IloExpr expr_obj(env);
         for(int item = 0; item < numitems; item++){
            x_vars[item] = IloNumVar(env, 0.0, 1.0, ILOBOOL); // add the binary variable for knapsack
            xc_vars[item] = IloNumVar(env, 0.0, sqrt(bs[item]), ILOFLOAT); // add the continuous variable for knapsack
            model.add(sqrt(bs[item]) * x_vars[item] <= xc_vars[item]); // add the coupling constraint for knapsack
            expr_sum_mu_x +=  mus[item] * x_vars[item];
            expr_sum_xc2 += xc_vars[item]*xc_vars[item];
            expr_obj += objs[item] * x_vars[item];      
         }

         // add the same and different constraints
         for(auto pair_item: items_diff){
            model.add(x_vars[pair_item.first] + x_vars[pair_item.second] <= 1);
         } 
         model.add(  expr_sum_mu_x +  Dalpha *z  <= capacity); 
         model.add( expr_sum_xc2 <= z*z ); // second order cone constraint
         model.add(  expr_obj   >= target_lb - tol); 
         model.add(IloMaximize(env, expr_obj)); // set the maximization objective

         // Extract model.

         cplex.extract(model);


         // set the time limit in CPU seconds
         cplex.setParam(IloCplex::Param::ClockType, 1);
         cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
         cplex.setParam(IloCplex::Param::Threads, algo_conf.is_parallelscplex ? 0 : 1);
         cplex.setOut(env.getNullStream());

        sol_type = Unknown;
         cplex.solve();
         if( cplex.getStatus() ==  IloAlgorithm::Infeasible){
            sol_type = Infeasible;
         }
         else if(cplex.getStatus() ==  IloAlgorithm::Optimal || cplex.getStatus() ==  IloAlgorithm::Feasible){
			   list<int> items_bin;
            for(int item = 0; item < numitems; item++){
                  if(IloRound(cplex.getValue(x_vars[item])) == 1){
                     items_bin.push_back(item);
                  }
            }
            sol_pool.push_back(items_bin);
            sol_val =  cplex.getObjValue();
            sol_ub = cplex.getBestObjValue();
            SCIPdebugMessage("%lf %lf\n", sol_val, sol_ub);
            if(cplex.getStatus() ==  IloAlgorithm::Optimal){
               sol_type = Optimal;
            }
            else{
               sol_type = Feasible_Exact;
            }
         }
         else{
            cout << "none";
         }
         env.end();
   } catch (IloException &e) {
      cerr << "IloException: " << e << endl;
      if (env.getImpl())
         env.end();
      ::abort();
   } catch (string& e) {
      cerr << e << endl;
      if (env.getImpl())
         env.end();
      ::abort();
   }  
}

/** solve linear knapsack problem by cplex before timie_limit,
*   find a solution with value greater than the tagret lower bound, and store the solution
*   possible return values for:
*  - Sol_Ge_Bd : a solution with value greater than the tagret lower bound is found
*  - No_Sol_Ge_Bd :  no solution with value greater than the tagret lower bound
*  - Unsolved: the problem greater than the tagret lower bound is unsolved.
*/
void linear_knap_cplex(
	const vector<SCIP_Real> & objs, /** objective coefficients of items */
	const vector<SCIP_Real> & mus, /** mu of items */
	SCIP_Real capacity, /** capacity */
	int numitems, /** number of items */
	const vector<pair<int,int>>& items_same, /** items in same constraints */
	const vector<pair<int,int>>& items_diff, /** items in same constraints */
   const list<int> & fixed, /** items fixed to zeros */
	vector<int> & items_bin, /* items in the maximal bin, assumed to be sorted */
   SCIP_Real & sol_val, /** solution value */
	SOLTYPE_CKNAP & sol_type, /** solution type */
   SCIP_Real time_limit,  /** solving time left */
   SCIP_Real target_lb /* the tagret lower bound bound  */
){
   IloEnv env;
   try {
         // Create the model, populate by row
         IloModel model(env);
         IloCplex cplex(env);
         IloObjective obj(env);
         IloNumVarArray x_vars(env, numitems);
         

         IloExpr expr_sum_mu_x(env);
         IloExpr expr_obj(env);
         for(int item = 0; item < numitems; item++){
            x_vars[item] = IloNumVar(env, 0.0, 1.0, ILOBOOL); // add the binary variable for knapsack
            expr_sum_mu_x +=  mus[item] * x_vars[item];
            expr_obj += objs[item] * x_vars[item];
         }
         // add the same and different constraints
         for(auto pair_item: items_same){
            model.add(x_vars[pair_item.first] == x_vars[pair_item.second]);
         }
         for(int item: fixed){
            //x_vars[item].setUb(0);//model.add(x_vars[item] == 0);
         }
         for(auto pair_item: items_diff){
            model.add(x_vars[pair_item.first] + x_vars[pair_item.second] <= 1);
         }
         model.add(  expr_sum_mu_x   <= capacity); 
         model.add(  expr_obj   >= target_lb - tol); 
         model.add(IloMaximize(env, expr_obj)); // set the maximization objective

         // Extract model.
         cplex.extract(model);

         // set the time limit in CPU seconds
         cplex.setParam(IloCplex::Param::ClockType, 1);
         cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
         cplex.setOut(env.getNullStream());

         cplex.solve();
        if( cplex.getStatus() ==  IloAlgorithm::Infeasible){
            sol_type = Infeasible;
         }
         else if(cplex.getStatus() ==  IloAlgorithm::Optimal || cplex.getStatus() ==  IloAlgorithm::Feasible){
            SCIP_Real objval = cplex.getObjValue();
            for(int item = 0; item < numitems; item++){
                  if(IloRound(cplex.getValue(x_vars[item])) == 1){
                     items_bin.push_back(item);
                  }
            }
            sol_val =  objval;
            if(cplex.getStatus() ==  IloAlgorithm::Optimal){
               sol_type = Optimal;
            }
            else{
               sol_type = Feasible_Exact;
            }
         }
         else{
            cout << "none";
         }
         env.end();
   } catch (IloException &e) {
      cerr << "IloException: " << e << endl;
      if (env.getImpl())
         env.end();
      ::abort();
   } catch (string& e) {
      cerr << e << endl;
      if (env.getImpl())
         env.end();
      ::abort();
   }
}

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
   SCIP_Real time_limit,  /** solving time left */
   SCIP_Real target_lb /* the tagret lower bound bound  */
){
   if(Dalpha > tol){
      SCIP_Real time_a, time_b;
      SCIP_Real summu = -1;
      SCIP_Real alpha = 0.3;
      time_a = SCIPgetSolvingTime(scip);
      //SCIPdebugMessage("is heur%d %d %d \n", int(algo_conf.is_heur), int(algo_conf.is_misocp), int(algo_conf.is_parallelscplex));
      for(int item = 0; item < numitems; item++){
         stable_center[item] = (1 - alpha) * objs[item] + alpha * stable_center[item];
      }
      if(algo_conf.is_heur)
      {
      // fix some variable with zero objectives
         SCIP_Real sol_val_heur = 0;
         vector<int> items_bin_heur;
         // solve heuristics
         solve_conic_knap_heur(objs, mus, bs, Dalpha, capacity, numitems, items_diff, conflict, items_bin_heur, sol_val_heur, time_limit, target_lb);
         for(int item :  items_bin_heur){
            summu += mus[item];
         }
         SCIP_Real  algo_ratio = (stat_prc.col_heur + 0.0) / stat_prc.col_exact;
         // no need for exact pricing
         if(sol_val_heur > target_lb + tol && sol_val_heur > stop_pricing_obj && algo_ratio < MAXFLOAT){
            // stablization
            if(algo_conf.is_stablize){ 
               // use stablized center to update duals
               SCIP_Real sol_val_heur_stab = 0;
               vector<int> items_bin_heur_stab;
               // run heursitcs again
               solve_conic_knap_heur(stable_center, mus, bs, Dalpha, capacity, numitems, items_diff, conflict, items_bin_heur_stab, sol_val_heur_stab, time_limit, target_lb);
               // objs in original duals
               sol_val_heur_stab = 0;
               for(int item: items_bin_heur_stab){
                  sol_val_heur_stab += objs[item];
               }
               // negative reduced cost
               if(sol_val_heur_stab > target_lb + tol){
                  items_bin_heur = items_bin_heur_stab;
                  sol_val_heur = sol_val_heur_stab;
               }
            }
            sol_type = Feasible_Heur;
            sol_val = sol_val_heur;
            sol_pool.push_back(list<int>(items_bin_heur.begin(), items_bin_heur.end()));
            stat_prc.col_heur++;
            //SCIPdebugMessage("%f/%f\n", summu, capacity);
            time_b = SCIPgetSolvingTime(scip);
            stat_prc.time_heur += time_b - time_a;
            return;
         }
      }


      // use quadratic estimator to solve MILP relaxation
      
      sol_val= 0;
      sol_ub = MAXFLOAT;
      sol_type = Unknown;
      time_a = SCIPgetSolvingTime(scip);
      SCIP_Real relative_gap = 100;
      if(!algo_conf.is_misocp){
         solve_conic_knap_approx(scip, objs,  mus,  bs,  Dalpha,  capacity,  numitems,   init_estimator,  algo_conf, summu, items_diff,  conflict,  sol_pool,  
         sol_val,  sol_ub, sol_type,  time_limit ,  target_lb);
         // update stablization center;
         relative_gap = fabs(sol_val- sol_ub)/fabs(max(sol_val, sol_ub)) * 100;
         stat_prc.shf_log_sum_gap += log(stat_prc.shf_param + relative_gap);
         stat_prc.col_exact++;
         time_b = SCIPgetSolvingTime(scip);    
         stat_prc.time_exact += time_b - time_a;
      }
      else{
         solve_micp(objs, mus, bs, Dalpha, capacity, numitems, algo_conf, items_diff,  sol_pool, sol_val, sol_ub,  sol_type, time_limit, target_lb);
         // update stablization center;      
         relative_gap = fabs(sol_val- sol_ub)/fabs(max(sol_val, sol_ub)) * 100;
         stat_prc.shf_log_sum_gap += log(stat_prc.shf_param + relative_gap);
         stat_prc.col_exact++;
         SCIPdebugMessage("%lf %lf\n", relative_gap, exp(stat_prc.shf_log_sum_gap / stat_prc.col_exact) - stat_prc.shf_param);
         time_b = SCIPgetSolvingTime(scip);    
         stat_prc.time_exact += time_b - time_a;
      }
   }
   else{
      //list<int> fixed;
      SCIPdebugMessage("linear binpack not callable\b");
      abort();
      //linear_knap_cplex(objs, mus, capacity, numitems, items_same, items_diff, fixed, items_bin, sol_val, sol_type, time_limit,  target_lb);
   }
}