/**@file   conic_knapsack_pwlbc.cpp
 * @brief  Conic Knapsack pwl bc algorithm
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include "conic_knapsack_pwlbc.h"
#include <limits>
#include <algorithm>
#include <functional>
#include <ilcplex/ilocplex.h>
#include <thread> 
#include <mutex>  




/** solve a conic knapsack problem to detemine the upper bound timie_limit,
 *   ub  = max y, mu* x == y,  D * sqrt(b * x) <= c -y. If summu > ub, no feaisbe solution
 *   If just a feasible solution is returned, it is not a valid upper bound.
 *   a binary search finds a valid upper bound y s.t. mu*x >=y, D * sqrt(b * x) <= c -y is infeasible.
*/
void ub_micp(
	const vector<SCIP_Real> & mus, /** mu of items */
	const vector<SCIP_Real> & bs, /** b of items */
	SCIP_Real Dalpha, /** Dalpha */
	SCIP_Real capacity, /** capacity */
	int numitems, /** number of items */
	const vector<pair<int,int>>& items_diff, /** items in same constraints */
   SCIP_Real & left_bd, /** left bound of ub */
   SCIP_Real & right_bd, /** right bound of ub */
	SOLTYPE_CKNAP & sol_type, /** solution type */
   SCIP_Real time_limit /** solving time left */
){
   IloEnv env;

   try {
		// Create the model, populate by row
		IloModel model(env);
		IloCplex cplex(env);
		IloNumVarArray x_vars(env, numitems);
		IloNumVarArray xc_vars(env, numitems);
		IloNumVar z(env, 0.0, IloInfinity, ILOFLOAT); 
		IloNumVar y(env, 0.0, capacity, ILOFLOAT); 

		IloExpr expr_sum_mu_x(env);
		IloExpr expr_sum_xc2(env);
		IloExpr expr_obj(env);
		for(int item = 0; item < numitems; item++){
		x_vars[item] = IloNumVar(env, 0.0, 1.0, ILOBOOL); // add the binary variable for knapsack
		xc_vars[item] = IloNumVar(env, 0.0, sqrt(bs[item]), ILOFLOAT); // add the continuous variable for knapsack
		model.add(sqrt(bs[item]) * x_vars[item] <= xc_vars[item]); // add the coupling constraint for knapsack
		expr_sum_mu_x +=  mus[item] * x_vars[item];
		expr_sum_xc2 += xc_vars[item]*xc_vars[item]; 
		}

		// add the same and different constraints
		for(auto pair_item: items_diff){
			model.add(x_vars[pair_item.first] + x_vars[pair_item.second] <= 1);
		} 
		model.add(  expr_sum_mu_x  == y);
		model.add( y +  Dalpha *z  <= capacity); 
		model.add( expr_sum_xc2 <= z*z ); // second order cone constraint
		model.add(IloMaximize(env, y)); // set the maximization objective

		// Extract model.
		cplex.extract(model);

		// set the time limit in CPU seconds
		cplex.setParam(IloCplex::Param::ClockType, 1);
		cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
		cplex.setOut(env.getNullStream());

		sol_type = Unknown;
		cplex.solve();
	if( cplex.getStatus() ==  IloAlgorithm::Infeasible){
		sol_type = Infeasible;
		}
		else if(cplex.getStatus() ==  IloAlgorithm::Optimal || cplex.getStatus() ==  IloAlgorithm::Feasible){
		left_bd = cplex.getObjValue();
		right_bd = cplex.getBestObjValue();
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


/** solve a mixed integer convex quadratically constrained problem to detemine the lower bound before timie_limit,
 *   lb = min y,  D*D*b*x >= (c-y)^2, mu * x<= y. If summu < lb, all mu*x <= summu, D*b*x <= (c -summu)*2
 *   If just a feasible solution is returned,  not a valid lower bound. 
 *    a binary search finds a valid upper bound y s.t. D*D*b*x >= (c-y)^2,  mu * x<= y is infeasibke
*/
void lb_miqcp(
	const vector<SCIP_Real> & mus, /** mu of items */
	const vector<SCIP_Real> & bs, /** b of items */
	SCIP_Real Dalpha, /** Dalpha */
	SCIP_Real capacity, /** capacity */
	int numitems, /** number of items */
	const vector<pair<int,int>>& items_diff, /** items in same constraints */
   SCIP_Real & left_bd, /** left bound of ub */
   SCIP_Real & right_bd, /** right bound of ub */
	SOLTYPE_CKNAP & sol_type, /** solution type */
   SCIP_Real time_limit /** solving time left */
){
   IloEnv env;

   try {
         // Create the model, populate by row
         IloModel model(env);
         IloCplex cplex(env);
         IloNumVarArray x_vars(env, numitems);
         IloNumVar z(env, 0.0, IloInfinity, ILOFLOAT); 
		 IloNumVar y(env, 0.0, capacity, ILOFLOAT); 
         

         IloExpr expr_sum_mu_x(env);
         IloExpr expr_sum_b_x(env);
         IloExpr expr_obj(env);
         for(int item = 0; item < numitems; item++){
            x_vars[item] = IloNumVar(env, 0.0, 1.0, ILOBOOL); 
            expr_sum_mu_x +=  mus[item] * x_vars[item];
			expr_sum_b_x += bs[item] * x_vars[item];
         }

         // add the same and different constraints
         for(auto pair_item: items_diff){
               model.add(x_vars[pair_item.first] + x_vars[pair_item.second] <= 1);
         } 
		 model.add(  expr_sum_mu_x  <= y);
         model.add(   Dalpha * Dalpha * expr_sum_b_x  >= (capacity - y)*(capacity - y)); 
         model.add(IloMinimize(env, y)); // set the minimization objective

         // Extract model.
         cplex.extract(model);

         // set the time limit in CPU seconds
         cplex.setParam(IloCplex::Param::ClockType, 1);
         cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
         cplex.setOut(env.getNullStream());

        sol_type = Unknown;
		cplex.solve();
        if( cplex.getStatus() ==  IloAlgorithm::Infeasible){
            sol_type = Infeasible;
         }
         else if(cplex.getStatus() ==  IloAlgorithm::Optimal || cplex.getStatus() ==  IloAlgorithm::Feasible){
            left_bd = cplex.getObjValue();
            right_bd = cplex.getBestObjValue();
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
){
	SCIP_Real valid_ub = capacity, valid_lb = 0; // valid upper/lower bound
	SCIP_Real search_left_ub = 0,  search_left_lb = 0, search_right_ub = capacity,  search_right_lb = capacity; // binary search left/right upper/lower bound
	SCIP_Real ub_time = time_limit, lb_time = time_limit; // set timie limit of search
	SOLTYPE_CKNAP sol_type_ub = Unknown, sol_type_lb = Unknown; // exact algorithm status
	SCIP_Real soltime_before,  soltime_after;
	SCIP_Real test_ub, test_lb;

	soltime_before	=SCIPgetSolvingTime(scip);
	ub_micp( mus,  bs,  Dalpha, capacity,  numitems, items_diff,   search_left_ub,  search_right_ub ,  sol_type_ub, ub_time  ); // find upper bound
	soltime_after	=SCIPgetSolvingTime(scip);
	ub_time -= (soltime_after - soltime_before);
	valid_ub = fmin(capacity, search_right_ub + tol);

	soltime_before	=SCIPgetSolvingTime(scip);
	lb_miqcp( mus,  bs,  Dalpha, capacity,  numitems, items_diff,  search_left_lb, search_right_lb,  sol_type_lb, lb_time);
	soltime_after	=SCIPgetSolvingTime(scip);
	lb_time -= (soltime_after - soltime_before);
	valid_lb = fmax(0, search_left_lb - tol);
	lb = valid_lb;
	ub = valid_ub;
}



// This is the class implementing the generic callback interface.
// It has a main function
//    - lazySOC: adds the second order cone separation cut as a lazy constraint.
//	
class SOCCallback: public IloCplex::Callback::Function {

   public:
	// Empty constructor is forbidden.
	SOCCallback ()
	{}

	// Copy constructor is forbidden.
	SOCCallback(const SOCCallback &tocopy);

    // Variables for items in {0,1}.
    IloNumVarArray x_s;

	// constraint data,
	vector<SCIP_Real> mus, bs;
	SCIP_Real capacity, Dalpha;
	int numitems;

	// Constructor with data.
	SOCCallback(const IloNumVarArray & x_s_,
					const vector<SCIP_Real>  &mus_,
					const vector<SCIP_Real>  &bs_,
					const SCIP_Real capacity_,
					const SCIP_Real Dalpha_,
					const int numitems_):
		x_s(x_s_), mus(mus_), bs(bs_), capacity(capacity_), Dalpha(Dalpha_), numitems(numitems_)
	{
	}

	// Lazy constraint callback to enforce the SOC constraints, i.e. mus^t x_s +Dalpha \sqrt{bs^t x_s} <= capacity
	// If used then the callback is invoked for every integer feasible
	// solution CPLEX finds. For solution x_s_
	// constraint
	//    mus^t x_s +Dalpha \sqrt{bs^t x_s} <= capacity
	// is satisfied. If not then it adds the violated constraint as lazy
	// constraint:
	//  Dalpha (( \sum_i  b_s[i] x_s_[i]  x_s[i] ) / sqrt{\sum_i bs[i]x_s_[i]}) <= capacity - \sum_i mus[i] x_s[i]
	inline void
	lazySOC (const IloCplex::Callback::Context &context)  {
		if ( !context.isCandidatePoint() )
			throw IloCplex::Exception(-1, "Unbounded solution");
		int coversize = 0;
		SCIP_Real  summu = 0, sumb = 0;
		for(int item = 0; item < numitems; item++){
			IloNum x_val = context.getCandidatePoint(x_s[item]);
			if(IloRound(x_val) == 1){ // find non zero items in the solution
				sumb += bs[item];
				summu += mus[item];
				coversize++;
			}
		}
		//SCIPdebugMessage("summu: %f\n", summu);
		if(summu + Dalpha*sqrt(sumb) < capacity + tol){ //  valid
			return;
		}
		SCIP_Real coef = Dalpha / sqrt(sumb);
		vector<SCIP_Real> coeffs(mus.begin(), mus.end());
		for(int item = 0; item < numitems; item++){  // construct outer approximation cut
			bool covered = IloRound(context.getCandidatePoint(x_s[item])) == 1;
			if(covered){
				coeffs[item] += bs[item]*coef;
			}
		}
		IloNumExpr coeffs_expr = IloExpr(context.getEnv());
		for(int item =0; item < numitems; item++){
			coeffs_expr += coeffs[item]*x_s[item];
		}
		SCIP_Real rhs = capacity + tol;
		context.rejectCandidate(coeffs_expr <= rhs );
		coeffs_expr.end();	
	}

	// This is the function that we have to implement and that CPLEX will call
	// during the solution process at the places that we asked for.
	virtual void invoke (const IloCplex::Callback::Context &context) ILO_OVERRIDE;

	/// Destructor
	virtual ~SOCCallback(){};
};


// This is the class implementing the generic callback interface.
// It has a main function
//    - nodeInfo: print the node information.
//	
class NodeCallback: public IloCplex::Callback::Function {

   public:
	// Empty constructor is forbidden.
	NodeCallback ()
	{}

	// Copy constructor is forbidden.
	NodeCallback(const NodeCallback &tocopy);

    // Variables for items in {0,1}.
    IloNumVarArray x_s;

	// constraint data,
	vector<SCIP_Real> mus, bs;
	SCIP_Real capacity, Dalpha;
	int numitems, max_record;
	SCIP_Real min_mu = MAXFLOAT, max_mu = 0;

	// Constructor with data.
	NodeCallback(const IloNumVarArray & x_s_,
					const vector<SCIP_Real>  &mus_,
					const vector<SCIP_Real>  &bs_,
					const SCIP_Real capacity_,
					const SCIP_Real Dalpha_,
					const int numitems_):
		x_s(x_s_), mus(mus_), bs(bs_), capacity(capacity_), Dalpha(Dalpha_), numitems(numitems_)
	{
	}

	// Lazy constraint callback to enforce the SOC constraints, i.e. mus^t x_s +Dalpha \sqrt{bs^t x_s} <= capacity
	// If used then the callback is invoked for every integer feasible
	// solution CPLEX finds. For solution x_s_
	// constraint
	//    mus^t x_s +Dalpha \sqrt{bs^t x_s} <= capacity
	// is satisfied. If not then it adds the violated constraint as lazy
	// constraint:
	//  Dalpha (( \sum_i  b_s[i] x_s_[i]  x_s[i] ) / sqrt{\sum_i bs[i]x_s_[i]}) <= capacity - \sum_i mus[i] x_s[i]
	inline void
	nodeInfo (const IloCplex::Callback::Context &context)  {
		if(context.getRelaxationStatus(0) == IloCplex::CplexStatus::Optimal || context.getRelaxationStatus(0) == IloCplex::CplexStatus::Infeasible){
			SCIP_Real  summu = 0, sumb = 0;
			for(int item = 0; item < numitems; item++){
				IloNum x_val = context.getRelaxationPoint(x_s[item]);
				sumb += bs[item] * x_val;
				summu += mus[item] * x_val;
			}
			min_mu =   summu < min_mu ? summu : min_mu;
			max_mu =   summu > max_mu ? summu : max_mu;
			//SCIPdebugMessage("branch: summu: %f in [%f,%f]\n", summu, min_mu, max_mu);
		}
	}

	// This is the function that we have to implement and that CPLEX will call
	// during the solution process at the places that we asked for.
	virtual void invoke (const IloCplex::Callback::Context &context) ILO_OVERRIDE;

	/// Destructor
	virtual ~NodeCallback(){};
};

// Implementation of the invoke method.
//
// This is the method that we have to implement to fulfill the
// generic callback contract. CPLEX will call this method during
// the solution process at the places that we asked for.
void
SOCCallback::invoke (const IloCplex::Callback::Context &context)
{
   if ( context.inCandidate() )
        lazySOC(context);
}

void
NodeCallback::invoke (const IloCplex::Callback::Context &context)
{
   if ( context.inRelaxation() )
        nodeInfo(context);
}


/** the method solves a relaxed MILP relaxation (2d piece-wise linear over Estimator knapsack problem) for the conic IP problem,
 * with separation cut
*/
void solve_conic_rel_milp_cut(
	const vector<SCIP_Real> & objs, /** objective coefficients of items */
	const vector<SCIP_Real> & mus, /** mu of items */
	const vector<SCIP_Real> & bs, /** b of items */
	SCIP_Real Dalpha, /** Dalpha */
	SCIP_Real capacity, /** capacity */
	int numitems, /** number of items */
	const vector<pair<int,int>>& items_diff, /** items in same constraints */
	BreakPoints & breakpoints, /* quadratic Estimator */
	const conf & algo_conf, /** algorithm configuration*/
	SCIP_Real time_limit ,  /** solving time left */
	SCIP_Real target_lb, /* the tagret lower bound bound  */
    list<list<int>> & sol_pool, /* solutions pools, solution items are assumed to be sorted */
	pair<SCIP_Real, SCIP_Real>& mubd, /* lower and upper bound of mu*/
	SCIP_Real & sol_val, /** solution value */
	SCIP_Real & sol_ub, /** solution value upper bound */
    SOLTYPE_CKNAP & sol_type, /** solution type */
	SCIP_Real & mu_val_rel, /** solution mu */
	SCIP_Real & sol_time /** solution time */
){
	const list<pt_info> & break_points = breakpoints.get_break_points();
	const SCIP_Real left_slope = breakpoints.get_left_slope();
	const SCIP_Real right_slope = breakpoints.get_right_slope();
	const SCIP_Real lb = breakpoints.get_lb();
	const SCIP_Real ub = breakpoints.get_ub();
	const int num_bps =  break_points.size();

    IloEnv env;
    try {
		// Create the model, populate by row
		IloModel model(env);
		IloNumVarArray x_vars(env, numitems);
		IloNumVar mu(env, 0.0, ub, ILOFLOAT); 
		IloNumVar b(env, 0.0, capacity* capacity, ILOFLOAT); 
		IloNumVar quad(env, 0.0, capacity*capacity, ILOFLOAT);
		IloNumArray sample_xs(env, break_points.size() );
		IloNumArray sample_fxs(env, break_points.size());

		IloExpr expr_sum_mu_x(env);
		IloExpr expr_sum_b_x(env);
		IloExpr expr_obj(env);
		IloExpr expr_pl(env); 
		for(int item = 0; item < numitems; item++){
			x_vars[item] = IloNumVar(env, 0.0, 1.0, ILOBOOL); // add the binary variable for knapsack
			expr_sum_mu_x +=  mus[item] * x_vars[item];
			expr_sum_b_x += bs[item] * x_vars[item];
			expr_obj += objs[item]  * x_vars[item];           
		}
		// add the same and different constraints
		for(auto pair_item: items_diff){
			model.add(x_vars[pair_item.first] + x_vars[pair_item.second] <= 1);
		} 

		int i = 0;
		for(auto it = break_points.begin(); it != break_points.end(); it++){
			sample_xs[i]= it->first; 
			sample_fxs[i] = it->second; 
			i++;
		}

		bool explicit_pwl = true;
		if(explicit_pwl)
		{
			model.add(b <= IloPiecewiseLinear(mu, left_slope, sample_xs,  sample_fxs, right_slope));
		}
		else{
			// multi choice model
			int num_pcs = num_bps - 1;
			IloExpr expr_mcobj(env); 
			IloExpr expr_mcvar(env); 
			IloNumVarArray y_vars(env, num_pcs);
			IloNumVarArray z_vars(env, num_pcs);
			expr_mcvar += lb;
			auto p = break_points.begin();
			for(int  i = 0; i < num_pcs; i++){
				SCIP_Real x0 = p->first;
				SCIP_Real fx0 = p->second;
				p++;
				SCIP_Real x1 = p->first;
				SCIP_Real fx1 = p->second;		
				y_vars[i] = IloNumVar(env, 0.0, 1.0, ILOBOOL);
				z_vars[i] = IloNumVar(env, 0.0, (x1 - lb), ILOFLOAT);		
				model.add((x0 - lb) * y_vars[i] <= z_vars[i]);
				model.add(z_vars[i] <= (x1 - lb) * y_vars[i]);	
				expr_mcvar += z_vars[i];
				expr_mcobj += fx0 * y_vars[i] + ((fx1 - fx0) / (x1 - x0)) * (z_vars[i] - (x0 - lb) * y_vars[i]);
			}
			model.add(expr_mcvar == mu);
			model.add(b <= expr_mcobj);
			model.add(b <= IloPiecewiseLinear(mu, left_slope, sample_xs,  sample_fxs, right_slope));
			IloSOS1(env, y_vars);
		}

		model.add(  expr_sum_mu_x   == mu); 
		model.add(  Dalpha*Dalpha * expr_sum_b_x    ==  b); 
        model.add(  expr_obj  >= (target_lb - tol)); 
		model.add(IloMaximize(env, expr_obj)); // set the maximization objective

		// set cplex
		IloCplex cplex(model);
		cplex.setParam(IloCplex::Param::ClockType, 1);
		cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
		cplex.setParam(IloCplex::Param::Emphasis::MIP, 3);
		cplex.setParam(IloCplex::Param::Threads,algo_conf.is_parallelscplex ? 0 : 1);
		cplex.setOut(env.getNullStream());
		
		
		SOCCallback cb( x_vars, mus, bs, capacity, Dalpha, numitems);
		NodeCallback ncb( x_vars, mus, bs, capacity, Dalpha, numitems);
		CPXLONG contextMask = 0;
		contextMask |= IloCplex::Callback::Context::Id::Candidate;
		contextMask |= IloCplex::Callback::Context::Id::Relaxation;

      	// If contextMask is not zero we add the callback.
      	if ( contextMask != 0 ){
         	cplex.use(&cb, contextMask);
			cplex.use(&ncb, contextMask);
		}

        sol_type = Unknown;
        cplex.solve();
		sol_time = cplex.getTime();
		if( cplex.getStatus() ==  IloAlgorithm::Infeasible){
			printf("infeas\n");
            sol_type = Infeasible;
        }
        else if(cplex.getStatus() ==  IloAlgorithm::Optimal || cplex.getStatus() ==  IloAlgorithm::Feasible){
			mubd.first = ncb.min_mu > ub ? lb: ncb.min_mu;
			mubd.second = ncb.max_mu < lb ? ub: ncb.max_mu;
			SCIP_Real summu = 0, sumb = 0;
			list<int> items_bin;
            for(int item = 0; item < numitems; item++){
                  if(IloRound(cplex.getValue(x_vars[item])) == 1){
                     items_bin.push_back(item);
					 summu += mus[item];
					 sumb += bs[item];
                  }
            }
			mu_val_rel = summu;
			//SCIPdebugMessage("final summu: %f %f %f %f\n", summu, mubd.first, mubd.second, cplex.getTime());
            sol_val =  cplex.getObjValue();
			sol_ub = cplex.getBestObjValue();
            if(cplex.getStatus() ==  IloAlgorithm::Optimal){
               sol_type = Optimal;
            }
            else{
               sol_type = Feasible_Exact;
            }
			if(summu + Dalpha*sqrt(sumb) < capacity + tol){
				sol_pool.push_back(items_bin);
			}
         }
         else{
            cout << "none";
         }
	} catch (IloException &e) {
		cerr << "IloException at relaxed MILP: " << e << endl;
		if (env.getImpl())
			env.end();
		::abort();
	} catch (string& e) {
		cerr << e << "unkonown error"<< endl;
		if (env.getImpl())
			env.end();
		::abort();
	} 
	env.end();
};




/** solve conic knapsack problem via PWLBC before timie_limit,
*   either find a high quality solution with the value greater than the tagret lower bound,  or conclude that there is no such solution
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
	SCIP_Real heur_mu,/* heur mu*/
	const vector<pair<int,int>>& items_diff, /** items in same constraints */
	const conflict_graph & conflict, /** the conflict graph*/
    list<list<int>> & sol_pool, /* solutions pools, solution items are assumed to be sorted */
	SCIP_Real & sol_val, /** solution value */
	SCIP_Real & sol_ub, /** solution value upper bound */
	SOLTYPE_CKNAP & sol_type, /** solution type */
    SCIP_Real time_limit,  /** solving time left */
    SCIP_Real target_lb /* the tagret lower bound bound  */
){
	SCIP_Real timelimit;
	SCIP_Real sol_time;
	BreakPoints breakpoints(init_estimator.breakpoints);
	auto & bb_breakpoints = breakpoints;
	// query knn results 
	SCIP_Real knn_time;
	pair<SCIP_Real,SCIP_Real> knn_mubd;
	if(algo_conf.knn_mode !=1 && !SCIPisEQ(scip, algo_conf.point_ratio,  1)){
		SCIP_Real knn_start = SCIPgetSolvingTime(scip);
		vector<pair<SCIP_Real,SCIP_Real>> knn_results;
		knn_mubd = init_estimator.knnregression(objs, algo_conf.kneighbors, knn_results, algo_conf.knn_mode);
		BreakPoints concen_breakpoints(breakpoints.getCap(), breakpoints.get_lb(), breakpoints.get_ub(), breakpoints.getNum(), knn_mubd.first, knn_mubd.second, algo_conf.point_ratio);	
		bb_breakpoints = concen_breakpoints;
		knn_time = SCIPgetSolvingTime(scip) - knn_start;
		stat_prc.cum_knn_time += knn_time;
	}
	
	int num_iter = 0;
	while(true){
		SCIPgetRealParam(scip, "limits/time", &timelimit) ;
		if( !SCIPisInfinity(scip, timelimit) )
			timelimit -= SCIPgetSolvingTime(scip);
		if(timelimit < 0){
			return;
		}
		timelimit  = timelimit < time_limit ? timelimit : time_limit;
		// solve the relxation
		SOLTYPE_CKNAP sol_type_rel = Unknown;
		SCIP_Real sol_val_rel = 0, sol_ub_rel = MAXFLOAT, mu_val_rel = 0;
		pair<SCIP_Real, SCIP_Real> mubd; /* lower and upper bound of mu*/
		solve_conic_rel_milp_cut(objs,  mus,  bs,  Dalpha,  capacity,  numitems,  items_diff, bb_breakpoints, algo_conf, timelimit,  target_lb, sol_pool, mubd, sol_val_rel,  sol_ub_rel,  sol_type_rel, mu_val_rel, sol_time); // solve the relaxation
		stat_prc.cum_sol_time += sol_time;
		// update estimator and knn
		init_estimator.add(objs, mubd);
		//SCIPdebugMessage("[%f, %f] [%f, %f] %f/%f %f/%f \n", knn_mubd.first, knn_mubd.second,  mubd.first, mubd.second,  sol_time, stat_prc.cum_sol_time, knn_time, stat_prc.cum_knn_time);
		if(sol_type_rel == Infeasible){ // conclude no feasible pricing solution
			sol_type = Infeasible;
			break;
		}
		else{ // check wether the solution for the relaxed problem is feasible in the conic knapsack problem
			if( !sol_pool.empty()){ // optimal solution 
				sol_type =sol_type_rel;
				sol_val = sol_val_rel;
				sol_ub = sol_ub_rel;
	    		SCIPgetRealParam(scip, "limits/time", &timelimit) ;
				if( !SCIPisInfinity(scip, timelimit) )
					timelimit -= SCIPgetSolvingTime(scip);
				break;
			}
			else{
				bb_breakpoints.insert_x(mu_val_rel, 1);
				SCIPdebugMessage("round %d\n", num_iter);
			}
		}
		num_iter++;
		//break; // debug
	}
};
