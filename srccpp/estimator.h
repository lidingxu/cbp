/**@file   BreakPoints.h
 * @brief  concave quadratic piece-wise linear BreakPoints
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#pragma once
#include "scip/scip_numerics.h"
#include <vector>
#include <list>
#include <tuple>
#include <utility>
using namespace std;

typedef pair<SCIP_Real,  SCIP_Real> pt_info; // x, fx

/** the over BreakPoints  of (capacity-z)^2*/
class BreakPoints{
	SCIP_Real capacity, lb, ub;
	list<pt_info> break_points;  //x increasing order , fx  decreasing order
public:
	BreakPoints(){
	};
	
	// BreakPoints constructor 
	BreakPoints(
		SCIP_Real capacity_, /**  the capacity*/
		int num_break_points_ /** the number of break points*/
	);

	// eststimator constructor with bound information
	BreakPoints(
		SCIP_Real capacity_, /**  the capacity*/
		SCIP_Real lb_, /** lower bound of z*/
		SCIP_Real ub_, /** upper bound of z*/
		int num_break_points_ /** the number of break points*/
	);



	// eststimator constructor with point vector (sorted x)
	BreakPoints(
		SCIP_Real capacity_, /**  the capacity*/
		SCIP_Real lb_, /** lower bound of z*/
		SCIP_Real ub_, /** upper bound of z*/
		vector<SCIP_Real> xs_ /** sorted x*/
	);

	// eststimator constructor with bound information and concentration information
	BreakPoints(
		SCIP_Real capacity_, /**  the capacity*/
		SCIP_Real lb_, /** lower bound of z*/
		SCIP_Real ub_, /** upper bound of z*/
		int num_break_points_, /** the number of break points*/
		SCIP_Real concen_lb, /** concentrate lower bound of z*/
		SCIP_Real concen_ub, /** concentrate upper bound of z*/
		SCIP_Real concen_ratio /** the ratio concentrate points*/
	);

	// copy constructor
	BreakPoints(
		BreakPoints & BreakPoints_ /**  the BreakPoints to copy*/
	);

	// get break points
	inline const list<pt_info> & get_break_points(){
		return 	break_points;
	} ;

	// compute pt_info from fx
	inline pt_info  compute_pt(SCIP_Real fx){
		SCIP_Real x = fx2x(fx);
		return make_pair(x, fx);
	};
	
	// print pt info	
	void print();

	// insert x and 2*num_points + 1 around x
	void insert_x(SCIP_Real x, int num_points);

	// insert x inner function
	void insert_x_inner(SCIP_Real x, int num_points, list<pt_info>::iterator iter);

	void insert_points(vector<SCIP_Real> points);

	// check function
	void check();


	// get the left slope at the left end
	SCIP_Real get_left_slope();

	// get the right slope at the right end
	SCIP_Real get_right_slope();

	// x to fx
	inline SCIP_Real x2fx(SCIP_Real x){
		return (capacity - x)*(capacity - x);
	}

	// x to f(x)
	inline SCIP_Real x2dfx(SCIP_Real x){
		return -2*(capacity-x);
	}

	// x to df(x)
	inline SCIP_Real fx2x(SCIP_Real fx){
		return capacity - sqrt(fx);
	}

	inline SCIP_Real getCap(){
		return capacity;
	}

	inline int getNum(){
		return break_points.size();
	}

	// get upper bound
	inline SCIP_Real get_ub(){
		return ub;
	}

	// get lower bound
	inline SCIP_Real get_lb(){
		return lb;
	}
};

// Estimator
class Estimator{
public:
	BreakPoints breakpoints;
	pair<SCIP_Real, SCIP_Real> default_bd;
	list<tuple<vector<double>, pair<SCIP_Real, SCIP_Real>, int>> estimations;

	Estimator(BreakPoints & breakpoints_){
		breakpoints = breakpoints_;
		default_bd = make_pair(breakpoints_.get_lb(), breakpoints_.get_ub());
	};

	Estimator(){

	};

	pair<SCIP_Real, SCIP_Real> knnregression(const vector<SCIP_Real> & norm_cvec, int k, vector<pair<SCIP_Real, SCIP_Real>> & k_queries,  int knn_mode);
	void add(const vector<SCIP_Real> & norm_cvec, pair<SCIP_Real, SCIP_Real> mubd);

};
