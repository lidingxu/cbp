/**@file   estimator.cpp
 * @brief  concave quadratic piece-wise linear estimator
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include "estimator.h"
#include <limits>
//#include <ilcplex/ilocplex.h>

// default estimator constructor 
estimator::estimator(
		SCIP_Real capacity_, /**  the capacity*/
		int num_break_points_ /** the number of break points*/
){	
	capacity = capacity_;
	lb = 0;
	ub = capacity;
	int num_break_points = num_break_points_ > 2 ? num_break_points_ : 2;
	SCIP_Real flb = capacity - ub;
	SCIP_Real fub = capacity - lb;
	SCIP_Real interval = (fub * fub - flb * flb) / (num_break_points - 1);
	for(int i = num_break_points - 1; i >= 0; i--){
		SCIP_Real fx = interval * i + flb*flb;
		SCIP_Real x = fx2x(fx);
		break_points.push_back(make_pair(x, fx));
	}
}

// eststimator constructor with bound information
estimator::estimator(
	SCIP_Real capacity_, /**  the capacity*/
	SCIP_Real lb_, /** lower bound of z*/
	SCIP_Real ub_, /** upper bound of z*/
	int num_break_points_ /** the number of break points*/
){
	capacity = capacity_;
	lb = lb_;
	ub = ub_;
	int num_break_points = num_break_points_ > 2 ? num_break_points_ : 2;
	SCIP_Real flb = capacity - ub;
	SCIP_Real fub = capacity - lb;
	SCIP_Real interval = (fub * fub - flb * flb) / (num_break_points - 1);
	for(int i = num_break_points - 1; i >= 0; i--){
		SCIP_Real fx = interval * i + flb*flb;
		SCIP_Real x = fx2x(fx);
		break_points.push_back(make_pair(x, fx));
	}
	if(lb > 1e-6){
		break_points.push_front(make_pair(0, capacity * capacity));
		num_break_points++;
	}
}

// copy constructor
estimator::estimator(
	estimator & estimator_ /**  the estimator to copy*/
){
	capacity = estimator_.capacity;
	lb = estimator_.lb;
	ub = estimator_.ub;
	break_points = list<pt_info>(estimator_.break_points);
}





// print pt info	
void estimator::print(){
	printf("%f, %f", lb, ub);
	int i = 0;
	for(auto iter = break_points.begin(); iter != break_points.end(); iter++){
		printf("%f %f\n", get<0>(*iter), get<1>(*iter));
	}
}

// insert x and 2*num_points + 1 around x
void estimator::insert_x(SCIP_Real x, int num_points){
	//print();
	auto iter = break_points.begin();		
	insert_x_inner(x, num_points, iter);
}

// insert x inner function
void estimator::insert_x_inner(SCIP_Real x, int num_points, list<pt_info>::iterator iter){
	//auto iter = break_points.begin();
	SCIP_Real fx = x2fx(x);
	list<pt_info>::iterator left_iter, right_iter;
	while(true){
		if(x >= get<0>(*iter) ){
			iter++;
		}
		else{    // *(iter-1) <= x < *(iter)
			if(num_points == 0){
				break_points.insert(iter, make_pair(x, fx));
				return;
			}
			SCIP_Real fx_l = get<1>(*iter);
			iter--;
			SCIP_Real fx_u = get<1>(*iter);
			fx = (fx_l + fx_u) / 2;
			x = fx2x(fx);
			left_iter = iter;
			iter++;
			break_points.insert(iter, make_pair(x, fx));
			//print();
			right_iter = iter;
			break;
		}
	}
	auto search_iter = left_iter;
	int left_num = num_points, right_num = num_points;
	while(left_num > 0){
		SCIP_Real fx = get<1>(*left_iter);
		auto next_iter = next(left_iter);
		SCIP_Real next_fx = get<1>(*next_iter);
		pt_info pt = compute_pt((fx + next_fx) / 2);
		break_points.insert(next_iter, pt);
		left_num--;
		if( left_iter == break_points.begin()){
			break;
		}
		left_iter--;
	}	
	while(right_num > 0){
		auto prev_iter = prev(right_iter);
		SCIP_Real fx = get<1>(*right_iter);
		SCIP_Real prev_fx = get<1>(*prev_iter);		
		pt_info pt = compute_pt((fx + prev_fx) / 2);
		break_points.insert(right_iter, pt);	
		right_num--;
		if(right_iter == break_points.end()){
			break;
		}
		right_iter++;
	}
	insert_x_inner(x, 0,  search_iter);
	//print();
	//check();
}

// check function
void estimator::check(){
	SCIP_Real tmpx = -1, tmpfx = capacity * capacity + 1;
	bool good = true;
	for(auto pt: break_points){
		SCIP_Real x, fx;
		x = get<0>(pt);
		fx = get<1>(pt);
		if(x >= tmpx + 1e-4 && fx <= tmpfx - 1e-4){
			tmpx =x;
			tmpfx = fx;
		}
		else{
			good = false;
			break;
		}
	}
	if(!good){
		print();
	}
	assert(good);
}


// get the left slope at the left end
SCIP_Real estimator::get_left_slope(){
	auto iter = break_points.begin();
	SCIP_Real x0 = get<0>(*iter) ;
	SCIP_Real fx0 = get<1>(*iter);
	iter++;
	SCIP_Real x1 = get<0>(*iter);
	SCIP_Real fx1 = get<1>(*iter); 
	return (fx1 - fx0)/(x1 - x0);
}

// get the right slope at the right end
SCIP_Real estimator::get_right_slope(){
	auto iter = break_points.rbegin();
	SCIP_Real x0 = get<0>(*iter) ;
	SCIP_Real fx0 = get<1>(*iter);
	iter++;
	SCIP_Real x1 = get<0>(*iter);
	SCIP_Real fx1 = get<1>(*iter); 
	return (fx1 - fx0)/(x1 - x0);
}

