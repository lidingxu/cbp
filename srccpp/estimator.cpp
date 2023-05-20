/**@file   BreakPoints.cpp
 * @brief  concave quadratic piece-wise linear BreakPoints
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include "estimator.h"
#include <limits>
//#include <ilcplex/ilocplex.h>

// default BreakPoints constructor 
BreakPoints::BreakPoints(
		SCIP_Real capacity_, /**  the capacity*/
		int num_break_points_ /** the number of break points*/
){	
	capacity = capacity_;
	lb = 0;
	ub = capacity;
	int num_break_points = num_break_points_ > 10 ? num_break_points_ : 10;
	SCIP_Real interval = (ub - lb) / (num_break_points - 1);
	for(int i = 0; i < num_break_points; i++){
		SCIP_Real x = interval * i + lb;
		SCIP_Real fx = x2fx(x);
		break_points.push_back(make_pair(x, fx));
	}
}

// eststimator constructor with bound information
BreakPoints::BreakPoints(
	SCIP_Real capacity_, /**  the capacity*/
	SCIP_Real lb_, /** lower bound of z*/
	SCIP_Real ub_, /** upper bound of z*/
	int num_break_points_ /** the number of break points*/
){
	capacity = capacity_;
	lb = lb_;
	ub = ub_;
	int num_break_points = num_break_points_ > 10 ? num_break_points_ : 10;
	SCIP_Real interval = (ub - lb) / (num_break_points - 1);
	for(int i = 0; i < num_break_points; i++){
		SCIP_Real x = interval * i + lb;
		SCIP_Real fx = x2fx(x);
		break_points.push_back(make_pair(x, fx));
	}
}

// eststimator constructor with point vector (sorted x)
BreakPoints::BreakPoints(
	SCIP_Real capacity_, /**  the capacity*/
	SCIP_Real lb_, /** lower bound of z*/
	SCIP_Real ub_, /** upper bound of z*/
	vector<SCIP_Real> xs_ /** sorted x*/
){
	capacity = capacity_;
	lb = lb_;
	ub = ub_;
	for(SCIP_Real x: xs_){
		SCIP_Real fx = x2fx(x);
		break_points.push_back(make_pair(x, fx));
	}
}

// eststimator constructor with bound information and concentration information
BreakPoints::BreakPoints(
	SCIP_Real capacity_, /**  the capacity*/
	SCIP_Real lb_, /** lower bound of z*/
	SCIP_Real ub_, /** upper bound of z*/
	int num_break_points_, /** the number of break points*/
	SCIP_Real concen_lb, /** concentrate lower bound of z*/
	SCIP_Real concen_ub, /** concentrate upper bound of z*/
	SCIP_Real ratio /** concentrate number*/
){
	capacity = capacity_;
	lb = lb_;
	ub = ub_;
	//printf("%f %f\n", concen_ub, concen_lb);
	concen_lb = concen_lb > lb ? concen_lb : lb;
	concen_ub = concen_ub < ub ? concen_ub : ub;
	if(fabs(concen_ub - concen_lb) < 1e-1){
		concen_lb = lb;
		concen_ub = ub;
	}
	assert(concen_ub >= concen_lb);
	assert(concen_ub <= ub);
	assert(concen_lb >= lb);
	

	SCIP_Real total_region = ub - lb;
	SCIP_Real concen_region = concen_ub - concen_lb;
	SCIP_Real nonconcen_region = total_region - concen_region;
	SCIP_Real left_nonconcen_region = concen_lb - lb;
	SCIP_Real right_nonconcen_region = ub - concen_ub;

	int num_break_points = num_break_points_ > 10 ? num_break_points_ : 10;
	SCIP_Real concen_ratio =  min(concen_region / total_region * ratio, 1.0);
	int num_concen = int(num_break_points * concen_region / total_region);
	int num_left_nonconcen = int(num_break_points * (1 - concen_ratio) * (left_nonconcen_region) / (nonconcen_region) + 1);
	int num_right_nonconcen = int(num_break_points * (1 - concen_ratio) *  (right_nonconcen_region) / (nonconcen_region) + 1);
	num_break_points = num_concen +  num_left_nonconcen + num_right_nonconcen;

	//printf("%d %d %d\n", num_concen, num_left_nonconcen, num_right_nonconcen);

	if(num_left_nonconcen > 1){
		SCIP_Real left_interval = left_nonconcen_region / (num_left_nonconcen - 1);
		for(int i = 0; i < num_left_nonconcen; i++){
			SCIP_Real x = left_interval * i + lb;
			SCIP_Real fx = x2fx(x);
			break_points.push_back(make_pair(x, fx));
		}	
	}

	SCIP_Real concen_interval = concen_region / (num_concen - 1);	
	for(int i = num_left_nonconcen > 1; i < num_concen; i++){
		SCIP_Real x = concen_interval * i + concen_lb;
		SCIP_Real fx = x2fx(x);
		break_points.push_back(make_pair(x, fx));
	}	

	if(num_right_nonconcen > 1){
		SCIP_Real right_interval = right_nonconcen_region / (num_right_nonconcen - 1);	
		for(int i = 1; i < num_right_nonconcen; i++){
			SCIP_Real x = right_interval * i + concen_ub;
			SCIP_Real fx = x2fx(x);
			break_points.push_back(make_pair(x, fx));
		}	
	}	
}

// copy constructor
BreakPoints::BreakPoints(
	BreakPoints & BreakPoints_ /**  the BreakPoints to copy*/
){
	capacity = BreakPoints_.capacity;
	lb = BreakPoints_.lb;
	ub = BreakPoints_.ub;
	break_points = list<pt_info>(BreakPoints_.break_points);
}







// print pt info	
void BreakPoints::print(){
	printf("BreakPoints info: %f, %f\n", lb, ub);
	int i = 0;
	for(auto iter = break_points.begin(); iter != break_points.end(); iter++){
		printf("%f %f\n", get<0>(*iter), get<1>(*iter));
	}
}

// insert x and 2*num_points + 1 around x
void BreakPoints::insert_x(SCIP_Real x, int num_points){
	//print();
	auto iter = break_points.begin();		
	insert_x_inner(x, num_points, iter);
}

// insert x inner function
void BreakPoints::insert_x_inner(SCIP_Real x, int num_points, list<pt_info>::iterator iter){
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
void BreakPoints::check(){
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
			printf("%f %f %f %f\n", tmpx, tmpfx, x, fx);
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
SCIP_Real BreakPoints::get_left_slope(){
	auto iter = break_points.begin();
	SCIP_Real x0 = get<0>(*iter) ;
	SCIP_Real fx0 = get<1>(*iter);
	iter++;
	SCIP_Real x1 = get<0>(*iter);
	SCIP_Real fx1 = get<1>(*iter); 
	return (fx1 - fx0)/(x1 - x0);
}

// get the right slope at the right end
SCIP_Real BreakPoints::get_right_slope(){
	auto iter = break_points.rbegin();
	SCIP_Real x0 = get<0>(*iter) ;
	SCIP_Real fx0 = get<1>(*iter);
	iter++;
	SCIP_Real x1 = get<0>(*iter);
	SCIP_Real fx1 = get<1>(*iter); 
	return (fx1 - fx0)/(x1 - x0);
}


pair<SCIP_Real, SCIP_Real> Estimator::knnregression(const vector<SCIP_Real> & cvec, int k, vector<pair<SCIP_Real, SCIP_Real>> & k_queries, int knn_mode){
	vector<SCIP_Real> norm_cvec(cvec);
	SCIP_Real sum = 0;
	for(int i = 0; i < cvec.size(); i++){
		sum += cvec[i];
	}
	for(int i = 0; i < cvec.size(); i++){
		norm_cvec[i] /= sum;
	}
	k_queries = vector<pair<SCIP_Real, SCIP_Real>>(k, default_bd);
	auto k_dists = vector<SCIP_Real>(k, MAXFLOAT);
	for(int i = 0; i < k; i++){	
		k_dists[i] -= (k - i);
	}
	for(auto & estimation: estimations){
		auto target_norm_cvec = get<0>(estimation);
		SCIP_Real dist = 0;
		for(int i = 0; i < target_norm_cvec.size(); i++){
			dist += (norm_cvec[i] - target_norm_cvec[i]) * (norm_cvec[i] - target_norm_cvec[i]);
			if(dist > k_dists[k-1]){
				break;
			}
		}
		if(dist > k_dists[k-1]){
			continue;
		}
		auto query =  get<1>(estimation);
		for(int i = 0; i < k; i++){
			if(dist < k_dists[i]){
				SCIP_Real tmp_dist = k_dists[i];
				auto tmp_query  = k_queries[i];
				k_dists[i] = dist;
				k_queries[i] = query;
				dist = tmp_dist;
				query = tmp_query;
			}
		}
	}

	
	pair<SCIP_Real, SCIP_Real>  knn_bd = make_pair(0,0);
	if(knn_mode == 2){ // uniformly weighted
		for(int i = 0; i < k; i++){
			knn_bd.first += k_queries[i].first;
			knn_bd.second += k_queries[i].second;
		}
		knn_bd.first /= k;
		knn_bd.second /= k;
	}
	else if(knn_mode == 3){ // distance weighted	
		SCIP_Real sum_dist = 1e-5;
		for(int i = 0; i < k; i++){
			k_dists[i] = sqrt(k_dists[i]);
			sum_dist += k_dists[i];
		}
		for(int i = 0; i < k; i++){
			SCIP_Real wt = (k_dists[i] / sum_dist);
			knn_bd.first += k_queries[i].first * wt;
			knn_bd.second += k_queries[i].second * wt;
		}
	}
	return knn_bd;
}


void Estimator::add(const vector<SCIP_Real> & cvec, pair<SCIP_Real, SCIP_Real> mubd){
	vector<SCIP_Real> norm_cvec(cvec);
	SCIP_Real sum = 0;
	for(int i = 0; i < cvec.size(); i++){
		sum += cvec[i];
	}
	for(int i = 0; i < cvec.size(); i++){
		norm_cvec[i] /= sum;
	}
	estimations.push_back(make_tuple(norm_cvec, mubd, 0));
}