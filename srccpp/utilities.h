/**@file   utilities.h
 * @brief  Conic Knapsack utility
 * @author Liding Xu
 *---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#pragma once
#include <vector>
#include <list> 
#include "scip/scip.h"


using namespace std;
#ifndef tol
#define tol 1e-6
#endif

/* type of solution of conic knapsack problem */
enum SolType
{ 
	Unknown = 0,
	Infeasible  = 1,
	Feasible_Heur  = 2, 
	Feasible_Exact  = 3,  
	Optimal = 4, 
	Aborted = 5
};
//typedef enum SolType SOLTYPE_CKNAP;
typedef SolType SOLTYPE_CKNAP;

/*stat*/
struct stat
{
   SCIP_Real time_heur = 0, time_exact = 0;
   int col_heur = 0,    col_exact = 0;
   SCIP_Real shf_param = 1;
   SCIP_Real shf_log_sum_gap = 0;
   SCIP_Real cum_sol_time = 0;
   SCIP_Real cum_knn_time = 0;
};



/*algorithm configuration*/
struct conf
{
   SCIP_Bool is_misocp = false; // use misocp or milp relaxation? default: false
   SCIP_Bool is_bd_tight = true; // use bound tightenning? default: true
   SCIP_Bool is_heur = true; // use  heuristics pricing first? default: true
   SCIP_Bool is_parallelscplex = false; // enbale cplex's parallelism
   int knn_mode = 1; // the mode of knn regression for learning breakpoints, 1: no knn search/learning, 2: uniformly weighted knn, 3: distance weighted knn. default: 1
   int kneighbors = 1; // the number neighbors of knn regression (k). default: 1
   SCIP_Real point_ratio = 1; // the concentration ratio of breakpoints. default: 1 (no knn search)
};

/* conflict graph */
class conflict_graph{
	int numitems;
	vector<list<int>> conflict_list;
public:
	explicit conflict_graph(
		int numitems_, /** number of items */
		const vector<pair<int,int>>& items_diff /** items in same constraints */
	): numitems(numitems_){
		conflict_list = vector<list<int>> (numitems_, list<int>());
		for(auto p: items_diff){
			conflict_list[p.first].push_back(p.second);
			conflict_list[p.second].push_back(p.first);
		}
	}
   explicit  conflict_graph(
	){

   };
	const list<int> & get_diffs(int item) const{
		return conflict_list[item];
	}
};

inline size_t C2( size_t n ) { return n * ( n - 1 ) / 2; }

/* template class 1d implementation of upper triangle matrix*/
class upper_triangle: private vector<SCIP_Real>
{
    size_t N, P;
public:
    explicit  upper_triangle( size_t n  = 0) : N( n - 1), P( n )
    {
        this->resize( C2( P + 1 ) );
    }

    void set( size_t i, size_t j, SCIP_Real k )
    {
        assert(N >= j and j >= i );  this->at( C2( P - i ) + N - j ) = k;
    }

   void add( size_t i, size_t j, SCIP_Real k )
    {
        assert(N >= j and j >= i);  this->at( C2( P - i ) + N - j ) += k;
    }

    SCIP_Real get( size_t i, size_t j ) const
    {
        return (N >= j and j >= i ) ? this->at( C2( P - i ) + N - j ) : 0;
    }

   void reset()
    {
        fill(this->begin(), this->end(), 0);
    }
};

