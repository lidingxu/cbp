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

/**@file   probdata_cbp.h
 * @brief  Problem data for conic binpacking problem
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_CONICBINPACKING__
#define __SCIP_PROBDATA_CONICBINPACKING__

#include "objscip/objscip.h"
#include "scip/cons_linear.h"
#include "estimator.h"
#include <map>
#include <list>
#include <vector>
#include <utility>

using namespace scip;
using namespace std;

/*stat*/
struct stat
{
   SCIP_Real time_heur = 0, time_exact = 0;
   int col_heur = 0,    col_exact = 0;
   SCIP_Real shf_param = 1;
   SCIP_Real shf_log_sum_gap = 0;
};


struct Pr_param
{
  SCIP_Real lb;
  SCIP_Real ub;
  int piece_sample_size;
};


/*algorithm configuration*/
struct conf
{
  SCIP_Bool is_misocp = false; // use misocp or milp relaxation? default: false
  SCIP_Bool is_bd_tight = true; // use bound tightenning? default: true
  SCIP_Bool is_heur = true; // use  heuristics pricing first? default: true
  SCIP_Bool is_cut_pool = false; // use cut pool
  SCIP_Bool is_stablize = false; // use stabilization for hybrid pricing
  SCIP_Bool is_parallelscplex = false; // enbale cplex's parallelism
  SCIP_Bool is_adapt_points = true; // enbale cplex's adaptive break points
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

/** Packing variable class */
class PackVar{
public:
	/** default constructor */
	PackVar(
		const vector<int>& item_array_, /**< items in the binpack */
		SCIP_VAR * p_var_ /**< pointer to the pack SCIP_VAR */
	): item_array(item_array_), p_var(p_var_) {};

   const vector <int> item_array; /**< items in the binpack, assumed to be ascending order */
	SCIP_VAR * p_var; /**< pointer to the packing var */
};

/** Variable data which is attached to packing variables.
 */
struct SCIP_VarData
{
   const PackVar * iter; // iterator to its PackVar object
};


extern SCIP_RETCODE vardataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata             /**< vardata to delete */
);

extern SCIP_RETCODE vardataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   const PackVar * iter  /**< iterator to its PackVar object */
);

/** Set covering constraint over item*/
class SC_Cons{
public:
	/** default constructor */
	SC_Cons(
		const int item_ind_, /**< item index */
      SCIP_CONS* sc_cons_  /**< pointer to hold the created constraint */
	): item_ind(item_ind_), sc_cons(sc_cons_){};

	SCIP_CONS* sc_cons;  /**< pointer to hold the created constraint */
   int item_ind; /**< item index */
};


/** SCIP user problem data for Conic BinPakcing */
class ProbDataCBP : public ObjProbData
{
public:

   /** default constructor */
   ProbDataCBP(
		const int numitems_,  /**< the number of items */
		const SCIP_Real capacity_, /**< capacity */
		const SCIP_Real Dalpha_, /**< Dalpha */
      const vector<SCIP_Real> & mus_, /**< mus */
		const vector<SCIP_Real> & bs_ /**< bs */
      )
      : numitems(numitems_), capacity(capacity_), Dalpha(Dalpha_), mus(mus_), bs(bs_){}

   /**< destructor */
   ~ProbDataCBP();

   /**< utility functions*/

   /** release scip reference in probelme data*/
   SCIP_RETCODE releaseAll(
	   SCIP*              scip                /**< SCIP data structure */
   );

   /** create constraints and initial columns */
   SCIP_RETCODE createConsInitialColumns(
	   SCIP*                 scip               /**< SCIP data structure */
   );

   /** add set partition constraint*/
   SCIP_RETCODE addSC_Cons(
	   SCIP * 	scip, /**< SCIP data structure */
	   const int item_ind_  /**< item index */
   );


   /** add packing variable and its data in scip and problem data*/
   SCIP_RETCODE addPackVar(
	   SCIP * 	scip, /**< SCIP data structure */
	   const vector<int>& item_array_, /**< items in the binpack */
	   const SCIP_Bool is_pricing /**< indicate pricing variable or not*/
   );

   /** return the number of packing variables */
   int getNumPackVars();

   /** destructor of user problem data to free original user data (called when original problem is freed)
    *
    *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual SCIP_RETCODE scip_delorig(
      SCIP*              scip                /**< SCIP data structure */
      );
   
   /** destructor of user problem data to free transformed user data (called when transformed problem is freed)
    *
    *  If the "*deleteobject" flag in the scip_trans() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "*deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual SCIP_RETCODE scip_deltrans(
      SCIP*              scip                /**< SCIP data structure */
      );

   /** creates user data of transformed problem by transforming the original user problem data
    *  (called after problem was transformed)
    *
    *  The user has two possibilities to implement this method:
    *   1. Return the pointer to the original problem data object (this) as pointer to the transformed problem data
    *      object. The user may modify some internal attributes, but he has to make sure, that these modifications are
    *      reversed in the scip_deltrans() method, such that the original problem data is restored. In this case,
    *      he should set *deleteobject to FALSE, because the problem data must not be destructed by SCIP after the
    *      solving process is terminated.
    *   2. Call the copy constructor of the problem data object and return the created copy as transformed problem
    *      data object. In this case, he probably wants to set *deleteobject to TRUE, thus letting SCIP call the
    *      destructor of the object if the transformed problem data is no longer needed.
    */
   virtual SCIP_RETCODE scip_trans(
      SCIP*              scip,               /**< SCIP data structure */
      ObjProbData**      objprobdata,        /**< pointer to store the transformed problem data object */
      SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
      );

   list<PackVar> p_vars; /**< list of packing variables */
   vector<SC_Cons> sc_conss; /**< set covering constraints indexed by items  */
   upper_triangle item_matrix; /**< item matrix is used for Ryan/Foster branching rule  */

   vector<pair<int, int>> items_same; // items in same constraints
	vector<pair<int, int>> items_differ; //items in differ constraints
   vector<int> belongs; // items corresponding to the to new items
   vector<vector<int>> new_items; // new items that contains the original items 
   vector<SCIP_Real> new_mus; // new items' mus
   vector<SCIP_Real> new_bs; // new items' bs
   vector<pair<int, int>> new_differ; // new different constraints 
   vector<SCIP_Real> stable_center; // stablization center 
   int num_new_items; // the number of new items
   int piece_size; // piece wise size
   estimator init_estimator; // initial quadratic estimator
   int piece_sample_size; //  sampled piece size
   SCIP_Real cbp_time; // pricing cbp time limit
   SCIP_Real global_lb; // global lower bound
   Pr_param pr_param; // pricing algorithm paramters */
   conflict_graph conflict; // conflict graph
   conf algo_conf; // algorithm configuration
   stat stat_pr; // statistics pricing


   long long currentnode; // current node, except for diving mode
   int numitems;
   SCIP_Real capacity;
   SCIP_Real Dalpha;
   const vector<SCIP_Real> mus;
   const vector<SCIP_Real> bs;
};/*lint !e1712*/


#endif
