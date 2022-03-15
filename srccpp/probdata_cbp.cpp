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

/**@file   probdata_cbp.cpp
 * @brief  Problem data for conic binpacking problem
 * @author Liding XU
 * /
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <string.h>
#include <utility>
#include <queue>        

#include "probdata_cbp.h"
#include "objscip/objscip.h"
#include "scip/struct_cons.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/var.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"

using namespace scip;
using namespace std;



SCIP_RETCODE vardataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata             /**< vardata to delete */
   )
{
   SCIPfreeBlockMemory(scip, vardata);
   return SCIP_OKAY;
};


extern SCIP_RETCODE vardataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   const PackVar * iter /**< iterator to its PackVar object */
)
{
   SCIP_CALL( SCIPallocBlockMemory(scip, vardata));
   (*vardata)->iter = iter;
   return SCIP_OKAY;
};


SCIP_DECL_VARTRANS(vardataTrans){
	assert(sourcedata != NULL);
	assert(sourcevar != NULL);
	SCIP_CALL(vardataCreate(scip, targetdata, sourcedata->iter));
	//SCIPdebugMessage("yes");
	return SCIP_OKAY;
};

SCIP_DECL_VARDELORIG(vardataDelOrig)
{
   SCIP_CALL( vardataDelete(scip, vardata) );

   return SCIP_OKAY;
};

SCIP_DECL_VARDELTRANS(vardataDelTrans)
{
   SCIP_CALL( vardataDelete(scip, vardata) );

   return SCIP_OKAY;
};

/** ProbDataCBP destructor */
ProbDataCBP::~ProbDataCBP()
{

}


/** release scip reference in probelme data*/
SCIP_RETCODE ProbDataCBP::releaseAll(
	SCIP*              scip                /**< SCIP data structure */
) {
	// release packing varibles
	int sizevar = 0, numvar =0 , numcons =0;
	for (auto it = p_vars.begin(); it != p_vars.end(); it++) {
		sizevar += sizeof(*it);
		numvar++;
		SCIP_CALL(SCIPreleaseVar(scip, &it->p_var));
	}

	// release set partition constraints
	for (int i = 0; i < sc_conss.size(); i++) {
		numcons++;
		SCIP_CALL(SCIPreleaseCons(scip, &sc_conss[i].sc_cons));
	}

	SCIPdebugMessage("freed %d %d %d\n ", int(sizevar), numvar, numcons);
	return SCIP_OKAY;
}

/** destructor of user problem data to free original user data (called when original problem is freed)
 *
 *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
 *  because all the work to delete the user problem data can be done in the destructor of the user problem
 *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
 *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
 *  longer needed.
 */
SCIP_RETCODE ProbDataCBP::scip_delorig(
   SCIP*              scip                /**< SCIP data structure */
   )
{
	SCIP_CALL(releaseAll(scip));
	return SCIP_OKAY;
}

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
SCIP_RETCODE ProbDataCBP::scip_trans(
   SCIP*              scip,               /**< SCIP data structure */
   ObjProbData**      objprobdata,        /**< pointer to store the transformed problem data object */
   SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
   )
{  /*lint --e{715}*/
   assert( objprobdata != NULL );
   assert( deleteobject != NULL );

	// collect algorithm configuration
	SCIP_CALL(SCIPgetBoolParam(scip,  "cbp/is_misocp", &algo_conf.is_misocp));
	SCIP_CALL(SCIPgetBoolParam(scip,  "cbp/is_bd_tight", &algo_conf.is_bd_tight));
	SCIP_CALL(SCIPgetBoolParam(scip,  "cbp/is_heur", &algo_conf.is_heur));
	SCIP_CALL(SCIPgetBoolParam(scip,  "cbp/is_parallelscplex", &algo_conf.is_parallelscplex));
	SCIP_CALL(SCIPgetBoolParam(scip,  "cbp/is_stablize", &algo_conf.is_stablize));
	SCIP_CALL(SCIPgetBoolParam(scip,  "cbp/is_adapt_points", &algo_conf.is_adapt_points));

   // create and cpature transformed path varibles 
   SCIPdebugMessage("start transform !!!!!!!!!!!\n");

    // allocate memory for target prob data

	ProbDataCBP * transprobdata = new ProbDataCBP(numitems, capacity, Dalpha,  mus, bs);

	transprobdata->item_matrix = item_matrix;
	transprobdata->belongs = belongs;
	transprobdata->currentnode = currentnode;
	transprobdata->global_lb = -SCIPinfinity(scip);
	transprobdata->stat_pr = stat_pr;
	transprobdata->algo_conf = algo_conf;
	SCIPdebugMessage("transformed data check!");
  // transform and cpature transformed set partition constraints
   for (int i = 0; i < sc_conss.size(); i++) {
	   transprobdata->sc_conss.push_back(SC_Cons(sc_conss[i]));
	   SCIP_CALL(SCIPtransformCons(scip,  sc_conss[i].sc_cons, &transprobdata->sc_conss[i].sc_cons));
   }


   transprobdata->p_vars = list<PackVar>();
	for (auto it = p_vars.begin(); it != p_vars.end(); it++) {
		transprobdata->p_vars.push_back(PackVar(*it));
		SCIP_CALL(SCIPtransformVar(scip, it->p_var, &(transprobdata->p_vars.back().p_var)));
	}

   SCIPdebugMessage("end transform \n");
   assert( transprobdata != NULL );
   *objprobdata = transprobdata;           
   
   *deleteobject = FALSE;

   return SCIP_OKAY;
}      

/** destructor of user problem data to free original user data (called when original problem is freed)
 *
 *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
 *  because all the work to delete the user problem data can be done in the destructor of the user problem
 *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
 *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
 *  longer needed.
 */
SCIP_RETCODE ProbDataCBP::scip_deltrans(
   SCIP*              scip                /**< SCIP data structure */
   )
{ 
   SCIP_CALL(releaseAll(scip));
   return SCIP_OKAY;
}

/** create constraints and initial columns */
SCIP_RETCODE ProbDataCBP::createConsInitialColumns(
	SCIP*                 scip               /**< SCIP data structure */
) {

   // create set cover constraint for each item
	for (int i = 0; i < numitems ; i++) {
		SCIP_CALL(addSC_Cons(scip, i));
	}
	
	SCIPdebugMessage("--set cover constraints added!\n");
   

   // create initial packing variables, a feasible solution is given by Best Fit for conic binpack problem
   queue<int> unpack_items;
   for(int i = 0; i < numitems; i++){
      unpack_items.push(i);
   }
   unpack_items.push(-1); // end flag
   int num_packs = 0; // number of packs in use
   SCIP_Real summu = 0;
   SCIP_Real sumb = 0;
   vector<int> assignment(numitems, -1);
   int best_fit = -1; 
   int skip = -2;
   SCIP_Real best_fit_cap_use = capacity + 1; // capacity used
   while(!unpack_items.empty()){
	  //SCIPdebugMessage("unpack_items %d\n", unpack_items.size());
      int item = unpack_items.front();
	  unpack_items.pop();
      if(item == -1){ // items in the queue are checked
         if(unpack_items.empty()){ // every item is assigned
		 	num_packs++;
            break;
         }
         else if ( best_fit != -1){ // there is an fitting item
            assignment[best_fit] = num_packs;
            summu += mus[best_fit];
            sumb += bs[best_fit];
            skip = best_fit;
         }
         else{ // need new pack
            num_packs++;
            summu = 0;
            sumb = 0;
            skip = -2;
         }
         unpack_items.push(-1);
         best_fit = -1;
         best_fit_cap_use =  capacity + 1;
         continue;
      }
      if(item == skip){
         continue;
      }
      SCIP_Real lhs = summu + mus[item] + Dalpha * sqrt(sumb + bs[item]); 
      if( lhs <= capacity){ // if fit into the bin
         SCIP_Real fit_cap_use = lhs - summu - Dalpha * sqrt(sumb); // use of the capacity
         if(fit_cap_use <  best_fit_cap_use){ // update the best fit
            best_fit = item;
            best_fit_cap_use = fit_cap_use;
         }
      } 
      // put the item in checked queue
      unpack_items.push(item);
   }
   
   SCIPdebugMessage("INTIAL FEAS!\n");
   // assign items to binpacks
   vector<vector<int>> binpacks(num_packs);
   for(int i = 0; i < numitems; i++){
      binpacks[assignment[i]].push_back(i);
   }


	for (int bin = 0; bin < num_packs; bin++) {
		SCIP_CALL(addPackVar(scip, binpacks[bin], FALSE));
	}

   // set the objective integer
   SCIP_CALL( SCIPsetObjIntegral(scip) );
   SCIPdebugMessage("--initial columns created!\n");
   return SCIP_OKAY;
}

/** add set cover constraint */
SCIP_RETCODE ProbDataCBP::addSC_Cons(
	SCIP * 	scip, /**< SCIP data structure */
	const int item_ /**< item index */
) {
	// create vertex constraint
	SCIP_CONS* cons;
	SCIP_CALL(SCIPcreateConsSetcover(
		scip, /**< SCIP data structure */
		&cons, /**< 	pointer to hold the created constraint */
		"set cover constraint", /**< 	name of constraint */
		0, /**< number of nonzeros in the constraint */
		NULL, /**< 	array with variables of constraint entries */
		TRUE,                   /* initial */
		FALSE,                  /* separate */
		TRUE,                   /* enforce */
		TRUE,                   /* check */
		TRUE,                   /* propagate */
		FALSE,                  /* local */
		TRUE,                   /* modifiable */
		FALSE,                  /* dynamic */
		FALSE,                  /* removable */
		FALSE));               /* stickingatnode */

	SCIP_CALL(SCIPaddCons(scip, cons));
	SCIP_CALL(SCIPcaptureCons(scip, cons));
	sc_conss.push_back(SC_Cons(item_ ,cons));
	SCIP_CALL(SCIPreleaseCons(scip, &cons));
	return SCIP_OKAY;
}

/** add binpack variable and its data in scip and problem data */
SCIP_RETCODE ProbDataCBP::addPackVar(
	SCIP * 	scip, /**< SCIP data structure */
	const vector<int>& item_array_, /**< items in the binpack */
	const SCIP_Bool is_pricing /**< indicate pricing variable or not*/
) {
	SCIP_VAR* p_var;
	if (is_pricing) {
		SCIP_CALL(SCIPcreateVar(
			scip, /**<	SCIP data structure*/
			&p_var, /**< 	pointer to variable object*/
			NULL, /**< name of variable, or NULL for automatic name creation*/
			0.0, /**<	lower bound of variable*/
			1.0, /**< 	upper bound of variable */
			1, /**<	objective function value */
			SCIP_VARTYPE_BINARY, /**< type of variable */
			FALSE, /**<	should var's column be present in the initial root LP?*/
			TRUE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
			NULL, NULL, NULL, NULL, NULL
		));
		SCIP_CALL(SCIPaddPricedVar(scip, p_var, 1.0));
	}
	else {
		SCIP_CALL(SCIPcreateVar(
			scip, /**<	SCIP data structure*/
			&p_var, /**< 	pointer to variable object*/
			NULL, /**< name of variable, or NULL for automatic name creation*/
			0.0, /**<	lower bound of variable*/
			1.0, /**< 	upper bound of variable */
			1, /**<	objective function value */
			SCIP_VARTYPE_BINARY, /**< type of variable */
			TRUE, /**<	should var's column be present in the initial root LP?*/
			TRUE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
			NULL, NULL, NULL, NULL, NULL
		));
		SCIP_CALL(SCIPaddVar(scip, p_var));
	}
	SCIP_CALL(SCIPchgVarUbLazy(scip, p_var, 1.0));
	// capture the variable
	SCIP_CALL(SCIPcaptureVar(scip, p_var));
	p_vars.push_back(PackVar(item_array_, p_var));

	SCIP_VARDATA * vardata;
	SCIP_CALL(vardataCreate(scip, &vardata, &p_vars.back() ) );
	/* set callback functions */
  	SCIPvarSetData(p_var, vardata);
	if(!is_pricing){
		SCIPvarSetDelorigData(p_var, vardataDelOrig);
		SCIPvarSetDeltransData(p_var, vardataDelTrans);
		SCIPvarSetTransData(p_var, vardataTrans);
	}
	else{
		SCIPvarSetDeltransData(p_var, vardataDelTrans);
	}
   // add the packing variable into the set cover constraints
	for (int item: item_array_) {
      SCIP_CALL(SCIPaddCoefSetppc(scip, sc_conss[item].sc_cons, p_var));
	}

	SCIP_CALL(SCIPreleaseVar(scip, &p_var));

	return SCIP_OKAY;
}


/** return the number of pack variables */
int ProbDataCBP::getNumPackVars() {
	return int(p_vars.size());
}



/**@} */
