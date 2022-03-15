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

/**@file   cons_samediff.cpp
 * @brief  Constraint handler stores the local branching decision data
 * @author Liding Xu
 * This constraint handler is used to store the branching decision of the \ref BINPACKING_BRANCHING "Ryan/Foster branching rule"
 * which is implemented in \ref branch_ryanfoster.c.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <algorithm>
#include <utility> 
#include <vector>  
#include <assert.h>
#include <string.h>

#include "cons_samediff.h"
#include "probdata_cbp.h"
#include "objscip/objscip.h"


using namespace scip;
using namespace std;


struct SCIP_ConsData
{
   int   itemid1;            /**< item id one */
   int   itemid2;            /**< item id two */
	CONSTYPE cons_type; /* type of constraint: differ or same */
	SCIP_NODE* node; /**< the node in the B&B-tree at which the cons is sticking */
	bool propagated; /**< is constraint already propagated? */
	int npropagatedvars; /**< number of variables that existed, the last time, the related node was
		*   propagated, used to determine whether the constraint should be
		*   repropagated*/
	int npropagations;  /**< stores the number propagations runs of this constraint */
};


/**@name Local methods
 *
 * @{
 */

/** create constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
   int                   itemid1,            /**< item id one */
   int                   itemid2,            /**< item id two */
   CONSTYPE              cons_type,               /**< stores whether the items have to be in the SAME or DIFFER packing */
   SCIP_NODE*            node                /**< the node in the B&B-tree at which the cons is sticking */
   )
{
   assert( scip != NULL );
   assert( consdata != NULL );
   assert( itemid1 >= 0 );
   assert( itemid2 >= 0 );
   assert( cons_type == DIFFER || cons_type == SAME );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->itemid1 = itemid1;
   (*consdata)->itemid2 = itemid2;
   (*consdata)->cons_type = cons_type;
   (*consdata)->npropagatedvars = 0;
   (*consdata)->npropagations = 0;
   (*consdata)->propagated = FALSE;
   (*consdata)->node = node;

   return SCIP_OKAY;
}


/** fixes a variable to zero if the corresponding packings are not valid for this constraint/node (due to branching) */
static
SCIP_RETCODE checkVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   PackVar *             p_var,                /**< variables to check  */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_Bool fixed;
   SCIP_Bool infeasible;


   assert(scip != NULL);
   assert(consdata != NULL);
   assert(p_var != NULL);
   assert(nfixedvars != NULL);
   assert(cutoff != NULL);

   SCIP_VAR * var = p_var->p_var;

   /* if variables is locally fixed to zero continue */
   if( SCIPvarGetUbLocal(var) < 0.5 )
      return SCIP_OKAY;

   /* check if the packing which corresponds to the variable is feasible for this constraint */

   bool existid1 = binary_search(p_var->item_array.begin(), p_var->item_array.end(), consdata->itemid1);
   bool existid2 = binary_search(p_var->item_array.begin(), p_var->item_array.end(), consdata->itemid2);
   CONSTYPE type = consdata->cons_type;

   if( (type == SAME && existid1 != existid2) || (type == DIFFER && existid1 && existid2) )
   {
      SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );

      if( infeasible )
      {
         assert( SCIPvarGetLbLocal(var) > 0.5 );
         //SCIPdebugMsg(scip, "-> cutoff\n");
         (*cutoff) = TRUE;
      }
      else
      {
         assert(fixed);
         (*nfixedvars)++;
      }
   }

   return SCIP_OKAY;
}

/** fixes variables to zero if the corresponding packings are not valid for this sonstraint/node (due to branching) */
static
SCIP_RETCODE consdataFixVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   list<PackVar>&        p_vars,              /**< list of pack variables */
   SCIP_RESULT*          result              /**< pointer to store the result of the fixing */
   )
{
   int nfixedvars;
   SCIP_Bool cutoff;

   nfixedvars = 0;
   cutoff = FALSE;

   //SCIPdebugMsg(scip, "check variables %d to %d\n", consdata->npropagatedvars, nvars);

	for (auto it = p_vars.begin(); it != p_vars.end() && !cutoff ; it++) {
		SCIP_CALL(checkVariable(scip, consdata, &(*it), &nfixedvars, &cutoff));
   }

   //SCIPdebugMsg(scip, "fixed %d variables locally\n", nfixedvars);

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** check if all variables are valid for the given consdata */
#ifndef NDEBUG
static
SCIP_Bool consdataCheck(
   SCIP*                 scip,               /**< SCIP data structure */
   ProbDataCBP *        probdata,           /**< problem data */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool             beforeprop          /**< is this check performed before propagation? */
   )
{
   int nvars = (beforeprop ? consdata->npropagatedvars : probdata->getNumPackVars());
   assert(nvars <= probdata->getNumPackVars());

   for(auto p_var: probdata->p_vars)
   {
      if(nvars == 0){
         break;
      }
      SCIP_VAR* var = p_var.p_var;

      /* if variables is locally fixed to zero continue */
      if( SCIPvarGetUbLocal(var) < 0.5 )
         continue;

      /* check if the packing which corresponds to the variable is feasible for this constraint */

      bool existid1 = binary_search(p_var.item_array.begin(), p_var.item_array.end(), consdata->itemid1);
      bool existid2 = binary_search(p_var.item_array.begin(), p_var.item_array.end(), consdata->itemid2);
      CONSTYPE type = consdata->cons_type;

      if( (type == SAME && existid1 != existid2) || (type == DIFFER && existid1 && existid2) )
      {
         //SCIPdebug( SCIPvardataPrint(scip, vardata, NULL) );
         //SCIPdebug( consdataPrint(scip, consdata, NULL) );
         //SCIPdebug( SCIPprintVar(scip, var, NULL) );
         return FALSE;
      }
      nvars--;
   }

   return TRUE;
}
#endif

/** frees samediff constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to the constraint data */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** frees specific constraint data
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 */
SCIP_DECL_CONSDELETE(ConshdlrSameDifferent::scip_delete)
{  /*lint --e{715}*/
   assert(consdata != NULL);
   assert(*consdata != NULL);
   /* free samediff constraint */
   SCIP_CALL( consdataFree(scip, consdata) );
   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(ConshdlrSameDifferent::scip_trans) 
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata,
         sourcedata->itemid1, sourcedata->itemid2, sourcedata->cons_type, sourcedata->node) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The propagation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_REDUCEDDOM : at least one domain reduction was found
 *  - SCIP_DIDNOTFIND : the propagator searched, but did not find any domain reductions
 *  - SCIP_DIDNOTRUN  : the propagator was skipped
 *  - SCIP_DELAYED    : the propagator was skipped, but should be called again
 */
SCIP_DECL_CONSPROP(ConshdlrSameDifferent::scip_prop) {
     /*lint --e{715}*/
   assert(scip != NULL);
   assert(result != NULL);

  // SCIPdebugMsg(scip, "propagation constraints of constraint handler <"CONSHDLR_NAME">\n");
	ProbDataCBP * probdata = NULL;
	probdata = dynamic_cast<ProbDataCBP *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);

   int nvars = probdata->getNumPackVars();
   *result = SCIP_DIDNOTFIND;

   for (int c_ind = 0; c_ind < nconss; c_ind++)
	{
		SCIP_CONSDATA* consdata = NULL;
		consdata = SCIPconsGetData(conss[c_ind]);
		assert(consdata != NULL);
      /* check if all previously generated variables are valid for this constraint */
      //assert( consdataCheck(scip, probdata, consdata, TRUE) );
      #ifndef NDEBUG
            {
               /* check if there are no equal consdatas */


               for(int i = c_ind +1; i < nconss; ++i )
               {
                  SCIP_CONSDATA* consdata2 = SCIPconsGetData(conss[i]);
                  assert( !(consdata->itemid1 == consdata2->itemid1
                        && consdata->itemid2 == consdata2->itemid2
                        && consdata->cons_type == consdata2->cons_type) );
                  assert( !(consdata->itemid1 == consdata2->itemid2
                        && consdata->itemid2 == consdata2->itemid1
                        && consdata->cons_type == consdata2->cons_type) );
               }
            }
      #endif
      if(!consdata->propagated)
      {
         //SCIPdebugMsg(scip, "propagate constraint <%s> ", SCIPconsGetName(conss[c]));
         //SCIPdebug( consdataPrint(scip, consdata, NULL) );
         SCIP_CALL( consdataFixVariables(scip, consdata, probdata->p_vars, result) );
         consdata->npropagations++;

         if( *result != SCIP_CUTOFF )
         {
            consdata->propagated = TRUE;
            consdata->npropagatedvars = nvars;
         }
         else{
            break;
         }
      }	
      /* check if constraint is completely propagated */
      //assert( consdataCheck(scip, probdata, consdata, FALSE) );
	}
	return SCIP_OKAY;
}

/** constraint activation notification method of constraint handler
 *
 *  @see SCIP_DECL_CONSACTIVE(x) in @ref type_cons.h
 */
SCIP_DECL_CONSACTIVE(ConshdlrSameDifferent::scip_active) 
{  /*lint --e{715}*/
	SCIP_CONSDATA* consdata = NULL;

	//SCIPdebugMessage("active debug\n");
	assert(scip != NULL);
	assert(cons != NULL);

	ProbDataCBP * probdata = NULL;
	probdata = dynamic_cast<ProbDataCBP *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   int nvars = probdata->getNumPackVars();
   assert(consdata->npropagatedvars <= nvars);

   //SCIPdebugMsg(scip, "activate constraint <%s> at node <%"SCIP_LONGINT_FORMAT"> in depth <%d>: ",
   //   SCIPconsGetName(cons), SCIPnodeGetNumber(consdata->node), SCIPnodeGetDepth(consdata->node));
   //SCIPdebug( consdataPrint(scip, consdata, NULL) );

   if( consdata->npropagatedvars != nvars )
   {
      //SCIPdebugMsg(scip, "-> mark constraint to be repropagated\n");
      consdata->propagated = FALSE;
      SCIP_CALL( SCIPrepropagateNode(scip, consdata->node) );
   }

   return SCIP_OKAY;
}

/** constraint deactivation notification method of constraint handler */
SCIP_DECL_CONSDEACTIVE(ConshdlrSameDifferent::scip_deactive) 
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->propagated || SCIPgetNChildren(scip) == 0);

	ProbDataCBP * probdata = NULL;
	probdata = dynamic_cast<ProbDataCBP *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);

   //SCIPdebugMsg(scip, "deactivate constraint <%s> at node <%"SCIP_LONGINT_FORMAT"> in depth <%d>: ",
   //   SCIPconsGetName(cons), SCIPnodeGetNumber(consdata->node), SCIPnodeGetDepth(consdata->node));
   //SCIPdebug( consdataPrint(scip, consdata, NULL) );

   /* set the number of propagated variables to current number of variables is SCIP */
   consdata->npropagatedvars = probdata->getNumPackVars();

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(ConshdlrSameDifferent::scip_check) {
	assert(scip != NULL);
	assert(conshdlr != NULL);
	assert(result != NULL);

	/* do nothing */
	*result = SCIP_FEASIBLE;

	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(ConshdlrSameDifferent::scip_enfolp) {
	assert(scip != NULL);
	assert(conshdlr != NULL);
	assert(result != NULL);

	/* do nothing */
	*result = SCIP_FEASIBLE;

	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(ConshdlrSameDifferent::scip_enfops) {
	assert(scip != NULL);
	assert(conshdlr != NULL);
	assert(result != NULL);

	/* do nothing */
	*result = SCIP_FEASIBLE;

	return SCIP_OKAY;
}
/** variable rounding lock method of constraint handler */
SCIP_DECL_CONSLOCK(ConshdlrSameDifferent::scip_lock) {
	assert(scip != NULL);
	assert(conshdlr != NULL);
	assert(cons != NULL);

	SCIPdebugMessage("Locking method for store graph constraint: <%s>.\n", SCIPconsGetName(cons));

	return SCIP_OKAY;
}

/**@} */



/** creates and captures a samediff constraint */
SCIP_RETCODE SCIPcreateConsSamediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   itemid1,            /**< item id one */
   int                   itemid2,            /**< item id two */
   CONSTYPE              cons_type,               /**< stores whether the items have to be in the SAME or DIFFER packing */
   SCIP_NODE*            node,               /**< the node in the B&B-tree at which the cons is sticking */
   SCIP_Bool             local               /**< is constraint only valid locally? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the samediff constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "cons_samediff");
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("samediff constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create the constraint specific data */
   SCIP_CALL( consdataCreate(scip, &consdata, itemid1, itemid2, cons_type, node) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
         local, FALSE, FALSE, FALSE, TRUE) );

   //SCIPdebugMsg(scip, "created constraint: ");
   //SCIPdebug( consdataPrint(scip, consdata, NULL) );

   return SCIP_OKAY;
}




/** return constraint type SAME or DIFFER */
CONSTYPE SCIPgetTypeSamediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->cons_type;
}

/** return the items pairs in the same constraints */
vector<pair<int,int>> getItemsInSame(
	SCIP*                 scip              /**< SCIP data structure */
) {
	SCIP_CONSHDLR* conshdlr = NULL;
	SCIP_CONS** conss = NULL;
	SCIP_CONS* cons = NULL;
	int nconss;
	/* find the forbidden edge constraint handler */
	conshdlr = SCIPfindConshdlr(scip, "cons_samediff");
	assert(scip != NULL);
	assert(conshdlr != NULL);

	/* collect all branching decision constraints */
	conss = SCIPconshdlrGetConss(conshdlr);
	nconss = SCIPconshdlrGetNConss(conshdlr);
	vector<pair<int, int>> items_same;
	for (int i = 0; i < nconss; i++) {
		if(!SCIPconsIsActive(conss[i])){
			continue;
		}
		SCIP_CONSDATA* consdata = NULL;
		consdata = SCIPconsGetData(conss[i]);
		assert(consdata != NULL);
		if(consdata->cons_type  == SAME ){
			items_same.push_back(make_pair(consdata->itemid1, consdata->itemid2));
		}
	}
   return items_same;
}

/** return the items pairs in the differ constraints */
vector<pair<int,int>>  getItemsInDiffer(
	SCIP*                 scip              /**< SCIP data structure */
) {
	SCIP_CONSHDLR* conshdlr = NULL;
	SCIP_CONS** conss = NULL;
	SCIP_CONS* cons = NULL;
	int nconss;
	/* find the forbidden edge constraint handler */
	conshdlr = SCIPfindConshdlr(scip, "cons_samediff");
	assert(scip != NULL);
	assert(conshdlr != NULL);

	/* collect all branching decision constraints */
	conss = SCIPconshdlrGetConss(conshdlr);
	nconss = SCIPconshdlrGetNConss(conshdlr);
	vector<pair<int, int>> items_differ;
	for (int i = 0; i < nconss; i++) {
		if(!SCIPconsIsActive(conss[i])){
			continue;
		}
		SCIP_CONSDATA* consdata = NULL;
		consdata = SCIPconsGetData(conss[i]);
		assert(consdata != NULL);
		if(consdata->cons_type  == DIFFER ){
			items_differ.push_back(make_pair(consdata->itemid1, consdata->itemid2));
		}
	}
   return items_differ;
}
	


/**@} */
