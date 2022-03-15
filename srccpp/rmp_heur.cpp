 /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
 /*                                                                           */
 /*                  This file is part of the program and library             */
 /*         SCIP --- Solving Constraint Integer Programs                      */
 /*                                                                           */
 /*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
 /*                            fuer Informationstechnik Berlin                */
 /*                                                                           */
 /*  SCIP is distributed under the terms of the ZIB Academic License.         */
 /*                                                                           */
 /*  You should have received a copy of the ZIB Academic License.             */
 /*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
 /*                                                                           */
 /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
 
 /**@file   rmp_heur.cpp
  * @brief  restricted master problem heursitic solves a sub MIP
  * @author Liding XU
  */
 
 /*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
 
 #define SCIP_DEBUG
 #include <iostream>
 #include <cassert>
 
 #include "objscip/objscip.h"
 #include "rmp_heur.h"
 #include "probdata_cbp.h"

 
 
 
 /** destructor of primal heuristic to free user data (called when SCIP is exiting) */
 SCIP_DECL_HEURFREE(HeurRMP::scip_free)
 {  /*lint --e{715}*/
    return SCIP_OKAY;
 }
 
 
 /** deinitialization method of primal heuristic (called before transformed problem is freed) */
 SCIP_DECL_HEUREXIT(HeurRMP::scip_exit)
 {  /*lint --e{715}*/
    SCIP_CALL( SCIPfreeSol(scip, &sol) );
    return SCIP_OKAY;
 }
 
 /** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
  *
  *  This method is called when the presolving was finished and the branch and bound process is about to begin.
  *  The primal heuristic may use this call to initialize its branch and bound specific data.
  *
  */
 SCIP_DECL_HEURINITSOL(HeurRMP::scip_initsol)
 {  /*lint --e{715}*/
    return SCIP_OKAY;
 }
 
 /** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
  *
  *  This method is called before the branch and bound process is freed.
  *  The primal heuristic should use this call to clean up its branch and bound data.
  */
 SCIP_DECL_HEUREXITSOL(HeurRMP::scip_exitsol)
 {  /*lint --e{715}*/
    return SCIP_OKAY;
 }
 

 /** perform solution lifting */
SCIP_RETCODE HeurRMP::performSolLifting(
    SCIP*                 scip,               /**< SCIP main data structure */
    ProbDataCBP * probdata,             /** < problem data */
    SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
    )
{
    SCIP_CALL(SCIPclearSol(scip, sol));
    //SCIPdebugMessage("STEP 1\n");
    for(auto p = probdata->p_vars.begin(); p!= probdata->p_vars.end(); p++){
        SCIP_CALL( SCIPsetSolVal(scip, sol, p->p_var, 0) );
    }
    auto & p_last = probdata->p_vars.back();
    vector<bool> packed(probdata->numitems, false);
    bool all_packed = false;
    for(int item: p_last.item_array){
        packed[item] = true;
    }
    /* store new solution value */
    SCIP_CALL( SCIPsetSolVal(scip, sol, p_last.p_var, 1) );
    //SCIPdebugMessage("STEP 2\n");
    int num_bins = 1;
    while(!all_packed){
        // search for the bin with the maximal pack
        int num_pack_max = 0;
        auto p_max = probdata->p_vars.end();
        for(auto p = probdata->p_vars.begin(); p!= probdata->p_vars.end(); p++){
            // may round up
            if(!SCIPvarMayRoundUp(p->p_var)){
                continue;
            }
            // compute the number of items (with conflics) packable
            int num_pack = 0;
            for(int item: p->item_array){
                num_pack += (!packed[item] ? 1 : 0);
            }
            if(num_pack > num_pack_max ){
                num_pack_max = num_pack;
                p_max = p;
            }
        }
        //SCIPdebugMessage("STEP 3 %d\n", num_pack_max);
        if(num_pack_max == 0){
            return SCIP_OKAY;
        }
        SCIP_CALL( SCIPsetSolVal(scip, sol, p_max->p_var, 1) );
        // try a new bin with the maximal pack
        num_bins++;
        // pack
        for(int item: p_max->item_array){
            packed[item] = TRUE;
        }
        // check all packed
        all_packed = true;
        for(int item = 0; item < probdata->numitems; item++){
            if(!packed[item]){
                all_packed = false;
                break;
            }
        }
    }
    
    if(num_bins <= SCIPgetPrimalbound(scip))
    {
        SCIP_Bool stored;
        SCIP_Bool checklprows;
        if( SCIPallColsInLP(scip) )
        {
            /* check solution for feasibility, and add it to solution store if possible
            * integrality need not be checked, because all fractional
            * variables were already moved in feasible direction to the next integer
            *
            * feasibility of LP rows must be checked again at the presence of
            * unroundable, implicit integer variables with fractional LP solution
            * value
            */
            SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, checklprows, &stored) );
        }
        else
        {
            /* if there are variables which are not present in the LP, e.g., for 
            * column generation, we need to check their bounds
            */
            SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, TRUE, FALSE, checklprows, &stored) );
        }


        
        if( stored )
        {
        //SCIPdebugMsg(scip, "found feasible rounded solution:\n");
    #ifdef SCIP_DEBUG__
            SCIPdebugMsg(scip, "found feasible rounded solution:\n");
            SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
    #endif
            *result = SCIP_FOUNDSOL;
        }

    }

    return SCIP_OKAY;
}

 /** execution method of primal heuristic
  *
  *  Searches for feasible primal solutions. The method is called in the node processing loop.
  *
  *  possible return values for *result:
  *  - SCIP_FOUNDSOL   : at least one feasible primal solution was found
  *  - SCIP_DIDNOTFIND : the heuristic searched, but did not find a feasible solution
  *  - SCIP_DIDNOTRUN  : the heuristic was skipped
  *  - SCIP_DELAYED    : the heuristic was skipped, but should be called again as soon as possible, disregarding
  *                      its frequency
  */
 SCIP_DECL_HEUREXEC(HeurRMP::scip_exec)
 {  /*lint --e{715}*/
    assert(scip != NULL);
    assert(!SCIPinDive(scip));
    //assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
    assert(result != NULL);
    assert(SCIPhasCurrentNodeLP(scip));
 
    *result = SCIP_DIDNOTRUN;
    ProbDataCBP * probdata = NULL;
	probdata = dynamic_cast<ProbDataCBP *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);
    //SCIPdebugMessage("in heur RMP\n");

    /* only call heuristic, if an optimal LP solution is at hand or if relaxation solution is available */
    //if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL && ! SCIPisRelaxSolValid(scip) )
    //   return SCIP_OKAY;
 

    if( heurtiming != SCIP_HEURTIMING_DURINGPRICINGLOOP &&  heurtiming != SCIP_HEURTIMING_DURINGLPLOOP ){
       return SCIP_OKAY;}
 
    *result = SCIP_DIDNOTFIND;
 
    /* try to round LP solution */
    //SCIP_CALL( performLPSimpleRounding(scip, probdata, heurtiming, result) );
 
    /* try to round relaxation solution */
    //SCIP_CALL( performRelaxSimpleRounding(scip, probdata, result) );
    //SCIPdebugMessage("perfrom heur RMP\n");

    SCIP_CALL(performSolLifting(scip, probdata, result));
 
    return SCIP_OKAY;
    //SCIP_CALL(SCIPinterruptSolve(scip));
    //return SCIP_OKAY;
 }
 
 /** clone method which will be used to copy a objective plugin */
 SCIP_DECL_HEURCLONE(scip::ObjCloneable* HeurRMP::clone) /*lint !e665*/
 {
    //SCIPdebugMessage("cloned\n");
    return new HeurRMP(scip);
 }

/** initialization method of primal heuristic (called after problem was transformed) */
 SCIP_DECL_HEURINIT(HeurRMP::scip_init) /*lint --e{715}*/
 {  /*lint --e{715}*/

    //assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

    /* create heuristic data */
    SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
    lastlp = -1;
    nroundablevars = -1;
 
    return SCIP_OKAY;
 }
 