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
 
 /**@file   rmp_heur.h
  * @brief  restricted master problem heursitic solves a sub MIP
  * @author Liding XU
  */
 
 /*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
 
 #ifndef __RMPHEUR_H__
 #define __RMPHEUR_H__
 
#include "objscip/objscip.h"
#include"probdata_cbp.h"
 
 /** C++ restriced master heuristics for covering program */
 class HeurRMP : public scip::ObjHeur
 {

 public:
   long long lastlp = -1;
   int nroundablevars = -1;
   SCIP_SOL* sol = NULL;
    /** default constructor */
    HeurRMP(
       SCIP* scip
       )
       : ObjHeur(scip, "rmp", "restriced master heuristics", 'R', 10, 1, 0, -1,
      SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_DURINGPRICINGLOOP , FALSE)
    {
    }
 
    /** destructor */
    virtual ~HeurRMP()
    {
    }
 
 
    /** destructor of primal heuristic to free user data (called when SCIP is exiting) */
    virtual SCIP_DECL_HEURFREE(scip_free);
 
    /** initialization method of primal heuristic (called after problem was transformed) */
    virtual SCIP_DECL_HEURINIT(scip_init);
 
    /** deinitialization method of primal heuristic (called before transformed problem is freed) */
    virtual SCIP_DECL_HEUREXIT(scip_exit);
 
    /** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
     *
     *  This method is called when the presolving was finished and the branch and bound process is about to begin.
     *  The primal heuristic may use this call to initialize its branch and bound specific data.
     *
     */
    virtual SCIP_DECL_HEURINITSOL(scip_initsol);
 
    /** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
     *
     *  This method is called before the branch and bound process is freed.
     *  The primal heuristic should use this call to clean up its branch and bound data.
     */
    virtual SCIP_DECL_HEUREXITSOL(scip_exitsol);
 
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
    virtual SCIP_DECL_HEUREXEC(scip_exec);
 
    /** clone method which will be used to copy a objective plugin */
    virtual SCIP_DECL_HEURCLONE(ObjCloneable* clone); /*lint !e665*/
 
    /** returns whether the objective plugin is copyable */
    virtual SCIP_DECL_HEURISCLONEABLE(iscloneable)
    {
       return false;
    }

   SCIP_RETCODE performSolLifting(
      SCIP*                 scip,               /**< SCIP main data structure */
      ProbDataCBP * probdata,             /** < problem data */
      SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
      );
 }; /*lint !e1712*/
 
 
 #endif