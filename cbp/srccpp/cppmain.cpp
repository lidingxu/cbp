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

/**@file   cppmain.cpp
 * @brief  Main file for branch and price conic binpacking problem
 * @author Liding Xu
 */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "branch_ryanfoster.h"
#include "cons_samediff.h"
#include "pricer_cknap.h"
#include "reader_cbp.h"
#include "reader_bp.h"
#include "statistics_table.h"
#include "rmp_heur.h"

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   
   /* we explicitly enable the use of a debug solution for this main SCIP instance */
   SCIPenableDebugSol(scip);

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include conic binpacking reader */
   //SCIP_CALL( SCIPincludeObjReader(scip, new ReaderBP(scip), TRUE));
   SCIP_CALL( SCIPincludeObjReader(scip, new ReaderCBP(scip), TRUE));

   /* include binpacking branching and branching data */
   SCIP_CALL(SCIPincludeObjBranchrule(scip, new BranchRyanFoster(scip), TRUE));
   SCIP_CALL( SCIPincludeObjConshdlr(scip, new ConshdlrSameDifferent(scip), TRUE) );
   
  /* include conic binpacking pricer  */
   static const char* CKNAP_PRICER_NAME = "CKNAP_Pricer";
   PricerConicKnap * pricer = new PricerConicKnap(scip, CKNAP_PRICER_NAME);
   SCIP_CALL( SCIPincludeObjPricer(scip, pricer, TRUE));

   /* include additional  statistics table */
   static const char* TABLE_NAME = "pricing_information";
   StatisticsTable * table = new  StatisticsTable(scip, TABLE_NAME);
   SCIP_CALL( SCIPincludeObjTable(scip, table, TRUE));   

   /* include pattern combination heuristics */
   SCIP_CALL(SCIPincludeObjHeur(scip, new HeurRMP(scip), TRUE));


    /* add cbp solver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "cbp/is_misocp","use cplex's misocp algorithm to solve pricing algorithm,otherwise use PLOA algorithm",  NULL,FALSE, FALSE, NULL,  NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "cbp/is_bd_tight","apply bound tightenning",  NULL, FALSE, TRUE,  NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "cbp/is_heur","run heuristic algorithm first (hybrid pricing)",  NULL, FALSE, TRUE, NULL,  NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "cbp/is_stablize","use stabilization for hybrid pricing",  NULL, FALSE, FALSE, NULL,  NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "cbp/is_parallelscplex","enbale cplex's parallelism",  NULL, FALSE, FALSE,  NULL, NULL) );


   /* for column generation instances, disable restarts */
   SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) );

   /* parameter setting */
   SCIP_CALL(SCIPsetRealParam(scip, "limits/gap", 1e-4));
   SCIP_CALL(SCIPsetRealParam(scip, "limits/absgap", 1e-6));

   /* turn off all separation algorithms */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   
   /**********************************
    * Process command line arguments *
    **********************************/
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runShell(argc, argv, "scip.set");
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}

