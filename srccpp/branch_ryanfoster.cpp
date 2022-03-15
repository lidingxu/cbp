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

/**@file   branch_ryanfoster.c
 * @brief  Ryan/Foster branching rule
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

#include <assert.h>
#include <string.h>
#include <vector>
#include "objscip/objscip.h"

#include "branch_ryanfoster.h"
#include "cons_samediff.h"
#include "probdata_cbp.h"



/** branching execution method for fractional LP solutions */
SCIP_DECL_BRANCHEXECLP(BranchRyanFoster::scip_execlp) 
{  /*lint --e{715}*/

   SCIP_Real** pairweights;

   assert(scip != NULL);


   //SCIPdebugMsg(scip, "start branching at node %"SCIP_LONGINT_FORMAT", depth %d\n", SCIPgetNNodes(scip), SCIPgetDepth(scip));

   *result = SCIP_DIDNOTRUN;

   ProbDataCBP * probdata = NULL;
	probdata = dynamic_cast<ProbDataCBP *>(SCIPgetObjProbData(scip));
   assert(probdata != NULL);


   /* get fractional LP candidates */
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   int nlpcands;
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands, NULL) );
   assert(nlpcands > 0);

   upper_triangle & item_mat = probdata->item_matrix;
   item_mat.reset();


   /* compute weights for each order pair */
   for(int v = 0; v < nlpcands; v++ )
   {

      assert(lpcands[v] != NULL);

      SCIP_Real solval = lpcandsfrac[v];

      /* get variable data which contains the information to which constraints/items the variable belongs */
      SCIP_VARDATA* vardata = SCIPvarGetData(lpcands[v]);

      const PackVar * packiter = vardata->iter;
      const vector<int>  & items_bin = packiter->item_array;
      int nitems_bin = items_bin.size();
      assert(nitems_bin > 0);

      /* loop over all constraints/items the variable belongs to */
      for(int i = 0; i < nitems_bin; i++ )
      {
         int id1 = items_bin[i];

         /* store the LP sum for single items in the diagonal */
         item_mat.add(id1,id1, solval);

         /* update LP sums for all pairs of items */
         for(int j = i+1; j < nitems_bin; j++ )
         {
            int id2 = items_bin[j];
            assert(id1 < id2);

            item_mat.add(id1,id2, solval);
         }
      }
   }

   SCIP_Real bestvalue = 0.0;
   int best_id1 = -1;
   int best_id2 = -1;

   int nitems = probdata->numitems;
   for(int i = 0; i < nitems; i++ )
   {
      SCIP_Real mat_ii = item_mat.get(i,i);
      for(int j = i + 1; j  < nitems; j++ )
      {
         SCIP_Real mat_ij = item_mat.get(i,j);
         SCIP_Real frac = mat_ij - SCIPfloor(scip, mat_ij);
         SCIP_Real value = MIN(frac, 1-frac);

         if( bestvalue < value )
         {
            /* there is no variable with (fractional) LP value > 0 that contains exactly one of the items */
            if( SCIPisEQ(scip, mat_ij, mat_ii) && SCIPisEQ(scip, mat_ij, item_mat.get(j,j)) )
               continue;
            bestvalue = value;
            best_id1 = i;
            best_id2 = j;
         }
      }
   }
   assert( bestvalue > 0.0 );
   assert( best_id1 >= -1 && best_id1 < nitems);
   assert( best_id2 >= -1 && best_id2 < nitems);


   SCIP_NODE* childsame;
   SCIP_NODE* childdiffer;
   SCIP_CONS* conssame;
   SCIP_CONS* consdiffer;

   /* create the branch-and-bound tree child nodes of the current node */
   SCIP_CALL( SCIPcreateChild(scip, &childsame, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &childdiffer, 0.0, SCIPgetLocalTransEstimate(scip)) );

   /* create corresponding constraints */
   SCIP_CALL( SCIPcreateConsSamediff(scip, &conssame, "same", best_id1, best_id2, SAME, childsame, TRUE) );
   SCIP_CALL( SCIPcreateConsSamediff(scip, &consdiffer, "differ", best_id1, best_id2, DIFFER, childdiffer, TRUE) );

  /* add constraints to nodes */
   SCIP_CALL( SCIPaddConsNode(scip, childsame, conssame, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdiffer, consdiffer, NULL) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &conssame) );
   SCIP_CALL( SCIPreleaseCons(scip, &consdiffer) );

  *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


