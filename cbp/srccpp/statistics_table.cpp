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

/**@file   statistics_table.cpp
 * @brief  additional statistics table with pricer information
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include "statistics_table.h"
#include <math.h> 

SCIP_DECL_TABLEOUTPUT(StatisticsTable::scip_output){
    assert(scip != NULL);
    ProbDataCBP * probdata = NULL;
	probdata = dynamic_cast<ProbDataCBP *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);
    int col_exact = probdata->stat_pr.col_exact;
    SCIP_Real shf_log_sum_gap = probdata->stat_pr.shf_log_sum_gap;
    SCIP_Real shf_avg_gap = exp(shf_log_sum_gap / col_exact) - probdata->stat_pr.shf_param;
    SCIPinfoMessage(scip, file, "pricing column exact: %d\n", col_exact);   
    SCIPinfoMessage(scip, file, "pricing log sum shifted gap: %lf\n", shf_log_sum_gap);     
    SCIPinfoMessage(scip, file, "pricing avg gap: %lf\n", shf_avg_gap);     
    return SCIP_OKAY;
}