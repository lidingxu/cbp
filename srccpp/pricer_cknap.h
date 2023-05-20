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

/**@file   pricer_cbp.h
 * @brief  Conic Binpacking variable pricer
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BIN_PRICER_CONICBINPACKING__
#define __BIN_PRICER_CONICBINPACKING__

#include "objscip/objscip.h"
#include "scip/pub_var.h"
#include "probdata_cbp.h"
#include "scip/scip_numerics.h"

#include <vector>
#include <list>
#include <utility>

using namespace std;
using namespace scip;


/** pricer class */
class PricerConicKnap : public ObjPricer
{
public:

	/** Constructs the pricer object with the data needed */
	PricerConicKnap(
		SCIP*                               scip,       /**< SCIP pointer */
		const char* CKNAP_PRICER_NAME /** < Pricer name */
	) : ObjPricer(scip, CKNAP_PRICER_NAME, "Finds pack with negative reduced cost.", 0, TRUE) 
	{
	}

	/** Destructs the pricer object. */
	virtual ~PricerConicKnap() {
	};

	/** reduced cost pricing method of variable pricer for feasible LPs */
	virtual SCIP_DECL_PRICERREDCOST(scip_redcost);

	
	/** farkas pricing method of variable pricer for infeasible LPs */
	virtual SCIP_DECL_PRICERFARKAS(scip_farkas);
};

/** solve conic knapsack problem
*   possible return values for:
*  - Negative_Sol: negative reduced cost solution is found, the correspoding primal feasible solution stored in items_bin
*  - Positive_Sol: no solution with negative reduces cost is found, empty items_bin
*  - Unsolved: the subproblem is unsolved.
*/

#endif
