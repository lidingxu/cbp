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

/**@file   branch_ryanfoster.h
 * @brief  Ryan/Foster branching rule
 * @author Liding Xu
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_RYANFOSTER_H__
#define __SCIP_BRANCH_RYANFOSTER_H__



#include "objscip/objscip.h"

using namespace scip;

/** C++ wrapper object for event handlers */
class BranchRyanFoster : public scip::ObjBranchrule
{
public:
	/** default constructor */
	BranchRyanFoster(
		SCIP* scip
	)
		: ObjBranchrule(scip, "branchRyanFoster", " branch rule for packing variables ", 50000, -1, 1.0)
	{

	}


	/** branching execution method for fractional LP solutions */
	virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);

}; /*lint !e1712*/


#endif
