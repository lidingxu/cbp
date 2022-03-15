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

/**@file   statistics_table.h
 * @brief  additional statistics table with pricer information
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STATISTICS_TABLE__
#define __STATISTICS_TABLE__

#include "objscip/objscip.h"
#include "objscip/objtable.h"
#include "scip/pub_var.h"
#include "probdata_cbp.h"
#include "scip/scip_numerics.h"
#include "objscip/objcloneable.h"

using namespace std;
using namespace scip;


/** StatisticsTable class */
class StatisticsTable : public ObjTable 
{
public:
	/** Constructs the table  object with the data needed */
	StatisticsTable(
		SCIP*                               scip,       /**< SCIP pointer */
		const char* TABLE_NAME /** < table name */
	) : ObjTable(scip, TABLE_NAME, "add additional statistics information", 20001,  SCIP_STAGE_SOLVING) 
	{
	};

	/** Destructs the table object. */
	virtual ~StatisticsTable() {
	};

    /** output method of statistics table to output file stream 'file'
     *
     *  @see SCIP_DECL_TABLEOUTPUT(x) in @ref type_table.h
     */
    virtual SCIP_DECL_TABLEOUTPUT(scip_output);

};


#endif
