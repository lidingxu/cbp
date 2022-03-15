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

/**@file   reader_bp.h
 * @brief  Binpacking problem reader file reader
 * @author Liding Xu
 *
 * This file implements the reader/parser used to read the binpacking input data. For more details see \ref BINPACKING_READER.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_BP_H__
#define __SCIP_READER_BP_H__


#define SCIP_DEBUG

#include <iostream>
#include <fstream>
#include <string>
#include "objscip/objscip.h"



/** SCIP file reader for NTRT data files */
class ReaderBP : public scip::ObjReader
{
public:

	/** default constructor */
	ReaderBP(SCIP* scip)
		: scip::ObjReader(scip, "reader", "file reader for Binpack files", "bp")
	{}


	/** destructor of file reader to free user data (called when SCIP is exiting) */
	virtual SCIP_DECL_READERFREE(scip_free);

	/** problem reading method of reader
	 *
	 *  possible return values for *result:
	 *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
	 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
	 *
	 *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
	 */
	virtual SCIP_DECL_READERREAD(scip_read);

	/** problem writing method of reader; NOTE: if the parameter "genericnames" is TRUE, then
	 *  SCIP already set all variable and constraint names to generic names; therefore, this
	 *  method should always use SCIPvarGetName() and SCIPconsGetName();
	 *
	 *  possible return values for *result:
	 *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
	 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
	 *
	 *  If the reader detected an error in the writing to the file stream, it should return
	 *  with RETCODE SCIP_WRITEERROR.
	 */
	virtual SCIP_DECL_READERWRITE(scip_write);



};/*lint !e1712*/



#endif
