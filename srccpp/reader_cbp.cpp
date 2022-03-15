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

/**@file   reader_cbp.cpp
 * @brief  Conic Binpacking problem reader
 * @author Liding XU
 *
 * This file implements the reader/parser used to read the cbp input data. For more details see \ref NETWORKROUTING_READER.
 *
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <utility>

#include "scip/cons_setppc.h"
#include "objscip/objscip.h"

#include "probdata_cbp.h"
#include "reader_cbp.h"

using namespace scip;
using namespace std;

/** destructor of file reader to free user data (called when SCIP is exiting) */
SCIP_DECL_READERFREE(ReaderCBP::scip_free)
{
	return SCIP_OKAY;
} /*lint !e715*/


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
SCIP_DECL_READERWRITE(ReaderCBP::scip_write)
{
	*result = SCIP_DIDNOTRUN;

	return SCIP_OKAY;
} /*lint !e715*/

/**@} */

/** problem reading method of reader
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
 *
 *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
 */
SCIP_DECL_READERREAD(ReaderCBP::scip_read) {
	*result = SCIP_DIDNOTRUN;

   SCIPdebugMessage("Start read!\n");
	// open the file
	ifstream filedata(filename);
	if (!filedata) {
		return SCIP_READERROR;
	}
	filedata.clear();

	// read parameters
    SCIP_Real capacity, alpha, Dalpha;
	string s;
	int numitems;
	filedata >> s >> capacity >> numitems >> alpha;
	filedata >> Dalpha;

	vector<SCIP_Real> mus(numitems);
	vector<SCIP_Real> bs(numitems);

	// read the mus
	int received_nummus = 0;
	while (!filedata.eof()) {
		filedata >> mus[received_nummus];
		received_nummus++;
		if (received_nummus == numitems) {
			break;
		}
	}
	assert(received_nummus == numitems);

	// read the bs
	int received_numbs = 0;
	while (!filedata.eof()) {
		filedata >> bs[received_numbs];
		received_numbs++;
		if (received_numbs == numitems) {
			break;
		}
	}
	assert(received_numbs == numitems);


	// create the problem's data structure
	ProbDataCBP * problemdata = NULL;
	problemdata = new ProbDataCBP(numitems, capacity, Dalpha, mus, bs);
	assert(problemdata != NULL);
	SCIPdebugMessage("--problem data completed!\n");
	// deletedobject
	SCIP_CALL(SCIPcreateObjProb(scip, filename, problemdata, FALSE));

	SCIPdebugMessage("objprob created and creating inital solutions!\n");

	// create constraints initial columns
	problemdata->item_matrix = upper_triangle(numitems);
	problemdata->item_matrix.reset();
	problemdata->currentnode = -1;
	problemdata->belongs = vector<int> (numitems);
	problemdata->global_lb = -SCIPinfinity(scip);


	SCIP_CALL(problemdata->createConsInitialColumns(scip));
 
	static const char* CKNAP_PRICER_NAME =  "CKNAP_Pricer";
   	SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, CKNAP_PRICER_NAME)) );
   
   *result = SCIP_SUCCESS;

	SCIPdebugMessage("--reader read completed!\n");
	return SCIP_OKAY;
}


/**@} */
