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

/**@file   cons_samediff.h
 * @brief  Constraint handler stores the local branching decision data
 * @author Liding Xu
 * This constraint handler is used to store the branching decision of the \ref BINPACKING_BRANCHING "Ryan/Foster branching rule"
 * which is implemented in \ref branch_ryanfoster.c.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_SAMEDIFF_H__
#define __SCIP_CONS_SAMEDIFF_H__


#include "objscip/objscip.h"
#include "probdata_cbp.h"

using namespace scip;
using namespace std;

/* type of constraint: differ or same */
enum ConsType
{
   DIFFER = 0,                               /**< constraint representing the branching decision differ(i,j) */
   SAME   = 1                                /**< constraint representing the branching decision same(i,j) */
};
typedef enum ConsType CONSTYPE;


/** C++ constraint handler for  constraints */
class ConshdlrSameDifferent : public ObjConshdlr
{
public:
	/** default constructor */
	ConshdlrSameDifferent(
		SCIP* scip /**< SCIP data structure */
	)
		: ObjConshdlr(scip, /**< SCIP data structure */
			"cons_samediff", /**< 	name of constraint handler */
			"stores the samedifferent decisions", /**< description of constraint handler*/
			0, /**< priority of the constraint handler for separation */
			0, /**< . priority of the constraint handler for constraint enforcing */
			9999999, /**< . priority of the constraint handler for checking infeasibility (and propagation) */
			-1,  /**< frequency for separating cuts; zero means to separate only in the root node */
			1,  /**< . frequency for propagating domains; zero means only preprocessing propagation */
			1, /**< . frequency for using all instead of only the useful constraints in separation,
                *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
			0,  /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
			FALSE, /**< should separation method be delayed, if other separators found cuts? */
			FALSE,   /**< should propagation method be delayed, if other propagators found reductions? */
			TRUE,  /**< should the constraint handler be skipped, if no constraints are available? */
			SCIP_PROPTIMING_BEFORELP, /**< positions in the node solving loop where propagation method of constraint handlers should be executed */
			SCIP_PRESOLTIMING_FAST  /**< timing mask of the constraint handler's presolving method */)
	{
	}

	/** feasibility check method of constraint handler for primal solutions */
	virtual SCIP_DECL_CONSCHECK(scip_check);

	/** constraint enforcing method of constraint handler for LP solutions */
	virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

	/** constraint enforcing method of constraint handler for pseudo solutions */
	virtual SCIP_DECL_CONSENFOPS(scip_enfops);

	/** variable rounding lock method of constraint handler */
	virtual SCIP_DECL_CONSLOCK(scip_lock);

	/** frees specific constraint data
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 */
	virtual SCIP_DECL_CONSDELETE(scip_delete);

	/** transforms constraint data into data belonging to the transformed problem */
	virtual SCIP_DECL_CONSTRANS(scip_trans);


	/** domain propagation method of constraint handler
	 *
	 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The propagation
	 *  method should process only the useful constraints in most runs, and only occasionally the remaining
	 *  nconss - nusefulconss constraints.
	 *
	 *  possible return values for *result:
	 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
	 *  - SCIP_REDUCEDDOM : at least one domain reduction was found
	 *  - SCIP_DIDNOTFIND : the propagator searched, but did not find any domain reductions
	 *  - SCIP_DIDNOTRUN  : the propagator was skipped
	 *  - SCIP_DELAYED    : the propagator was skipped, but should be called again
	 */
	virtual SCIP_DECL_CONSPROP(scip_prop);

	/** constraint activation notification method of constraint handler
	 *
	 *  @see SCIP_DECL_CONSACTIVE(x) in @ref type_cons.h
	 */
	virtual SCIP_DECL_CONSACTIVE(scip_active);

	/** constraint deactivation notification method of constraint handler
	 *
	 *  @see SCIP_DECL_CONSDEACTIVE(x) in @ref type_cons.h
	 */
	virtual SCIP_DECL_CONSDEACTIVE(scip_deactive);
	
}; /*lint !e1712*/

/** creates and captures a samediff constraint */
SCIP_RETCODE SCIPcreateConsSamediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   itemid1,            /**< item id one */
   int                   itemid2,            /**< item id two */
   CONSTYPE              type,               /**< stores whether the items have to be in the SAME or DIFFER packing */
   SCIP_NODE*            node,               /**< the node in the B&B-tree at which the cons is sticking */
   SCIP_Bool             local               /**< is constraint only valid locally? */
   );

/** return constraint type SAME or DIFFER */
CONSTYPE SCIPgetTypeSamediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

/** return the items pairs in the same constraints */
vector<pair<int,int>>  getItemsInSame(
	SCIP*                 scip              /**< SCIP data structure */
);

/** return the items pairs in the differ constraints */
vector<pair<int,int>> getItemsInDiffer(
	SCIP*                 scip              /**< SCIP data structure */
);
#endif
