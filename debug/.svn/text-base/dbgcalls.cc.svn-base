/*
 * =====================================================================================
 *
 *       Filename:  dbgcalls.cc
 *
 *    Description:  
 *
 *        Created:  03/03/2011 03:57:30 PM EST
 *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * =====================================================================================
 * $Id:$
 */

#include <cstdio>
using namespace std;

#include "dbgcalls.h"


void init_debug_vars (int myid, DbgVars *dbg_struct, int write2f)
{
#ifdef DEBUG
  if ( write2f != 0 )
  {
    sprintf(dbg_struct.filename, "debug%03d.log", myid);
    dbg_struct.fileptr = fopen(dbg_struct.filename, "w+");
    dbg_struct.write2file = true;
  }
  else
  {
    dbg_struct.fileptr = stdout;
    dbg_struct.write2file = false;
  }
#else
  dbg_struct.trace_enabled = false;
#endif
}

void dbg_printf(const char *format, ...)

