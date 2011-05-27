/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: BSFC_update_element_proc.C,v 1.1.1.1 2003/08/13 19:26:11 sorokine Exp $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include <hashtab.h>
#include <bucket.h>

#include "exvar.h"
#include "repartition_BSFC.h"
 
// this routine updates the processor that an element is assigned to

void BSFC_update_element_proc(int myid, int numprocs, 
			      HashTable* BG_mesh,
			      BSFC_VERTEX_PTR sfc_vert_ptr)
{
  int j;
  HTIterator *itr = new HTIterator(BG_mesh);
  Bucket *buck;
  while ((buck = (Bucket *) itr->next()))
    if ( buck->get_bucket_type () != GHOST )
      if(buck->get_myprocess() != sfc_vert_ptr[j].destination_proc) 
      // this element will get moved to a new processor
      { 
        buck->put_myprocess(sfc_vert_ptr[j].destination_proc);
      }

  delete itr;
  return;
}

