
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
 * $Id:$ 
 */

#ifndef BUCKET_HEAD__H
#  define BUCKET_HEAD__H

#  include <hashtab.h>

class BucketHead
{
protected:
  Key head;
  double coord[2];

public:
  //! constructor
  BucketHead(Key keyi, double * crd)
  {
    head = keyi;
    for (int i = 0; i < 2; i++)
      coord[i] = crd[i];
  }

  //! get key of first bucket
  Key get_head () const
  {
    return head;
  }

  //! get cooord
  const double * get_coord () const
  {
    return coord;
  }

  //! comparison operator for sorting
  bool operator < (const BucketHead & rhs) const
  {
    return compare_keys(head, rhs.get_head());
  }
};

#endif //BUCKET_HEAD__H
