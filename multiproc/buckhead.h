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

#ifndef BUCKET__HEAD_
#define BUCKET__HEAD_

#include <hashtab.h>

class BucketHead {
  protected:
    Key bucket;
    double coord;

  public:
    //! constructor
    BucketHead(Key k, double x)
    {
      bucket = k;
      coord  = x;
    }

    //! get key of first bucket
    Key get_bucket() 
    {
      return bucket;
    }

    //! get cooord
    double get_coord() const
    {
      return coord;
    }

    //! comparison operator for sorting
    bool operator< (BucketHead rhs) const
    {
      return ( coord < rhs.get_coord());
    }
};

#endif //BUCKET__HEAD_
