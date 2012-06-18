
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
 * $Id: hashtab.h 206 2009-01-26 17:32:10Z dkumar $ 
 */

/* Hash table */

/* Every table can process NBUCKETS entries */

/* All buckets with the same entry are linked together*/

/* The pointer to the first bucket is stored in the table */

/*---------------------------------------------------------*/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;

#include <constants.h>

typedef struct
{
  unsigned key[KEYLENGTH];
  void copy_key (unsigned k[])
  {
    int i;
    for (i = 0; i < KEYLENGTH; i++)
        key[i] = k[i];
  };

  void fill_key (unsigned k[])
  {
    for (int i = 0; i < KEYLENGTH; i++)
      k[i] = key[i];
  };
} Key;

inline bool
compare_keys (const Key & K1, const Key & K2)
{

  if (K1.key[0] < K2.key[0])
    return true;
  else if (K1.key[0] > K2.key[0])
    return false;
  else if (K1.key[1] < K2.key[1])
    return true;
  else if (K1.key[1] > K2.key[1])
    return false;
  else
    return false;
}

struct HashEntry
{
  unsigned key[KEYLENGTH];      // key: object key word
  void * value;                 // value: poiter to record
  HashEntry * pre;              // pre, next: objects with same entry 
  HashEntry * next;             // will be stored in a two-way linked-list

  HashEntry (unsigned * keyi)
  {
    int i;
    for (i = 0; i < KEYLENGTH; i++)
        key[i] = keyi[i];
      next = NULL;
      pre = NULL;
  }

  HashEntry ()
  {
    value = NULL;
    next = NULL;
    pre = NULL;
  }

  ~HashEntry ()                 //keep the follower when deleting an object
  {
    if (next)
      next->pre = pre;
    if (pre)
      pre->next = next;
  }
};

typedef HashEntry *HashEntryPtr;

class HashTable
{

protected:
  unsigned MinKey[KEYLENGTH];
  unsigned MaxKey[KEYLENGTH];
  unsigned Range;
  double umax;
  double minDom[DIMENSION];
  double maxDom[DIMENSION];

  HashEntryPtr *bucket;
  int NBUCKETS;
  int PRIME;
  int ENTRIES;
  int SIZE01;
  int SIZE02;

  HashEntryPtr addElement (int entry, unsigned *key);
  HashEntryPtr searchBucket (HashEntryPtr p, unsigned *key);

public:
  HashTable (unsigned *, unsigned *, int, int);
  HashTable (unsigned *, unsigned *,  int,  int, double *, double *);
  HashTable (int,  int, double *, double *);
  ~HashTable ();

  void add (unsigned * key, void * value);
  void *lookup (unsigned * key);
  void *lookup (Key);
  void remove (unsigned * key);
  int hash (unsigned * key);

  // remove entery
  void remove (Key k)
  {
    unsigned key[KEYLENGTH];
      k.fill_key (key);
      remove (key);
  }

  //! add entry to the hash-table
  void add (Key inkey, void * value)
  {
    unsigned key[KEYLENGTH];
    inkey.fill_key (key);
    add (key, value);
  }

  //! get size of the hash-table
  int get_no_of_buckets ()
  {
    return NBUCKETS;
  }

  //! get first entry
  HashEntryPtr * getbucketptr ()
  {
    return bucket;
  }

  //! get i-th bucket
  HashEntryPtr getBucket (int entry)
  {
    return bucket[entry];
  }

  //! get min of domain
  double *get_minDom ()
  {
    return minDom;
  }

  //! get max of domain
  double *get_maxDom ()
  {
    return maxDom;
  }

  //! get min of all keys
  unsigned *get_MinKey ()
  {
    return MinKey;
  }

  //! get max of all keys
  unsigned *get_MaxKey ()
  {
    return MaxKey;
  }
};


/*********************************
 *  Hash-Table iterator
 ********************************/
class HTIterator
{
private:
  HashTable * table;
  HashEntryPtr current;
  int index;
  int size;

public:
   ~HTIterator ()
  {
  };

  HTIterator (HashTable * ht);
  HashEntryPtr getNextBucket ();
  void *next ();
  void reset ()
  {
    current = table->getBucket (0);
    index = 0;
  };
};

#endif // HASHTABLE_H
