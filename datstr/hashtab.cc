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
 * $Id: hashtab2.C 206 2009-01-26 17:32:10Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#ifdef HAVE_MPI_H
#include <mpi.h>
#endif

#include <cassert>
#include <climits>
using namespace std;

#include "hashtab.h"

#define HASHTABLE_EXTENDER 1000000
#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

HashTable:: HashTable(unsigned* min, unsigned* max, int size, int prime)   
{
  MinKey[0] = *min;
  MinKey[1] = *(min+1); 
  MaxKey[0] = *max;
  MaxKey[1] = *(max+1);

  // extend the hashtable bounds a little bit to make it more efficient for adaptive meshes
  unsigned hashtable_extender = HASHTABLE_EXTENDER;
  if(MinKey[0] >= hashtable_extender)
    MinKey[0] -= hashtable_extender;
  else
    MinKey[0] = 0;
  unsigned umax = IScale;
  if((hashtable_extender/2 + MaxKey[0]/2) <= (umax/2))
    MaxKey[0] += hashtable_extender;
  else 
    MaxKey[0] = umax;
  
  NBUCKETS = size;
  PRIME  = prime;
  //Range  = *(MaxKey);
  Range  = *(MaxKey)-*(MinKey); //Keith Made this change 20061109

  bucket = new HashEntryPtr[NBUCKETS];  

  for(int i = 0; i< NBUCKETS;i++)
    bucket[i] = 0;

}


HashTable::HashTable(double *keyrangein, int size, int prime, 
                     double  minR[], double maxR[])
{
  int i;
  NBUCKETS = size;
  PRIME  = prime;
  
  for(i=0;i<KEYLENGTH;i++)
    keyrange[i]=keyrangein[i];
  
  hashconstant=8.0*NBUCKETS/
    (keyrange[0]*keyrange[1]+keyrange[1]);

  bucket = new HashEntryPtr[NBUCKETS];
  for(i = 0; i< NBUCKETS; i++)
    bucket[i] = NULL;

  for(i=0; i<DIMENSION; i++)
   {
     minDom[i]=minR[i];
     maxDom[i]=maxR[i];
     invrange[i]=1.0/(maxR[i]-minR[i]);
   }
}


HashTable:: ~HashTable()              //evacuate the table
{
  for (int i = 0; i<NBUCKETS; i++)
  {
    HashEntryPtr p = bucket[i]; 
    while (p)
    { 
       HashEntryPtr p_next = p->next; 
       delete p; 
       p = p_next;  
    } 
  }
  delete [] bucket;
}


HashEntryPtr HashTable::searchBucket(HashEntryPtr p, unsigned* keyi)
{
  int i;
  while(p) 
  {
    for(i=0;i<KEYLENGTH;i++) 
    {
      if(p->key[i] != *(keyi+i))
      { 
        p = p->next;                     
        break;                          //not found, check next element
      } 
      else if(i == KEYLENGTH -1)            
        return p;                       //found, return element pointer
    }   
  }
  return p;
}

HashEntryPtr HashTable:: addElement(int entry, unsigned key[])
{    
  HashEntryPtr p = new HashEntry(key);  
  if((bucket[entry])) //this place is already occupied
  {
    HashEntryPtr currentPtr = bucket[entry];
    while(currentPtr!=0 && (key[0] > currentPtr->key[0] ))
    {
      p->pre=currentPtr;
      currentPtr=currentPtr->next;
    }
    if(currentPtr!=0 && key[0] == currentPtr->key[0])
    {
      while(currentPtr!=0 && (key[1] > currentPtr->key[1] ))
      {
        p->pre=currentPtr;
        currentPtr=currentPtr->next;
      }  
    }
	 
    if(currentPtr) currentPtr->pre=p;
    p->next = currentPtr;
    currentPtr = p->pre;
    if(currentPtr) currentPtr->next=p;
    else bucket[entry]=p;
  }			

  //  p->next = *(bucket+entry);        //add the bucket to the head
  else 
    bucket[entry] = p;                //else eliminate it
  return p;
}

void * HashTable::lookup(unsigned* key)
{
  int entry = hash (key);
  HashEntryPtr p = searchBucket( bucket[entry], key );
  if(!p) 
    return NULL;                      //if not found, return 0
  return p->value;                    //if found return a pointer  
}

// lookup using key structure
void * HashTable::lookup(Key kstr)
{
  unsigned key[KEYLENGTH];

  kstr.fill_key(key);
  int entry = hash (key);
  HashEntryPtr p = searchBucket( bucket[entry], key );
  if (!p)
    return NULL;
  return p->value;
}

void HashTable:: add (unsigned* key, void* value)
{
  int entry = hash(key);
  HashEntryPtr p = searchBucket (bucket[entry], key);
  if (p == NULL)  //was (!p)
  {
     p = addElement(entry, key);
     p->value = value;
  }
  return;
}

void HashTable:: remove (unsigned* key)
{
  int entry = hash(key);
  HashEntryPtr p = searchBucket(bucket[entry], key);
  if(!p)
    return;
  if(bucket[entry] == p) 
  {
    bucket[entry] = p->next;
    delete p;
  }
  else
  {
    if(!(p->next))
       delete p;    
    else
    {
      (p->pre)->next = p->next;
      (p->next)->pre = p->pre;
      delete p;
    }
  }
}

int HashTable::hash(unsigned* key)
{ 
  int which_entry;
  which_entry = (((int) ((key[0]*keyrange[1]+key[1])*hashconstant+0.5))%NBUCKETS);
  return which_entry;
}



// Hashtable iterator
HTIterator::HTIterator(HashTable *ht)
{
  table = ht;
  current = table->getBucket(0);
  size  = table->get_no_of_buckets();
  index = 0;
}   

HashEntryPtr HTIterator::getNextBucket()
{
  while (table->getBucket(++index) == NULL && index < size );
  if (index == size ) return NULL;
  return table->getBucket(index);
}

void * HTIterator::next()
{
  void * value;
  if (current)  // if current pointer has more links
  {
    value = current->value;
    current = current->next;  // move to next link, after getting the value
    return value;
  }
  current = getNextBucket(); // get next valid Hashtable entry
  if (current == NULL) return NULL; 
  value = current->value;
  current = current->next;   // move to next link, after getting the value
  return value;
}
