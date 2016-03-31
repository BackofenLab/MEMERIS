/*
 * $Id: hash.c,v 1.2 2005/10/25 19:06:39 nadya Exp $
 * 
 * $Log: hash.c,v $
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 00:20:29  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/**********************************************************************/
/*
	hashing functions
*/
/**********************************************************************/
/* hash.c */
/* 5-24-00 tlb; fix hashing function and make generic */

#include "macros.h"
#include "hash.h"

#define two24 	8388610			/* 2**24 */

/**********************************************************************/
/*
	hash
*/
/**********************************************************************/
extern int hash(
  char *key1, 				/* character key */
  int key2,				/* integer key */
  int n					/* modulo */
)
{
  int i, p, d;

  for (i=0, d=1, p=key2; key1[i] != 0; i++, d *= 256) {
    if (d >= two24) d = 1;		/* d<2**23 insures d*key1<2**31 */
    p += (d * key1[i]) % n;
  }
  return p % n;
}
  

/**********************************************************************/
/*
	hash_create
*/
/**********************************************************************/
extern HASH_TABLE hash_create(
  int n
)
{
  int i;
  HASH_TABLE ht = (HASH_TABLE) mymalloc(sizeof(struct hash_table));

  /* initialize hash table */
  ht->n = n;
  ht->table = (HASH_TABLE_ENTRY **) mymalloc(n * sizeof(HASH_TABLE_ENTRY *));
  for (i=0; i<n; i++) ht->table[i] = NULL;

  return ht;
}

/**********************************************************************/
/*
	hash_destroy
*/
/**********************************************************************/
extern void hash_destroy(
  HASH_TABLE ht 		/* hash table to destroy */
)
{
  int i;
  for (i=0; i < ht->n; i++) {
    HASH_TABLE_ENTRY *hte = ht->table[i];
    while (hte != NULL) {
      HASH_TABLE_ENTRY *next = hte->next;
      myfree(hte->key1);
      myfree(hte);
      hte = next;
    }
  }
  myfree(ht->table);
  myfree(ht);
}

/**********************************************************************/
/*
	hash_insert
*/
/**********************************************************************/
extern void hash_insert(
  char *key1,			/* character key */
  int key2,			/* integer key */
  HASH_TABLE ht			/* the hash table */
)
{
  HASH_TABLE_ENTRY *hte;

  /* compute the position in hash table of entry */
  int p = hash(key1, key2, ht->n);

  /* see if already present */
  hte = ht->table[p];
  while (hte != NULL) {
    if ((hte->key2 == key2) && (!strcmp(hte->key1, key1))) break;
    hte = hte->next;
  }
  
  if (!hte) { 			/* not already present--insert it */
    /* make a hash-table entry */
    hte = (HASH_TABLE_ENTRY *) mymalloc(sizeof(HASH_TABLE_ENTRY));
    hte->key1 = (char *) mymalloc (strlen(key1) + 1);
    strcpy(hte->key1, key1);
    hte->key2 = key2;

    /* insert the entry into the table */
    hte->next = ht->table[p];
    ht->table[p] = hte;
  }
}

/**********************************************************************/
/*
	hash_lookup
*/
/**********************************************************************/
extern BOOLEAN hash_lookup(
  char *key1,			/* character key */
  int key2,			/* integer key */ 
  HASH_TABLE ht			/* the hash table */
)
{
  /* compute the position in hash table of entry */
  int p = hash(key1, key2, ht->n);

  /* get collision list for this hash index */
  HASH_TABLE_ENTRY *hte = ht->table[p];
  while (hte != NULL) {
    if ((hte->key2 == key2) && (!strcmp(hte->key1, key1))) return TRUE; 
    hte = hte->next;
  }
  return FALSE;
}

