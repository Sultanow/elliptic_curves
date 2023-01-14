#ifndef _SMALLJACTAB_INCLUDE_
#define _SMALLJACTAB_INCLUDE_

/*
    Copyright (c) 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <stdint.h>
#include "ff_poly.h"
#include "jac.h"
#include "hecurve.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
	Simple table lookup/insert code for smalljac.  This is designed for small tables
        that fit in cache memory, using 64 bits per entry (128 bits per overflow entry).
	This requires that the caller to either store the group entries in a seperate list
	or reconstruct them as required.
	
        The table size is always a power of two and hashing is done by simply bit-masking
	the group element for speed.  This can occasionally cause performance problems
	but works well on average.
	
	Note that the table size used in smalljac_init may be any power of 2 smaller than
	the allocated table size - this is desirable as it reduces the cose of initialization
	and improves locality when this is as small as possible without creating too many
	collisions.  A load factor of around 0.5 is works well.
*/

#if HECURVE_GENUS == 1
#define SMALLJAC_VBITS			17			// must be < 32
#define SMALLJAC_ITEM_BITS		21			// ITEM_BITS+1+VBITS < 64
#endif
#if HECURVE_GENUS == 2
#define SMALLJAC_VBITS			24			// must be < 32
#define SMALLJAC_ITEM_BITS		19			// 2*ITEM_BITS+1+VBITS < 64
#endif
#if HECURVE_GENUS == 3
#define SMALLJAC_VBITS			20			// must be < 32
#define SMALLJAC_ITEM_BITS		14			// 3*ITEM_BITS+1+VBITS < 64
#endif
#define SMALLJAC_VMASK		((1UL<<SMALLJAC_VBITS)-1)
#define SMALLJAC_MAX_MATCHES	128


struct smalljac_htab_list_item {
	unsigned long item;
	struct smalljac_htab_list_item *pNext;
};
unsigned long *smalljac_htab;
struct smalljac_htab_list_item *smalljac_htab_lists, *smalljac_htab_next, *smalljac_htab_end;
ff_t smalljac_tabmask;
int smalljac_table_bits;

void smalljac_table_alloc (int bits);
void smalljac_table_init (int bits);

void _smalljac_table_collision_insert (int index, unsigned long item);
int _smalljac_table_get_matches (uint32_t *matches, jac_t a[1], int index);

static inline void smalljac_table_insert (jac_t a[1], uint32_t value)
{
	register ff_t i;
	register unsigned long item;
	
#if HECURVE_GENUS == 1
	i = a[0].u[0] & smalljac_tabmask;
	item = a[0].u[0];
#else  // genus >= 2
	i = a[0].u[1] & smalljac_tabmask;
	item = (_ff_one(a[0].u[HECURVE_GENUS]) ? (1UL<<SMALLJAC_ITEM_BITS) : 0) + a[0].u[HECURVE_GENUS-1];
	item <<= SMALLJAC_ITEM_BITS;
	item += a[0].u[HECURVE_GENUS-2];
#if HECURVE_GENUS == 3
	item <<= SMALLJAC_ITEM_BITS;
	item += a[0].u[HECURVE_GENUS-3];
#endif	
#endif
	item <<= SMALLJAC_VBITS;
	item += value;
	if ( ! smalljac_htab[i] ) { smalljac_htab[i] = item;  smalljac_htab_lists[i].item = 0;  return; }
	_smalljac_table_collision_insert ((int)i, item);
}


static inline int smalljac_table_lookup (uint32_t matches[SMALLJAC_MAX_MATCHES], jac_t a[1])
{
	register ff_t i;
	
#if HECURVE_GENUS == 1
	i = a[0].u[0] & smalljac_tabmask;
#else
	i = a[0].u[1] & smalljac_tabmask;
#endif
	if ( ! smalljac_htab[i] ) return 0;		// usual case
	return _smalljac_table_get_matches (matches, a, (int)i);
}

#ifdef __cplusplus
}
#endif

#endif
