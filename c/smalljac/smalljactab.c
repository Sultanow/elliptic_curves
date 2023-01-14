#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "cstd.h"
#include "smalljactab.h"
#include "ff_poly.h"


/*
    Copyright (c) 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/


/*
	We allocate three sections of memory, all of which will be small in practice (a few megabytes).
	
	(1) the hash table - 64 bits per entry
	(2) the lookaside table - 128 bits per entry
	(3) overflow lists - 128 bits per entry
	
	In most cases, only (1) gets used, and this is what we want to keep in cache.   We don't really
	care about the size of the rest.
	
	Currently the number of overflow entries is equal the size of the hash table.  If we exceed this we
	are getting a whole lot of collisions and should increase the table size anway.
*/

void smalljac_table_alloc (int bits)
{
	assert (bits>0);
	if ( bits > 31 ) { err_printf ("smalljac_table_alloc: bits=%d too large!\n", bits);  exit (0); }
	smalljac_htab = malloc(sizeof(*smalljac_htab)*(1<<bits));				// use malloc rather than mem_alloc - we don't need the memory initialized
	// allocate space equal to table size for overflow - this could be made bigger
	smalljac_htab_lists = malloc(sizeof(*smalljac_htab_lists)*2*(1<<bits));		// ditto
	smalljac_htab_end = smalljac_htab_lists + 2*(1<<bits);
	smalljac_table_bits = bits;
}


void smalljac_table_init (int bits)
{
	static int warn;
	
	if ( bits > smalljac_table_bits )  {
		if ( ! warn ) { printf ("%lu: smalljac_table_init: bits=%d exceeds allocated table size of %d bits\n", _ff_p, bits, smalljac_table_bits);  warn = 1; }
		if ( bits-smalljac_table_bits > 4 ) { printf ("you must increase the allocated table size");  exit (0); }
		bits = smalljac_table_bits;
	}
	smalljac_tabmask = (1UL<<bits)-1;
	memset (smalljac_htab, 0, (1<<bits)*sizeof(unsigned long));
	smalljac_htab_next = smalljac_htab_lists + (1<<bits);
}


void _smalljac_table_collision_insert (int index, unsigned long item)
{
//	printf ("collision at index %d\n", index);
	if ( ! smalljac_htab_lists[index].item ) { smalljac_htab_lists[index].item = item;  smalljac_htab_lists[index].pNext = 0;  return; }
	if ( smalljac_htab_next >= smalljac_htab_end ) {
		err_printf ("%lu: Ran out of htab overflow entries smalljac_htab_next-smalljac_htab_lists = %ld - increase table size.\n", _ff_p, smalljac_htab_next-smalljac_htab_lists);  exit (0); }
	smalljac_htab_next->item = item;
	smalljac_htab_next->pNext = smalljac_htab_lists[index].pNext;
	smalljac_htab_lists[index].pNext = smalljac_htab_next++;
}


int _smalljac_table_get_matches (uint32_t matches[SMALLJAC_MAX_MATCHES], jac_t a[1], int index)
{
	register unsigned long item;
	register struct smalljac_htab_list_item *pItem;
	register int i;
	
#if HECURVE_GENUS == 1
	item = a[0].u[0];
#else
	item = (_ff_one(a[0].u[HECURVE_GENUS]) ? (1UL<<SMALLJAC_ITEM_BITS) : 0) + a[0].u[HECURVE_GENUS-1];
	item <<= SMALLJAC_ITEM_BITS;
	item += a[0].u[HECURVE_GENUS-2];
#endif
#if HECURVE_GENUS == 3
	item <<= SMALLJAC_ITEM_BITS;
	item += a[0].u[HECURVE_GENUS-3];
#endif	
	item <<= SMALLJAC_VBITS;
	i = 0;
	if ( (smalljac_htab[index] ^ item) <= SMALLJAC_VMASK ) matches[i++] = smalljac_htab[index]&SMALLJAC_VMASK;
	if ( ! smalljac_htab_lists[index].item ) return i;
	if ( (smalljac_htab_lists[index].item ^ item) <= SMALLJAC_VMASK ) matches[i++] = smalljac_htab_lists[index].item&SMALLJAC_VMASK;
	if ( ! smalljac_htab_lists[index].pNext ) return i;
	pItem = smalljac_htab_lists[index].pNext;
	do {
		if ( (pItem->item ^ item) <= SMALLJAC_VMASK ) {
			if ( i >= SMALLJAC_MAX_MATCHES ) { err_printf ("Exceeded SMALLJAC_MAX_MATCHES=%d with p=%ld!\n", i, _ff_p);  exit (0); }
			matches[i++] = pItem->item&SMALLJAC_VMASK;
		}
		pItem = pItem->pNext;
	} while ( pItem );
	return i;
}


