TODO:
Replace all instances of getbase() with a header defining it (still as a
static inline function)
Fix the locate() function, it probably doesn't work right now (in particular
the problem is with retrieving the full SA from the FM-index; this is
theoretically doable in constant time (for any given element) but I've messed
something up)
Write some basic test for RNA-seq (should really fix fmitest.h first, I suspect
that I'll need to include it in some things)
Ideas: stitching, semi-backtracking search? (for when we run out of seeds
of sufficient length)

Things to consider:
Is it advantageous for us to keep the suffix array, the original sequence, and
the FM-index all in memory at once? (~18.25x the size of the sequence I think)
It dramatically reduces the cache eviction rate to not have to retrieve the
original sequence (notice that the BWT is essentially (that is, with regards
to caches) a random permutation of the original sequence and therefore causes
massive numbers of cache misses (the percentage of cache *hits* scales
inversely with the length of the sequence!); on the other hand we should
not have to retrieve that much of the original sequence at a time; on the
third hand, the retrieval costs us on the order of 10ns, which is 30+ cycles,
per nucleotide, so retrieving a sequence of ~50 nucleotides could easily
cost us 2000+ cycles instead of 70!
