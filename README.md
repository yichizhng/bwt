TODO:

Proper scoring (also, increase affine gap penalty so it doesn't do stupid
things)
Maybe actually write alignment data out in SAM format

Things to consider:

LF(x) is a fairly inexpensive function (although often resulting in cache
misses probably...); how much, on average, does storing the entire SA help
(basically, we expect 15 calls to LF() on average)?
As written, the FM-index currently takes up 37.5% of the space of an
uncompressed chromosome (meaning that it would take up 150% that of
a compressed one, since we have the 1:4 byte packing, or 75% of a .2bit file
(since those have masks)). Do we get significant advantages by storing more
data? (For example, we could try storing seq_index, although that one is quite
a pain in the ass, considering that it's actually an array of pointers, or
more of the suffix array (Bowtie stores the same ratio I do, but its is
adjustable))

Methods:

The basic idea of most modern alignment tools is to use the FM-index (Full-text
minute space), which consists of the Burrows-Wheeler transform of the sequence
along with some information which allows us to do a backward search in linear
time (in particular, we need to be able to do a rank lookup in constant time,
and also know the count of each symbol).

It is very inconvenient to construct the Burrows-Wheeler transform without
calculating the entire suffix array (i.e. the sorted list of all suffixes of the
sequence), and since we need at least a partial suffix array to implement
locate() in reasonable (ish) time we may as well calculate it.

The most obvious and naive way to do this is sorting the suffixes. Obviously we
don't want to store them all, so we can represent them as offsets into the
sequence. This is demonstrated in bwt.c, which uses the library quicksort. Since
strcmp is O(n) (even though the library version is probably accelerated somehow)
this gives us an average runtime of O(n^2 log n), which is quite undesirable.
(As well as very memory expensive, since we sort an array of pointers)

Try: ./bwt mississippi

A more advanced idea is to use a linear time sort, such as bucketsort or
histogram sort. Each bucketing step can be done in linear time (there are no
comparisons, and the bucket is determined simply by indexing into the sequence),
which gives us an average runtime of O(n log n). Histogram sort, in particular,
has the advantages of being easily parallelized and having good cache coherence,
allowing it to outperform the more "optimal" algorithms at the cost of higher
memory usage. histsort.c contains the "naive" implementation of a bucket sort
on an alphabet of size 4 (0, 1, 2, and 3 refer to A, C, G, and T respectively)

It is, of course, useful to note that a nucleotide can be stored in only two
bits, which allows us a compression ratio of 25%, which may even speed up
operations (despite the additional calls to getbase() needed - see any file
where getbase() is defined) - this is done in histsortcomp.c. This speeds up the
computation of the BWT somewhat, due mainly to minor improvements in cache
coherence (which is a large factor in this kind of task).

fmitest.c contains a "full" implementation of an FM-index. The suffix array
is built using either the histogram sort in histsortcomp.c or SACA-K (in
csacak.c; this is an adaptation of Ge Nong's sacak.cpp), then used to build the
Burrows-Wheeler transformed string. The idxs member of a fm_index holds a
partial suffix array, allowing SA[i] to be computed in constant time by using
the LF-mapping (note, however, that this is a bit slow). Lookup is an elegantly
generated array which holds counts of each nucleotide for any given 4-nucleotide
sequence (see the lookup_table() function in seq_index.c to see how it's
generated and how the data is organized). This, along with the array rank_index,
allows us to calculate the rank (often referred to as the Occ() function in
the context of FM-indexes) of a given nucleotide in the Burrows-Wheeler
transformed string.

More "optimal" algorithms for constructing the suffix array usually aim to
achieve O(n) runtime while using minimal auxiliary space. SACA-K is one of these
algorithms, and is an improvement on SA-IS; it uses practically constant space
(to be precise, it's O(log n)) and is essentially linear in runtime, but is
slower than the multithreaded histogram sort in most practical cases (a highly
repetitive sequence is the pathological case for histogram sort). SACA-K was
chosen mainly to minimize memory usage; it has constant auxiliary space
requirement.

Backward search can be done in O(m) time (i.e. constant in sequence length), but
the locate() function (i.e. associating a particular match with its position
on the genome) requires O(m + log n) time (in particular the association
requires O(log n) time and is essentially the calculation of SA[i]; this is
the unc_sa() function - unc is short for "uncompress"). (A more naive approach,
though sometimes useful, is to do a binary search over the index; however, this
is generally a worse idea than doing the backwards search on the reverse
FM-index, since the binary search is O(m log n) or slower)

One important note is that many operations involving the Burrows-Wheeler
transform and the original sequence are extremely cache unfriendly; for example,
the main loop of the sprintcbwt() function (which takes the original sequence
and the suffix array and prints the BWT) results in a cache miss above 99.9% of
the time for longer sequences. This is due to the fact that the access pattern
required to reconstruct the original sequence from the BWT or vice versa is
essentially (i.e. for the purposes of cache line eviction strategy) random
access over a large array, which is essentially an impossible case (even
sequential access with large stride can (theoretically) be dealt with via
prefetching, although few architectures actually do this). This suggests
that we should avoid things like retrieving the original sequence from the
Burrows-Wheeler transform if possible (or if we do so, only do it once.
1 billion bps * 100 cycles/bp (and that's a rather pessimistic estimate) means
it will take less than a minute anyway. We want O(n) of these poor accesses,
not O(m), since m is going to be a lot bigger).

The main problem we have to cope with here is mismatches, which can be caused
by transcription errors or by differences between the reference genome and
the individual's genome. Different alignment tools deal with this in different
ways; however, it is certainly a computationally expensive problem to deal with.
Bowtie uses greedy backtracking reverse search (I haven't looked at the source
code, but I imagine it tries replacing a character in the sequence with one
that results in more matches at that particular point). Seed-and-extend methods
(which originated when the main method of sequence alignment was hashing) are
also viable, with seeds being generated via various types of search (STAR, for
example, uses binary search to find what it refers to as a maximum mappable
prefix)

rnaseqtest.c contains an implementation which does gapped alignment using
maximum mappable suffixes, augmented by several tricks. The main idea behind
using suffixes is that we can achieve O(m + log(n)) performance by doing a
backwards search on the FM-index rather than a binary search, which is
O(m log (n)) (we could use prefixes by reversing the corpus, but this seems
immaterial). Instead of using backtracking search, we instead simply skip
ahead a few nucleotides (with the intention of later using the Needleman-Wunsch
algorithm (implemented in smw.c) to stitch together these matches. Under
reasonable assumptions regarding the number and kind of transcription errors
we should be able to achieve O(m + log (n)) speed

On input format:

The genome is expected to be given as a single text file; the characters
A, C, T, and G will be treated as their corresponding nucleotides; all others
will be treated as A. There is a utility (filread.cc) which turns FastA genomes
into the expected format (it also changes all unrecognized characters into G)

Reads are expected to be given one per line. Unrecognized characters will be
treated as if they are 'N'. 'N' causes some odd behavior in alignment algorithms
(some slowness, mostly, and possible lack of sensitivity if too many occur)

The algorithm has trouble detecting certain patterns of errors (two errors
within about 10 nts of each other will cause the gap calculator to fail).
On the other hand, not backtracking is fast and it's something that can
be dealt with by the rest of the alignment algorithm.