CC = gcc
CFLAGS = -pthread -std=gnu99 -O3 -m64

# Requirements: Some sort of reasonable x86 or x86-64 system (for the former,
# compile with -m32 and edit rdtscll.h to use the 32-bit version), some sort
# of C compiler that speaks C99 (in particular initial loop declarations),
# Posix threads (although it's relatively simple to remove that requirement)
# Memory usage is rather high, especially if you use the histogram sort

# Has no compiler warnings, unless you're the kind of person who likes turning
# on extra warnings and reading through them ("of *course* I'm indexing this
# array with a char, that's why it has size 1024")

# With -Wall -pedantic -Wextra gcc will complain about things like
# char indexing into arrays (intended and correct behavior), comparison
# between unsigned and signed ints (intended and correct behavior),
# values possibly being used uninitialized (they aren't), argc being
# ignored (it is, but the program will segfault quickly with wrong arguments),
# unused variables (true enough, but harmless), passing pointer targets with
# different signedness (that's a bug, although one I haven't gotten around to
# caring about), and some questioning of my knowledge of operator precedence
# (unwarranted)

# TODO: Remove some of the functions from seqindex.c so that it can be
# used without including histsortcomp and csacak; the program which
# aligns reads should not need to 

TESTS =  bwt histtest histcomptest fmitest searchtest rnaseqtest filetest gaptest build_index index_test search_reads single_align

all: $(TESTS)

single_align: histsortcomp.o csacak.o single_align.o fileio.o seqindex.o smw.o stack.o
	gcc -o $@ $^ $(CFLAGS)

search_reads: histsortcomp.o seqindex.o csacak.o search_reads.o fileio.o
	gcc -o $@ $^ $(CFLAGS)

rnaseqtest: rnaseqtest.o histsortcomp.o seqindex.o csacak.o smw.o stack.o
	gcc -o $@ $^ $(CFLAGS)

#smw: smw.o
#	gcc -o $@ $^ $(CFLAGS)

index_test: index_test.o fileio.o seqindex.o csacak.o histsortcomp.o
	gcc -o $@ $^ $(CFLAGS)

build_index: build_index.o histsortcomp.o csacak.o fileio.o seqindex.o
	gcc -o $@ $^ $(CFLAGS)

gaptest: gaptest.o histsortcomp.o seqindex.o csacak.o 
	gcc -o $@ $^ $(CFLAGS)

filetest: filetest.o histsortcomp.o seqindex.o csacak.o fileio.o
	gcc -o $@ $^ $(CFLAGS)

searchtest: searchtest.o histsortcomp.o seqindex.o csacak.o
	gcc -o $@ $^ $(CFLAGS)

fmitest: histsortcomp.o fmitest.o seqindex.o csacak.o
	gcc -o $@ $^ $(CFLAGS)

histcomptest: histsortcomp.o histsortcomptest.o csacak.o
	gcc -o $@ $^ $(CFLAGS)

histtest: histtest.o histsort.o sacak.o
	gcc -o $@ $^ $(CFLAGS)

bwt: bwt.o
	gcc -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o $(TESTS) *~

.PHONY: clean all
