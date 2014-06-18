CC = gcc
CFLAGS = -pthread -std=gnu99 -O3 -m64 -fomit-frame-pointer

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

all: bwt histtest histcomptest fmitest searchtest rnaseqtest filetest

rnaseqtest: rnaseqtest.o histsortcomp.o seqindex.o csacak.o smw.o
	gcc -o $@ $^ $(CFLAGS)

#smw: smw.o
#	gcc -o $@ $^ $(CFLAGS)

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
	rm -f *.o histtest bwt histcomptest fmitest searchtest rnaseqtest smw filetest *~

.PHONY: clean all
