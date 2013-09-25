CC = gcc
CFLAGS = -pthread -std=gnu99 -m64 -fomit-frame-pointer -g
# Has no compiler warnings, unless you're the kind of silly person who
# likes turning on extra warnings and reading through them ("of *course* I'm
# indexing this array with a char, that's why it has size 256")
# At the very least, valgrind has nothing to complain about

# With -Wall -pedantic -Wextra gcc will complain about things like
# char indexing into arrays (intended and correct behavior), comparison
# between unsigned and signed ints (intended and correct behavior),
# values possibly being used uninitialized (they aren't), argc being
# ignored (it is, but the program will segfault quickly with wrong arguments),
# passing pointer targets with different signedness (that's a bug, it should
# be obvious that it should be unsigned), and some questioning of my knowledge
# of operator precedence

all: bwt histtest histcomptest fmitest searchtest rnaseqtest

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

rnaseqtest: rnaseqtest.o histsortcomp.o seqindex.o csacak.o smw.o
	gcc -o $@ $^ $(CFLAGS)

#smw: smw.o
#	gcc -o $@ $^ $(CFLAGS)

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
	rm -f *.o histtest bwt histcomptest fmitest *~

.PHONY: clean all
