CC = gcc
CFLAGS = -pthread -O3 -std=gnu99 -m64 -fomit-frame-pointer -g

all: bwt histtest histcomptest fmitest searchtest

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

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
