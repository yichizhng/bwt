CC = gcc
CFLAGS = -pthread -O3 -std=gnu99 -m64 -fomit-frame-pointer
OBJS = histsort.o bwt.o histsortcomp.o fmitest.o seqindex.o

all: bwt histtest histcomptest fmitest

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

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
