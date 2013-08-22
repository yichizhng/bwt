// Generates a randomish genome for us to test

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  int i, x;
  char *arr;
  srand(time(0));
  if(argc > 1) {
    x = atoi(argv[1]);
    if (!x) {
      printf("Usage: gen_seq [x]\n");
    }
  }
  else {
    x = 100000;
  }
  arr = malloc(4);
  arr[0] = 'A';
  arr[1] = 'C';
  arr[2] = 'G';
  arr[3] = 'T';
  for (i = 0; i < x; ++i) {
    putchar(arr[rand() & 3]);
  }
}
