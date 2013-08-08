#ifndef _HISTSORTCOMP_H
#define _HISTSORTCOMP_H

void histhelper(const char *, int, int *, int *, int, int, int);

int * histsort(const char *, int);

void putsg(const char *, int, int);

void putbwt(const char *, const int *, int);

void sprintbwt(char *, const char *, const int *, int);

char * makebwt(const char *, int);

int makecbwt(const char *, int, char *);

int saca_makecbwt(const char *, int, char *);

int sprintcbwt(const char *, int *, int, char *);

#endif // _HISTSORTCOMP_H
