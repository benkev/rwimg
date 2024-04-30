/*
 * dummy.c
 * Insert argc line with a string in double quotes
 */

/* #ifndef ARRAYOBJECT */
/* #define ARRAYOBJECT  /usr/lib/python2.6/site-packages/numpy/core/include/numpy/arrayobject.h */
/* #endif */
#include "arrayobject.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  //  FILE *inpf, *outf;
  size_t n = 1024;
  char *line = (char *) malloc(n*sizeof(char));
  
  if (argc != 2) {
    printf("Must be one parameter: the line to insert\n");
    return -1;
  }

  fprintf(stdout, "#define ARRAYOBJECT \"%s\"\n", argv[1]);
  
  while (-1 != getline (&line, &n, stdin)) printf("%s\n", line);
    
  return 0;
}
