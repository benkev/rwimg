/*
 * insert_define_arrayobject.c
 *
 * Insert at the beginning of the stdin a line 
 * #define ARRAYOBJECT "path"
 * with a path string in double quotes.
 * The path must be the only command line parameter.
 * Example:
 * $ insert_define_arrayobject /usr/lib/python2.6/site-packages \
 *       /numpy/core/include/numpy/arrayobject.h
 */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  size_t n = 256;
  char *line = (char *) malloc(n*sizeof(char));
  
  if (argc != 2) {
    printf("Must be one parameter: the line to insert\n");
    return -1;
  }

  fprintf(stdout, "#define ARRAYOBJECT \"%s\"\n", argv[1]);
  
  while (-1 != getline (&line, &n, stdin)) printf("%s\n", line);
    
  return 0;
}
