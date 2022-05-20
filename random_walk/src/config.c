#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

#include "random_walk.h"
#include "data.h"



int write_config(Par *par, int *walk_buff, char *fname)
{
  create_if_not_exists("conf/");
  char filename[FNAMESIZE + 5] = "conf/";
  strcat(filename, fname);
  int fdesc = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (fdesc < 0) {
    printf("Can't open %s\n", filename);
    return 0;
  }
 
  int nw = write(fdesc, walk_buff, par->N * 2 * sizeof(int));
  close(fdesc);

  return 1;
}


// Return the size in bytes of the file for the given file descriptor.
int filebytes(int fdesc)
{
  struct stat buf;

  buf.st_size = 0;
  fstat(fdesc, &buf);
  return buf.st_size;
}


int read_config(Par *par, int *walk_buff, char *fname)
{
  int bytes, st, fdesc;
  char filename[FNAMESIZE + 5] = "conf/";

  strcat(filename, fname);
  bytes = par->N * 2 * sizeof(int);

  printf("Read config file %s", filename);

  // Try to open file
  fdesc = open(filename, O_RDONLY);
  if (fdesc < 0) {
    printf("...not successful.\n");
    return 0;
  }

  // Check file size
  if (bytes == filebytes(fdesc)) {
    int nr = read(fdesc, walk_buff, bytes);
    printf(".\n");
    st = 1;
  }
  else {
    printf("...bad file size: %d bytes expected.\n", bytes);
    st = 0;
  }

  close(fdesc);

  return st;
}
