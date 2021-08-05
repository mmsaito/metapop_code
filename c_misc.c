#include <stdlib.h>

/* qsort wrapper for fortran */

/* sharing box between sorting func. and comparator */
static struct {
  int (*comp)(const void *x, const void *y);
  int sizeOf;
  void *items;
  int idxbase;
} share__idx_qsort;

void *items_share;

/* qsort's comparator */
int comp__idx_qsort(const void *i, const void *j)
{
  void *p, *q;
  int ii, jj;

  ii = *(int *)i - share__idx_qsort.idxbase;
  jj = *(int *)j - share__idx_qsort.idxbase;
  //printf("%d,%d,%d,%d\n", ii, jj, *i, *j);
  p = (unsigned char*)share__idx_qsort.items + share__idx_qsort.sizeOf*ii;
  q = (unsigned char*)share__idx_qsort.items + share__idx_qsort.sizeOf*jj;
  return (*share__idx_qsort.comp)(p,q);
}

/* A wrapper of qsort for index sorting and called from fortran */
void idx32_qsort_(void *items, int *n, int *sizeOf, int indexes[], int (*comp)(const void *x, const void *y), int *idxbase)
{
   share__idx_qsort.comp   = comp;
   share__idx_qsort.sizeOf = *sizeOf;
   share__idx_qsort.items  = items;

   if (idxbase) share__idx_qsort.idxbase = *idxbase;
   else         share__idx_qsort.idxbase = 0;

   qsort(indexes, *n, sizeof(int), comp__idx_qsort);
}
