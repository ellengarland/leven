#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include "hashmap.h"

#define INTERN_STRINGS 1

#define MIN(a,b,c) (a < b ? (a < c ? a : c) : (b < c ? b : c))

int intern(map_t map, int* next, const char* data)
{
   uintptr_t index;
   if (hashmap_get(map, (char*)data, (void**)&index) == MAP_OK)
      return index;
   index = (*next)++;
   hashmap_put(map, (char*)data, (void*)index);
   return index;
}

SEXP cleven(SEXP x, SEXP y)
{
   int * rows[2];
   SEXP result;
   int i;
   int write_row = 0;
#ifdef INTERN_STRINGS
   map_t hashmap;
   int* xx;
   int* yy;
#endif

   if ((TYPEOF(x) != STRSXP && TYPEOF(x) != NILSXP) || (TYPEOF(y) != STRSXP && TYPEOF(y) != NILSXP))
      error("invalid input");

   result = PROTECT(allocVector(REALSXP, 1));

   if (TYPEOF(x) == NILSXP && TYPEOF(y) == NILSXP)
   {
      REAL(result)[0] = 0;
      UNPROTECT(1);
      return result;
   }
   if (TYPEOF(x) == NILSXP)
   {
      REAL(result)[0] = length(y);
      UNPROTECT(1);
      return result;
   }

   if (TYPEOF(y) == NILSXP)
   {
      REAL(result)[0] = length(x);
      UNPROTECT(1);
      return result;
   }

#ifdef INTERN_STRINGS
   // This converts each unique input string to a number, which may make comparisons faster
   hashmap = hashmap_new();
   int next = 0;
   xx = malloc(length(x) * sizeof(int));
   yy = malloc(length(y) * sizeof(int));
   for (int i = 0; i < length(x); i++)
      xx[i] = intern(hashmap, &next, CHAR(STRING_ELT(x, i)));
   for (int i = 0; i < length(y); i++)
      yy[i] = intern(hashmap, &next, CHAR(STRING_ELT(y, i)));
   hashmap_free(hashmap);
#endif

   rows[0] = malloc((length(y) + 1) * sizeof(int));
   rows[1] = malloc((length(y) + 1) * sizeof(int));
   memset(rows[0], 0, sizeof(int) * (length(y) + 1));
   memset(rows[1], 0, sizeof(int) * (length(y) + 1));
   for (i = 0; i < length(x); i++)
   {
      rows[write_row & 1][0] = i+1;
      for (int j = 0; j < length(y); j++)
      {
#ifdef INTERN_STRINGS
         int cost = (xx[i] == yy[j])?0:1;
#else
         int cost = strcmp(CHAR(STRING_ELT(x, i)), CHAR(STRING_ELT(y, j))) == 0?0:1;
#endif
         int add = rows[write_row & 1][j] + 1;
         int del = rows[(write_row + 1) & 1][j+1] + 1;
         int subs = rows[(write_row + 1) & 1][j] + cost;
         rows[write_row & 1][j+1] = MIN(add, del, subs);
      }
      write_row++;
   }
   REAL(result)[0] = rows[(write_row + 1) & 1][length(y)];
   UNPROTECT(1);
   return result;
}
