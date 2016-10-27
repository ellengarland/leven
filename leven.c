#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include "hashmap.h"

// This file contains two R-callable functions:
//  * cleven(Vector1, Vector2, CostMap): A very fast implementation of the levenshtein distance.
//    Vector1 and Vector2 are assumed to be vectors of strings. First, the strings are all interned
//    into an atom table, essentially numbering all the possible input values and replacing them with
//    those numbers. This means comparisons between values is very fast (and independent of the length
//    of the values). Then we do the 2-row, iterativel leven calculation, with substitution costs from
//    CostMap. To get a CostMap from a cost matrix, use optimise_cost_matrix()
//  * optimise_cost_matrix(): Converts a labelled matrix of costs into a O(1) lookup array. This
//    operation is relatively slow, which is why it cannot be done for every call to cleven().
//    Note that we must share the label indices between this map and the interned values used in the
//    cleven() function.


#define MIN(a,b,c) (a < b ? (a < c ? a : c) : (b < c ? b : c))

int intern(map_t map, int* next, const char* data)
{
   uintptr_t index;
   if (hashmap_get(map, (char*)data, (void**)&index) == MAP_OK)
      return index;
   index = (*next)++;
   hashmap_put(map, strdup(data), (void*)index);
   return index;
}

typedef struct
{
   double** costs;  // Actual transition costs
   map_t hashmap;   // Map of strings -> atoms
   int next;        // The next index in the atom table
   int cost_atoms;  // The number of atoms for which there is a cost value
   int rows;        // The number of rows in the table
} costmap_t;

static int _free_intern(any_t key, any_t data)
{
   free(key);
   return MAP_OK;
}

static void _finalizer(SEXP ext)
{
   if (NULL == R_ExternalPtrAddr(ext))
      return;
   costmap_t *ptr = (costmap_t *) R_ExternalPtrAddr(ext);
   for (int i = 0; i < ptr->rows; i++)
      free(ptr->costs[i]);
   free(ptr->costs);
   hashmap_iterate(ptr->hashmap, _free_intern);
   hashmap_free(ptr->hashmap);
   free(ptr);
   R_ClearExternalPtr(ext);
}

// optimise_cost_matrix() returns a pointer to a C structure that we can efficiently search for costs
// If we use such an object, it is imperative that we use the same dictionary for later runs. This means
// the dictionary will grow unbounded (well, bounded only by the number of atoms in the system)
// R should free the resources associated with the costmap when the object is finalized by the GC
SEXP optimise_cost_matrix(SEXP cost_matrix)
{
   if (TYPEOF(cost_matrix) == REALSXP)
   {
      costmap_t* costmap = malloc(sizeof(costmap_t));
      costmap->hashmap = hashmap_new();
      costmap->next = 0;
      SEXP names = getAttrib(cost_matrix, R_DimNamesSymbol);
      SEXP colnames = VECTOR_ELT(names, 0);
      SEXP rownames = VECTOR_ELT(names, 1);
      // We must also intern all the colnames and rownames in case they do not overlap exactly
      for (int i = 0; i < length(colnames); i++)
         intern(costmap->hashmap, &costmap->next, CHAR(STRING_ELT(colnames, i)));
      for (int j = 0; j < length(rownames); j++)
         intern(costmap->hashmap, &costmap->next, CHAR(STRING_ELT(rownames, j)));

      // Now allocate the cost array
      costmap->rows = costmap->next;
      costmap->costs = malloc(sizeof(double*) * costmap->next);
      for (int i = 0; i < costmap->next; i++)
         costmap->costs[i] = malloc(sizeof(double*) * costmap->next);

      // And then for each element in the cost matrix, put that value into our cost array
      for (int i = 0; i < length(colnames); i++)
      {
         int ix = intern(costmap->hashmap, &costmap->next, CHAR(STRING_ELT(colnames, i)));
         for (int j = 0; j < length(rownames); j++)
         {
            int jx = intern(costmap->hashmap, &costmap->next, CHAR(STRING_ELT(rownames, j)));
            double value = REAL(cost_matrix)[j * length(rownames) + i];
            costmap->costs[ix][jx] = value;
         }
      }
      SEXP ptr = PROTECT(R_MakeExternalPtr(costmap, R_NilValue, R_NilValue));
      R_RegisterCFinalizerEx(ptr, _finalizer, TRUE);
      UNPROTECT(1);
      costmap->cost_atoms = costmap->next;
      return ptr;
   }
   else
      error("Cost matrix has the wrong type");
}

SEXP cleven(SEXP x, SEXP y, SEXP cost_matrix)
{
   double * rows[2];
   SEXP result;
   int i;
   int write_row = 0;
   map_t hashmap;
   int* xx;
   int* yy;
   costmap_t* costmap = NULL;
   int* next;

   // First, do some type checks and check for degenerate cases
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

   // If we are using a costmap, we MUST use its dictionary. This implies the dictionary will mutate as
   // more calls to cleven() occur. Atoms are never released (until the costmap is freed)
   if (TYPEOF(cost_matrix) == EXTPTRSXP)
   {
      costmap = (costmap_t*)R_ExternalPtrAddr(cost_matrix);
      next = &costmap->next;
      hashmap = costmap->hashmap;
   }
   else if (TYPEOF(cost_matrix) == NILSXP)
   {
      int tmp = 0;
      costmap = NULL;
      hashmap = hashmap_new();
      next = &tmp;
   }
   else
      error("cost_matrix has the wrong type. Make sure you have run it through optimise_cost_matrix() before passing it to this function");

   // Intern all the input values, building xx and yy as a C array of ints representing the atom values corresponding to x and y
   xx = malloc(length(x) * sizeof(int));
   yy = malloc(length(y) * sizeof(int));
   for (int i = 0; i < length(x); i++)
      xx[i] = intern(hashmap, next, CHAR(STRING_ELT(x, i)));
   for (int i = 0; i < length(y); i++)
      yy[i] = intern(hashmap, next, CHAR(STRING_ELT(y, i)));

   // If we have no costmap, we can free the dictionary now
   if (costmap == NULL)
   {
      hashmap_iterate(hashmap, _free_intern);
      hashmap_free(hashmap);
   }

   // Allocate space for the calculation
   rows[0] = malloc((length(y) + 1) * sizeof(double));
   rows[1] = malloc((length(y) + 1) * sizeof(double));

   // Populate the first row
   for (int i = 0; i < length(y)+1; i++)
      rows[0][i] = i;
   // Clear the second row for safety
   memset(rows[1], 0, sizeof(double) * (length(y) + 1));
   // Start by writing to the second row
   write_row = 1;

   // Iterate, building the tables
   for (i = 0; i < length(x); i++)
   {
      rows[write_row & 1][0] = i+1;
      for (int j = 0; j < length(y); j++)
      {
         double cost;
         if (costmap != NULL && xx[i] < costmap->cost_atoms && yy[j] < costmap->cost_atoms)
            cost = (xx[i] == yy[j])?0:costmap->costs[xx[i]][yy[j]];
         else
            cost = (xx[i] == yy[j])?0:1;
         double add = rows[write_row & 1][j] + 1;
         double del = rows[(write_row + 1) & 1][j+1] + 1;
         double subs = rows[(write_row + 1) & 1][j] + cost;
         rows[write_row & 1][j+1] = MIN(add, del, subs);
      }
      write_row++;
   }
   // Return the value in the bottom-right corner
   REAL(result)[0] = rows[(write_row + 1) & 1][length(y)];
   free(xx);
   free(yy);
   free(rows[0]);
   free(rows[1]);
   UNPROTECT(1);
   return result;
}
