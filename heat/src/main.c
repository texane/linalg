/*
** Made by fabien le mentec <texane@gmail.com>
** 
** Started on  Mon May 31 16:51:53 2010 texane
** Last update Mon May 31 21:32:53 2010 texane
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>



/* ij pair */

typedef struct ij_pair
{
  size_t i;
  size_t j;
} ij_pair_t;

static inline ij_pair_t* make_ij_pair
(ij_pair_t* ijp, size_t i, size_t j)
{
  ijp->i = i;
  ijp->j = j;
  return ijp;
}

static inline size_t ij_to_index
(size_t dim, const ij_pair_t* ijp)
{
  return ijp->i * dim + ijp->j;
}

static inline ij_pair_t* ij_to_grid
(ij_pair_t* ijp)
{
  ijp->i += 1;
  ijp->j += 1;
  return ijp;
}

static inline ij_pair_t* ij_to_matrix
(ij_pair_t* ijp)
{
  ijp->i -= 1;
  ijp->j -= 1;
  return ijp;
}


/* grid structure */

typedef struct grid
{
  size_t dim;
  double side_temps[4];
} grid_t;

static void grid_init_once(grid_t* g)
{
  g->dim = 7;
}

static void grid_load_string(grid_t* g, const char* s)
{
  /* north:east:south:west */

  double n[4];

#define GRID_NORTH_SIDE 0
#define GRID_EAST_SIDE 1
#define GRID_SOUTH_SIDE 2
#define GRID_WEST_SIDE 3

  sscanf(s, "%lf:%lf:%lf:%lf", &n[0], &n[1], &n[2], &n[3]);

#if 0 /* initial values */
  g->side_temps[GRID_NORTH_SIDE] = 20.f;
  g->side_temps[GRID_EAST_SIDE] = 20.f;
  g->side_temps[GRID_SOUTH_SIDE] = 30.f;
  g->side_temps[GRID_WEST_SIDE] = 25.f;
#endif

  g->side_temps[GRID_NORTH_SIDE] = n[0];
  g->side_temps[GRID_EAST_SIDE] = n[1];
  g->side_temps[GRID_SOUTH_SIDE] = n[2];
  g->side_temps[GRID_WEST_SIDE] = n[3];
}

static inline double grid_get_side_temp
(const grid_t* g, unsigned int i)
{
  /* north, clockwise */
  return g->side_temps[i];
}

static void get_neigh_indices
(size_t adim, size_t index, size_t* nis)
{
#define NEIGH_NORTH_SIDE 0
#define NEIGH_EAST_SIDE 1
#define NEIGH_SOUTH_SIDE 2
#define NEIGH_WEST_SIDE 3

  nis[NEIGH_NORTH_SIDE] = index - adim;
  nis[NEIGH_EAST_SIDE] = index + 1;
  nis[NEIGH_SOUTH_SIDE] = index + adim;
  nis[NEIGH_WEST_SIDE] = index - 1;
}

static void generate_a(const grid_t* g, gsl_matrix* a)
{
  /* generate the coefficient matrix */

  const unsigned int adim = g->dim - 2;
  const unsigned int nrows = g->dim - 2;
  const unsigned int ncols = g->dim - 2;

  size_t index;
  size_t i;
  size_t j;
  size_t k;
  ij_pair_t ijp;

  size_t nis[4];

  gsl_matrix_set_zero(a);

  /* border columns */

  for (i = 1; i < (nrows - 1); ++i)
  {
    /* first col, north, east, south */
    index = ij_to_index(adim, make_ij_pair(&ijp, i, 0));
    gsl_matrix_set(a, index, index, 4.f);
    get_neigh_indices(adim, index, nis);
    gsl_matrix_set(a, index, nis[NEIGH_NORTH_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_EAST_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_SOUTH_SIDE], -1.f);

    /* last col, north, west, south */
    index = ij_to_index(adim, make_ij_pair(&ijp, i, ncols - 1));
    gsl_matrix_set(a, index, index, 4.f);
    get_neigh_indices(adim, index, nis);
    gsl_matrix_set(a, index, nis[NEIGH_NORTH_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_WEST_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_SOUTH_SIDE], -1.f);
  }

  /* border rows */

  for (j = 1; j < (ncols - 1); ++j)
  {
    /* first row, east, south, west */
    index = ij_to_index(adim, make_ij_pair(&ijp, 0, j));
    gsl_matrix_set(a, index, index, 4.f);
    get_neigh_indices(adim, index, nis);
    gsl_matrix_set(a, index, nis[NEIGH_EAST_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_SOUTH_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_WEST_SIDE], -1.f);

    /* last row, north, east, west */
    index = ij_to_index(adim, make_ij_pair(&ijp, nrows - 1, j));
    gsl_matrix_set(a, index, index, 4.f);
    get_neigh_indices(adim, index, nis);
    gsl_matrix_set(a, index, nis[NEIGH_NORTH_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_EAST_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_WEST_SIDE], -1.f);
  }

  /* corner cases, north west first, clockwise */

  index = ij_to_index(adim, make_ij_pair(&ijp, 0, 0));
  get_neigh_indices(adim, index, nis);
  gsl_matrix_set(a, index, index, 4.f);
  gsl_matrix_set(a, index, nis[NEIGH_EAST_SIDE], -1.f);
  gsl_matrix_set(a, index, nis[NEIGH_SOUTH_SIDE], -1.f);

  index = ij_to_index(adim, make_ij_pair(&ijp, 0, ncols - 1));
  get_neigh_indices(adim, index, nis);
  gsl_matrix_set(a, index, index, 4.f);
  gsl_matrix_set(a, index, nis[NEIGH_WEST_SIDE], -1.f);
  gsl_matrix_set(a, index, nis[NEIGH_SOUTH_SIDE], -1.f);

  index = ij_to_index(adim, make_ij_pair(&ijp, nrows - 1, ncols - 1));
  get_neigh_indices(adim, index, nis);
  gsl_matrix_set(a, index, index, 4.f);
  gsl_matrix_set(a, index, nis[NEIGH_WEST_SIDE], -1.f);
  gsl_matrix_set(a, index, nis[NEIGH_NORTH_SIDE], -1.f);

  index = ij_to_index(adim, make_ij_pair(&ijp, nrows - 1, 0));
  get_neigh_indices(adim, index, nis);
  gsl_matrix_set(a, index, index, 4.f);
  gsl_matrix_set(a, index, nis[NEIGH_EAST_SIDE], -1.f);
  gsl_matrix_set(a, index, nis[NEIGH_NORTH_SIDE], -1.f);

  /* inner points */

  for (i = 1; i < nrows - 1; ++i)
  {
    for (j = 1; j < ncols - 1; ++j)
    {
      index = ij_to_index(adim, make_ij_pair(&ijp, i, j));
      gsl_matrix_set(a, index, index, 4.f);

      get_neigh_indices(adim, ij_to_index(adim, &ijp), nis);
      for (k = 0; k < 4; ++k)
	gsl_matrix_set(a, index, nis[k], -1.f);
    }
  }
}

static void generate_b(const grid_t* g, gsl_vector* b)
{
  /* generate the result matrix */

  const unsigned int adim = g->dim - 2;
  const unsigned int nrows = g->dim - 2;
  const unsigned int ncols = g->dim - 2;

  size_t index;
  size_t i;
  size_t j;
  ij_pair_t ijp;
  double temp;

  gsl_vector_set_zero(b);

  /* border columns */

  for (i = 1; i < (nrows - 1); ++i)
  {
    /* first col, north, east, south */
    index = ij_to_index(adim, make_ij_pair(&ijp, i, 0));
    gsl_vector_set(b, index, grid_get_side_temp(g, GRID_WEST_SIDE));

    /* last col, north, west, south */
    index = ij_to_index(adim, make_ij_pair(&ijp, i, ncols - 1));
    gsl_vector_set(b, index, grid_get_side_temp(g, GRID_EAST_SIDE));
  }

  /* border rows */

  for (j = 1; j < (ncols - 1); ++j)
  {
    /* first row, east, south, west */
    index = ij_to_index(adim, make_ij_pair(&ijp, 0, j));
    gsl_vector_set(b, index, grid_get_side_temp(g, GRID_NORTH_SIDE));

    /* last row, north, east, west */
    index = ij_to_index(adim, make_ij_pair(&ijp, nrows - 1, j));
    gsl_vector_set(b, index, grid_get_side_temp(g, GRID_SOUTH_SIDE));
  }

  /* corner cases, north west first, clockwise */

  temp = grid_get_side_temp(g, GRID_NORTH_SIDE);
  temp += grid_get_side_temp(g, GRID_WEST_SIDE);
  index = ij_to_index(adim, make_ij_pair(&ijp, 0, 0));
  gsl_vector_set(b, index, temp);

  temp = grid_get_side_temp(g, GRID_NORTH_SIDE);
  temp += grid_get_side_temp(g, GRID_EAST_SIDE);
  index = ij_to_index(adim, make_ij_pair(&ijp, 0, ncols - 1));
  gsl_vector_set(b, index, temp);

  temp = grid_get_side_temp(g, GRID_SOUTH_SIDE);
  temp += grid_get_side_temp(g, GRID_EAST_SIDE);
  index = ij_to_index(adim, make_ij_pair(&ijp, nrows - 1, ncols - 1));
  gsl_vector_set(b, index, temp);

  temp = grid_get_side_temp(g, GRID_SOUTH_SIDE);
  temp += grid_get_side_temp(g, GRID_WEST_SIDE);
  index = ij_to_index(adim, make_ij_pair(&ijp, nrows - 1, 0));
  gsl_vector_set(b, index, temp);
}

static void __attribute__((unused)) generate_ab
(const grid_t* g, gsl_matrix* a, gsl_vector* b)
{
  /* generate both a and b matrices */

  generate_a(g, a);
  generate_b(g, b);
}


static void __attribute__((unused)) solve_linear_system
(gsl_matrix* a, gsl_vector* x, const gsl_vector* b)
{
  gsl_permutation* const p = gsl_permutation_alloc(a->size1);

  int s;

  gsl_linalg_LU_decomp(a, p, &s);
  gsl_linalg_LU_solve(a, p, b, x);

  gsl_permutation_free(p);
}


static void alloc_linear_system
(const grid_t* g, gsl_matrix** a, gsl_vector** x, gsl_vector** b)
{
  /* there are dim * dim unknowns */

  const size_t dim = g->dim - 2;
  const size_t unknown_count = dim * dim;

  *a = gsl_matrix_alloc(unknown_count, unknown_count);
  *x = gsl_vector_alloc(unknown_count);
  *b = gsl_vector_alloc(unknown_count);
}


static void free_linear_system
(gsl_matrix* a, gsl_vector* x, gsl_vector* b)
{
  gsl_matrix_free(a);
  gsl_vector_free(x);
  gsl_vector_free(b);
}

static inline void print_double(double value)
{
  if (value >= 0.f)
    printf(" ");
  printf("%g ", value);
}

static void __attribute__((unused)) print_linear_system
(const gsl_matrix* a, const gsl_vector* x, const gsl_vector* b)
{
  size_t i;
  size_t j;

  for (i = 0; i < a->size1; ++i)
  {
    for (j = 0; j < a->size2; ++j)
      print_double(gsl_matrix_get(a, i, j));

    printf("\t");
    print_double(gsl_vector_get(x, i));

    printf("\t");
    print_double(gsl_vector_get(b, i));

    printf("\n");
  }
}

static void print_vector(const gsl_vector* v)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
  {
    print_double(gsl_vector_get(v, i));
    printf("\n");
  }
}

int main(int ac, char** av)
{
  grid_t g;
  gsl_matrix* a;
  gsl_vector* x;
  gsl_vector* b;

  gsl_permutation* p;
  int s;
  int i;

  grid_init_once(&g);

  alloc_linear_system(&g, &a, &x, &b);
  generate_a(&g, a);

  p = gsl_permutation_alloc(a->size1);
  gsl_linalg_LU_decomp(a, p, &s);

  for (i = 1; i < ac; ++i)
  {
    grid_load_string(&g, av[i]);
    generate_b(&g, b);
    gsl_linalg_LU_solve(a, p, b, x);
    print_vector(x);
    printf("---\n");
  }

  gsl_permutation_free(p);
  free_linear_system(a, x, b);

  return 0;
}
