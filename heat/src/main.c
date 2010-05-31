/*
** Made by fabien le mentec <texane@gmail.com>
** 
** Started on  Mon May 31 16:51:53 2010 texane
** Last update Mon May 31 18:25:55 2010 texane
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>



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
  /* M[i,j] -> Xn */
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

static void init_grid(grid_t* g)
{
  g->dim = 7;

#define GRID_NORTH_SIDE 0
#define GRID_EAST_SIDE 1
#define GRID_SOUTH_SIDE 2
#define GRID_WEST_SIDE 3

  g->side_temps[GRID_NORTH_SIDE] = 20.f;
  g->side_temps[GRID_EAST_SIDE] = 20.f;
  g->side_temps[GRID_SOUTH_SIDE] = 30.f;
  g->side_temps[GRID_WEST_SIDE] = 25.f;
}

static inline double grid_get_side_temp
(const grid_t* g, unsigned int i)
{
  /* north, clockwise */
  return g->side_temps[i];
}

static void grid_get_neigh_indices
(size_t adim, const ij_pair_t* ijp, size_t* nis)
{
  /* ijp in the grid coord, nis in the matrix coord */

  ij_pair_t tmp = *ijp;
  const size_t index = ij_to_index(adim, ij_to_matrix(&tmp));

#define NEIGH_NORTH_SIDE 0
#define NEIGH_EAST_SIDE 1
#define NEIGH_SOUTH_SIDE 2
#define NEIGH_WEST_SIDE 3

  nis[NEIGH_NORTH_SIDE] = index - adim;
  nis[NEIGH_EAST_SIDE] = index + 1;
  nis[NEIGH_SOUTH_SIDE] = index + adim;
  nis[NEIGH_WEST_SIDE] = index - 1;
}

#if 0
static void get_neighbor_coords
(const grid_t* g, neigh_coord_t* c, unsigned int i, unsigned int j)
{
  /* assume i, j valid */

  static const int vals[4][2] =
    {
      /* x, y */
      { 0, -1},
      { 1,  0},
      { 0,  1},
      {-1,  0}
    };

  unsigned int k;

  for (k = 0; k < 4; ++k)
  {
    c[k].i = i + vals[k][1];
    c[k].j = j + vals[k][0];
  }
}

static void get_neighbor_indices
(const grid_t* g, neigh_index_t* n, unsigned int i, unsigned int j)
{
}
#endif

static void generate_linear_system
(const grid_t* g, gsl_matrix* a, gsl_vector* b)
{
  const unsigned int adim = g->dim - 2;
  const unsigned int nrows = g->dim - 2;
  const unsigned int ncols = g->dim - 2;

  size_t index;
  size_t i;
  size_t j;
  size_t k;
  ij_pair_t ijp;
  double temp;

  size_t nis[4];

  /* zero a and b */

  gsl_matrix_set_zero(a);
  gsl_vector_set_zero(b);

  /* border columns */

  for (i = 1; i < (nrows - 1); ++i)
  {
    /* first col, north, east, south */
    index = ij_to_index(adim, make_ij_pair(&ijp, i, 0));
    gsl_matrix_set(a, index, index, 4.f);
    grid_get_neigh_indices(adim, ij_to_grid(&ijp), nis);
    gsl_matrix_set(a, index, nis[NEIGH_NORTH_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_EAST_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_SOUTH_SIDE], -1.f);
    gsl_vector_set(b, index, grid_get_side_temp(g, GRID_WEST_SIDE));

    /* last col, north, west, south */
    index = ij_to_index(adim, make_ij_pair(&ijp, i, ncols - 1));
    gsl_matrix_set(a, index, index, 4.f);
    grid_get_neigh_indices(adim, ij_to_grid(&ijp), nis);
    gsl_matrix_set(a, index, nis[NEIGH_NORTH_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_WEST_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_SOUTH_SIDE], -1.f);
    gsl_vector_set(b, index, grid_get_side_temp(g, GRID_EAST_SIDE));
  }

  /* border rows */

  for (j = 1; j < (ncols - 1); ++j)
  {
    /* first row, east, south, west */
    index = ij_to_index(adim, make_ij_pair(&ijp, 0, j));
    gsl_matrix_set(a, index, index, 4.f);
    grid_get_neigh_indices(adim, ij_to_grid(&ijp), nis);
    gsl_matrix_set(a, index, nis[NEIGH_EAST_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_SOUTH_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_WEST_SIDE], -1.f);
    gsl_vector_set(b, index, grid_get_side_temp(g, GRID_NORTH_SIDE));

    /* last row, north, east, west */
    index = ij_to_index(adim, make_ij_pair(&ijp, nrows - 1, j));
    gsl_matrix_set(a, index, index, 4.f);
    grid_get_neigh_indices(adim, ij_to_grid(&ijp), nis);
    gsl_matrix_set(a, index, nis[NEIGH_NORTH_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_EAST_SIDE], -1.f);
    gsl_matrix_set(a, index, nis[NEIGH_WEST_SIDE], -1.f);
    gsl_vector_set(b, index, grid_get_side_temp(g, GRID_SOUTH_SIDE));
  }

  /* corner cases, north west first, clockwise */

  temp = grid_get_side_temp(g, GRID_NORTH_SIDE);
  temp += grid_get_side_temp(g, GRID_WEST_SIDE);
  index = ij_to_index(adim, make_ij_pair(&ijp, 0, 0));
  grid_get_neigh_indices(adim, ij_to_grid(&ijp), nis);
  gsl_matrix_set(a, index, nis[NEIGH_EAST_SIDE], -1.f);
  gsl_matrix_set(a, index, nis[NEIGH_SOUTH_SIDE], -1.f);
  gsl_vector_set(b, index, temp);

  temp = grid_get_side_temp(g, GRID_NORTH_SIDE);
  temp += grid_get_side_temp(g, GRID_EAST_SIDE);
  index = ij_to_index(adim, make_ij_pair(&ijp, 0, ncols - 1));
  grid_get_neigh_indices(adim, ij_to_grid(&ijp), nis);
  gsl_matrix_set(a, index, nis[NEIGH_WEST_SIDE], -1.f);
  gsl_matrix_set(a, index, nis[NEIGH_SOUTH_SIDE], -1.f);
  gsl_vector_set(b, index, temp);

  temp = grid_get_side_temp(g, GRID_SOUTH_SIDE);
  temp += grid_get_side_temp(g, GRID_EAST_SIDE);
  index = ij_to_index(adim, make_ij_pair(&ijp, nrows - 1, ncols - 1));
  grid_get_neigh_indices(adim, ij_to_grid(&ijp), nis);
  gsl_matrix_set(a, index, nis[NEIGH_WEST_SIDE], -1.f);
  gsl_matrix_set(a, index, nis[NEIGH_NORTH_SIDE], -1.f);
  gsl_vector_set(b, index, temp);

  temp = grid_get_side_temp(g, GRID_SOUTH_SIDE);
  temp += grid_get_side_temp(g, GRID_WEST_SIDE);
  index = ij_to_index(adim, make_ij_pair(&ijp, nrows - 1, 0));
  grid_get_neigh_indices(adim, ij_to_grid(&ijp), nis);
  gsl_matrix_set(a, index, nis[NEIGH_EAST_SIDE], -1.f);
  gsl_matrix_set(a, index, nis[NEIGH_NORTH_SIDE], -1.f);
  gsl_vector_set(b, index, temp);

  /* inner points */

  for (i = 1; i < nrows - 1; ++i)
  {
    gsl_matrix_set(a, i, i, 4.f);

    for (j = 1; j < ncols - 1; ++j)
    {
      make_ij_pair(&ijp, i, j);
      grid_get_neigh_indices(adim, ij_to_grid(&ijp), nis);

      make_ij_pair(&ijp, i, j);
      for (k = 0; k < 4; ++k)
	gsl_matrix_set(a, ij_to_index(adim, &ijp), nis[k], -1.f);

      gsl_vector_set(b, ij_to_index(adim, &ijp), 0.f);
    }
  }
}


static void solve_linear_system
(const gsl_matrix* a, const gsl_vector* b)
{
  /* todo */
}


static void alloc_linear_system
(const grid_t* g, gsl_matrix** a, gsl_vector** b)
{
  /* there are dim * dim unknowns */

  const size_t dim = g->dim - 2;
  const size_t unknown_count = dim * dim;

  *a = gsl_matrix_alloc(unknown_count, unknown_count);
  *b = gsl_vector_alloc(unknown_count);
}


static void free_linear_system
(gsl_matrix* a, gsl_vector* b)
{
  gsl_matrix_free(a);
  gsl_vector_free(b);
}


static void print_linear_system
(const gsl_matrix* a, const gsl_vector* b)
{
/*   gsl_matrix_fprintf(stdout, a, "%g"); */
  gsl_vector_fprintf(stdout, b, "%g");
}


int main(int ac, char** av)
{
  grid_t g;
  gsl_matrix* a;
  gsl_vector* b;

  init_grid(&g);
  alloc_linear_system(&g, &a, &b);
  generate_linear_system(&g, a, b);

  print_linear_system(a, b);
  solve_linear_system(a, b);
  free_linear_system(a, b);

  return 0;
}
