#include "BSP.h"
#include "conforming_mesh.h"
#include "extended_predicates.h"
#include "string.h"
#include <stdarg.h>
#include <time.h>

int vertex_compare(const void *void_v1, const void *void_v2) {
  const vertex_t *v1 = (vertex_t *)void_v1;
  const vertex_t *v2 = (vertex_t *)void_v2;
  const double dx = v1->coord[0] - v2->coord[0];
  const double dy = v1->coord[1] - v2->coord[1];
  const double dz = v1->coord[2] - v2->coord[2];
  return (4 * ((dx > 0) - (dx < 0)) + 2 * ((dy > 0) - (dy < 0)) +
          ((dz > 0) - (dz < 0)));
}

inline bool coincident_points(const vertex_t *a, const vertex_t *b) {
  return !vertex_compare(a->coord, b->coord);
}

void remove_duplicated_points(vertex_t **vertices_p, uint32_t *npts,
                              vertex_t *vrts_copy, uint32_t *map,
                              uint32_t *diff) {

  // Sorting vertices by coordinates lexicographic order (x,y,z)
  qsort(vrts_copy, *npts, sizeof(vertex_t), vertex_compare);

  // Count and memory position of duplicated vertices.
  uint32_t vrts_counter = 0;
  for (uint32_t i = 1; i < (*npts); i++) {
    if (coincident_points(vrts_copy + i - 1, vrts_copy + i))
      vrts_counter++;
    diff[i] = vrts_counter;
  }

  // Set original_index to follow vertices permutation.
  for (uint32_t i = 0; i < (*npts); i++)
    map[vrts_copy[i].original_index] = i;

  // Allocating memory to store uinque mesh vertices (vertices_p).
  *vertices_p = (vertex_t *)malloc(sizeof(vertex_t) * (*npts - vrts_counter));

  // Fill mesh vertices (vertices_p)
  memcpy(*vertices_p, vrts_copy, sizeof(vertex_t));
  for (uint32_t i = vrts_counter = 1; i < (*npts); i++)
    if (!coincident_points(vrts_copy + i - 1, vrts_copy + i))
      memcpy(*vertices_p + (vrts_counter++), vrts_copy + i, sizeof(vertex_t));

  // Update uinque vertices number.
  (*npts) = vrts_counter;
}

/// //////////////////////////////////////////////////////////////////////////////////////////

void read_nodes_and_constraints(double *coords_A, uint32_t npts_A,
                                uint32_t *tri_idx_A, uint32_t ntri_A,
                                vertex_t **vertices_p, uint32_t *npts,
                                uint32_t **tri_vertices_p, uint32_t *ntri,
                                bool verbose) {

  // Reading points coordinates.
  *npts = npts_A;
  *ntri = ntri_A;
  vertex_t *tmp = (vertex_t *)malloc(*npts * sizeof(vertex_t));
  uint32_t *diff = (uint32_t *)calloc(*npts, sizeof(uint32_t));
  uint32_t *map = (uint32_t *)malloc(*npts * sizeof(uint32_t));
  *tri_vertices_p = (uint32_t *)malloc(sizeof(uint32_t) * 3 * (*ntri));

  for (uint32_t i = 0; i < (*npts); i++) {
    tmp[i].coord[0] = coords_A[i * 3];
    tmp[i].coord[1] = coords_A[i * 3 + 1];
    tmp[i].coord[2] = coords_A[i * 3 + 2];
    tmp[i].original_index = i;
  }

  if (*npts > 1)
    remove_duplicated_points(vertices_p, npts, tmp, map, diff);
  if (verbose)
    printf("Using %u unique vertices\n", *npts);

  free(tmp);

  for (uint32_t i = 0, j = 0; i < (*ntri); j++) {

    const uint32_t i1 = tri_idx_A[j * 3];
    const uint32_t i2 = tri_idx_A[j * 3 + 1];
    const uint32_t i3 = tri_idx_A[j * 3 + 2];
    (*tri_vertices_p)[3 * i] = map[i1] - diff[map[i1]];
    (*tri_vertices_p)[3 * i + 1] = map[i2] - diff[map[i2]];
    (*tri_vertices_p)[3 * i + 2] = map[i3] - diff[map[i3]];

    const double *v1c = ((*vertices_p) + (*tri_vertices_p)[3 * i])->coord;
    const double *v2c = ((*vertices_p) + (*tri_vertices_p)[3 * i + 1])->coord;
    const double *v3c = ((*vertices_p) + (*tri_vertices_p)[3 * i + 2])->coord;

    if (misAlignment(v1c, v2c, v3c) == 0)
      (*ntri)--;
    else
      i++;
  }
  free(map);
  free(diff);

  if (verbose)
    printf("Using %u non-degenerate constraints\n", *ntri);
}

void read_nodes_and_constraints_twoInput(
    double *coords_A, uint32_t npts_A, uint32_t *tri_idx_A, uint32_t ntri_A,
    double *coords_B, uint32_t npts_B, uint32_t *tri_idx_B, uint32_t ntri_B,
    vertex_t **vertices_p, uint32_t *npts, uint32_t **tri_vertices_p,
    uint32_t *ntri, uint32_t **tri_group, bool verbose) {

  // Global numbers of points and triangles.
  *npts = npts_A + npts_B;
  *ntri = ntri_A + ntri_B;

  // Reading points coordinates.
  vertex_t *tmp = (vertex_t *)malloc(*npts * sizeof(vertex_t));
  uint32_t *diff = (uint32_t *)calloc(*npts, sizeof(uint32_t));
  uint32_t *map = (uint32_t *)malloc(*npts * sizeof(uint32_t));

  for (uint32_t i = 0; i < npts_A; i++) {
    tmp[i].coord[0] = coords_A[i * 3];
    tmp[i].coord[1] = coords_A[i * 3 + 1];
    tmp[i].coord[2] = coords_A[i * 3 + 2];
    tmp[i].original_index = i;
  }
  for (uint32_t i = 0; i < npts_B; i++) {
    tmp[i + npts_A].coord[0] = coords_B[i * 3];
    tmp[i + npts_A].coord[1] = coords_B[i * 3 + 1];
    tmp[i + npts_A].coord[2] = coords_B[i * 3 + 2];
    tmp[i + npts_A].original_index = i + npts_A;
  }

  if (*npts > 1)
    remove_duplicated_points(vertices_p, npts, tmp, map, diff);
  if (verbose)
    printf("Using %u unique vertices\n", *npts);

  free(tmp);

  // Reading triangle vertices indices.
  *tri_vertices_p = (uint32_t *)malloc(*ntri * 3 * sizeof(uint32_t));
  *tri_group = (uint32_t *)malloc(*ntri * sizeof(uint32_t));

  uint32_t i1, i2, i3, base_ntri_A = ntri_A;
  for (uint32_t ti = 0, j = 0; ti < (*ntri); j++) {
    if (ti < ntri_A) {
      i1 = tri_idx_A[3 * j];
      i2 = tri_idx_A[3 * j + 1];
      i3 = tri_idx_A[3 * j + 2];

      (*tri_vertices_p)[3 * ti] = map[i1] - diff[map[i1]];
      (*tri_vertices_p)[3 * ti + 1] = map[i2] - diff[map[i2]];
      (*tri_vertices_p)[3 * ti + 2] = map[i3] - diff[map[i3]];
      (*tri_group)[ti] = CONSTR_A;

      const double *v1c = (*vertices_p + (*tri_vertices_p)[3 * ti])->coord;
      const double *v2c = (*vertices_p + (*tri_vertices_p)[3 * ti + 1])->coord;
      const double *v3c = (*vertices_p + (*tri_vertices_p)[3 * ti + 2])->coord;
      if (misAlignment(v1c, v2c, v3c) == 0) {
        (*ntri)--;
        ntri_A--;
      } else
        ti++;
    } else {
      if (ti == ntri_A)
        j = 0;
      i1 = tri_idx_B[3 * j];
      i2 = tri_idx_B[3 * j + 1];
      i3 = tri_idx_B[3 * j + 2];

      i1 += npts_A;
      i2 += npts_A;
      i3 += npts_A;
      (*tri_vertices_p)[3 * ti] = map[i1] - diff[map[i1]];
      (*tri_vertices_p)[3 * ti + 1] = map[i2] - diff[map[i2]];
      (*tri_vertices_p)[3 * ti + 2] = map[i3] - diff[map[i3]];
      (*tri_group)[ti] = CONSTR_B;

      const double *v1c = (*vertices_p + (*tri_vertices_p)[3 * ti])->coord;
      const double *v2c = (*vertices_p + (*tri_vertices_p)[3 * ti + 1])->coord;
      const double *v3c = (*vertices_p + (*tri_vertices_p)[3 * ti + 2])->coord;
      if (misAlignment(v1c, v2c, v3c) == 0)
        (*ntri)--;
      else
        ti++;
    }
  }

  *tri_vertices_p =
      (uint32_t *)realloc(*tri_vertices_p, *ntri * 3 * sizeof(uint32_t));
  *tri_group = (uint32_t *)realloc(*tri_group, *ntri * sizeof(uint32_t));
  free(map);
  free(diff);

  if (verbose)
    printf("Using %u non-degenerate constraints\n", *ntri);
}

/// <summary>
/// Main function - Create a polyhedral mesh out of the input
/// Input may be made of either one or two models to be combined into a boolean
/// composition
/// </summary>
/// <param name="coords_A">Serialized coordinates of first model
/// vertices</param> <param name="npts_A">Number of first model vertices</param>
/// <param name="tri_idx_A">Serialized indexes of first model triangles</param>
/// <param name="ntri_A">Number of first model triangles</param>
/// <param name="coords_B">Serialized coordinates of second model
/// vertices</param> <param name="npts_B">Number of second model
/// vertices</param> <param name="tri_idx_B">Serialized indexes of second model
/// triangles</param> <param name="ntri_B">Number of second model
/// triangles</param> <param name="bool_opcode">Boolean operation (0 = no op, U
/// = union, D = difference, I = intersection</param> <param
/// name="verbose">Print useful info during the process</param> <param
/// resulting BSPcomplex structure</returns>
BSPcomplex *makePolyhedralMesh(double *coords_A, uint32_t npts_A,
                               uint32_t *tri_idx_A, uint32_t ntri_A,
                               double *coords_B, uint32_t npts_B,
                               uint32_t *tri_idx_B, uint32_t ntri_B,
                               char bool_opcode, bool verbose) {
  bool two_input = (bool_opcode != '0');

  if (verbose) {
    if (coords_B == NULL) {
      printf("\nResolve auto-intersections and/or repair.\n\n");
    } else {
      printf("\nBoolean operator: ");
      if (bool_opcode == 'U')
        printf(" union.\n\n");
      else if (bool_opcode == 'I')
        printf(" intersection.\n\n");
      else if (bool_opcode == 'D')
        printf(" difference.\n\n");
      else {
        ip_error("INVALID\n\n");
      }
    }
  }

  //--Initialization-----------------
  TetMesh *mesh = new TetMesh;
  Constraint *constraints = new Constraint;

  if (!two_input) {
    read_nodes_and_constraints(coords_A, npts_A, tri_idx_A, ntri_A,
                               &mesh->vertices, &mesh->num_vertices,
                               &constraints->tri_vertices,
                               &constraints->num_triangles, verbose);
    constraints->constr_group =
        (uint32_t *)calloc(constraints->num_triangles, sizeof(uint32_t));
  } else { // two input
    read_nodes_and_constraints_twoInput(
        coords_A, npts_A, tri_idx_A, ntri_A, coords_B, npts_B, tri_idx_B,
        ntri_B, &mesh->vertices, &mesh->num_vertices,
        &constraints->tri_vertices, &constraints->num_triangles,
        &constraints->constr_group, verbose);
  }

  if (mesh->num_vertices < 4)
    ip_error("Cannot mesh less than 4 vertices.");
  if (constraints->num_triangles < 1)
    ip_error("No non-degenerate constraints loaded.");

  // (free_mem)
  {
    free(coords_A);
    free(tri_idx_A);
    if (two_input) {
      free(coords_B);
      free(tri_idx_B);
    }
  }

  clock_t time0 = clock();
  clock_t time1 = time0;

  //--Delaunay-Insertion-----------------
  // Setup for following vertex permutation in tetrahedrize
  for (uint32_t i = 0; i < mesh->num_vertices; i++)
    mesh->vertices[i].original_index = i;

  // Create Delaunay tetrahedrization of the vertices
  mesh->tetrahedrize();

  // Align constraint vertices after vertex permutation
  if (mesh->vertices[2].original_index != 3) {
    for (uint32_t k = 0; k < 3 * constraints->num_triangles; k++)
      constraints->tri_vertices[k] =
          mesh->vertices[constraints->tri_vertices[k]].original_index;
  } else {
    for (uint32_t k = 0; k < 3 * constraints->num_triangles; k++) {
      uint32_t l = mesh->vertices[3].original_index;
      uint32_t c = constraints->tri_vertices[k];
      if (c == 2)
        constraints->tri_vertices[k] = l;
      else if (c == 3)
        constraints->tri_vertices[k] = 2;
      else if (c == l)
        constraints->tri_vertices[k] = 3;
    }
  }
  clock_t time2 = clock();
  if (verbose)
    printf("\tDelaunay insertion: %f s\n",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  //--Half-Edges-and-Virtual-Constraint----------------------
  half_edge_t *half_edges = (half_edge_t *)calloc(
      3 * constraints->num_triangles, sizeof(half_edge_t));
  fill_half_edges(constraints, half_edges);
  sort_half_edges(half_edges, 3 * constraints->num_triangles);

  uint32_t nvc = place_virtual_constraints(mesh, constraints, half_edges);
  if (verbose)
    printf("\t%u virtual constraints added\n", nvc);

  free(half_edges);

  clock_t time3 = clock();
  if (verbose)
    printf("\tHalf-edges: %f s\n", (double)(time3 - time2) / CLOCKS_PER_SEC);

  //--Map-Tetrahedra-Constraint-Intersections----------------
  // Initialize an array to map the tetrahedra improperly intersecated
  // by the constraints:
  //    the i-th element is a pointer (of tetrahedra index type) used to store
  //    the indices of the constraints that improper intersect the i-th
  //    tetrahedron,
  //    it points to an array (map[i] of lenght num_map[i]) whose elements
  //    are the indices of the constraints which improperly intersect the
  //    i-th tetrahedron.
  uint32_t *num_map = (uint32_t *)calloc(mesh->tet_num, sizeof(uint32_t));
  uint32_t **map = (uint32_t **)calloc(mesh->tet_num, sizeof(uint32_t *));

  // Initalize other 4 maps to memory the tet-faces that are partially or
  // completely overlap with a constraint.
  uint32_t *num_map_f0 = (uint32_t *)calloc(mesh->tet_num, sizeof(uint32_t));
  uint32_t **map_f0 = (uint32_t **)calloc(mesh->tet_num, sizeof(uint32_t *));
  uint32_t *num_map_f1 = (uint32_t *)calloc(mesh->tet_num, sizeof(uint32_t));
  uint32_t **map_f1 = (uint32_t **)calloc(mesh->tet_num, sizeof(uint32_t *));
  uint32_t *num_map_f2 = (uint32_t *)calloc(mesh->tet_num, sizeof(uint32_t));
  uint32_t **map_f2 = (uint32_t **)calloc(mesh->tet_num, sizeof(uint32_t *));
  uint32_t *num_map_f3 = (uint32_t *)calloc(mesh->tet_num, sizeof(uint32_t));
  uint32_t **map_f3 = (uint32_t **)calloc(mesh->tet_num, sizeof(uint32_t *));

  insert_constraints(mesh, constraints, 
                     num_map, map, 
                     num_map_f0, map_f0,
                     num_map_f1, map_f1,
                     num_map_f2, map_f2,
                     num_map_f3, map_f3);

  clock_t time4 = clock();
  if (verbose)
    printf("\tMap creation: %f s\n", (double)(time4 - time3) / CLOCKS_PER_SEC);

  double DEL_time = (double)(time4 - time0) / CLOCKS_PER_SEC;
  if (verbose)
    printf("TOTAL Delaunay + map: %f s\n\n", DEL_time);

  initFPU(); // From here on we need indirect predicates

  //-Init BSP with mesh and constraints---------------------------------------
  BSPcomplex *complex = new BSPcomplex(
      mesh, constraints,
      (const uint32_t **)map, num_map,
      (const uint32_t **)map_f0, num_map_f0,
      (const uint32_t **)map_f1, num_map_f1, 
      (const uint32_t **)map_f2, num_map_f2,
      (const uint32_t **)map_f3, num_map_f3);

  // Free the memory used by map(s) and num_map(s)
  for (uint64_t i = 0; i < mesh->tet_num; i++) {
    if (num_map[i])
      free(map[i]);
    if (num_map_f0[i])
      free(map_f0[i]);
    if (num_map_f1[i])
      free(map_f1[i]);
    if (num_map_f2[i])
      free(map_f2[i]);
    if (num_map_f3[i])
      free(map_f3[i]);
  }
  free(num_map);
  free(map);
  free(num_map_f0);
  free(map_f0);
  free(num_map_f1);
  free(map_f1);
  free(num_map_f2);
  free(map_f2);
  free(num_map_f3);
  free(map_f3);

  delete mesh;
  delete constraints;

  clock_t time5 = clock();
  if (verbose)
    printf("\tDelaunay -> Complex %lf s\n",
           (double)(time5 - time4) / CLOCKS_PER_SEC);
  if (verbose)
    printf("\tInitial cells: %lu\n", complex->cells.size());

  //-Subdivision----------------------------------------------------------------
  for (size_t i = 0; i < complex->cells.size(); /*ignore */) {
    if (complex->cells[i].constraints.size() > 0)
      complex->splitCell(i);
    else
      i++;
  }
  clock_t time6 = clock();
  if (verbose)
    printf("\tCell subdivision %f s\n",
           (double)(time6 - time5) / CLOCKS_PER_SEC);
  if (verbose)
    printf("\tFinal cells: %lu\n", complex->cells.size());

  //--Decide colour of GREY faces-----------------------------------------------
  for (size_t i = 0; i < complex->faces.size(); i++) {
    BSPface &face = complex->faces[i];
    if (face.colour == GREY)
      face.colour = complex->blackAB_or_white(i, bool_opcode != '0');
  }

  clock_t time7 = clock();
  if (verbose)
    printf("\tFind black faces %f s\n",
           (double)(time7 - time6) / CLOCKS_PER_SEC);

  //-Classification:intrenal/external cells-------------------------------------
  complex->constraintsSurface_complexPartition(bool_opcode != '0');

  clock_t time8 = clock();
  if (verbose)
    printf("\tInt-ext class. %f s\n", (double)(time8 - time7) / CLOCKS_PER_SEC);

  clock_t time9 = clock();
  double BSP_time = (double)(time9 - time4) / CLOCKS_PER_SEC;
  if (verbose) {
    printf("TOTAL BSP: %f s\n\n", BSP_time);
    printf("TOTAL time: %f s\n\n", BSP_time + DEL_time);
  }

  return complex;
}
