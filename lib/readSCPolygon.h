#ifndef __READ_SC_POLYGON_HEADER__
#define __READ_SC_POLYGON_HEADER__
#include<stdio.h>
#include<stdlib.h>

#define MALLOC(p, b, s, t) {if ((b) > 0) { \
                            p= (t*)malloc(b); if (!(p)) { \
                            fprintf(stderr, "gpc malloc failure: %s\n", s); \
                            exit(0);}} else p= NULL;}

#define FREE(p)            {if (p) {free(p); (p)= NULL;}}

struct gpc_vertex                     /* Polygon vertex structure          */
{
  double              x;            /* Vertex x component                */
  double              y;            /* vertex y component                */
};

struct gpc_vertex_list                      /* Vertex list structure             */
{
  int                 num_vertices; /* Number of vertices in list        */
  gpc_vertex         *vertex;       /* Vertex array pointer              */
};

struct gpc_polygon                      /* Polygon set structure             */
{
  int                 num_contours; /* Number of contours in polygon     */
  int                *hole;         /* Hole / external contour flags     */
  gpc_vertex_list    *contour;      /* Contour array pointer             */
};

void gpc_read_polygon(FILE *fp, gpc_polygon *p);

#endif