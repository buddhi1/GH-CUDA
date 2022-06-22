#include "readSCPolygon.h"
#include <iostream>

using namespace std;

// compile: g++ readPolygon.cpp
// run: ./a.out

void gpc_read_polygon(FILE *fp, gpc_polygon *p)
{
    int v;

    //fscanf(fp, "%d", &(p->num_contours));
  
    MALLOC(p->contour, 1
         * sizeof(gpc_vertex_list), "contour creation", gpc_vertex_list);
  
    fscanf(fp, "%d", &(p->contour[0].num_vertices));

    MALLOC(p->contour[0].vertex, p->contour[0].num_vertices
           * sizeof(gpc_vertex), "vertex creation", gpc_vertex);
           
    for (v= 0; v < p->contour[0].num_vertices; v++)
    {
      fscanf(fp, "%lf,%lf", &(p->contour[0].vertex[v].x),
                            &(p->contour[0].vertex[v].y));
     // printf(" %f %f, ", p->contour[0].vertex[v].x,
     //                    p->contour[0].vertex[v].y );
    }
}

int main()
{

    gpc_polygon subject, clip;
    FILE *sfp, *cfp;

    //uncomment for clipper comparison
    // Paths SClipper(1), CClipper(1), solution;
    // Clipper c;
    
    sfp= fopen("s.txt", "r");
    cfp= fopen("c.txt", "r");
    gpc_read_polygon(sfp, &subject);
    gpc_read_polygon(cfp, &clip);

    cout<<"Reading polygons done "<<endl;
    
    cout<<"Number of vertices in Subject polygon "<<subject.contour->num_vertices<<endl;
    cout<<"Number of vertices in Clip polygon "<<clip.contour->num_vertices<<endl;
        
// More detail here on intersection: 
// Clipper library: http://www.angusj.com/delphi/clipper.php   

    //uncomment for clipper comparison
     // gpc_polygon_to_path(&subject, &SClipper);
//      gpc_polygon_to_path(&clip, &CClipper);
//     
// 
//      c.AddPaths(SClipper, ptSubject, true);
//      c.AddPaths(CClipper, ptClip, true);
//      
//      c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);

}