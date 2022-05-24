#include <iostream>
#include<bits/stdc++.h>

#include "lib/polyclip.cpp"
#include "ghcuda.h"


// CUDA kernel for intersection point calculation


int main(int argc, char* argv[])
{
  // check input parameters
  if (argc < 4) {
    cout << "insufficient number of parameters" << endl;
    exit(0);
  }
  
  int argn = 1;
  if (string(argv[1]) == "-union") {
    cout << "\n!!! computing UNION instead of INTERSECTION !!!\n";
    UNION = true;
    argn++;
  }

  // read input polygons
  cout << "\nP "; loadPolygon(PP,string(argv[argn++]));
  cout <<   "Q "; loadPolygon(QQ,string(argv[argn++]));
  
  // arrays for polygon P and polygon Q
  // array format polyPX=[x_1, x_2, ...]
  // array format polyPY=[y_1, y_2, ...]
  double polyPX[PP[0].numVertices];
  double polyPY[PP[0].numVertices];
  double polyQX[QQ[0].numVertices];
  double polyQY[QQ[0].numVertices];

  int i=0;
  // copy polygon P values
  for (vertex* V : PP[0].vertices(ALL)){
    polyPX[i] = V->p.x;
    polyPY[i++] = V->p.y;
	  // cout << "--- " << setprecision (15) << V->p.x << endl;
	}

  i=0;
  // copy polygon Q values
  for (vertex* V : QQ[0].vertices(ALL)){
    polyQX[i] = V->p.x;
    polyQY[i++] = V->p.y;
	  // cout << "--- " << setprecision (15) << V->p.x << endl;
	}

  // how to traverse through edges
  // for (edge E : PP[0].edges(SOURCE)){
  //   // polyP[i] = 
	//   cout << "--- " << setprecision (15) << E.one->p  << " " << E.two->p << endl;
	// }
  // cout << "**my count " << QQ[0].numVertices << endl;

  // testhello();
  countIntersections(polyPX, polyPY, polyQX, polyQY, PP[0].numVertices, QQ[0].numVertices);
  // calculateIntersections(polyPX, polyPY, polyQX, polyQY, PP[0].numVertices, QQ[0].numVertices);

  // phase 1
  // computeIntersections();  
  
  //*************** requires a barrier before starting labling *****************

  // // phase 2
  // labelIntersections();

  // // phase 3
  // createResult();
  
  // // post-processing
  // cleanUpResult();
  
  // // write output polygon
  // cout << "R "; savePolygon(RR,string(argv[argn]));
}