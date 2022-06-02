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
  double *intersectionsP, *intersectionsQ;
  int countNonDegenIntP, countNonDegenIntQ, *initLabelsP, *initLabelsQ;
  vertex *tmpVertex, *current;
  countIntersections(polyPX, polyPY, polyQX, polyQY, PP[0].numVertices, QQ[0].numVertices, &countNonDegenIntP, &countNonDegenIntQ, &intersectionsP, &intersectionsQ, &initLabelsP, &initLabelsQ);
  // calculateIntersections(polyPX, polyPY, polyQX, polyQY, PP[0].numVertices, QQ[0].numVertices);
  
  // Polygon P: (PP)insert intersection vertices and change alpha value in the degenerate cases
  i=0;
  cout << "&& " << PP[0].root->next->p.x << "," << PP[0].root->prev->p.y << endl;
  // for (vertex* V : PP[0].vertices(ALL)){
  vertex* V=PP[0].root;
  do{
    current=V;
    while(*(intersectionsP+i)!=V->p.x || *(intersectionsP+i+1)!=V->p.y){
      tmpVertex=new vertex(*(intersectionsP+i), *(intersectionsP+i+1));
      tmpVertex->alpha=*(intersectionsP+i+2);
      tmpVertex->label=(IntersectionLabel)(*(initLabelsP+(i/3)));
      tmpVertex->next=current;
      current->prev->next=tmpVertex;
      tmpVertex->prev=current->prev;
      current->prev=tmpVertex;
      cout << i << " " << tmpVertex->p.x << " // " << tmpVertex->p.y << " " << tmpVertex->alpha << endl; 
      i+=3;
    }
    V->alpha=*(intersectionsP+i+2);
    V->label=(IntersectionLabel)(*(initLabelsP+(i/3)));
    cout << i << " " << V->p.x << " ** " << V->p.y << " " << V->alpha << endl;
    i+=3;
    V=current->next;
  }while(V->p.x!=PP[0].root->p.x || V->p.y!=PP[0].root->p.y);

  current=current->next;
  for(; i<countNonDegenIntP*3; i+=3){
    tmpVertex=new vertex(*(intersectionsP+i), *(intersectionsP+i+1));
    tmpVertex->alpha=*(intersectionsP+i+2);
    tmpVertex->label=(IntersectionLabel)(*(initLabelsP+(i/3)));
    tmpVertex->next=current;
    current->prev->next=tmpVertex;
    tmpVertex->prev=current->prev;
    current->prev=tmpVertex;
    cout << tmpVertex->p.x << " >> " << tmpVertex->p.y << " " << tmpVertex->alpha << endl;
  }

// Polygon Q: (QQ)insert intersection vertices and change alpha value in the degenerate cases
  i=0;
  cout << "&&&& " << QQ[0].root->next->p.x << "," << QQ[0].root->prev->p.y << endl;
  V=QQ[0].root;
  do{
    current=V;
    while(*(intersectionsQ+i)!=V->p.x || *(intersectionsQ+i+1)!=V->p.y){
      tmpVertex=new vertex(*(intersectionsQ+i), *(intersectionsQ+i+1));
      tmpVertex->alpha=*(intersectionsQ+i+2);
      tmpVertex->label=(IntersectionLabel)(*(initLabelsQ+(i/3)));
      tmpVertex->next=current;
      current->prev->next=tmpVertex;
      tmpVertex->prev=current->prev;
      current->prev=tmpVertex;
      cout << i << " " << tmpVertex->p.x << " // " << tmpVertex->p.y << " " << tmpVertex->alpha << endl; 
      i+=3;
    }
    V->alpha=*(intersectionsQ+i+2);
    V->label=(IntersectionLabel)(*(initLabelsQ+(i/3)));
    cout << i << " " << V->p.x << " ** " << V->p.y << " " << V->alpha << endl;
    i+=3;
    V=current->next;
  }while(V->p.x!=QQ[0].root->p.x || V->p.y!=QQ[0].root->p.y);

  current=current->next;
  for(; i<countNonDegenIntQ*3; i+=3){
    tmpVertex=new vertex(*(intersectionsQ+i), *(intersectionsQ+i+1));
    tmpVertex->alpha=*(intersectionsQ+i+2);
    tmpVertex->label=(IntersectionLabel)(*(initLabelsQ+(i/3)));
    tmpVertex->next=current;
    current->prev->next=tmpVertex;
    tmpVertex->prev=current->prev;
    current->prev=tmpVertex;
    cout << tmpVertex->p.x << " >> " << tmpVertex->p.y << " " << tmpVertex->alpha << endl;
  }
  
  cout << "\ncount degen " << countNonDegenIntP << endl;
  for(i=0; i<countNonDegenIntP*3; ++i){
      if(i%3==0)
        cout << "\n" << i/3;
    cout << " " << *(intersectionsP+i) << " ";
  }
  cout << "\nprint from PP" << endl;
  for (vertex* V : PP[0].vertices(ALL)){
    cout << V->p.x << ", " << V->p.y << " " << V->alpha << " " << V->label << endl;
  }

  cout << "\ncount degen " << countNonDegenIntQ << endl;
  for(i=0; i<countNonDegenIntQ*3; ++i){
    if(i%3==0)
      cout << "\n" << i/3;
    cout << " " << *(intersectionsQ+i) << " ";
  }
  cout << "\nprint from QQ" << endl;
  for (vertex* V : QQ[0].vertices(ALL)){
    cout << V->p.x << ", " << V->p.y << " " << V->alpha << " " << V->label << endl;
  }

  // phase 1
  // computeIntersections();  
  
  //*************** requires a barrier before starting labling *****************

  // // phase 2
  labelIntersections();

  // // phase 3
  createResult();
  
  // // post-processing
  cleanUpResult();
  
  // // write output polygon
  cout << "R "; savePolygon(RR,string(argv[argn]));
}