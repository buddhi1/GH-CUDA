#include <iostream>
#include<bits/stdc++.h>

#include "lib/polyclip.cpp"
#include "ghcuda.h"

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

  // -------------------------------------------------------------------------------------------
  // PHASE:1 read input polygons
  // -------------------------------------------------------------------------------------------
  cout << "\nP "; loadPolygon(PP,string(argv[argn++]));
  cout <<   "Q "; loadPolygon(QQ,string(argv[argn++]));
  // -------------------------------------------------------------------------------------------
  
  // -------------------------------------------------------------------------------------------
  // arrays for polygon P and polygon Q
  // array format polyPX=[x_1, x_2, ...]
  // array format polyPY=[y_1, y_2, ...]
  // -------------------------------------------------------------------------------------------
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
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // how to traverse through edges
  // -------------------------------------------------------------------------------------------
  // for (edge E : PP[0].edges(SOURCE)){
  //   // polyP[i] = 
	//   cout << "--- " << setprecision (15) << E.one->p  << " " << E.two->p << endl;
	// }
  // cout << "**my count " << QQ[0].numVertices << endl;
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // PHASE:2 calculate intersections using GPU acceleration
  // -------------------------------------------------------------------------------------------
  double *intersectionsP, *intersectionsQ;
  int countNonDegenIntP, countNonDegenIntQ, *initLabelsP, *initLabelsQ, *neighborP, *neighborQ, *neighborMapP, *neighborMapQ;
  int *alphaValuesP, *alphaValuesQ;
  vertex *tmpVertex, *current;
  calculateIntersections(
      polyPX, polyPY, 
      polyQX, polyQY, 
      PP[0].numVertices, QQ[0].numVertices, 
      &countNonDegenIntP, &countNonDegenIntQ, 
      &intersectionsP, &intersectionsQ, &alphaValuesP, &alphaValuesQ,
      &initLabelsP, &initLabelsQ, 
      &neighborMapP, &neighborMapQ, &neighborP, &neighborQ);
  // -------------------------------------------------------------------------------------------
// return 0;
  // -------------------------------------------------------------------------------------------
  // Polygon P: (PP)insert intersection vertices and change alpha value in the degenerate cases
  // -------------------------------------------------------------------------------------------
  i=0;
  // cout << "&& " << PP[0].root->next->p.x << "," << PP[0].root->prev->p.y << endl;
  // for (vertex* V : PP[0].vertices(ALL)){
  vertex* V=PP[0].root;
  do{
    current=V;
    while(*(intersectionsP+i)!=V->p.x || *(intersectionsP+i+1)!=V->p.y){
      tmpVertex=new vertex(*(intersectionsP+i), *(intersectionsP+i+1));
      // tmpVertex->alpha=*(intersectionsP+i+2);
      tmpVertex->label=(IntersectionLabel)(*(initLabelsP+(i/2)));
      tmpVertex->source=false;
      tmpVertex->intersection=true;
      tmpVertex->next=current;
      current->prev->next=tmpVertex;
      tmpVertex->prev=current->prev;
      current->prev=tmpVertex;
      // cout << i << " " << tmpVertex->p.x << " // " << tmpVertex->p.y << " " << tmpVertex->intersection << endl; 
      i+=2;
    }
    // V->alpha=*(intersectionsP+i+2);
    V->label=(IntersectionLabel)(*(initLabelsP+(i/2)));
    // if(*(intersectionsP+i+2)!=-100){
    if(*(alphaValuesP+(i/2))!=-100){
      V->intersection=true;
    }
    // cout << i << " " << V->p.x << " ** " << V->p.y << " " << V->intersection << endl;
    i+=2;
    V=current->next;
  }while(V->p.x!=PP[0].root->p.x || V->p.y!=PP[0].root->p.y);

  current=current->next;
  for(; i<countNonDegenIntP*2; i+=2){
    tmpVertex=new vertex(*(intersectionsP+i), *(intersectionsP+i+1));
    // tmpVertex->alpha=*(intersectionsP+i+2);
    tmpVertex->label=(IntersectionLabel)(*(initLabelsP+(i/2)));
    tmpVertex->source=false;
    tmpVertex->intersection=true;
    tmpVertex->next=current;
    current->prev->next=tmpVertex;
    tmpVertex->prev=current->prev;
    current->prev=tmpVertex;
    // cout << tmpVertex->p.x << " >> " << tmpVertex->p.y << " " << tmpVertex->alpha << endl;
  }
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // Polygon Q: (QQ)insert intersection vertices and change alpha value in the degenerate cases
  // -------------------------------------------------------------------------------------------
  i=0;
  // cout << "&&&& " << QQ[0].root->next->p.x << "," << QQ[0].root->prev->p.y << endl;
  V=QQ[0].root;
  do{
    current=V;
    while(*(intersectionsQ+i)!=V->p.x || *(intersectionsQ+i+1)!=V->p.y){
      tmpVertex=new vertex(*(intersectionsQ+i), *(intersectionsQ+i+1));
      // tmpVertex->alpha=*(intersectionsQ+i+2);
      tmpVertex->label=(IntersectionLabel)(*(initLabelsQ+(i/2)));
      tmpVertex->source=false;
      tmpVertex->intersection=true;
      tmpVertex->next=current;
      current->prev->next=tmpVertex;
      tmpVertex->prev=current->prev;
      current->prev=tmpVertex;
      // cout << i << " " << tmpVertex->p.x << " // " << tmpVertex->p.y << " " << tmpVertex->alpha << endl; 
      i+=2;
    }
    // V->alpha=*(intersectionsQ+i+2);
    V->label=(IntersectionLabel)(*(initLabelsQ+(i/2)));
    // if(*(intersectionsQ+i+2)!=-100){ 
    if(*(alphaValuesQ+(i/2))!=-100){ 
      V->intersection=true;
    }
    // cout << i << " " << V->p.x << " ** " << V->p.y << " " << V->alpha << endl;
    i+=2;
    V=current->next;
  }while(V->p.x!=QQ[0].root->p.x || V->p.y!=QQ[0].root->p.y);

  current=current->next;
  for(; i<countNonDegenIntQ*2; i+=2){
    tmpVertex=new vertex(*(intersectionsQ+i), *(intersectionsQ+i+1));
    // tmpVertex->alpha=*(intersectionsQ+i+2);
    tmpVertex->label=(IntersectionLabel)(*(initLabelsQ+(i/2)));
    tmpVertex->source=false;
    tmpVertex->intersection=true;
    tmpVertex->next=current;
    current->prev->next=tmpVertex;
    tmpVertex->prev=current->prev;
    current->prev=tmpVertex;
    // cout << tmpVertex->p.x << " >> " << tmpVertex->p.y << " " << tmpVertex->alpha << endl;
  }
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // linking polygon P and Polygon Q with neighbor property
  // ******RULE: Each vertex will only have ONE NEIGHBOR
  // -------------------------------------------------------------------------------------------
  // cout << "\n\n-----------\n";
  V=PP[0].root;
  vertex *VQ=QQ[0].root;
  int j=0;
  i=0;
  do{
    if(*(neighborP+i)!=0){
      VQ=QQ[0].root;
      for(j=0; j<(*(neighborP+i)-1); ++j){
        VQ=VQ->next;
      }
      V->neighbour=VQ;
      VQ->neighbour=V;
    }
    V=V->next;
    ++i;
  }while(V->p.x!=PP[0].root->p.x || V->p.y!=PP[0].root->p.y);
  // cout << "\n-----------\n";
  // -------------------------------------------------------------------------------------------
  
  // -------------------------------------------------------------------------------------------
  // Test print to check PP and QQ updated with intersection points
  // -------------------------------------------------------------------------------------------
  // cout << "\ncount degen " << countNonDegenIntP << endl;
  // for(i=0; i<countNonDegenIntP*2; ++i){
  //     if(i%2==0)
  //       cout << "\n" << i/2;
  //   cout << " " << *(intersectionsP+i) << " ";
  // }
  // cout << "\nprint from PP" << endl;
  // for (vertex* V : PP[0].vertices(ALL)){
  //   if(V->intersection)
  //     cout << V->p.x << ", " << V->p.y << /*" " << V->alpha << */" **" << V->label << "** -> " << V->neighbour->p.x << ", " << V->neighbour->p.y << endl;
  //   else
  //     cout << V->p.x << ", " << V->p.y << /*" " << V->alpha << */" " << V->label << endl;
  // }

  // cout << "\ncount degen " << countNonDegenIntQ << endl;
  // for(i=0; i<countNonDegenIntQ*2; ++i){
  //   if(i%2==0)
  //     cout << "\n" << i/2;
  //   cout << " " << *(intersectionsQ+i) << " ";
  // }
  // cout << "\nprint from QQ" << endl;
  // for (vertex* V : QQ[0].vertices(ALL)){
  //   if(V->intersection)
  //     cout << V->p.x << ", " << V->p.y << " " << /*V->alpha << */" **" << V->label << "** -> " << V->neighbour->p.x << ", " << V->neighbour->p.y << endl;
  //   else
  //     cout << V->p.x << ", " << V->p.y << " " << /*V->alpha <<*/ " " << V->label << " " << V->intersection << endl;
  // }
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // PHASE: 3
  // -------------------------------------------------------------------------------------------
  labelIntersections();
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // PHASE: 4
  // -------------------------------------------------------------------------------------------
  createResult();
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // post-processing
  // -------------------------------------------------------------------------------------------
  cleanUpResult();
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // write output polygon
  // -------------------------------------------------------------------------------------------
  cout << "R "; savePolygon(RR,string(argv[argn]));
  // -------------------------------------------------------------------------------------------
}