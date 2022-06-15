#include <iostream>
#include<bits/stdc++.h>

#include "lib/polyclip.cpp"
#include "ghcuda.h"

int argn;

// handles polygons without holes
void regularPolygonHandler(int argc, char* argv[]){
  // check input parameters
  if (argc < 4) {
    cout << "insufficient number of parameters" << endl;
    exit(0);
  }
  
  argn = 1;
  if (string(argv[1]) == "-union") {
    cout << "\n!!! computing UNION instead of INTERSECTION !!!\n";
    UNION = true;
    argn++;
  }
  int i=0, j;

  int sizeP=0, sizeQ=0;

  // -------------------------------------------------------------------------------------------
  // Alternate 1 -> PHASE:1 read input polygons
  // -------------------------------------------------------------------------------------------
  double *polyPX;
  double *polyPY;
  double *polyQX;
  double *polyQY;

  FILE *pfile, *qfile;
  pfile=fopen(argv[argn++], "r");
  qfile=fopen(argv[argn++], "r");
  gpc_read_polygon(pfile, &polyPX, &polyPY, &sizeP, "PP");
  gpc_read_polygon(qfile, &polyQX, &polyQY, &sizeQ, "QQ");

  // -------------------------------------------------------------------------------------------
  // Alternate 2 -> PHASE:1 read input polygons
  // -------------------------------------------------------------------------------------------
  /*cout << "\nP "; loadPolygon(PP,string(argv[argn]));
  cout <<   "Q "; loadPolygon(QQ,string(argv[argn+1]));
  // -------------------------------------------------------------------------------------------
  
  // -------------------------------------------------------------------------------------------
  // Handle multiple components in PP and QQ
  // -------------------------------------------------------------------------------------------
  // Read input polygons 
  // -------------------------------------------------------------------------------------------


  // cout << "size " << sizeP << endl;
  // for(i=0; i<sizeP; ++i){
  //   cout << polyPX[i] << "," << polyPY[i] << endl;
  // }

  // return;

  // -------------------------------------------------------------------------------------------
  // arrays for polygon P and polygon Q
  // array format polyPX=[x_1, x_2, ...]
  // array format polyPY=[y_1, y_2, ...]
  // -------------------------------------------------------------------------------------------
  double polyPX[PP[0].numVertices];
  double polyPY[PP[0].numVertices];
  double polyQX[QQ[0].numVertices];
  double polyQY[QQ[0].numVertices];

  i=0;
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
	}*/
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
  j=0;
  i=0;
  do{
    if(*(neighborP+i)!=0){
      VQ=QQ[0].root;
      for(j=0; j<(*(neighborP+i)-1); ++j){
        VQ=VQ->next;
      }
      V->neighbour=VQ;
      VQ->neighbour=V;
      // if(i<35) cout << "neigh " << i << " " << j << " (" << V->p.x << "," << V->p.y << "-" << VQ->p.x << "," << VQ->p.y << ") " << V->label << endl;
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
  // // for(i=0; i<countNonDegenIntP*2; ++i){
  // for(i=0; i<12*2; ++i){
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
}

// handles polygons with holes
void multiComponentPolygonHandler(int argc, char* argv[]){
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
  int i=0, j;

  // -------------------------------------------------------------------------------------------
  // PHASE:1 read input polygons
  // -------------------------------------------------------------------------------------------
  cout << "\nP "; loadPolygon(PP,string(argv[argn++]));
  cout <<   "Q "; loadPolygon(QQ,string(argv[argn++]));
  // -------------------------------------------------------------------------------------------
  
  // -------------------------------------------------------------------------------------------
  // Handle multiple components in PP and QQ
  // -------------------------------------------------------------------------------------------
  // Read input polygons 
  // -------------------------------------------------------------------------------------------
  int sizeP=0, sizeQ=0, sizePP=PP.size(), sizeQQ=QQ.size();
  for(i=0; i<sizePP; ++i){
    sizeP+=PP[i].numVertices;
  }
  for(i=0; i<sizeQQ; ++i){
    sizeQ+=QQ[i].numVertices;
  }
  double polyPX[sizeP];
  double polyPY[sizeP];
  double polyQX[sizeQ];
  double polyQY[sizeQ];

  int sizesPP[sizePP+1], sizesQQ[sizeQQ+1];

  cout << sizeP << " sizes " << sizeQ << endl;
  i=0;
  sizesPP[0]=0;
  for(j=0; j<sizePP; ++j){
    sizesPP[j+1]=sizesPP[j]+PP[j].numVertices;
    // cout << j << " +++ " << sizesPP[j] << endl;
    // copy polygon P values
    for (vertex* V : PP[j].vertices(ALL)){
      polyPX[i] = V->p.x;
      polyPY[i++] = V->p.y;
      // cout << "--- " << setprecision (15) << polyPX[i-1] << "," << polyPY[i-1] << endl;
    }
  }
  i=0;
  sizesQQ[0]=0;
  for(j=0; j<sizeQQ; ++j){
    sizesQQ[j+1]=sizesQQ[j]+QQ[j].numVertices;
    // cout << j << " --- " << sizesQQ[j] << endl;
    // copy polygon Q values
    for (vertex* V : QQ[j].vertices(ALL)){
      polyQX[i] = V->p.x;
      polyQY[i++] = V->p.y;
      // cout << "+++ " << setprecision (15) << polyQX[i-1] << "," << polyQY[i-1] << endl;
    }
  }

  for(j=0; j<sizePP+1; ++j)
    cout << " PP " << sizesPP[j];
  cout << endl;
  for(j=0; j<sizeQQ+1; ++j)
    cout << " QQ " << sizesQQ[j];
  cout << endl;

  double *intersectionsP, *intersectionsQ;
  int countNonDegenIntArrayP[sizePP], countNonDegenIntArrayQ[sizeQQ], *initLabelsP, *initLabelsQ, *neighborP, *neighborQ, *neighborMapP, *neighborMapQ;
  int *alphaValuesP, *alphaValuesQ;
  vertex *tmpVertex, *current;
  vertex* V;

  calculateIntersectionsMultipleComponents(
      polyPX, polyPY, 
      polyQX, polyQY, 
      sizeP, sizeQ, sizesPP, sizesQQ, PP.size(), QQ.size(),
      countNonDegenIntArrayP, countNonDegenIntArrayQ, 
      &intersectionsP, &intersectionsQ, &alphaValuesP, &alphaValuesQ,
      &initLabelsP, &initLabelsQ, 
      &neighborMapP, &neighborMapQ, &neighborP, &neighborQ);
  // -------------------------------------------------------------------------------------------
// return 0;  
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // Polygon P: (PP)insert intersection vertices and change alpha value in the degenerate cases
  // -------------------------------------------------------------------------------------------
  i=0;
  // cout << "&& " << PP[0].root->next->p.x << "," << PP[0].root->prev->p.y << endl;
  // for (vertex* V : PP[0].vertices(ALL)){
  for(int polyId=0; polyId<PP.size(); ++polyId){
    V=PP[polyId].root;
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
    }while(V->p.x!=PP[polyId].root->p.x || V->p.y!=PP[polyId].root->p.y);

    current=current->next;
    for(; i<countNonDegenIntArrayP[polyId]*2; i+=2){
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
  }
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // Polygon Q: (QQ)insert intersection vertices and change alpha value in the degenerate cases
  // -------------------------------------------------------------------------------------------
  i=0;
  // cout << "&&&& " << QQ[0].root->next->p.x << "," << QQ[0].root->prev->p.y << endl;
  for(int polyId=0; polyId<QQ.size(); ++polyId){
    V=QQ[polyId].root;
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
      cout << i << "insert " << V->p.x << ", " << V->p.y << " " << V->alpha << endl;
      i+=2;
      V=current->next;
    }while(V->p.x!=QQ[polyId].root->p.x || V->p.y!=QQ[polyId].root->p.y);

    current=current->next;
    for(; i<countNonDegenIntArrayQ[polyId]*2; i+=2){
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
  }
  // -------------------------------------------------------------------------------------------
  cout << "\nprint from QQ" << endl;
  for(int polyId=0; polyId<sizeQQ; ++polyId){
    cout << "poly " << polyId << endl;
    for (vertex* V : QQ[polyId].vertices(ALL)){
      // if(V->intersection)
      //   cout << V->p.x << ", " << V->p.y << " " << /*V->alpha << */" **" << V->label << "** -> " << V->neighbour->p.x << ", " << V->neighbour->p.y << endl;
      // else
        cout << V->p.x << ", " << V->p.y << " " << /*V->alpha <<*/ " " << V->label << " " << V->intersection << endl;
    }
  }
  cout <<"-------end " << endl;
  // -------------------------------------------------------------------------------------------
  // linking polygon P and Polygon Q with neighbor property
  // ******RULE: Each vertex will only have ONE NEIGHBOR
  // -------------------------------------------------------------------------------------------
  // cout << "\n\n-----------\n";
  int polyId2, QQStart;
  vertex *VQ;
  i=0;
  for(int polyId=0; polyId<PP.size(); ++polyId){
    V=PP[polyId].root;
    do{
      if(*(neighborP+i)!=0){
        for(polyId2=0; countNonDegenIntArrayQ[polyId2]<(*(neighborP+i)-1); ++polyId2);
        // polyId2=0;
        VQ=QQ[polyId2].root;
        QQStart=0;
        if(polyId2>0) QQStart=countNonDegenIntArrayQ[polyId2-1];
        // cout << "=== " << (*(neighborP+i)-1) << " " << countNonDegenIntArrayQ[polyId2] << " " << QQStart << endl;
        for(j=QQStart; j<(*(neighborP+i)-1); ++j){
          VQ=VQ->next;
          if(VQ->p.x==QQ[polyId2].root->p.x && VQ->p.y==QQ[polyId2].root->p.y){
            VQ=QQ[++polyId2].root;
          }
        }
        V->neighbour=VQ;
        VQ->neighbour=V;
        if(V->p.x!=VQ->p.x || V->p.y!=VQ->p.y){
          cout << "neigh " << i << " " << j << " (" << V->p.x << "," << V->p.y << " <> " << VQ->p.x << "," << VQ->p.y << ") " << V->label << endl;

          //  cout << (V->p.x!=VQ->p.x) << " wrong Neigh " << i << " " << j << " "<<  (V->p.y!=VQ->p.y) << endl; 
        }
        // if(i<35) cout << "neigh " << i << " " << j << " (" << V->p.x << "," << V->p.y << " - " << VQ->p.x << "," << VQ->p.y << ") " << V->label << endl;
      }
      V=V->next;
      ++i;
    }while(V->p.x!=PP[polyId].root->p.x || V->p.y!=PP[polyId].root->p.y);
  }

  // -------------------------------------------------------------------------------------------
  // Test print to check PP and QQ updated with intersection points
  // -------------------------------------------------------------------------------------------
  // cout << "\ncount degen " << countNonDegenIntP << endl;
  // // for(i=0; i<countNonDegenIntP*2; ++i){
  // for(i=0; i<12*2; ++i){
  //     if(i%2==0)
  //       cout << "\n" << i/2;
  //   cout << " " << *(intersectionsP+i) << " ";
  // }
  cout << "\nprint from PP" << endl;
  for(int polyId=0; polyId<sizePP; ++polyId){
    for (vertex* V : PP[polyId].vertices(ALL)){
      if(V->intersection)
        cout << V->p.x << ", " << V->p.y << /*" " << V->alpha << */" **" << V->label << "** -> " << V->neighbour->p.x << ", " << V->neighbour->p.y << endl;
      else
        cout << V->p.x << ", " << V->p.y << /*" " << V->alpha << */" " << V->label << endl;
    }
  }

  // cout << "\ncount degen " << countNonDegenIntQ << endl;
  // for(i=0; i<countNonDegenIntQ*2; ++i){
  //   if(i%2==0)
  //     cout << "\n" << i/2;
  //   cout << " " << *(intersectionsQ+i) << " ";
  // }
  cout << "\nprint from QQ" << endl;
  for(int polyId=0; polyId<sizeQQ; ++polyId){
    for (vertex* V : QQ[polyId].vertices(ALL)){
      if(V->intersection)
        cout << V->p.x << ", " << V->p.y << " " << /*V->alpha << */" **" << V->label << "** -> " << V->neighbour->p.x << ", " << V->neighbour->p.y << endl;
      else
        cout << V->p.x << ", " << V->p.y << " " << /*V->alpha <<*/ " " << V->label << " " << V->intersection << endl;
    }
  }
  // -------------------------------------------------------------------------------------------
}

int main(int argc, char* argv[])
{
  regularPolygonHandler(argc, argv);
  // multiComponentPolygonHandler(argc, argv);

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