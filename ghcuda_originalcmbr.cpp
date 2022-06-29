#include <iostream>
#include<bits/stdc++.h>
#include <vector>
#include <chrono>

#include "lib/polyclip.cpp"
#include "ghcuda.h"

using namespace std::chrono;

#include "lib/readShapefile.cpp"


int argn;
string outputFile=string("results/Fig-def-R.poly");

// read from shape files
void readInputFromShapeFiles(double **polyPX, double **polyPY, double **polyQX, double **polyQY){
	int PPID=36; //ne_10m_ocean
  string inputShp1=string("../datasets/ne_10m_ocean.csv");
  loadPolygonFromShapeFile2(PPTmp, inputShp1, PPID+1);

  // string inputShp2=string("../datasets/continents.csv");
  // int QQID=521; //continents
	// loadPolygonFromShapeFile2(QQTmp, inputShp2, QQID+1);

  string inputShp2=string("../datasets/ne_10m_land.csv");
  int QQID=1; //ne_10m_land
	loadPolygonFromShapeFile2(QQTmp, inputShp2, QQID+1);

  if(DEBUG_INFO_PRINT){
    cout<<"Shape file1: "<<inputShp1<<" PPID: "<<PPID<<endl;;
    cout<<"Shape file2: "<<inputShp2<<" QQID: "<<QQID<<endl;
    cout << "PP Polygon size " << PPTmp[PPID].size;
    cout << " QQ Polygon size " << QQTmp[QQID].size << endl;
  }

  PP.push_back(PPTmp[PPID]);
  QQ.push_back(QQTmp[QQID]);

  *polyPX=(double *)malloc(PPTmp[PPID].size*sizeof(double));
  *polyPY=(double *)malloc(PPTmp[PPID].size*sizeof(double));
  *polyQX=(double *)malloc(QQTmp[QQID].size*sizeof(double));
  *polyQY=(double *)malloc(QQTmp[QQID].size*sizeof(double));
  // PPVertexPointers=new vertex*[PPTmp[PPID].size];


  int i=0;
  // copy polygon P values
  for (vertex* V : PP[0].vertices(ALL)){
    *(*polyPX+i) = V->p.x;
    *(*polyPY+i++) = V->p.y;
    // *(PPVertexPointers+i++)=V;
	  // cout << "--- " << setprecision (15) << polyPX[i-1] << ", " << polyPY[i-1] << endl;
	}
  if(DEBUG_INFO_PRINT) cout<<"PP Count "<<i;

  i=0;
  // copy polygon Q values
  for (vertex* V : QQ[0].vertices(ALL)){
    *(*polyQX+i) = V->p.x;
    *(*polyQY+i++) = V->p.y;

	  // cout << "--- " << setprecision (15) << V->p.x << endl;
	}
  if(DEBUG_INFO_PRINT) cout<<" QQ Count "<<i<<endl;
  cout<<PP[0].root->p.x<<" val \n";
  // vertex *xx=*(PPVertexPointers+0);
  // cout<<" v** "<<xx->p.x<<endl;
  // cout<<" v "<<(*PPVertexPointers+0)->p.x<<endl;
}

// get CMBR for PP and QQ
void getCMBR(double *cmbr){
  vector<double> PPMBR;
  vector<double> QQMBR;
  // double *cmbr; //minx, miny, maxx, maxy

  PPMBR=getMBR(PP[0]);
  QQMBR=getMBR(QQ[0]);
  cmbr[0]=max(PPMBR[0], QQMBR[0]);
  cmbr[1]=max(PPMBR[1], QQMBR[1]);
  cmbr[2]=min(PPMBR[2], QQMBR[2]);
  cmbr[3]=min(PPMBR[3], QQMBR[3]);

  cout<<"mbr1 "<<PPMBR[0]<<", "<<PPMBR[1]<<", "<<PPMBR[2]<<", "<<PPMBR[3]<<endl;
  cout<<"mbr2 "<<QQMBR[0]<<", "<<QQMBR[1]<<", "<<QQMBR[2]<<", "<<QQMBR[3]<<endl;
  cout<<"cmbr "<<cmbr[0]<<", "<<cmbr[1]<<", "<<cmbr[2]<<", "<<cmbr[3]<<endl;
}

// read polygons into vectors
void readPolygons(int argc, char* argv[], double **polyPX, double **polyPY, double **polyQX, double **polyQY){
  // check input parameters
  if (argc < 4) {
    readInputFromShapeFiles(polyPX, polyPY, polyQX, polyQY);
  }else{
    argn = 1;
    if (string(argv[1]) == "-union") {
      cout << "\n!!! computing UNION instead of INTERSECTION !!!\n";
      UNION = true;
      argn++;
    }
    // ******** MAKE SURE sizeP>sizeQ
    int sizeP=0, sizeQ=0, i=0;

    // -------------------------------------------------------------------------------------------
    // Alternate 1 -> PHASE:1 read input polygons from polygon XY format file
    // -------------------------------------------------------------------------------------------
    // /*
    FILE *pfile, *qfile;
    pfile=fopen(argv[argn++], "r");
    qfile=fopen(argv[argn++], "r");
    gpc_read_polygon(pfile, polyPX, polyPY, &sizeP, "PP");
    gpc_read_polygon(qfile, polyQX, polyQY, &sizeQ, "QQ");
    // */

    // -------------------------------------------------------------------------------------------
    // Alternate 2 -> PHASE:1 read input polygons from given XY format with , and ; seperators
    // -------------------------------------------------------------------------------------------
    /*
    cout << "\nP "; loadPolygon(PP,string(argv[argn++]));
    cout <<   "Q "; loadPolygon(QQ,string(argv[argn++]));
    int i=0;
    // copy polygon P values
    for (vertex* V : PP[0].vertices(ALL)){
      *(*polyPX+i) = V->p.x;
      *(*polyPY+i++) = V->p.y;
      // cout << "--- " << setprecision (15) << polyPX[i-1] << ", " << polyPY[i-1] << endl;
    }
    if(DEBUG_INFO_PRINT) cout<<"PP Count "<<i;

    i=0;
    // copy polygon Q values
    for (vertex* V : QQ[0].vertices(ALL)){
      *(*polyQX+i) = V->p.x;
      *(*polyQY+i++) = V->p.y;
      // cout << "--- " << setprecision (15) << V->p.x << endl;
    }
    if(DEBUG_INFO_PRINT) cout<<" QQ Count "<<i<<endl;
    // -------------------------------------------------------------------------------------------
  */
  outputFile=string(argv[argn]);
  }
}

// handles polygons without holes
void regularPolygonHandler(double *polyPX, double *polyPY, double *polyQX, double *polyQY){
  // PHASE:1 calculate CMBR
  // -------------------------------------------------------------------------------------------
  

  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // PHASE:2 calculate intersections using GPU acceleration
  // -------------------------------------------------------------------------------------------
  double *intersectionsP, *intersectionsQ;
  int countNonDegenIntP, countNonDegenIntQ, *initLabelsP, *initLabelsQ, *neighborP, *neighborQ;
  int *alphaValuesP, *alphaValuesQ;
  vertex *tmpVertex, *current;
  double *cmbr;
  cmbr=(double *)malloc(4*sizeof(double));
  getCMBR(cmbr);

  calculateIntersections(
      polyPX, polyPY, 
      polyQX, polyQY, 
      PP[0].size, QQ[0].size, cmbr, 
      &countNonDegenIntP, &countNonDegenIntQ, 
      &intersectionsP, &intersectionsQ, &alphaValuesP, &alphaValuesQ,
      &initLabelsP, &initLabelsQ, 
      &neighborP, &neighborQ);
  // -------------------------------------------------------------------------------------------
// return 0;
  // -------------------------------------------------------------------------------------------
  // Polygon P: (PP)insert intersection vertices and change alpha value in the degenerate cases
  // -------------------------------------------------------------------------------------------
  int i=0, j=0, pi=0;
  int intersectionPArrayMax=countNonDegenIntP*2;

  PPVertexPointers=new vertex*[countNonDegenIntP];
  // QQVertexPointers=new vertex*[countNonDegenIntQ];
// cout <<"conondegen "<< countNonDegenIntP<< endl;
  // cout << "&& " << PP[0].root->next->p.x << "," << PP[0].root->prev->p.y << endl;
  // for (vertex* V : PP[0].vertices(ALL)){
  vertex* V=PP[0].root;
  // do{
  for(int ii=0; ii<=PP[0].size; ii++){
    current=V;

    while(*(intersectionsP+(i%intersectionPArrayMax))!=V->p.x || *(intersectionsP+((i+1)%intersectionPArrayMax))!=V->p.y){
      tmpVertex=new vertex(*(intersectionsP+i), *(intersectionsP+i+1));
      // tmpVertex->alpha=*(intersectionsP+i+2);
      tmpVertex->label=(IntersectionLabel)(*(initLabelsP+(i/2)));
      tmpVertex->source=false;
      tmpVertex->intersection=true;
      tmpVertex->next=current;
      current->prev->next=tmpVertex;
      tmpVertex->prev=current->prev;
      current->prev=tmpVertex;
      PPVertexPointers[pi++]=tmpVertex;
      // cout << i << " " << tmpVertex->p.x << " // " << tmpVertex->p.y << " " << tmpVertex->intersection << endl; 
      i+=2;
    }
    if(ii<PP[0].size){ 
      PPVertexPointers[pi]=V;
      // if(i==12381*2) cout<<PPVertexPointers[pi]->p.x<<", "<<PPVertexPointers[pi]->p.y<<" >> "<<V->p.x<<", "<<V->p.y<<endl;
      pi++;
      // V->alpha=*(intersectionsP+i+2);
      V->label=(IntersectionLabel)(*(initLabelsP+(i/2)));
      // if(*(intersectionsP+i+2)!=-100){
      if(*(alphaValuesP+(i/2))!=-100){
        V->intersection=true;
      }
    }
    // cout << i << " " << V->p.x << " ** " << V->p.y << " " << V->intersection << endl;
    i+=2;
    V=current->next;
  }
  // }while(V->p.x!=PP[0].root->p.x || V->p.y!=PP[0].root->p.y);
  // current=current->next;
  // for(; i<countNonDegenIntP*2; i+=2){
  //   tmpVertex=new vertex(*(intersectionsP+i), *(intersectionsP+i+1));
  //   // tmpVertex->alpha=*(intersectionsP+i+2);
  //   tmpVertex->label=(IntersectionLabel)(*(initLabelsP+(i/2)));
  //   tmpVertex->source=false;
  //   tmpVertex->intersection=true;
  //   tmpVertex->next=current;
  //   current->prev->next=tmpVertex;
  //   tmpVertex->prev=current->prev;
  //   current->prev=tmpVertex;
  //   *(PPVertexPointers+pi++)=tmpVertex;
  //   // cout << tmpVertex->p.x << " >> " << tmpVertex->p.y << " " << tmpVertex->alpha << endl;
  // }
  // -------------------------------------------------------------------------------------------
// cout<<*(intersectionsP+12381*2)<<", "<<*(intersectionsP+12381*2+1)<<"\n\n*****************\n";
// vertex* cc=PP[0].root;
//   for(int c=0; c<12509; c++){
//     if(PPVertexPointers[c]->p.x!=cc->p.x && PPVertexPointers[c]->p.y!=cc->p.y)
//       cout<<c<<" > "<<PPVertexPointers[c]->p.x<<", "<<PPVertexPointers[c]->p.y<<" - "<<cc->p.x<<", "<<cc->p.y<<endl;
//     cc=cc->next;
//   }
// cout<<"\n*****************\n\n";

  // -------------------------------------------------------------------------------------------
  // Polygon Q: (QQ)insert intersection vertices and change alpha value in the degenerate cases
  // -------------------------------------------------------------------------------------------
  i=0;
  // cout << "&&&& " << QQ[0].root->next->p.x << "," << QQ[0].root->prev->p.y << endl;
  V=QQ[0].root;
  int intersectionQArrayMax=countNonDegenIntQ*2;
  // do{
  for(int ii=0; ii<=QQ[0].size; ++ii){
    current=V;
    while(*(intersectionsQ+(i%intersectionQArrayMax))!=V->p.x || *(intersectionsQ+((i+1)%intersectionQArrayMax))!=V->p.y){
      tmpVertex=new vertex(*(intersectionsQ+i), *(intersectionsQ+i+1));
      // tmpVertex->alpha=*(intersectionsQ+i+2);
      tmpVertex->label=(IntersectionLabel)(*(initLabelsQ+(i/2)));
      tmpVertex->source=false;
      tmpVertex->intersection=true;
      tmpVertex->next=current;
      current->prev->next=tmpVertex;
      tmpVertex->prev=current->prev;
      current->prev=tmpVertex;
      tmpVertex->neighbour=PPVertexPointers[(*(neighborQ+(i/2)))-1];
      PPVertexPointers[(*(neighborQ+(i/2)))-1]->neighbour=tmpVertex;
      if(ii<10) cout<<(*(neighborQ+(i/2)))-1<<" *** "<<PPVertexPointers[(*(neighborQ+(i/2)))-1]->p.x<<", "<<PPVertexPointers[(*(neighborQ+(i/2)))-1]->p.y<<" >> "<<tmpVertex->p.x<<", "<<tmpVertex->p.y<<endl;
    
      // cout << i << " " << tmpVertex->p.x << " // " << tmpVertex->p.y << " " << tmpVertex->alpha << endl; 
      i+=2;
    }
    if(ii<QQ[0].size){
      // V->alpha=*(intersectionsQ+i+2);
      V->label=(IntersectionLabel)(*(initLabelsQ+(i/2)));
      // if(*(intersectionsQ+i+2)!=-100){ 
      if(*(alphaValuesQ+(i/2))!=-100){ 
        V->intersection=true;
        V->neighbour=PPVertexPointers[(*(neighborQ+(i/2)))-1];
        PPVertexPointers[(*(neighborQ+(i/2)))-1]->neighbour=V;
        // if(ii<10) cout<<(*(neighborQ+(i/2)))-1<<" -- "<<PPVertexPointers[(*(neighborQ+(i/2)))-1]->p.x<<", "<<PPVertexPointers[(*(neighborQ+(i/2)))-1]->p.y<<" >> "<<V->p.x<<", "<<V->p.y<<endl;
      }
    }
    // cout << i << " " << V->p.x << " ** " << V->p.y << " " << V->alpha << endl;
    i+=2;
    V=current->next;
  }
  printf("Copying completed");
  // }while(V->p.x!=QQ[0].root->p.x || V->p.y!=QQ[0].root->p.y);

  // current=current->next;
  // for(; i<countNonDegenIntQ*2; i+=2){
  //   tmpVertex=new vertex(*(intersectionsQ+i), *(intersectionsQ+i+1));
  //   // tmpVertex->alpha=*(intersectionsQ+i+2);
  //   tmpVertex->label=(IntersectionLabel)(*(initLabelsQ+(i/2)));
  //   tmpVertex->source=false;
  //   tmpVertex->intersection=true;
  //   tmpVertex->next=current;
  //   current->prev->next=tmpVertex;
  //   tmpVertex->prev=current->prev;
  //   current->prev=tmpVertex;
  //   tmpVertex->neighbour=(*PPVertexPointers+ *(neighborQ+(i/2)-1));
  //   (*PPVertexPointers+ *(neighborQ+(i/2)-1))->neighbour=tmpVertex;
  //   count2++;
  //   // cout << tmpVertex->p.x << " >> " << tmpVertex->p.y << " " << tmpVertex->alpha << endl;
  // }

  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // linking polygon P and Polygon Q with neighbor property
  // ******RULE: Each vertex will only have ONE NEIGHBOR
  // -------------------------------------------------------------------------------------------
  //  cout << "\n\n-----------\n";
  // int count=0;
  // V=QQ[0].root;
  // j=0;
  // i=0;
  // do{
  //   if(*(neighborQ+i)!=0){
  //     j++;
  //     cout<<V->p.x<<", "<<V->p.y<<" "<< *(neighborQ+i)<<" i "<<i<<endl;
  //     if(j<10) cout<<PPVertexPointers[(*(neighborQ+(i/2)))-1]->p.x<<", "<<PPVertexPointers[(*(neighborQ+(i/2)))-1]->p.y<<" >> "<<V->p.x<<", "<<V->p.y<<endl;
  //     V->neighbour=PPVertexPointers[(*(neighborQ+(i/2)))-1];
  //     PPVertexPointers[(*(neighborQ+(i/2)))-1]->neighbour=V;
  //     // if(V->p.x != VQ->p.x &&  V->p.y != VQ->p.y )
  //       // cout << count++ <<" neigh " << i << " " << j << " (" << V->p.x << "," << V->p.y << " | " << VQ->p.x << "," << VQ->p.y << ") " << V->label << endl;
  //   }
  //   V=V->next;
  //   ++i;
  // }while(V->p.x!=QQ[0].root->p.x || V->p.y!=QQ[0].root->p.y);
  // cout << "\n-----------\n";
 
 
 
  // cout << "\n\n-----------\n";
  // int count=0;
  // V=QQ[0].root;
  // vertex *VQ=PP[0].root;
  // j=0;
  // i=0;
  // do{
  //   if(*(neighborQ+i)!=0){
  //     VQ=PP[0].root;
  //     for(j=0; j<(*(neighborQ+i)-1); ++j){
  //       VQ=VQ->next;
  //     }
  //     V->neighbour=VQ;
  //     VQ->neighbour=V;
  //     // if(V->p.x != VQ->p.x &&  V->p.y != VQ->p.y )
  //       // cout << count++ <<" neigh " << i << " " << j << " (" << V->p.x << "," << V->p.y << " | " << VQ->p.x << "," << VQ->p.y << ") " << V->label << endl;
  //   }
  //   V=V->next;
  //   ++i;
  // }while(V->p.x!=QQ[0].root->p.x || V->p.y!=QQ[0].root->p.y);
  // cout << "\n-----------\n";




  // cout << "\n\n-----------\n";
  // int count=0;
  // V=PP[0].root;
  // vertex *VQ=QQ[0].root;
  // j=0;
  // i=0;
  // do{
  //   if(*(neighborP+i)!=0){
  //     VQ=QQ[0].root;
  //     for(j=0; j<(*(neighborP+i)-1); ++j){
  //       VQ=VQ->next;
  //     }
  //     V->neighbour=VQ;
  //     VQ->neighbour=V;
  //     // if(V->p.x != VQ->p.x &&  V->p.y != VQ->p.y )
  //       // cout << count++ <<" neigh " << i << " " << j << " (" << V->p.x << "," << V->p.y << " | " << VQ->p.x << "," << VQ->p.y << ") " << V->label << endl;
  //   }
  //   V=V->next;
  //   ++i;
  // }while(V->p.x!=PP[0].root->p.x || V->p.y!=PP[0].root->p.y);
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

int main(int argc, char* argv[]){
  double *polyPX, *polyPY, *polyQX, *polyQY;

  readPolygons(argc, argv, &polyPX, &polyPY, &polyQX, &polyQY);
  high_resolution_clock::time_point start, end;

  if(DEBUG_TIMING) start = high_resolution_clock::now();

  regularPolygonHandler(polyPX, polyPY, polyQX, polyQY);
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
  if(DEBUG_TIMING) end = high_resolution_clock::now();
  // -------------------------------------------------------------------------------------------
  // post-processing
  // -------------------------------------------------------------------------------------------
  cleanUpResult();
  // -------------------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------------------
  // write output polygon
  // -------------------------------------------------------------------------------------------
  if(DEBUG_INFO_PRINT) {
    cout << "R ";
    savePolygon(RR, outputFile);
  }
  // -------------------------------------------------------------------------------------------

  if(DEBUG_TIMING){
    auto duration = duration_cast<microseconds>(end - start);
    cout << "All time in microseconds\nTime: Total : " << fixed
    << duration.count() << setprecision(10) << endl;
  }
}