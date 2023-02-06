#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/scan.h>
#include <cooperative_groups.h>

#include "../lib/constants.h"

typedef struct{
  double x, y;
} point;

__device__ double A(const point& P, const point& Q, const point& R){
	return (Q.x-P.x) * (R.y-P.y) - (Q.y-P.y) * (R.x-P.x);
}

// difference of two 2D points
__device__ point sub(const point& a, const point& b){ 
  point r;
  r.x=a.x-b.x;
  r.y=a.y-b.y;
  return r; 
}

// add two 2D points
__device__ point add(const point& a, const point& b){ 
  point r;
  r.x=a.x+b.x;
  r.y=a.y+b.y;
  return r; 
}

// multiply two 2D points
__device__ double mul(const point& a, const point& b){ 
  point r;
  r.x=a.x*b.x;
  r.y=a.y*b.y;
  return (r.x+r.y); 
}

// multiply scalar with 2D points
__device__ point mulScalar(const double c, const point& b){ 
  point r;
  r.x=c*b.x;
  r.y=c*b.y;
  return r; 
}

// find min
__device__ double getMin(double a, double b){
  if(a<b) return a;
  return b;
}

// find max
__device__ double getMax(double a, double b){
  if(a<b) return b;
  return a;
}

/*
-----------------------------------------------------------------
Function to returns the start index of the current id's 
intersections
Returns the intersection starting index
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
__device__ int getIntersectionStartIndex(int id, int *ps1){
  if(id==0) return 0;
  else return ps1[id];
}

/*
-----------------------------------------------------------------
Function to return intersection  type
Returns the type of the intersection
Runs in GPU
Called from Device
  NO_INTERSECTION, //0
  X_INTERSECTION,  //1
  T_INTERSECTION_Q, //2
  T_INTERSECTION_P, //3
  V_INTERSECTION, //4
  X_OVERLAP,      //5
  T_OVERLAP_Q,    //6
  T_OVERLAP_P,    //7
  V_OVERLAP       //8
-------------------------------------------------------------------
*/
__device__ int getIntersectType(
            const point& P1, const point& P2, 
            const point& Q1, const point& Q2,  
            double& alpha, double& beta){
	double AP1 = A(P1,Q1,Q2);
	double AP2 = A(P2,Q1,Q2);

	if (fabs(AP1-AP2) > EPSILON){
		// from here: [P1,P2] and [Q1,Q2] are not parallel
		// analyse potential intersection
		double AQ1 = A(Q1,P1,P2);
		double AQ2 = A(Q2,P1,P2);
		// compute alpha and beta
		alpha = AP1 / (AP1-AP2);
		beta  = AQ1 / (AQ1-AQ2);
		// classify alpha
		bool alpha_is_0 = false;
		bool alpha_in_0_1 = false;
		if ( (alpha > EPSILON) && (alpha < 1.0-EPSILON) )
			alpha_in_0_1 = true;
		else
			if (fabs(alpha) <= EPSILON)
				alpha_is_0 = true;
		// classify beta
		bool beta_is_0 = false;
		bool beta_in_0_1 = false;
		if ( (beta > EPSILON) && (beta < 1.0-EPSILON) )
			beta_in_0_1 = true;
		else
			if (fabs(beta) <= EPSILON)
				beta_is_0 = true;
		// distinguish intersection types
		if (alpha_in_0_1 && beta_in_0_1) return (1);  // return (X_INTERSECTION);
		if (alpha_is_0 && beta_in_0_1) return (2);    // return (T_INTERSECTION_Q);
		if (beta_is_0 && alpha_in_0_1) return (3);    // return (T_INTERSECTION_P);
		if (alpha_is_0 && beta_is_0) return (4);      // return (V_INTERSECTION);
	}else if (fabs(AP1) < EPSILON){
			// from here: [P1,P2] and [Q1,Q2] are collinear
			// analyse potential overlap			
      point dP = sub(P2, P1);
			point dQ = sub(Q2, Q1);
			point PQ = sub(Q1, P1);
			alpha = mul(PQ,dP) / mul(dP,dP);
			beta = -mul(PQ,dQ) / mul(dQ,dQ);
			// classify alpha
			bool alpha_is_0 = false;
			bool alpha_in_0_1 = false;
			bool alpha_not_in_0_1 = false;
			if ((alpha > EPSILON) && (alpha < 1.0-EPSILON))
				alpha_in_0_1 = true;
			else
				if (fabs(alpha) <= EPSILON)
					alpha_is_0 = true;
				else
					alpha_not_in_0_1 = true;
			// classify beta
			bool beta_is_0 = false;
			bool beta_in_0_1 = false;
			bool beta_not_in_0_1 = false;
			if ((beta > EPSILON) && (beta < 1.0-EPSILON))
				beta_in_0_1 = true;
			else
				if (fabs(alpha) <= EPSILON)
					beta_is_0 = true;
				else
					beta_not_in_0_1 = true;

			// distinguish intersection types
			if (alpha_in_0_1 && beta_in_0_1) return (5);      // return (X_OVERLAP);
			if (alpha_not_in_0_1 && beta_in_0_1) return (6);  // return (T_OVERLAP_Q);
			if (beta_not_in_0_1 && alpha_in_0_1) return (7);  // return (T_OVERLAP_P);
			if (alpha_is_0 && beta_is_0) return (8);          // return (V_OVERLAP);
		}
  return (0);	// return (NO_INTERSECTION); 
}

/*
-----------------------------------------------------------------
Function to get circular id of a given id 
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
__device__ int getCircularId(int id, int maxCount){
  if(maxCount==id) return 0;
  else if(id==-1) return maxCount-1;
  else return id;
}

/*
-----------------------------------------------------------------
Function to get relative position type
Runs in GPU
Called from Device
0 -> LEFT,
1 -> RIGHT,
2 -> IS_P_m,
3 -> IS_P_p
-------------------------------------------------------------------
*/
__device__ int oracle(int pMNId, int pPNId, int qId, const point& Q, const point& P1, const point& P2, const point& P3) {
  // is Q linked to P1 ?
  if(pMNId!=-100 && pMNId==qId) return 2;
  // is Q linked to P2 ?
  else if(pPNId!=-100 && pPNId==qId) return 3;
  // check relative position of Q with respect to chain (P1,P2,P3)
  double s1 = A(Q, P1, P2);
  double s2 = A(Q, P2, P3);
  double s3 = A(P1, P2, P3);
  if(s3>0){ 
    // chain makes a left turn
    if (s1>0 && s2>0)
      return 0;
    else
      return 1;
  }else{
    // chain makes a right turn (or is straight)
    if(s1<0 && s2<0)
      return 1;
    else
      return 0;
  }
}

/*
-----------------------------------------------------------------
Function to get initial classification label
Runs in GPU
Called from Device
Intersection Labels
0  NONE,
1  CROSSING,
2  BOUNCING,
3  LEFT_ON,
4  RIGHT_ON,
5  ON_ON,
6  ON_LEFT,
7  ON_RIGHT,
8  DELAYED_CROSSING,
9  DELAYED_BOUNCING
-------------------------------------------------------------------
*/
__device__ int getInitialLabel(int qMType, int qPType){
  // check non-overlapping cases
  if((qMType==0  && qPType==1)||(qMType==1 && qPType==0)){
    return 1;
  }
  if((qMType==0  && qPType==0)||(qMType==1 && qPType==1)){
    return 2;
  }
  // check overlapping cases
  if(((qPType==3) && (qMType==1))||((qMType==3) && (qPType==1))) return 3;
  if(((qPType==3) && (qMType==0))||((qMType==3) && (qPType==0))) return 4;
  if(((qPType==3) && (qMType==2))||((qMType==3) && (qPType==2))) return 5;
  if(((qMType==2) && (qPType==1))||((qPType==2) && (qMType==1))) return 6;
  if(((qMType==2) && (qPType==0))||((qPType==2) && (qMType==0))) return 7;
  else return -102;
}

/*
-----------------------------------------------------------------
Function to do counting sort of arr[] according to
  the digit represented by exp.
Returns sorted by single base digit
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
__device__ void gpuCountSort(int arr[], int tmpBucket[], int sortedIndicies[], int start, int end, int exp){
  int *output=tmpBucket; // used to track indices w.r.t original araay values
  int i, count[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  // Store count of occurrences in count[]
  for(i=start; i<end; i++){
    *(output+i)=sortedIndicies[i];
    count[(arr[*(output+i)] / exp) % 10]++;
  }
  // count prefix sum contains actual positions
  for(i=1; i<10; i++){
    count[i] += count[i - 1];
  }
  // Build the output array indices
  for(i=end-1; i>=start; i--){
    sortedIndicies[start+(count[(arr[*(output+i)] / exp) % 10]-1)]=*(output+i);
    count[(arr[*(output+i)] / exp) % 10]--;
  }
}

/*
-----------------------------------------------------------------
Function that sorts arr[] of size n using Radix Sort
Returns sorted array
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
__device__ void gpuRadixsort(int arr[], int tmpBucket[], int alphaSortedIndicies[], int start, int end){  
  // Do counting sort for every digit. Note that instead
  // of passing digit number, exp is passed. exp is 10^i
  // where i is current digit number
  int i, exp=1;
  for(i=start; i<end; i++){
      alphaSortedIndicies[i]=i;
  }
  for (i=1; i<=EPSILON_POSITIONS; i++){
    gpuCountSort(arr, tmpBucket, alphaSortedIndicies, start, end, exp);
    exp*=10;
  }
  // record sorted alpha values in tmpBucket
  for(i=start; i<end; ++i)
    tmpBucket[i]=arr[alphaSortedIndicies[i]];
}

/*
-----------------------------------------------------------------
Function to return vertex 2 of a given vertex 1
Returns index of vertex 2 
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
__device__ int gpuGetVertex2Index(int vertex1Index, int polySize[], int polyId){
  if(vertex1Index<polySize[polyId+1]-1) return vertex1Index+1;
  else if(vertex1Index=polySize[polyId+1]-1) return polySize[polyId];
}

/*
-----------------------------------------------------------------
Function: iterative search 
Returns location of x in given array arr[l..r] if present,
  otherwise -1
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
__device__ int gpuSearchPolygonId(int arr[], int numPol, int x){
  for(int i=0; i<numPol; ++i){
    if(arr[i]<=x && arr[i+1]>x)
      return i;
  }
  return -1;
}

/*
-----------------------------------------------------------------
Function to check if there is a overlap between given 2 edges 
Returns 1 if there is a overlap; else 0
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
__device__ int gpuLSMF(point P1, point P2, point Q1, point Q2){
  double minPX=P1.x, minPY=P1.y;
  double maxPX=P2.x, maxPY=P2.y;
  double minQX=Q1.x, minQY=Q1.y;
  double maxQX=Q2.x, maxQY=Q2.y;
  // this toggle way optimizes this computation well compared to using 8 min max calls seperately
  if(minPX>P2.x){
    minPX=P2.x;
    maxPX=P1.x;
  }
  if(minPY>P2.y){
    minPY=P2.y;
    maxPY=P1.y;
  }
  if(minQX>Q2.x){
    minQX=Q2.x;
    maxQX=Q1.x;
  }
  if(minQY>Q2.y){
    minQY=Q2.y;
    maxQY=Q1.y;
  }
  // check intersection between MBRs
  if(minPX>maxQX || maxPX<minQX) return 0;
  if(minPY>maxQY || maxPY<minQY) return 0;
  return 1;
}

/*
-----------------------------------------------------------------
Function to check if edegs are intersecting with the CMBR
Return prefix sum arrays.
  if a marked boolean array if the edges are intersecting with it
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
__global__ void gpuCMBRFilter(
                double *polyX, double *polyY, 
                double cmbrMinX, double cmbrMinY, double cmbrMaxX, double cmbrMaxY,
                int size, int *boolPs, int *ps1, int *ps2){
  int id=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
  
  if(id>size) return;
  
  point P1, P2;
  P1.x=polyX[id];
  P1.y=polyY[id];
  P2.x=polyX[(id+1)%size];
  P2.y=polyY[(id+1)%size];

  double minX=getMin(P1.x, P2.x), minY=getMin(P1.y, P2.y);
  double maxX=getMax(P1.x, P2.x), maxY=getMax(P1.y, P2.y);

  boolPs[id]=1;
  ps1[id]=0;
  ps2[id]=1; //by default paren is in the list. Hence the initial value
  if(minX>cmbrMaxX || maxX<cmbrMinX) boolPs[id]=0;
  if(minY>cmbrMaxY || maxY<cmbrMinY) boolPs[id]=0;
  // if(boolPs[id]!=1) printf("/// %d\n", id);
}

/*
-----------------------------------------------------------------
Function to record all indicies which intersects with CMBR 
Return prefix sum arrays.
  index arrays
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
__global__ void gpuSaveCMBRIntersectedIndicies(
                double *polyX, double *polyY, 
                double cmbrMinX, double cmbrMinY, double cmbrMaxX, double cmbrMaxY,
                int size, int *boolPol, int *boolPs){
  int id=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
  
  if(id>size) return;
  
  point P1, P2;
  P1.x=polyX[id];
  P1.y=polyY[id];
  P2.x=polyX[(id+1)%size];
  P2.y=polyY[(id+1)%size];

  double minX=getMin(P1.x, P2.x), minY=getMin(P1.y, P2.y);
  double maxX=getMax(P1.x, P2.x), maxY=getMax(P1.y, P2.y);

  int intersect=1;
  if(minX>cmbrMaxX || maxX<cmbrMinX) intersect=0;
  if(minY>cmbrMaxY || maxY<cmbrMinY) intersect=0;
  if(intersect){
    boolPol[boolPs[id]]=id;
    // if(boolPs[id]!=id) printf("Error %d %d \n", id, boolPs[id]);
  }
}

/*
-----------------------------------------------------------------
Function to count all intersections. Simple bool check CMBR filter
Return prefix sum arrays.
  *prefix sum of count of all intersection vertices x2 (P and Q)
  *prefix sum of count of all intersection vertices excluding 
   degenerate cases x2 (P and Q)
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
__global__ void gpuCountIntersections(
                  double *polyPX, double *polyPY, 
                  double *polyQX,  double *polyQY, 
                  int sizeP, int sizeQ,
                  int *psP1, int *psP2, int *boolPIndex){
  int id=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
  int idx=threadIdx.x;  
  __shared__ double poly2X_shared[MAX_POLY2_SIZE+1], poly2Y_shared[MAX_POLY2_SIZE+1] /*+1 for halo next*/;
  double alpha;
  double beta;
  point I;
  int count1=0, count2=0, size=0, qid;
  point P1, P2, Q1, Q2;

  int tiles=(sizeQ+MAX_POLY2_SIZE-1)/MAX_POLY2_SIZE;
  int tileCellsPerThread=MAX_POLY2_SIZE/blockDim.x;
  if(id<sizeP){
    P1.x = polyPX[id];
    P1.y = polyPY[id];
    P2.x = polyPX[(id+1)%sizeP];
    P2.y = polyPY[(id+1)%sizeP];
  }
  for(int tileId=0; tileId<tiles; tileId++){
    size=MAX_POLY2_SIZE;
    qid=idx*SHARED_MEMORY_PADDING;
    if(tileId==tiles-1 && sizeQ%MAX_POLY2_SIZE!=0){
      size=sizeQ%MAX_POLY2_SIZE;
      qid=0;
    }
    for(int localId=0; localId<tileCellsPerThread; ++localId){
      if(tileId!=tiles-1 || (tileId==tiles-1 && idx<size)){
        // load data into shared memory collaboratively
        poly2X_shared[idx+(blockDim.x*localId)]=polyQX[idx+(blockDim.x*localId)+(tileId*MAX_POLY2_SIZE)];
        poly2Y_shared[idx+(blockDim.x*localId)]=polyQY[idx+(blockDim.x*localId)+(tileId*MAX_POLY2_SIZE)];
        if(tileId!=tiles-1 && idx==blockDim.x-1 && localId==tileCellsPerThread-1){
          poly2X_shared[idx+(blockDim.x*localId)+1]=polyQX[idx+(blockDim.x*localId)+1+(tileId*MAX_POLY2_SIZE)];
          poly2Y_shared[idx+(blockDim.x*localId)+1]=polyQY[idx+(blockDim.x*localId)+1+(tileId*MAX_POLY2_SIZE)];
        }
      }
    } 
    __syncthreads();
    if(boolPIndex[id]) 
    {
      for(int qCount=0; qCount<size; qid=((qid+1)%size), ++qCount){   
      // for(int qid=0; qid<size; qid++){  
        Q1.x = poly2X_shared[qid];
        Q1.y = poly2Y_shared[qid];
        // reset P2 vertex of last edge to first vertex
        if(tileId==tiles-1 && qid==size-1){
          Q2.x=polyQX[0];
          Q2.y=polyQY[0];
        }else{
          Q2.x=poly2X_shared[qid+1];
          Q2.y=poly2Y_shared[qid+1];
        }      
        // if MBRs of two edges does not have a CMBR, there cannot be any intersection at all
        // if(gpuLSMF(P1, P2, Q1, Q2))
        {
          // determine intersection or overlap type
          int i = getIntersectType(P1, P2, Q1, Q2, alpha, beta);
          if(i!=0){
            count1++;
            if(i==1 || i==3 || i==5 || i==7)
              count2++;
          }        
        }
      } 
    }
    __syncthreads();
  }
  if(id<sizeP){
    count2++; //represent the parent vertex 
    psP1[id]=count1;
    psP2[id]=count2; 
  }
}
__global__ void gpuNeighborMap(
                  double *polyPX, double *polyPY, 
                  double *polyQX,  double *polyQY, 
                  int sizeP, int sizeQ, 
                  int *psP1, int *psQ1, int *psQ2,
                  int *neighborMapQ){
  int id=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
  double alpha;
  double beta;
  point I;
  int count1=0, count2=0, nonDegenCount=0;

  if(id>=sizeQ) return;

  neighborMapQ[psQ2[id]+count2]=-100;  
  // check if the current edge has any intersections. If not return
  // printf("id %d %d %d \n", id, psQ1[id], psQ1[id+1]);
  // CMBR filter: check if the edge intersect with CMBR (from boolPIndex)
  // prefix sum filter: check if the current edge has any intersection count
  if(psQ1[id+1]!=psQ1[id])
  {
    point P1, P2, Q1, Q2;

    P1.x = polyQX[id];
    P1.y = polyQY[id];
    P2.x = polyQX[(id+1)%sizeQ];
    P2.y = polyQY[(id+1)%sizeQ];

    for(int qid=0; qid<sizeP; qid++){        
      // prefix sum filter: check if the current edge has any intersection count      
      if(psP1[qid+1]!=psP1[qid])
      {
        Q1.x = polyPX[qid];
        Q1.y = polyPY[qid];
        Q2.x = polyPX[(qid+1)%sizeP];
        Q2.y = polyPY[(qid+1)%sizeP];

        // if(gpuLSMF(P1, P2, Q1, Q2))
        {
          // determine intersection or overlap type
          int i = getIntersectType(P1, P2, Q1, Q2, alpha, beta);
          if(i!=0){
            count1++;
            if((id<sizeP && (i==1 || i==3 || i==5 || i==7)) || (id>=sizeP && (i==1 || i==3 || i==5 || i==7))){
              nonDegenCount++;
              count2=nonDegenCount;
            }
            else if((id<sizeP && (i==2 || i==4 || i==6 || i==8)) || (id>=sizeP && (i==2 || i==4 || i==6 || i==8)))
              count2=0;

            neighborMapQ[psQ2[id]+count2]=qid;      
          }
        }
      }
    }
  } 
}

/*
-----------------------------------------------------------------
Function to calculate all intersections save them in the correct 
location using prefixsum arrays and make neighbor connections
Returns 
  *intersection arrays with orginal vertices in them x2 (P and Q)
  *neighbor arrays x2 (P and q)
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
__global__ void gpuCalculateIntersections(
                  double *polyPX, double *polyPY, 
                  double *polyQX, double *polyQY, 
                  int sizeP, int sizeQ, 
                  int *psP1, int *psP2, int *psQ1, int *psQ2, 
                  double *intersectionsP, double *intersectionsQ, double *intersectionsP2, double *intersectionsQ2,
                  int *alphaValuesP, int *alphaValuesQ, int *tmpBucketP, int *alphaSortedIndiciesP,
                  int *neighborP, int *neighborQ, int *neighborP2, int *neighborQ2,
                  int *neighborMapQ /*, int *boolPIndex, int *boolQIndex*/){
  int id=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
  double alpha;
  double beta;
  point I;
  int count1=0, count2=0, nonDegenCount=0, start, end, localI, neighborQId;

  if(id>=sizeP) return;

  point P1, P2, Q1, Q2;
  int pid=id;

  intersectionsP[psP2[pid]*2]=polyPX[pid];       //consider edge for the intersection array
  intersectionsP[psP2[pid]*2+1]=polyPY[pid];
  intersectionsP2[psP2[pid]*2]=polyPX[pid];       //consider edge for the intersection array
  intersectionsP2[psP2[pid]*2+1]=polyPY[pid];
  alphaValuesP[psP2[pid]]=-100;

  if(id<sizeQ){
    intersectionsQ[psQ2[pid]*2]=polyQX[pid];       //consider edge for the intersection array
    intersectionsQ[psQ2[pid]*2+1]=polyQY[pid];
    intersectionsQ2[psQ2[pid]*2]=polyQX[pid];       //consider edge for the intersection array
    intersectionsQ2[psQ2[pid]*2+1]=polyQY[pid];
  }

  // prefix sum filter: check if the current edge has any intersection count      
  if(psP1[id+1]!=psP1[id])
  // CMBR filter followed by prefix sum filter
  // if(boolPIndex[id] && psP1[id+1]!=psP1[id])
  {
    P1.x = polyPX[pid];
    P1.y = polyPY[pid];
    P2.x = polyPX[(pid+1)%sizeP];
    P2.y = polyPY[(pid+1)%sizeP];

    for(int qid=0; qid<sizeQ; qid++){
      // prefix sum filter: check if the current edge has any intersection count      
      if(psQ1[qid+1]!=psQ1[qid])
      // CMBR filter followed by prefix sum filter
      // if(boolQIndex[qid] && psQ1[qid+1]!=psQ1[qid])
      {
        Q1.x = polyQX[qid];
        Q1.y = polyQY[qid];
        Q2.x = polyQX[(qid+1)%sizeQ];
        Q2.y = polyQY[(qid+1)%sizeQ];

        // if(gpuLSMF(P1, P2, Q1, Q2))
        {
          // determine intersection or overlap type
          int i = getIntersectType(P1, P2, Q1, Q2, alpha, beta);
          if(i){
            count1++;
            if(i==1 || i==3 || i==5 || i==7){
              nonDegenCount++;
              count2=nonDegenCount;
            }
            else if(i==2 || i==4 || i==6 || i==8)
              count2=0;
            start=psQ2[qid];
            end=psQ2[qid+1];

            if(i!=5){
              // local search to find the index of qid
              for(localI=start; localI<end; ++localI){
                if(pid==neighborMapQ[localI]){
                  neighborQId=localI;
                  neighborP[psP2[pid]+count2]=neighborQId+1;   //+1 acting as a padding and helps to identify 0 being empty 
                  neighborP2[psP2[pid]+count2]=neighborQId+1;   //+1 acting as a padding and helps to identify 0 being empty 
                  neighborQ[neighborQId]=psP2[pid]+count2+1;   //+1 acting as a padding and helps to identify 0 being empty 
                  neighborQ2[neighborQId]=psP2[pid]+count2+1;   //+1 acting as a padding and helps to identify 0 being empty 
                  localI=end+2; // break; 
                }
              }
            }else{
              neighborQId=start;
              neighborP[psP2[pid]+count2]=neighborQId+1;   //+1 acting as a padding and helps to identify 0 being empty 
              neighborP2[psP2[pid]+count2]=neighborQId+1;   //+1 acting as a padding and helps to identify 0 being empty 
              neighborQ[neighborQId]=psP2[pid]+count2+1;   //+1 acting as a padding and helps to identify 0 being empty 
              neighborQ2[neighborQId]=psP2[pid]+count2+1;
              
              for(localI=start; localI<end; ++localI){
                if(pid==neighborMapQ[localI]){
                  neighborQId=localI;
                  neighborP[psP2[pid]]=neighborQId+1;   //+1 acting as a padding and helps to identify 0 being empty 
                  neighborP2[psP2[pid]]=neighborQId+1;   //+1 acting as a padding and helps to identify 0 being empty 
                  neighborQ[neighborQId]=psP2[pid]+1;   //+1 acting as a padding and helps to identify 0 being empty 
                  neighborQ2[neighborQId]=psP2[pid]+1;   //+1 acting as a padding and helps to identify 0 being empty 
                  localI=end+2; // break; 
                }
              }
            }
            switch(i) {
              // case X_INTERSECTION:
              // I and I
              case 1:
                I = add(mulScalar((1.0-alpha), P1), mulScalar(alpha, P2));
                intersectionsP[(psP2[pid]+count2)*2]=I.x;       //consider edge for the intersection array
                intersectionsP[(psP2[pid]+count2)*2+1]=I.y;
                intersectionsP2[(psP2[pid]+count2)*2]=I.x;       //consider edge for the intersection array
                intersectionsP2[(psP2[pid]+count2)*2+1]=I.y;
                alphaValuesP[psP2[pid]+count2]=(int)pow(10, EPSILON_POSITIONS)*alpha;

                intersectionsQ[neighborQId*2]=I.x;       //consider edge for the intersection array
                intersectionsQ[neighborQId*2+1]=I.y;
                intersectionsQ2[neighborQId*2]=I.x;       //consider edge for the intersection array
                intersectionsQ2[neighborQId*2+1]=I.y;
                alphaValuesQ[neighborQId]=(int)pow(10, EPSILON_POSITIONS)*beta;
                break;
              // X-overlap
              // P1 and I(=P1 I is in Q)
              // I(=Q1 I is in P) and Q1
              case 5:
                intersectionsP[(psP2[pid]+count2)*2]=Q1.x;
                intersectionsP[(psP2[pid]+count2)*2+1]=Q1.y;
                intersectionsP2[(psP2[pid]+count2)*2]=Q1.x;
                intersectionsP2[(psP2[pid]+count2)*2+1]=Q1.y;
                alphaValuesP[psP2[pid]+count2]=(int)pow(10, EPSILON_POSITIONS)*alpha;

                intersectionsQ[neighborQId*2]=P1.x;    
                intersectionsQ[neighborQId*2+1]=P1.y;
                intersectionsQ2[neighborQId*2]=P1.x;    
                intersectionsQ2[neighborQId*2+1]=P1.y;
                alphaValuesQ[neighborQId]=(int)pow(10, EPSILON_POSITIONS)*beta;
                break;
              // case T_INTERSECTION_Q:
              // case T_OVERLAP_Q:
              // P1 and I(=P1 is in Q)
              case 2:
              case 6:
                alphaValuesP[psP2[pid]]=(int)pow(10, EPSILON_POSITIONS)*alpha;

                intersectionsQ[neighborQId*2]=P1.x;
                intersectionsQ[neighborQId*2+1]=P1.y;
                intersectionsQ2[neighborQId*2]=P1.x;
                intersectionsQ2[neighborQId*2+1]=P1.y;
                alphaValuesQ[neighborQId]=(int)pow(10, EPSILON_POSITIONS)*beta;
              break;
              // case T_INTERSECTION_P:
              // case T_OVERLAP_P:
              // I(=Q1 is in P) and Q1
              case 3:
              case 7:
                intersectionsP[(psP2[pid]+count2)*2]=Q1.x;
                intersectionsP[(psP2[pid]+count2)*2+1]=Q1.y;
                intersectionsP2[(psP2[pid]+count2)*2]=Q1.x;
                intersectionsP2[(psP2[pid]+count2)*2+1]=Q1.y;
                alphaValuesP[psP2[pid]+count2]=(int)pow(10, EPSILON_POSITIONS)*alpha;

                alphaValuesQ[psQ2[qid]]=(int)pow(10, EPSILON_POSITIONS)*beta;
                break;
              // case V_INTERSECTION:
              // case V_OVERLAP:
              // P1 and Q1
              case 4:
              case 8:
                alphaValuesP[psP2[pid]]=(int)pow(10, EPSILON_POSITIONS)*alpha;
                
                alphaValuesQ[psQ2[qid]]=(int)pow(10, EPSILON_POSITIONS)*beta;
                break;
            } 
          } 
        }
      }
    }
  
    // --------------------------------------------------------------------------------------------
    // local sort for each edge, start to end
    // --------------------------------------------------------------------------------------------
    start=psP2[pid];
    end=psP2[pid+1];
    // sort intersection vertices in this edge locally
    if((end-start)>2){
      gpuRadixsort(alphaValuesP, tmpBucketP, alphaSortedIndiciesP, start+1, end);
      // using sorted index array, change intersection locations in the array and neighbors
      // decending order JUST FOR TESING
      // for(int i=start+1, j=end-1; i<end; ++i, j--){
      // acending order of alpha values 
      for(int i=start+1, j=start+1; i<end; i++, j++){
        alphaValuesP[i]=tmpBucketP[j];
        intersectionsP[i*2]=intersectionsP2[alphaSortedIndiciesP[j]*2];
        intersectionsP[i*2+1]=intersectionsP2[alphaSortedIndiciesP[j]*2+1];
        neighborP[i]=neighborP2[alphaSortedIndiciesP[j]];
        neighborQ[neighborP2[alphaSortedIndiciesP[j]]-1]=i+1; //+1 is the padding. When reading do -1
        neighborQ2[neighborP2[alphaSortedIndiciesP[j]]-1]=i+1; //updates neighborQ2 as the new original to be used with sorted Q array
      } 
    } 
  // --------------------------------------------------------------------------------------------
  }
}

/*
-----------------------------------------------------------------
Function to save vertices of Q in edge wise sorted order
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
__global__ void gpuSortPolyQ(
                  int sizeQ, 
                  int *psQ2, 
                  double *intersectionsQ, double *intersectionsQ2,
                  int *alphaValuesQ, int *tmpBucketQ,  int *alphaSortedIndiciesQ,
                  int *neighborP, int *neighborQ, int *neighborQ2){
  int id=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;

  if(id<sizeQ){
    int start=psQ2[id], end=psQ2[id+1];
    // sort intersection vertices in this edge locally
    if((end-start)>2){
      gpuRadixsort(alphaValuesQ, tmpBucketQ, alphaSortedIndiciesQ, start+1, end);
      // using sorted index array, change intersection locations in the array and neighbors
      // decending order JUST FOR TESING
      // for(int i=start+1, j=end-1; i<end; ++i, j--){
      // acending order of alpha values 
      for(int i=start+1, j=start+1; i<end; i++, j++){
        alphaValuesQ[i]=tmpBucketQ[j];////////////////?????????????????????? need to swap alpha too!!!
        // (x,y,alpha) tuple change in sorted order
        intersectionsQ[i*2]=intersectionsQ2[alphaSortedIndiciesQ[j]*2];
        intersectionsQ[i*2+1]=intersectionsQ2[alphaSortedIndiciesQ[j]*2+1];
        //neighbor array update
        neighborQ[i]=neighborQ2[alphaSortedIndiciesQ[j]];
        neighborP[neighborQ2[alphaSortedIndiciesQ[j]]-1]=i+1; //+1 is the padding. When reading do -1 //[]= i+1
      } 
    } 
  }
}

/*
-----------------------------------------------------------------
Function to calculate initial label
Returns 
  *initial labels x2 (P and Q)
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
__global__ void gpuCalculateInitLabel(
                int sizeP, int *psP2,
                double *intersectionsP, double *intersectionsQ, int *alphaValuesP, 
                int *neighborP,
                int sizeNP, int sizeNQ, int *initLabelsP, int *initLabelsQ){
  int id=(blockIdx.y*gridDim.x+blockIdx.x)*blockDim.x+threadIdx.x;
  int pid=id;
  if(id>=sizeP) return;
  int start=psP2[pid], end=psP2[pid+1];
  // int start=psP2[id], end=psP2[id+1];
  int tmpId, nId, pMNId, pPNId;
  point pM, pP, qM, qP, current;
  int qMType, qPType, tmpIniLabel;
  int i;
  for(i=start; i<end; i++){
    initLabelsP[i]=-100;
    if(alphaValuesP[i]!=-100){    //consider intersections only
      current.x=intersectionsP[i*2]; 
      current.y=intersectionsP[i*2+1]; 
      tmpId=getCircularId(i-1, sizeNP);
      // determine local configuration at this intersection vertex
      pM.x=intersectionsP[tmpId*2];                // P-, predecessor of I on P
      pM.y=intersectionsP[tmpId*2+1];                // P-, predecessor of I on P
      // if(intersectionsP[tmpId*2+2]!=-100)
      if(alphaValuesP[tmpId]!=-100)
        pMNId=neighborP[tmpId]-1; //get neighbor id of P_m vertex
      else pMNId=-100;

      tmpId=getCircularId(i+1, sizeNP);
      pP.x=intersectionsP[tmpId*2];                // P+, successor of I on P
      pP.y=intersectionsP[tmpId*2+1];                // P+, successor of I on P
      if(alphaValuesP[tmpId]!=-100)
        pPNId=neighborP[tmpId]-1; //get neighbor id of P_p vertex
      else pPNId=-100;

      // nId=getNeighborIndex(i, neighborMapP, neighborQ);
      nId=neighborP[i]-1;
      tmpId=getCircularId(nId-1, sizeNQ);
      qM.x=intersectionsQ[tmpId*2];     // Q-, predecessor of I on Q
      qM.y=intersectionsQ[tmpId*2+1];     // Q-, predecessor of I on Q
      qMType=oracle(pMNId, pPNId, tmpId, qM, pM, current, pP);

      tmpId=getCircularId(nId+1, sizeNQ);
      qP.x=intersectionsQ[tmpId*2];     // Q+, successor of I on P
      qP.y=intersectionsQ[tmpId*2+1];     // Q+, successor of I on P
      qPType=oracle(pMNId, pPNId, tmpId, qP, pM, current, pP);

      tmpIniLabel=getInitialLabel(qMType, qPType);
      initLabelsP[i]=tmpIniLabel;
      initLabelsQ[nId]=tmpIniLabel;
    }
  }
}

/*
-----------------------------------------------------------------
Function to count how many intersection points and prefix sums
Returns 
  *count of non degenerate vertices x2 (P and Q)
  *intersection points with non degenrate vertices included x2
  *neighbor map x2
  *neighbor arrays x2
  *initial labels x2
Neighbor of a vertex (assume index i) in P can be read in O(1) time using
  neighborQ[neighborMapP[i]]
  for Q
    neighborP[neighborMapQ[i]]
Runs in CPU
Called from Host
-------------------------------------------------------------------
*/
void calculateIntersections(
                  double *polyPX, double *polyPY, 
                  double *polyQX,  double *polyQY, 
                  int sizeP, int sizeQ, double *cmbr,
                  int *countNonDegenIntP, int *countNonDegenIntQ, 
                  double **intersectionsP, double **intersectionsQ, int **alphaValuesP, int **alphaValuesQ,
                  int **initLabelsP, int **initLabelsQ,
                  int **neighborP, int **neighborQ){
    double *dev_polyPX, *dev_polyPY, *dev_polyQX, *dev_polyQY;
    int *dev_psP1, *dev_psP2, *dev_psQ1, *dev_psQ2, *dev_boolPsPX, *dev_boolPsQX, *dev_boolPX, *dev_boolQX;
    int psP1[sizeP+1], psP2[sizeP+1], psQ1[sizeQ+1], psQ2[sizeQ+1];
    int boolPsPX[sizeP+1], boolPsQX[sizeQ+1];
    cudaEvent_t kernelStart0, kernelStart1, kernelStart12, kernelStart2, kernelStart3, kernelStart4, kernelStart5, kernelStart6, kernelStart7, kernelStart8;
    cudaEvent_t kernelStop0, kernelStop1, kernelStop12, kernelStop2, kernelStop3, kernelStop4, kernelStop5, kernelStop6, kernelStop7, kernelStop8;
    int countCMBRP,countCMBRQ, sum;

    // printf("cmbr %f %f %f %f\n",*(cmbr+0), *(cmbr+1), *(cmbr+2), *(cmbr+3));
    
    // Phase1: Count intersections in each block. Create prefix sums to find local locations in each thread 
    // Allocate memory in device 
    if(DEBUG_TIMING){
        cudaEventCreate(&kernelStart0);
        cudaEventCreate(&kernelStop0);
    }
    cudaMalloc((void **) &dev_polyPX, sizeP*sizeof(double));
    cudaMalloc((void **) &dev_polyPY, sizeP*sizeof(double));
    cudaMalloc((void **) &dev_polyQX, sizeQ*sizeof(double));
    cudaMalloc((void **) &dev_polyQY, sizeQ*sizeof(double));
    cudaMalloc((void **) &dev_psP1, (sizeP+1)*sizeof(int));
    cudaMalloc((void **) &dev_psP2, (sizeP+1)*sizeof(int));
    cudaMalloc((void **) &dev_psQ1, (sizeQ+1)*sizeof(int));
    cudaMalloc((void **) &dev_psQ2, (sizeQ+1)*sizeof(int));

    // cudaMalloc((void **) &dev_boolPX, sizeP*sizeof(int));
    // cudaMalloc((void **) &dev_boolQX, sizeQ*sizeof(int));
    cudaMalloc((void **) &dev_boolPsPX, (sizeP+1)*sizeof(int));
    cudaMalloc((void **) &dev_boolPsQX, (sizeQ+1)*sizeof(int));

    // Copy input vectors from host memory to GPU buffers.
    cudaMemcpy(dev_polyPX, polyPX, sizeP*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_polyPY, polyPY, sizeP*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_polyQX, polyQX, sizeQ*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_polyQY, polyQY, sizeQ*sizeof(double), cudaMemcpyHostToDevice);

    int blocksPerGrid=((sizeP+sizeQ) + xThreadPerBlock - 1) / xThreadPerBlock;
    int xBlocksPerGrid=(blocksPerGrid + yBlockPerGrid - 1) / yBlockPerGrid;

    int blocksPerGridQ=(sizeQ + xThreadPerBlock - 1) / xThreadPerBlock;
    int xBlocksPerGridQ=(blocksPerGridQ + yBlockPerGrid - 1) / yBlockPerGrid;
    int blocksPerGridP=(sizeP + xThreadPerBlock - 1) / xThreadPerBlock;
    int xBlocksPerGridP=(blocksPerGridP + yBlockPerGrid - 1) / yBlockPerGrid;
    
    // ******size_t number_of_blocks = N/threads_per_block + (size_t)(N % threads_per_block != 0);
    dim3 dimBlock(xThreadPerBlock, yThreadPerBlock, 1);
    dim3 dimGridP(xBlocksPerGridP, yBlockPerGrid, 1); 
    dim3 dimGridQ(xBlocksPerGridQ, yBlockPerGrid, 1); 


    // CMBR filter 

    if(DEBUG_TIMING) cudaEventRecord(kernelStart0);
    gpuCMBRFilter<<<dimGridP, dimBlock>>>(
                dev_polyPX, dev_polyPY, 
                cmbr[0], cmbr[1], cmbr[2], cmbr[3],
                sizeP, dev_boolPsPX, dev_psP1, dev_psP2);
    gpuCMBRFilter<<<dimGridQ, dimBlock>>>(
                dev_polyQX, dev_polyQY, 
                cmbr[0], cmbr[1], cmbr[2], cmbr[3],
                sizeQ, dev_boolPsQX, dev_psQ1, dev_psQ2);

    if(DEBUG_TIMING) cudaEventRecord(kernelStop0);

    if(DEBUG_TIMING) cudaEventSynchronize(kernelStop0);

    cudaDeviceSynchronize();
  
    if(DEBUG_INFO_PRINT){
      cudaMemcpy(&boolPsPX, dev_boolPsPX, (sizeP+1)*sizeof(int), cudaMemcpyDeviceToHost);
      cudaMemcpy(&boolPsQX, dev_boolPsQX, (sizeQ+1)*sizeof(int), cudaMemcpyDeviceToHost);
      // count how many edges overlap with CMBRs
      countCMBRP=0;
      for(int x=0; x<sizeP; ++x) if(boolPsPX[x]) countCMBRP++;
      printf("\nP overlap count with CMBR %d ",countCMBRP);
      countCMBRQ=0;
      for(int x=0; x<sizeQ; ++x) if(boolPsQX[x]) countCMBRQ++;
      printf("Q overlap count with CMBR %d \n\n",countCMBRQ);
    }

    if(DEBUG_TIMING){
        cudaEventCreate(&kernelStart1);
        cudaEventCreate(&kernelStop1);
    }

    if(DEBUG_TIMING) cudaEventRecord(kernelStart1);
    gpuCountIntersections<<<dimGridQ, dimBlock>>>(
          dev_polyQX, dev_polyQY, 
          dev_polyPX, dev_polyPY, 
          sizeQ, sizeP,
          dev_psQ1, dev_psQ2, dev_boolPsQX);
    
    if(DEBUG_TIMING) cudaEventRecord(kernelStop1);
    if(DEBUG_TIMING) cudaEventSynchronize(kernelStop1);


    if(DEBUG_TIMING){
        cudaEventCreate(&kernelStart12);
        cudaEventCreate(&kernelStop12);
    }
    if(DEBUG_TIMING) cudaEventRecord(kernelStart12);

    gpuCountIntersections<<<dimGridP, dimBlock>>>(
          dev_polyPX, dev_polyPY, 
          dev_polyQX, dev_polyQY, 
          sizeP, sizeQ,
          dev_psP1, dev_psP2, dev_boolPsPX);

    if(DEBUG_TIMING) cudaEventRecord(kernelStop12);

    cudaDeviceSynchronize();

    cudaFree(dev_boolPsPX);
    cudaFree(dev_boolPsQX);

    dim3 dimGrid2(xBlocksPerGrid, yBlockPerGrid, 1);

    cudaMemcpy(&psP1, dev_psP1, (sizeP+1)*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&psP2, dev_psP2, (sizeP+1)*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&psQ1, dev_psQ1, (sizeQ+1)*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&psQ2, dev_psQ2, (sizeQ+1)*sizeof(int), cudaMemcpyDeviceToHost);

    if(DEBUG_TIMING) cudaEventSynchronize(kernelStop12);
    cudaDeviceSynchronize();

    if(DEBUG_TIMING){
        cudaEventCreate(&kernelStart2);
        cudaEventCreate(&kernelStop2);
    }
    if(DEBUG_TIMING) cudaEventRecord(kernelStart2);
    thrust::exclusive_scan(thrust::host, psP1, psP1 + sizeP+1, psP1);   //sizeP location contains the total size of the count1
    thrust::exclusive_scan(thrust::host, psP2, psP2 + sizeP+1, psP2);
    thrust::exclusive_scan(thrust::host, psQ1, psQ1 + sizeQ+1, psQ1);   //sizeQ location contains the total size of the count1
    thrust::exclusive_scan(thrust::host, psQ2, psQ2 + sizeQ+1, psQ2);
    if(DEBUG_TIMING) cudaEventRecord(kernelStop2);

    if(DEBUG_TIMING) cudaEventSynchronize(kernelStop2);

    cudaDeviceSynchronize();

    //Phase2: NEW- Fill neighborMap
    int *dev_neighborMapQ;
    int *neighborMapQ;
    *countNonDegenIntP=psP2[sizeP];
    *countNonDegenIntQ=psQ2[sizeQ];

    if(DEBUG_INFO_PRINT){
      printf("Non-degen count P %d *****--- Q %d\n", *countNonDegenIntP-sizeP, *countNonDegenIntQ-sizeQ);
      printf("Intersection count P %d *****--- Q %d\n", psP1[sizeP], psQ1[sizeQ]);
    }

    dim3 dimGrid(xBlocksPerGrid, yBlockPerGrid, 1);

    neighborMapQ=(int *)malloc(*countNonDegenIntQ*sizeof(int));

    cudaMalloc((void **) &dev_neighborMapQ, *countNonDegenIntQ*sizeof(int));

    if(DEBUG_TIMING){
        cudaEventCreate(&kernelStart3);
        cudaEventCreate(&kernelStop3);
    }
    cudaMemcpy(dev_psP1, psP1, (sizeP+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_psP2, psP2, (sizeP+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_psQ1, psQ1, (sizeQ+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_psQ2, psQ2, (sizeQ+1)*sizeof(int), cudaMemcpyHostToDevice);

    if(DEBUG_TIMING) cudaEventRecord(kernelStart3);

    gpuNeighborMap<<<dimGridQ, dimBlock>>>(
            dev_polyPX, dev_polyPY, 
            dev_polyQX, dev_polyQY, 
            sizeP, sizeQ,  
            dev_psP1, dev_psQ1, dev_psQ2,
            dev_neighborMapQ);
    if(DEBUG_TIMING) cudaEventRecord(kernelStop3);
  
  if(DEBUG_TIMING) cudaEventSynchronize(kernelStop3);

  // Phase 3: Calcualte intersections and save them in the arrays. Make neighbor connections
  int countIntersections=psP1[sizeP];

  int *alphaSortedIndiciesP, *alphaSortedIndiciesQ;
  double *dev_intersectionsP, *dev_intersectionsQ, *dev_intersectionsP2, *dev_intersectionsQ2;
  int *dev_neighborP, *dev_neighborQ, *dev_neighborP2, *dev_neighborQ2;
  int *dev_initLabelsP, *dev_initLabelsQ;
  int *dev_alphaValuesP, *dev_alphaValuesQ, *dev_tmpBucketP, *dev_tmpBucketQ, *dev_alphaSortedIndiciesP, *dev_alphaSortedIndiciesQ;

  *intersectionsP=(double *)malloc(*countNonDegenIntP*2*sizeof(double));
  *intersectionsQ=(double *)malloc(*countNonDegenIntQ*2*sizeof(double));
  *alphaValuesP=(int *)malloc(*countNonDegenIntP*sizeof(int));
  *alphaValuesQ=(int *)malloc(*countNonDegenIntQ*sizeof(int));
  alphaSortedIndiciesP=(int *)malloc(*countNonDegenIntP*sizeof(int));
  alphaSortedIndiciesQ=(int *)malloc(*countNonDegenIntQ*sizeof(int));
  *initLabelsP=(int *)malloc(*countNonDegenIntP*sizeof(int));
  *initLabelsQ=(int *)malloc(*countNonDegenIntQ*sizeof(int));
  *neighborP=(int *)malloc(*countNonDegenIntP*sizeof(int));
  *neighborQ=(int *)malloc(*countNonDegenIntQ*sizeof(int));

  for(int i=0; i<*countNonDegenIntQ; ++i){
    *(*initLabelsQ+i)=-100;
    *(*alphaValuesQ+i)=-100;
  }

  cudaDeviceSynchronize();

  // Allocate memory in device 
  cudaMalloc((void **) &dev_intersectionsP, *countNonDegenIntP*2*sizeof(double));
  cudaMalloc((void **) &dev_intersectionsP2, *countNonDegenIntP*2*sizeof(double));
  cudaMalloc((void **) &dev_intersectionsQ, *countNonDegenIntQ*2*sizeof(double));
  cudaMalloc((void **) &dev_intersectionsQ2, *countNonDegenIntQ*2*sizeof(double));
  cudaMalloc((void **) &dev_alphaValuesP, *countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_alphaValuesQ, *countNonDegenIntQ*sizeof(int));

  cudaMalloc((void **) &dev_tmpBucketP, *countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_tmpBucketQ, *countNonDegenIntQ*sizeof(int));
  cudaMalloc((void **) &dev_alphaSortedIndiciesP, *countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_alphaSortedIndiciesQ, *countNonDegenIntQ*sizeof(int));

  cudaMalloc((void **) &dev_neighborP, *countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_neighborP2, *countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_neighborQ, *countNonDegenIntQ*sizeof(int));
  cudaMalloc((void **) &dev_neighborQ2, *countNonDegenIntQ*sizeof(int));

  cudaMemcpy(dev_alphaValuesQ, *alphaValuesQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyHostToDevice);

  
  if(DEBUG_TIMING){
    cudaEventCreate(&kernelStart4);
    cudaEventCreate(&kernelStop4);
  }

  if(DEBUG_TIMING) cudaEventRecord(kernelStart4);
  gpuCalculateIntersections<<<dimGridP, dimBlock>>>(
        dev_polyPX, dev_polyPY, 
        dev_polyQX, dev_polyQY, 
        sizeP, sizeQ, 
        dev_psP1, dev_psP2, dev_psQ1, dev_psQ2, 
        dev_intersectionsP, dev_intersectionsQ, dev_intersectionsP2, dev_intersectionsQ2,
        dev_alphaValuesP, dev_alphaValuesQ, dev_tmpBucketP, dev_alphaSortedIndiciesP,
        dev_neighborP, dev_neighborQ, dev_neighborP2, dev_neighborQ2,
        dev_neighborMapQ);
  if(DEBUG_TIMING) cudaEventRecord(kernelStop4);
  if(DEBUG_TIMING) cudaEventSynchronize(kernelStop4);

  cudaDeviceSynchronize();

  cudaFree(dev_polyPX);
  cudaFree(dev_polyPY);
  cudaFree(dev_polyQX);
  cudaFree(dev_polyQY);
  cudaFree(dev_neighborMapQ);
  cudaFree(dev_intersectionsP2);
  cudaFree(dev_tmpBucketP);
  cudaFree(dev_alphaSortedIndiciesP);
  cudaFree(dev_neighborP2);
  cudaFree(dev_psP1);
  cudaFree(dev_psQ1);

  if(DEBUG_TIMING){
    cudaEventCreate(&kernelStart5);
    cudaEventCreate(&kernelStop5);
  }
  if(DEBUG_TIMING) cudaEventRecord(kernelStart5);
  gpuSortPolyQ<<<dimGridQ, dimBlock>>>(
        sizeQ, 
        dev_psQ2, 
        dev_intersectionsQ, dev_intersectionsQ2,
        dev_alphaValuesQ, dev_tmpBucketQ,  dev_alphaSortedIndiciesQ,
        dev_neighborP, dev_neighborQ, dev_neighborQ2);
  if(DEBUG_TIMING) cudaEventRecord(kernelStop5);
  if(DEBUG_TIMING) cudaEventSynchronize(kernelStop5);

  cudaDeviceSynchronize();

  cudaFree(dev_psQ2);
  cudaFree(dev_intersectionsQ2);
  cudaFree(dev_tmpBucketQ);
  cudaFree(dev_alphaSortedIndiciesQ);
  cudaFree(dev_neighborQ2);

  // Phase4: Inital label classificaiton
  cudaMalloc((void **) &dev_initLabelsP, *countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_initLabelsQ, *countNonDegenIntQ*sizeof(int));

  cudaMemcpy(dev_initLabelsQ, *initLabelsQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyHostToDevice);
 
  // negative alpha values are not handled explicitly since they are original vertices
  // ******No need to copy alpha values since they are only used to sort edge wise******
  // cudaMemcpy(alphaSortedIndicies, dev_alphaSortedIndicies, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  if(DEBUG_TIMING){
    cudaEventCreate(&kernelStart6);
    cudaEventCreate(&kernelStop6);
  }

  if(DEBUG_TIMING) cudaEventRecord(kernelStart6);
  gpuCalculateInitLabel<<<dimGridP, dimBlock>>>(
      sizeP,  dev_psP2,
      dev_intersectionsP, dev_intersectionsQ, dev_alphaValuesP,
      dev_neighborP,
      *countNonDegenIntP, *countNonDegenIntQ, dev_initLabelsP, dev_initLabelsQ);
  if(DEBUG_TIMING) cudaEventRecord(kernelStop6);

  cudaMemcpy(*intersectionsP, dev_intersectionsP, *countNonDegenIntP*2*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(*intersectionsQ, dev_intersectionsQ, *countNonDegenIntQ*2*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(*neighborP, dev_neighborP, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*neighborQ, dev_neighborQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
 
  cudaMemcpy(*initLabelsP, dev_initLabelsP, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*initLabelsQ, dev_initLabelsQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*alphaValuesP, dev_alphaValuesP, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*alphaValuesQ, dev_alphaValuesQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  
  if(DEBUG_TIMING) cudaEventSynchronize(kernelStop6);
  
  cudaDeviceSynchronize();

  float kernelTiming0=0, kernelTiming1=0, kernelTiming12=0, kernelTiming2=0, kernelTiming3=0, kernelTiming4=0, kernelTiming5=0, kernelTiming6=0;
  if(DEBUG_TIMING){
    cudaEventElapsedTime(&kernelTiming0, kernelStart0, kernelStop0);
    cudaEventElapsedTime(&kernelTiming1, kernelStart1, kernelStop1);
    cudaEventElapsedTime(&kernelTiming12, kernelStart12, kernelStop12);
    cudaEventElapsedTime(&kernelTiming2, kernelStart2, kernelStop2);
    cudaEventElapsedTime(&kernelTiming3, kernelStart3, kernelStop3);
    cudaEventElapsedTime(&kernelTiming4, kernelStart4, kernelStop4);
    cudaEventElapsedTime(&kernelTiming5, kernelStart5, kernelStop5);
    cudaEventElapsedTime(&kernelTiming6, kernelStart6, kernelStop6);
    // printf("gpuCMBR kernel exe time(microsecond) %f\n", kernelTiming0*1000);
    // printf("gpuCountIntersections kernel exe time(microsecond) %f\n", kernelTiming1*1000);
    // printf("gpuCountIntersections2 kernel exe time(microsecond) %f\n", kernelTiming12*1000);
    // printf("prefixsum kernels exe time(microsecond) %f\n", kernelTiming2*1000);
    // printf("gpuNeighborMap kernel exe time(microsecond) %f\n", kernelTiming3*1000);
    // printf("gpuCalculateIntersections kernel exe time(microsecond) %f\n", kernelTiming4*1000);
    // printf("gpuSortPolyQ kernel exe time(microsecond) %f\n", kernelTiming5*1000);
    // printf("gpuCalculateInitLabel kernel exe time(microsecond) %f\n\n", kernelTiming6*1000);
    
    printf("%f, %f, %f, %f, %f, %f, ", (kernelTiming1*1000 + kernelTiming12*1000), 
          kernelTiming2*1000, kernelTiming3*1000, kernelTiming4*1000, 
          kernelTiming5*1000, kernelTiming6*1000);
  }

  int limitP=10;
  int limitQ=10;

  cudaFree(dev_psP2);
  cudaFree(dev_intersectionsP);
  cudaFree(dev_intersectionsQ);
  cudaFree(dev_alphaValuesP);
  cudaFree(dev_alphaValuesQ);
  cudaFree(dev_neighborP);
  cudaFree(dev_neighborQ);
  cudaFree(countNonDegenIntP);
  cudaFree(countNonDegenIntQ);
  cudaFree(dev_initLabelsP);
  cudaFree(dev_initLabelsQ);

  // cudaFree(dev_polyPX);
  // cudaFree(dev_polyPY);
  // cudaFree(dev_polyQX);
  // cudaFree(dev_polyQY);
}