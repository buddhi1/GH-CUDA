#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/scan.h>
#include <cooperative_groups.h>

#include "constants.h"

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
Function to serach neighbor of a given vertex
Returns the index of the neighbor
Runs in GPU
Called from Device
id: id of the vertex that need to get the neighbor
neighborMap: Own neighborMap. If vertex is from P, Map also from P
neighbor: Other polygon's neighor array. If Id from P, neighbor from Q
-------------------------------------------------------------------
*/
// __device__ int getNeighborIndex(int id, int *neighborMap, int *neighbor){
//   int neighborId=neighborMap[id];
//   if(neighborId!=-100) return neighbor[neighborId];
//   return -1; //no neighbor for this vertex
// }

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
Function to get a given double value within tolerance 
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
__device__ double getValueTolarence(double val){
  if(val<EPSILON)
    return 0.0;
  return val;
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
  if(vertex1Index<polySize[polyId]-1) return vertex1Index+1;
  else if(polyId==0) return 0; 
  else return polyId-1;
}

/*
-----------------------------------------------------------------
Function to count all intersections. 
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
                  int *psP1, int *psP2, int *psQ1, int *psQ2){
  int id=blockDim.x*blockIdx.x+threadIdx.x;
  double alpha;
  double beta;
  point I;
  int count1=0, count2=0, size=sizeQ;
  double *poly1X=polyPX, *poly1Y=polyPY, *poly2X=polyQX, *poly2Y=polyQY;

  if(id>=sizeP+sizeQ) return;
  point P1, P2, Q1, Q2;
  int pid=id;
  if(id>=sizeP){
    size=sizeP;
    poly1X=polyQX; 
    poly1Y=polyQY; 
    poly2X=polyPX;
    poly2Y=polyPY;
    pid=id-sizeP;
  }
  for(int qid=0; qid<size; qid++){
    P1.x = poly1X[pid];
    P1.y = poly1Y[pid];

    Q1.x = poly2X[qid];
    Q1.y = poly2Y[qid];
    Q2.x = poly2X[qid+1];
    Q2.y = poly2Y[qid+1];

    // reset P2 vertex of last edge to first vertex
    if(qid == size-1){
      Q2.x = poly2X[0];
      Q2.y = poly2Y[0];
    }
    //polygon1 is P and polygon2 is Q
    if(pid==id && pid==sizeP-1){
      P2.x = poly1X[0];
      P2.y = poly1Y[0];
      // printf("sp %d\n", pid);
    }else if(pid!=id && pid == sizeQ-1){ //polygon2 is P and polygon1 is Q
      P2.x = poly1X[0];
      P2.y = poly1Y[0];
      // printf("sp %d\n", pid);
    } else { //no need reset. Normal case
      P2.x = poly1X[pid+1];
      P2.y = poly1Y[pid+1];
    }

    //
    // determine intersection or overlap type
    //
    // IntersectionType i = intersect(edgeP, edgeQ, alpha, beta);
    int i = getIntersectType(P1, P2, Q1, Q2, alpha, beta);
    if(i!=0){
      count1++;
      // if((id<sizeP && (i==1 || i==3 || i==5 || i==7)) || (id>=sizeP && (i==1 || i==2 || i==5 || i==6)))
      if((id<sizeP && (i==1 || i==3 || i==5 || i==7)) || (id>=sizeP && (i==1 || i==3 || i==5 || i==7)))
        count2++;
    }
    // if(id==9 || id==0)
      printf("id %d count1 %d count2  %d P(%f,%f) Q(%f,%f) %d\n", id, count1, count2, P1.x, P1.y, Q1.x, Q1.y, i);

  }
  count2++; //represent the parent vertex 
  if(id<sizeP){
    psP1[pid]=count1;
    psP2[pid]=count2;
  } else{
    psQ1[pid]=count1;
    psQ2[pid]=count2;
  }

  // __syncthreads();
  // thrust::exclusive_scan(thrust::device, psP1, psP1 + sizeP+1, psP1);   //sizeP location contains the total size of the count1
  // thrust::exclusive_scan(thrust::device, psP2, psP2 + sizeP+1, psP2);
  // thrust::exclusive_scan(thrust::device, psQ1, psQ1 + sizeQ+1, psQ1);   //sizeQ location contains the total size of the count1
  // thrust::exclusive_scan(thrust::device, psQ2, psQ2 + sizeQ+1, psQ2);
  // printf("id %d count1 %d count2  %d (%f,%f) (%f,%f)\n", id, count1, count2, P1.x, P1.y, P2.x, P2.y);
  // __syncthreads();
  // if(id==0){
  //   printf("%d \n", sizeP);
  //   for(int ii=0; ii<sizeP; ++ii){
  //     // printf("%d *%d ", psP1[ii], psP2[ii]);
  //     printf("%d ", psP1[ii]);
  //   }
  //   printf("\nend\n");
  //   for(int ii=0; ii<sizeQ; ++ii){
  //     // printf("%d *%d ", psQ1[ii], psQ2[ii]);
  //     printf("%d ", psQ1[ii]);
  //   }
  //   printf("\nend\n");
  // }
}


/*
-----------------------------------------------------------------
Function to neighbor map intersections. 
Return prefix sum arrays.
  *neighbor map all intersection vertices x2 (P and Q)
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
__global__ void gpuNeighborMap(
                  double *polyPX, double *polyPY, 
                  double *polyQX,  double *polyQY, 
                  int sizeP, int sizeQ, 
                  int *psP2, int *psQ2,
                  int *neighborMapP, int *neighborMapQ){
  int id=blockDim.x*blockIdx.x+threadIdx.x;
  double alpha;
  double beta;
  point I;
  int count1=0, count2=0, size=sizeQ, nonDegenCount=0;
  double *poly1X=polyPX, *poly1Y=polyPY, *poly2X=polyQX, *poly2Y=polyQY;

  if(id>=sizeP+sizeQ) return;
  point P1, P2, Q1, Q2;
  int pid=id;
  if(id>=sizeP){
    size=sizeP;
    poly1X=polyQX; 
    poly1Y=polyQY; 
    poly2X=polyPX;
    poly2Y=polyPY;
    pid=id-sizeP;
    neighborMapQ[psQ2[pid]+count2]=-100;   
  }else{
    neighborMapP[psP2[pid]+count2]=-100;
  }
  for(int qid=0; qid<size; qid++){
    P1.x = poly1X[pid];
    P1.y = poly1Y[pid];

    Q1.x = poly2X[qid];
    Q1.y = poly2Y[qid];
    Q2.x = poly2X[qid+1];
    Q2.y = poly2Y[qid+1];

    // reset P2 vertex of last edge to first vertex
    if(qid == size-1){
      Q2.x = poly2X[0];
      Q2.y = poly2Y[0];
    }
    //polygon1 is P and polygon2 is Q
    if(pid==id && pid==sizeP-1){
      P2.x = poly1X[0];
      P2.y = poly1Y[0];
      // printf("sp %d\n", pid);
    }else if(pid!=id && pid == sizeQ-1){ //polygon2 is P and polygon1 is Q
      P2.x = poly1X[0];
      P2.y = poly1Y[0];
      // printf("sp %d\n", pid);
    } else { //no need reset. Normal case
      P2.x = poly1X[pid+1];
      P2.y = poly1Y[pid+1];
    }
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
      if(id<sizeP){
        // if(pid<35) printf("#PPPP# %d %d %d (%f,%f - %f,%f)\n", pid, psP2[pid]+count2, qid, P1.x, P1.y, Q1.x, Q1.y);
        neighborMapP[psP2[pid]+count2]=qid;
      }else{
        // if(pid<35) printf("#qqqq# %d %d %d (%f,%f - %f,%f)\n", pid, psQ2[pid]+count2, qid, P1.x, P1.y, Q1.x, Q1.y);
        neighborMapQ[psQ2[pid]+count2]=qid;
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
                  int *alphaValuesP, int *alphaValuesQ, int *tmpBucketP, int *tmpBucketQ, int *alphaSortedIndiciesP, int *alphaSortedIndiciesQ,
                  int *neighborP, int *neighborQ, int *neighborP2, int *neighborQ2,
                  int *neighborMapP, int *neighborMapQ, int *neighborMapP2, int *neighborMapQ2,
                  int *initLabelsQ){
  int id=blockDim.x*blockIdx.x+threadIdx.x;
  double alpha;
  double beta;
  point I;
  int count1=0, count2=0, nonDegenCount=0, size=sizeQ, indexIntP, indexIntQ, start, end, localI, neighborQId;
  double *poly1X=polyPX, *poly1Y=polyPY, *poly2X=polyQX, *poly2Y=polyQY;

  if(id>=sizeP+sizeQ) return;

  point P1, P2, Q1, Q2;
  int pid=id;
  if(id>=sizeP){
    size=sizeP;
    poly1X=polyQX; 
    poly1Y=polyQY; 
    poly2X=polyPX;
    poly2Y=polyPY;
    pid=id-sizeP;
    intersectionsQ[psQ2[pid]*2]=poly1X[pid];       //consider edge for the intersection array
    intersectionsQ[psQ2[pid]*2+1]=poly1Y[pid];
    // intersectionsQ[psQ2[pid]*3+2]=-100;     //alpha value define it is a parent, not intersection
    intersectionsQ2[psQ2[pid]*2]=poly1X[pid];       //consider edge for the intersection array
    intersectionsQ2[psQ2[pid]*2+1]=poly1Y[pid];
    alphaValuesQ[psQ2[pid]]=-100;
    // neighborMapQ[psQ2[pid]+count2]=-100;    //default neighbor value. No neighbor
    // neighborMapQ2[psQ2[pid]+count2]=-100;    //default neighbor value. No neighbor
    // printf("id %d loc %d x:%f y:%f\n", id, psQ2[pid], intersectionsQ[psQ2[pid]*3], intersectionsQ[psQ2[pid]*3+1]);
    indexIntQ=getIntersectionStartIndex(pid, psQ1);
  } else {
    intersectionsP[psP2[pid]*2]=poly1X[pid];       //consider edge for the intersection array
    intersectionsP[psP2[pid]*2+1]=poly1Y[pid];
    // intersectionsP[psP2[pid]*3+2]=-100;
    
    intersectionsP2[psP2[pid]*2]=poly1X[pid];       //consider edge for the intersection array
    intersectionsP2[psP2[pid]*2+1]=poly1Y[pid];
    alphaValuesP[psP2[pid]]=-100;
    // neighborMapP[psP2[pid]+count2]=-100;    //default neighbor value. No neighbor
    // neighborMapP2[psP2[pid]+count2]=-100;    //default neighbor value. No neighbor
    // printf("id %d loc %d x:%f y:%f\n", id, psP2[pid], intersectionsP[psP2[pid]*3], intersectionsP[psP2[pid]*3+1]);
    indexIntP=getIntersectionStartIndex(pid, psP1);
  }

  for(int qid=0; qid<size; qid++){
    P1.x = poly1X[pid];
    P1.y = poly1Y[pid];

    Q1.x = poly2X[qid];
    Q1.y = poly2Y[qid];
    Q2.x = poly2X[qid+1];
    Q2.y = poly2Y[qid+1];

    // reset P2 vertex of last edge to first vertex
    if(qid == size-1){
      Q2.x = poly2X[0];
      Q2.y = poly2Y[0];
    }
    //polygon1 is P and polygon2 is Q
    if(pid==id && pid==sizeP-1){
      P2.x = poly1X[0];
      P2.y = poly1Y[0];
      // printf("sp %d\n", pid);
    }else if(pid!=id && pid == sizeQ-1){ //polygon2 is P and polygon1 is Q
      P2.x = poly1X[0];
      P2.y = poly1Y[0];
      // printf("sp %d\n", pid);
    } else { //no need reset. Normal case
      P2.x = poly1X[pid+1];
      P2.y = poly1Y[pid+1];
    }

    //
    // determine intersection or overlap type
    //
    // IntersectionType i = intersect(edgeP, edgeQ, alpha, beta);
    int i = getIntersectType(P1, P2, Q1, Q2, alpha, beta);
    // if(i!=0){
    //   printf("%d %d (%f,%f) %d(%f,%f) %d\n", id, pid, P1.x, P1.y, qid, Q1.x, Q1.y, i);
    // }
    if(i && id<sizeP){
      count1++;
      if(i==1 || i==3 || i==5 || i==7){
        nonDegenCount++;
        count2=nonDegenCount;
      }
      else if(i==2 || i==4 || i==6 || i==8)
        count2=0;
      // neighborMapP[psP2[pid]+count2]=indexIntP+count1-1;
      // neighborMapP2[psP2[pid]+count2]=indexIntP+count1-1;
      // neighborP[indexIntP+count1-1]=psP2[pid]+count2;                    //neighbor of new vertex
      // neighborP2[indexIntP+count1-1]=psP2[pid]+count2; 

      start=psQ2[neighborMapP[psP2[pid]+count2]];
      end=psQ2[neighborMapP[psP2[pid]+count2]+1];
      // printf("***-----***** %d %d %d (%d %d)\n", i, neighborMapP[psP2[pid]+count2], pid, start, end);
      // local search to find the index of qid
      for(localI=start; localI<end; ++localI){
        if(pid==neighborMapQ[localI]){
          neighborQId=localI;
          // if(pid<35) printf("&&& %d %d (%d %d) %d %d\n", id, pid, psP2[pid]+count2, neighborQId, start, neighborMapP[psP2[pid]+count2]);
          neighborP[psP2[pid]+count2]=neighborQId+1;   //+1 acting as a padding and helps to identify 0 being empty 
          neighborP2[psP2[pid]+count2]=neighborQId+1;   //+1 acting as a padding and helps to identify 0 being empty 
          neighborQ[neighborQId]=psP2[pid]+count2+1;   //+1 acting as a padding and helps to identify 0 being empty 
          neighborQ2[neighborQId]=psP2[pid]+count2+1;   //+1 acting as a padding and helps to identify 0 being empty 
          localI=end+2; // break;
        }
      }

      switch(i) {
        // case X_INTERSECTION:
        case 1:
          I = add(mulScalar((1.0-alpha), P1), mulScalar(alpha, P2));
          I.x=getValueTolarence(I.x);
          I.y=getValueTolarence(I.y);
          // printf("* %d %d %d %d %f %f\n", (psP2[pid]+count2), indexIntP, count1-1, psP2[pid]+count2, I.x, I.y);
          intersectionsP[(psP2[pid]+count2)*2]=I.x;       //consider edge for the intersection array
          intersectionsP[(psP2[pid]+count2)*2+1]=I.y;
          intersectionsP2[(psP2[pid]+count2)*2]=I.x;       //consider edge for the intersection array
          intersectionsP2[(psP2[pid]+count2)*2+1]=I.y;
          alphaValuesP[psP2[pid]+count2]=(int)pow(10, EPSILON_POSITIONS)*alpha;
          break;
        // X-overlap
        case 5:
          // printf("** %d %d %d %d\n", (psP2[pid]+count2), indexIntP, count1-1, psP2[pid]+count2);
          intersectionsP[(psP2[pid]+count2)*2]=Q1.x;
          intersectionsP[(psP2[pid]+count2)*2+1]=Q1.y;
          intersectionsP2[(psP2[pid]+count2)*2]=Q1.x;
          intersectionsP2[(psP2[pid]+count2)*2+1]=Q1.y;
          alphaValuesP[psP2[pid]+count2]=(int)pow(10, EPSILON_POSITIONS)*alpha;
          break;
        // case T_INTERSECTION_Q:
        // case T_OVERLAP_Q:
        case 2:
        case 6:
          // intersectionsP[psP2[pid]*3+2]=alpha;          //***** error prone. Did not checked in depth
          alphaValuesP[psP2[pid]]=(int)pow(10, EPSILON_POSITIONS)*alpha;
        break;
        // case T_INTERSECTION_P:
        // case T_OVERLAP_P:
        case 3:
        case 7:
          // printf("*** %d %d %d %d\n", (psP2[pid]+count2), indexIntP, count1-1, psP2[pid]+count2);
          intersectionsP[(psP2[pid]+count2)*2]=Q1.x;
          intersectionsP[(psP2[pid]+count2)*2+1]=Q1.y;
          intersectionsP2[(psP2[pid]+count2)*2]=Q1.x;
          intersectionsP2[(psP2[pid]+count2)*2+1]=Q1.y;
          alphaValuesP[psP2[pid]+count2]=(int)pow(10, EPSILON_POSITIONS)*alpha;
          break;
        // case V_INTERSECTION:
        // case V_OVERLAP:
        case 4:
        case 8:
          // printf("%d %d (%f,%f) %d(%f,%f) [%d %d %d] %d\n", id, pid, P1.x, P1.y, qid, Q1.x, Q1.y, (psP2[pid]+count2), indexIntP, count1-1, i);
          alphaValuesP[psP2[pid]]=(int)pow(10, EPSILON_POSITIONS)*alpha;
          break;
      } 
    } else if(i && id>=sizeP){
      initLabelsQ[(psQ2[pid]+count2)]=-100;    //make init label to default -100 
      count1++;
      if(i==1 || i==3 || i==5 || i==7){
      // if(i==1 || i==2 || i==5 || i==6){
        nonDegenCount++;
        count2=nonDegenCount;
      }
      else if(i==2 || i==4 || i==6 || i==8)
        count2=0;

        // if(i==1 || i==3 || i==5 || i==7)
        // // if(i==1 || i==2 || i==5 || i==6)
        //   count2++;
        // neighborMapQ[psQ2[pid]+count2]=indexIntQ+count1-1;
        // neighborMapQ2[psQ2[pid]+count2]=indexIntQ+count1-1;
        // neighborQ[indexIntQ+count1-1]=psQ2[pid]+count2;                    //neighbor of new vertex
        // neighborQ2[indexIntQ+count1-1]=psQ2[pid]+count2;  
      switch(i) {
        // case X_INTERSECTION:
        case 1:
          I = add(mulScalar((1.0-alpha), P1), mulScalar(alpha, P2));
          I.x=getValueTolarence(I.x);
          I.y=getValueTolarence(I.y);
          // printf("/* %d %d %d %d %f %f\n", (psQ2[pid]+count2), indexIntQ, count1-1, psQ2[pid]+count2, I.x, I.y);
          intersectionsQ[(psQ2[pid]+count2)*2]=I.x;       //consider edge for the intersection array
          intersectionsQ[(psQ2[pid]+count2)*2+1]=I.y;
          intersectionsQ2[(psQ2[pid]+count2)*2]=I.x;       //consider edge for the intersection array
          intersectionsQ2[(psQ2[pid]+count2)*2+1]=I.y;
          alphaValuesQ[psQ2[pid]+count2]=(int)pow(10, EPSILON_POSITIONS)*alpha;
          break;
        // case X_OVERLAP:
        case 5:
          // printf("/** %d %d %d %d\n", (psQ2[pid]+count2), indexIntQ, count1-1, psQ2[pid]+count2);
          intersectionsQ[(psQ2[pid]+count2)*2]=Q1.x;    
          intersectionsQ[(psQ2[pid]+count2)*2+1]=Q1.y;
          intersectionsQ2[(psQ2[pid]+count2)*2]=Q1.x;    
          intersectionsQ2[(psQ2[pid]+count2)*2+1]=Q1.y;
          alphaValuesQ[psQ2[pid]+count2]=(int)pow(10, EPSILON_POSITIONS)*beta;
          break;
        // case T_INTERSECTION_Q:
        // case T_OVERLAP_Q: 
        // was 2, 6
        case 3:
        case 7:
          // printf("/*** %d %d %d %d\n", (psQ2[pid]+count2), indexIntQ, count1-1, psQ2[pid]+count2);
          intersectionsQ[(psQ2[pid]+count2)*2]=Q1.x;
          intersectionsQ[(psQ2[pid]+count2)*2+1]=Q1.y;
          intersectionsQ2[(psQ2[pid]+count2)*2]=Q1.x;
          intersectionsQ2[(psQ2[pid]+count2)*2+1]=Q1.y;
          alphaValuesQ[psQ2[pid]+count2]=(int)pow(10, EPSILON_POSITIONS)*alpha;
          break;
        // case T_INTERSECTION_P:
        // case T_OVERLAP_P:
        // was 3, 7
        case 2:
        case 6:
          alphaValuesQ[psQ2[pid]]=(int)pow(10, EPSILON_POSITIONS)*alpha;
        break;
        // case V_INTERSECTION:
        // case V_OVERLAP:
        case 4:
        case 8:
          // printf("/**** %d %d %d %d\n", (psQ2[pid]+count2), indexIntQ, count1-1, psQ2[pid]+count2);
          // printf("%d %d (%f,%f) %d(%f,%f) [%d %d %d] %d\n", id, pid, P1.x, P1.y, qid, Q1.x, Q1.y, (psQ2[pid]+count2), indexIntQ, count1-1, i);
          alphaValuesQ[psQ2[pid]]=(int)pow(10, EPSILON_POSITIONS)*alpha;
          break;
      } 
    }
  }
  // --------------------------------------------------------------------------------------------
  // local sort for each edge, start to end
  // --------------------------------------------------------------------------------------------
  if(id<sizeP){
    int start=psP2[pid], end=psP2[pid+1];
    // printf(".. %d %d %d\n", id, start+1, end);
    // for(i=start; i<end; ++i){
    //   printf("(%d %f %f %f %d) \n", id, intersectionsP[i*3], intersectionsP[i*3+1], intersectionsP[i*3+2], alphaValuesP[i]);
    // }
    // sort intersection vertices in this edge locally
    if((end-start)>2){
      gpuRadixsort(alphaValuesP, tmpBucketP, alphaSortedIndiciesP, start+1, end);
      // using sorted index array, change intersection locations in the array and neighbors
      // printf("ssss %d %d %d\n", id, start, end);
      // decending order JUST FOR TESING
      // for(int i=start+1, j=end-1; i<end; ++i, j--){
      // acending order of alpha values 
      for(int i=start+1, j=start+1; i<end; ++i, j++){
        alphaValuesP[i]=tmpBucketP[j];////////////////?????????????????????? need to swap alpha too!!!
        // (x,y,alpha) tuple change in sorted order
        intersectionsP[alphaSortedIndiciesP[j]*2]=intersectionsP2[i*2];
        intersectionsP[alphaSortedIndiciesP[j]*2+1]=intersectionsP2[i*2+1];
        //neighborMap update
        // neighborMapP[alphaSortedIndiciesP[j]]=neighborMapP2[i];
        //neighbor array update
        neighborP[alphaSortedIndiciesP[j]]=neighborP2[i];
        neighborQ[neighborP2[i]-1]=alphaSortedIndiciesP[j]+1; //+1 is the padding. When reading do -1
        neighborQ2[neighborP2[i]-1]=neighborQ[neighborP2[i]-1]; //updates neighborQ2 as the new originla to be used with sorted Q array
      } 
      // for(int i=start, j=end-1; i<end; ++i, --j){
      //   printf("*(%d %d %f %f %d) reverse->%d \n", id, i, intersectionsP[i*2], intersectionsP[i*2+1], alphaValuesP[i], alphaSortedIndiciesP[j]);
      // }
    } 
  }
  // else{
  //   int start=psQ2[pid], end=psQ2[pid+1];
  //   // printf(".. %d %d %d\n", id, start+1, end);
  //   // for(i=start; i<end; ++i){
  //   //   printf("(%d %f %f %f %d) \n", id, intersectionsP[i*3], intersectionsP[i*3+1], intersectionsP[i*3+2], alphaValuesP[i]);
  //   // }
  //   // sort intersection vertices in this edge locally
  //   if((end-start)>2){
  //     gpuRadixsort(alphaValuesQ, tmpBucketQ, alphaSortedIndiciesQ, start+1, end);
  //     // using sorted index array, change intersection locations in the array and neighbors
  //     // printf("ssss %d %d %d\n", id, start, end);
  //     // decending order JUST FOR TESING
  //     // for(int i=start+1, j=end-1; i<end; ++i, j--){
  //     // acending order of alpha values 
  //     for(int i=start+1, j=start+1; i<end; ++i, j++){
  //       alphaValuesQ[i]=tmpBucketQ[j];////////////////?????????????????????? need to swap alpha too!!!
  //       // (x,y,alpha) tuple change in sorted order
  //       intersectionsQ[alphaSortedIndiciesQ[j]*2]=intersectionsQ2[i*2];
  //       intersectionsQ[alphaSortedIndiciesQ[j]*2+1]=intersectionsQ2[i*2+1];
  //       //neighborMap update
  //       // neighborMapQ[alphaSortedIndiciesQ[j]]=neighborMapQ2[i];
  //       //neighbor array update
  //       neighborQ[alphaSortedIndiciesQ[j]]=neighborQ2[i];
  //       neighborP[neighborQ2[i]-1]=alphaSortedIndiciesQ[j]+1; //+1 is the padding. When reading do -1
  //     } 
  //     // for(int i=start, j=end-1; i<end; ++i, --j){
  //     //   printf("****(%d %d %f %f %d) reverse->%d \n", id, i, intersectionsQ[i*2], intersectionsQ[i*2+1], alphaValuesQ[i], alphaSortedIndiciesQ[j]);
  //     // }
  //   } 
  // }
  // --------------------------------------------------------------------------------------------
}

/*
-----------------------------------------------------------------
Function to save vertices of Q in edge wise sorted order
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
__global__ void gpuSortPolyQ(
                  int sizeP, int sizeQ, 
                  int *psQ2, 
                  double *intersectionsQ, double *intersectionsQ2,
                  int *alphaValuesQ, int *tmpBucketQ,  int *alphaSortedIndiciesQ,
                  int *neighborP, int *neighborQ, int *neighborQ2){
  int id=blockDim.x*blockIdx.x+threadIdx.x;
  if(id>=sizeP && id<(sizeP+sizeQ)){
    int pid=id-sizeP;
    int start=psQ2[pid], end=psQ2[pid+1];
    // sort intersection vertices in this edge locally
    if((end-start)>2){
      gpuRadixsort(alphaValuesQ, tmpBucketQ, alphaSortedIndiciesQ, start+1, end);
      // using sorted index array, change intersection locations in the array and neighbors
      // printf("ssss %d %d %d\n", id, start, end);
      // decending order JUST FOR TESING
      // for(int i=start+1, j=end-1; i<end; ++i, j--){
      // acending order of alpha values 
      for(int i=start+1, j=start+1; i<end; ++i, j++){
        alphaValuesQ[i]=tmpBucketQ[j];////////////////?????????????????????? need to swap alpha too!!!
        // (x,y,alpha) tuple change in sorted order
        intersectionsQ[alphaSortedIndiciesQ[j]*2]=intersectionsQ2[i*2];
        intersectionsQ[alphaSortedIndiciesQ[j]*2+1]=intersectionsQ2[i*2+1];
        //neighbor array update
        neighborQ[alphaSortedIndiciesQ[j]]=neighborQ2[i];
        neighborP[neighborQ2[i]-1]=alphaSortedIndiciesQ[j]+1; //+1 is the padding. When reading do -1
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
  int id=blockDim.x*blockIdx.x+threadIdx.x;
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
    // if(intersectionsP[i*2+2]!=-100){    //consider intersections only
    if(alphaValuesP[i]!=-100){    //consider intersections only
      current.x=intersectionsP[i*2]; 
      current.y=intersectionsP[i*2+1]; 
      tmpId=getCircularId(i-1, sizeNP);
      // determine local configuration at this intersection vertex
      pM.x=intersectionsP[tmpId*2];                // P-, predecessor of I on P
      pM.y=intersectionsP[tmpId*2+1];                // P-, predecessor of I on P
      // if(intersectionsP[tmpId*2+2]!=-100)
      if(alphaValuesP[tmpId]!=-100)
        // pMNId=getNeighborIndex(tmpId, neighborMapP, neighborQ); //get neighbor id of P_m vertex
        pMNId=neighborP[tmpId]-1; //get neighbor id of P_m vertex
      else pMNId=-100;

      tmpId=getCircularId(i+1, sizeNP);
      pP.x=intersectionsP[tmpId*2];                // P+, successor of I on P
      pP.y=intersectionsP[tmpId*2+1];                // P+, successor of I on P
      // if(intersectionsP[tmpId*2+2]!=-100)
      if(alphaValuesP[tmpId]!=-100)
        // pPNId=getNeighborIndex(tmpId, neighborMapP, neighborQ); //get neighbor id of P_p vertex
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
      // printf("%d %d (%f, %f) (%f, %f) (%f, %f) (%f, %f)\n", i, nId, pM.x, pM.y, pP.x, pP.y, qM.x, qM.y, qP.x, qP.y);
      // printf(">>> %d %d %d %d\n", i, qMType, qPType, getInitialLabel(qMType, qPType));
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
                  int sizeP, int sizeQ, 
                  int *countNonDegenIntP, int *countNonDegenIntQ, 
                  double **intersectionsP, double **intersectionsQ, int **alphaValuesP, int **alphaValuesQ,
                  int **initLabelsP, int **initLabelsQ,
                  int **neighborMapP, int **neighborMapQ, int **neighborP, int **neighborQ){
  double *dev_polyPX, *dev_polyPY, *dev_polyQX, *dev_polyQY;
  int *dev_psP1, *dev_psP2, *dev_psQ1, *dev_psQ2;
  int psP1[sizeP+1], psP2[sizeP+1], psQ1[sizeQ+1], psQ2[sizeQ+1];

  // Phase1: Count intersections in each block. Create prefix sums to find local locations in each thread 
  // Allocate memory in device 
  cudaMalloc((void **) &dev_polyPX, sizeP*sizeof(double));
  cudaMalloc((void **) &dev_polyPY, sizeP*sizeof(double));
  cudaMalloc((void **) &dev_polyQX, sizeQ*sizeof(double));
  cudaMalloc((void **) &dev_polyQY, sizeQ*sizeof(double));
  cudaMalloc((void **) &dev_psP1, (sizeP+1)*sizeof(int));
  cudaMalloc((void **) &dev_psP2, (sizeP+1)*sizeof(int));
  cudaMalloc((void **) &dev_psQ1, (sizeQ+1)*sizeof(int));
  cudaMalloc((void **) &dev_psQ2, (sizeQ+1)*sizeof(int));

  // Copy input vectors from host memory to GPU buffers.
  cudaMemcpy(dev_polyPX, polyPX, sizeP*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_polyPY, polyPY, sizeP*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_polyQX, polyQX, sizeQ*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_polyQY, polyQY, sizeQ*sizeof(double), cudaMemcpyHostToDevice);

  // int threadsPerBlock = 2;
  // int blocksPerGrid = (2*(sizeP+sizeQ) + threadsPerBlock - 1) / threadsPerBlock;
  // dim3 dimBlock(threadsPerBlock, 1, 1), dimGrid(blocksPerGrid, 1, 1);   
  int xblocksPerGrid = (2*(sizeP+sizeQ) + xThreadPerBlock - 1) / xThreadPerBlock;
  int yblocksPerGrid = (2*(sizeP+sizeQ) + yThreadPerBlock - 1) / yThreadPerBlock;
  dim3 dimBlock(xThreadPerBlock, yThreadPerBlock, 1), dimGrid(xblocksPerGrid, 1, 1); 
  printf("blockDim %d gridDim %d\n", dimBlock.x, dimGrid.x);

  gpuCountIntersections<<<dimGrid, dimBlock>>>(
        dev_polyPX, dev_polyPY, 
        dev_polyQX, dev_polyQY, 
        sizeP, sizeQ, 
        dev_psP1, dev_psP2, dev_psQ1, dev_psQ2);

  cudaMemcpy(&psP1, dev_psP1, (sizeP+1)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psP2, dev_psP2, (sizeP+1)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psQ1, dev_psQ1, (sizeQ+1)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psQ2, dev_psQ2, (sizeQ+1)*sizeof(int), cudaMemcpyDeviceToHost);

  thrust::exclusive_scan(thrust::host, psP1, psP1 + sizeP+1, psP1);   //sizeP location contains the total size of the count1
  thrust::exclusive_scan(thrust::host, psP2, psP2 + sizeP+1, psP2);
  thrust::exclusive_scan(thrust::host, psQ1, psQ1 + sizeQ+1, psQ1);   //sizeQ location contains the total size of the count1
  thrust::exclusive_scan(thrust::host, psQ2, psQ2 + sizeQ+1, psQ2);

  // for (int i = 0; i < sizeQ+1; ++i){
  // for (int i = 0; i < 15+1; ++i){
  //   printf(" %d-%d ", i, psP2[i]);
  // }
  // printf("--- \n");

  // // for (int i = 0; i < sizeQ+1; ++i){
  // for (int i = 0; i < 15+1; ++i){
  //   printf(" %d-%d ", i, psQ2[i]);
  // }
  // printf("--- \n");
  cudaDeviceSynchronize();

  //Phase2: NEW- Fill neighborMap
  int *dev_neighborMapP, *dev_neighborMapQ;
  *countNonDegenIntP=psP2[sizeP];
  *countNonDegenIntQ=psQ2[sizeQ];

  *neighborMapP=(int *)malloc(*countNonDegenIntP*sizeof(int));
  *neighborMapQ=(int *)malloc(*countNonDegenIntQ*sizeof(int));
  
  cudaMalloc((void **) &dev_neighborMapP, *countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_neighborMapQ, *countNonDegenIntQ*sizeof(int));

  cudaMemcpy(dev_psP1, psP1, (sizeP+1)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_psP2, psP2, (sizeP+1)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_psQ1, psQ1, (sizeQ+1)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_psQ2, psQ2, (sizeQ+1)*sizeof(int), cudaMemcpyHostToDevice);

  gpuNeighborMap<<<dimGrid, dimBlock>>>(
        dev_polyPX, dev_polyPY, 
        dev_polyQX, dev_polyQY, 
        sizeP, sizeQ,  
        dev_psP2, dev_psQ2,
        dev_neighborMapP, dev_neighborMapQ);
  
// -----------------------------------------------------------------------------------------------------
  // remove after kernel testing
  // cudaMemcpy(*neighborMapP, dev_neighborMapP, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  // cudaMemcpy(*neighborMapQ, dev_neighborMapQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
// -----------------------------------------------------------------------------------------------------


  // Phase 3: Calcualte intersections and save them in the arrays. Make neighbor connections
  int countIntersections=psP1[sizeP];

  int *alphaSortedIndiciesP, *alphaSortedIndiciesQ;
  double *dev_intersectionsP, *dev_intersectionsQ, *dev_intersectionsP2, *dev_intersectionsQ2;
  int *dev_neighborP, *dev_neighborQ, *dev_neighborP2, *dev_neighborQ2;
  int *dev_neighborMapP2, *dev_neighborMapQ2, *dev_initLabelsP, *dev_initLabelsQ;
  int  *dev_alphaValuesP, *dev_alphaValuesQ, *dev_tmpBucketP, *dev_tmpBucketQ, *dev_alphaSortedIndiciesP, *dev_alphaSortedIndiciesQ;

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
  cudaMalloc((void **) &dev_neighborMapP2, *countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_neighborMapQ2, *countNonDegenIntQ*sizeof(int));
  cudaMalloc((void **) &dev_initLabelsQ, *countNonDegenIntQ*sizeof(int));

  gpuCalculateIntersections<<<dimGrid, dimBlock>>>(
        dev_polyPX, dev_polyPY, 
        dev_polyQX, dev_polyQY, 
        sizeP, sizeQ, 
        dev_psP1, dev_psP2, dev_psQ1, dev_psQ2, 
        dev_intersectionsP, dev_intersectionsQ, dev_intersectionsP2, dev_intersectionsQ2,
        dev_alphaValuesP, dev_alphaValuesQ, dev_tmpBucketP, dev_tmpBucketQ, dev_alphaSortedIndiciesP, dev_alphaSortedIndiciesQ,
        dev_neighborP, dev_neighborQ, dev_neighborP2, dev_neighborQ2,
        dev_neighborMapP, dev_neighborMapQ, dev_neighborMapP2, dev_neighborMapQ2,
        dev_initLabelsQ);

  gpuSortPolyQ<<<dimGrid, dimBlock>>>(
        sizeP, sizeQ, 
        dev_psQ2, 
        dev_intersectionsQ, dev_intersectionsQ2,
        dev_alphaValuesQ, dev_tmpBucketQ,  dev_alphaSortedIndiciesQ,
        dev_neighborP, dev_neighborQ, dev_neighborQ2);

  // Phase4: Inital label classificaiton
  // cudaMemcpy(*initLabelsQ, dev_initLabelsQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMalloc((void **) &dev_initLabelsP, *countNonDegenIntP*sizeof(int));
  // cudaMemcpy(dev_initLabelsQ, *initLabelsQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyHostToDevice);
 
  // negative alpha values are not handled explicitly since they are original vertices
  // ******No need to copy alpha values since they are only used to sort edge wise******
  // cudaMemcpy(alphaSortedIndicies, dev_alphaSortedIndicies, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);

  gpuCalculateInitLabel<<<dimGrid, dimBlock>>>(
      sizeP,  dev_psP2,
      dev_intersectionsP, dev_intersectionsQ, dev_alphaValuesP,
      dev_neighborP,
      *countNonDegenIntP, *countNonDegenIntQ, dev_initLabelsP, dev_initLabelsQ);

  cudaMemcpy(*intersectionsP, dev_intersectionsP, *countNonDegenIntP*2*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(*intersectionsQ, dev_intersectionsQ, *countNonDegenIntQ*2*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(*neighborP, dev_neighborP, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*neighborQ, dev_neighborQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*neighborMapP, dev_neighborMapP, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*neighborMapQ, dev_neighborMapQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);

  cudaMemcpy(*initLabelsP, dev_initLabelsP, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*initLabelsQ, dev_initLabelsQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*alphaValuesP, dev_alphaValuesP, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*alphaValuesQ, dev_alphaValuesQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  
  cudaDeviceSynchronize();

  // int limitP=*countNonDegenIntP;
  // int limitQ=*countNonDegenIntQ;
  int limitP=10;
  int limitQ=10;

  printf("intersectionP");
  for (int i = 0; i < limitP*2; ++i){
    if(i%2==0) 
      printf("\n%d %d ", i/2, *(*alphaValuesP+(i/2)));
    // printf(" %f ", intersectionsP[i]);
    printf(" %f ", *(*intersectionsP+i));
  }
  printf("\n\nintersectionQ");
  for (int i = 0; i < limitQ*2; ++i){
    if(i%2==0)
      printf("\n%d %d ", i/2, *(*alphaValuesQ+(i/2)));
    printf(" %f ", *(*intersectionsQ+i));
  }
  // printf("\n\nalpha P\n");
  // for (int i = 0; i < *countNonDegenIntP; ++i){
  //   printf(" %d>%d ", i, alphaValuesP[i]);
  // }
  // printf("\n\nalpha Q\n");
  // for (int i = 0; i < *countNonDegenIntQ; ++i){
  //   printf(" %d>%d ", i, alphaValuesQ[i]);
  // }
  // printf("\n");
  printf("\nneighbor P\n");
  for (int i = 0; i < limitP; ++i){
    printf(" %d-%d ", i, *(*neighborP+i));
  }
  printf("\nnneighbor Q\n");
  for (int i = 0; i < limitQ; ++i){
    printf(" %d-%d ", i, *(*neighborQ+i));
  }
  // printf("\n");
  // for (int i = 0; i < *countNonDegenIntP; ++i)
  // {
  //   printf(" %d-%d ", i, *(*neighborMapP+i));
  // }
  // printf("\n");
  // for (int i = 0; i < *countNonDegenIntQ; ++i)
  // {
  //   printf(" %d-%d ", i, *(*neighborMapQ+i));
  // }
  printf("\nLabel P\n");
  for (int i = 0; i < limitP; ++i){
    printf(" %d>%d ", i, *(*initLabelsP+i));
  }
  printf("\nLable Q\n");
  for (int i = 0; i < limitQ; ++i){
    printf(" %d>%d ", i, *(*initLabelsQ+i));
  }
  printf("\n");


  cudaFree(dev_polyPX);
  cudaFree(dev_polyPY);
  cudaFree(dev_polyQX);
  cudaFree(dev_polyQY);
}

/*
-----------------------------------------------------------------
Function to count how many intersection points and prefix sums
  Works with multiple components in PP and QQ
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
void calculateIntersectionsMultipleComponents(
                  double *polyPX, double *polyPY, 
                  double *polyQX,  double *polyQY, 
                  int sizeP, int sizeQ, int *sizesPP, int *sizesQQ, int sizePP, int sizeQQ,
                  int *countNonDegenIntArrayP, int *countNonDegenIntArrayQ, 
                  double **intersectionsP, double **intersectionsQ, int **alphaValuesP, int **alphaValuesQ,
                  int **initLabelsP, int **initLabelsQ,
                  int **neighborMapP, int **neighborMapQ, int **neighborP, int **neighborQ){
  double *dev_polyPX, *dev_polyPY, *dev_polyQX, *dev_polyQY;
  int *dev_psP1, *dev_psP2, *dev_psQ1, *dev_psQ2;
  int psP1[sizeP+1], psP2[sizeP+1], psQ1[sizeQ+1], psQ2[sizeQ+1];

  
  printf("\nP ");
  for (int i = 0; i < sizeP; ++i){
  // for (int i = 0; i < 15+1; ++i){
    printf(" %d-%f,%f \n", i, polyPX[i], polyPY[i]);
  }
  printf("--- \nQ ");
  for (int i = 0; i < sizeQ; ++i){
  // for (int i = 0; i < 15+1; ++i){
    printf(" %d-%f,%f\n ", i, polyQX[i], polyQY[i]);
  }
  printf("\n");
  // Phase1: Count intersections in each block. Create prefix sums to find local locations in each thread 
  // Allocate memory in device 
  cudaMalloc((void **) &dev_polyPX, sizeP*sizeof(double));
  cudaMalloc((void **) &dev_polyPY, sizeP*sizeof(double));
  cudaMalloc((void **) &dev_polyQX, sizeQ*sizeof(double));
  cudaMalloc((void **) &dev_polyQY, sizeQ*sizeof(double));
  cudaMalloc((void **) &dev_psP1, (sizeP+1)*sizeof(int));
  cudaMalloc((void **) &dev_psP2, (sizeP+1)*sizeof(int));
  cudaMalloc((void **) &dev_psQ1, (sizeQ+1)*sizeof(int));
  cudaMalloc((void **) &dev_psQ2, (sizeQ+1)*sizeof(int));

  // Copy input vectors from host memory to GPU buffers.
  cudaMemcpy(dev_polyPX, polyPX, sizeP*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_polyPY, polyPY, sizeP*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_polyQX, polyQX, sizeQ*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_polyQY, polyQY, sizeQ*sizeof(double), cudaMemcpyHostToDevice);

  // int threadsPerBlock = 2;
  // int blocksPerGrid = (2*(sizeP+sizeQ) + threadsPerBlock - 1) / threadsPerBlock;
  // dim3 dimBlock(threadsPerBlock, 1, 1), dimGrid(blocksPerGrid, 1, 1);   
  int xblocksPerGrid = (2*(sizeP+sizeQ) + xThreadPerBlock - 1) / xThreadPerBlock;
  int yblocksPerGrid = (2*(sizeP+sizeQ) + yThreadPerBlock - 1) / yThreadPerBlock;
  dim3 dimBlock(xThreadPerBlock, yThreadPerBlock, 1), dimGrid(xblocksPerGrid, 1, 1); 
  printf("blockDim %d gridDim %d\n", dimBlock.x, dimGrid.x);

  gpuCountIntersections<<<dimGrid, dimBlock>>>(
        dev_polyPX, dev_polyPY, 
        dev_polyQX, dev_polyQY, 
        sizeP, sizeQ, 
        dev_psP1, dev_psP2, dev_psQ1, dev_psQ2);

  cudaMemcpy(&psP1, dev_psP1, (sizeP+1)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psP2, dev_psP2, (sizeP+1)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psQ1, dev_psQ1, (sizeQ+1)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psQ2, dev_psQ2, (sizeQ+1)*sizeof(int), cudaMemcpyDeviceToHost);

  thrust::exclusive_scan(thrust::host, psP1, psP1 + sizeP+1, psP1);   //sizeP location contains the total size of the count1
  thrust::exclusive_scan(thrust::host, psP2, psP2 + sizeP+1, psP2);
  thrust::exclusive_scan(thrust::host, psQ1, psQ1 + sizeQ+1, psQ1);   //sizeQ location contains the total size of the count1
  thrust::exclusive_scan(thrust::host, psQ2, psQ2 + sizeQ+1, psQ2);

  printf("\nsizesPP ");
  for (int i = 0; i < sizePP; ++i){
  // for (int i = 0; i < 15+1; ++i){
    printf(" %d-%d ", i, sizesPP[i]);
  }
  printf("--- \nsizesQQ ");
  for (int i = 0; i < sizeQQ; ++i){
  // for (int i = 0; i < 15+1; ++i){
    printf(" %d-%d ", i, sizesQQ[i]);
  }
  printf("--- \npsP2 ");
  for (int i = 0; i < sizeP+1; ++i){
  // for (int i = 0; i < 15+1; ++i){
    printf(" %d-%d ", i, psP2[i]);
  }
  printf("--- \npsQ2 ");

  for (int i = 0; i < sizeQ+1; ++i){
  // for (int i = 0; i < 15+1; ++i){
    printf(" %d-%d ", i, psQ2[i]);
  }
  printf("--- \n");
  cudaDeviceSynchronize();

  //Phase2: NEW- Fill neighborMap
  int *dev_neighborMapP, *dev_neighborMapQ;
  int countNonDegenIntP, countNonDegenIntQ;
  int count=0;
  for(int i=0; i<sizePP; ++i){
    count+=sizesPP[i];
    countNonDegenIntArrayP[i]=psP2[count];
  }
  count=0;
  for(int i=0; i<sizeQQ; ++i){
    count+=sizesQQ[i];
    countNonDegenIntArrayQ[i]=psQ2[count];
  }
  countNonDegenIntP=psP2[sizeP];
  countNonDegenIntQ=psQ2[sizeQ];

  *neighborMapP=(int *)malloc(countNonDegenIntP*sizeof(int));
  *neighborMapQ=(int *)malloc(countNonDegenIntQ*sizeof(int));
  
  cudaMalloc((void **) &dev_neighborMapP, countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_neighborMapQ, countNonDegenIntQ*sizeof(int));

  cudaMemcpy(dev_psP1, psP1, (sizeP+1)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_psP2, psP2, (sizeP+1)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_psQ1, psQ1, (sizeQ+1)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_psQ2, psQ2, (sizeQ+1)*sizeof(int), cudaMemcpyHostToDevice);

  gpuNeighborMap<<<dimGrid, dimBlock>>>(
        dev_polyPX, dev_polyPY, 
        dev_polyQX, dev_polyQY, 
        sizeP, sizeQ,  
        dev_psP2, dev_psQ2,
        dev_neighborMapP, dev_neighborMapQ);
  
// -----------------------------------------------------------------------------------------------------
  // remove after kernel testing
  // cudaMemcpy(*neighborMapP, dev_neighborMapP, countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  // cudaMemcpy(*neighborMapQ, dev_neighborMapQ, countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
// -----------------------------------------------------------------------------------------------------


  // Phase 3: Calcualte intersections and save them in the arrays. Make neighbor connections
  int countIntersections=psP1[sizeP];

  int *alphaSortedIndiciesP, *alphaSortedIndiciesQ;
  double *dev_intersectionsP, *dev_intersectionsQ, *dev_intersectionsP2, *dev_intersectionsQ2;
  int *dev_neighborP, *dev_neighborQ, *dev_neighborP2, *dev_neighborQ2;
  int *dev_neighborMapP2, *dev_neighborMapQ2, *dev_initLabelsP, *dev_initLabelsQ;
  int  *dev_alphaValuesP, *dev_alphaValuesQ, *dev_tmpBucketP, *dev_tmpBucketQ, *dev_alphaSortedIndiciesP, *dev_alphaSortedIndiciesQ;

  *intersectionsP=(double *)malloc(countNonDegenIntP*2*sizeof(double));
  *intersectionsQ=(double *)malloc(countNonDegenIntQ*2*sizeof(double));
  *alphaValuesP=(int *)malloc(countNonDegenIntP*sizeof(int));
  *alphaValuesQ=(int *)malloc(countNonDegenIntQ*sizeof(int));
  alphaSortedIndiciesP=(int *)malloc(countNonDegenIntP*sizeof(int));
  alphaSortedIndiciesQ=(int *)malloc(countNonDegenIntQ*sizeof(int));
  *initLabelsP=(int *)malloc(countNonDegenIntP*sizeof(int));
  *initLabelsQ=(int *)malloc(countNonDegenIntQ*sizeof(int));
  *neighborP=(int *)malloc(countNonDegenIntP*sizeof(int));
  *neighborQ=(int *)malloc(countNonDegenIntQ*sizeof(int));

  // Allocate memory in device 
  cudaMalloc((void **) &dev_intersectionsP, countNonDegenIntP*2*sizeof(double));
  cudaMalloc((void **) &dev_intersectionsP2, countNonDegenIntP*2*sizeof(double));
  cudaMalloc((void **) &dev_intersectionsQ, countNonDegenIntQ*2*sizeof(double));
  cudaMalloc((void **) &dev_intersectionsQ2, countNonDegenIntQ*2*sizeof(double));
  cudaMalloc((void **) &dev_alphaValuesP, countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_alphaValuesQ, countNonDegenIntQ*sizeof(int));

  cudaMalloc((void **) &dev_tmpBucketP, countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_tmpBucketQ, countNonDegenIntQ*sizeof(int));
  cudaMalloc((void **) &dev_alphaSortedIndiciesP, countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_alphaSortedIndiciesQ, countNonDegenIntQ*sizeof(int));

  cudaMalloc((void **) &dev_neighborP, countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_neighborP2, countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_neighborQ, countNonDegenIntQ*sizeof(int));
  cudaMalloc((void **) &dev_neighborQ2, countNonDegenIntQ*sizeof(int));
  cudaMalloc((void **) &dev_neighborMapP2, countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_neighborMapQ2, countNonDegenIntQ*sizeof(int));
  cudaMalloc((void **) &dev_initLabelsQ, countNonDegenIntQ*sizeof(int));

  gpuCalculateIntersections<<<dimGrid, dimBlock>>>(
        dev_polyPX, dev_polyPY, 
        dev_polyQX, dev_polyQY, 
        sizeP, sizeQ, 
        dev_psP1, dev_psP2, dev_psQ1, dev_psQ2, 
        dev_intersectionsP, dev_intersectionsQ, dev_intersectionsP2, dev_intersectionsQ2,
        dev_alphaValuesP, dev_alphaValuesQ, dev_tmpBucketP, dev_tmpBucketQ, dev_alphaSortedIndiciesP, dev_alphaSortedIndiciesQ,
        dev_neighborP, dev_neighborQ, dev_neighborP2, dev_neighborQ2,
        dev_neighborMapP, dev_neighborMapQ, dev_neighborMapP2, dev_neighborMapQ2,
        dev_initLabelsQ);

  gpuSortPolyQ<<<dimGrid, dimBlock>>>(
        sizeP, sizeQ, 
        dev_psQ2, 
        dev_intersectionsQ, dev_intersectionsQ2,
        dev_alphaValuesQ, dev_tmpBucketQ,  dev_alphaSortedIndiciesQ,
        dev_neighborP, dev_neighborQ, dev_neighborQ2);

  // Phase4: Inital label classificaiton
  // cudaMemcpy(*initLabelsQ, dev_initLabelsQ, countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMalloc((void **) &dev_initLabelsP, countNonDegenIntP*sizeof(int));
  // cudaMemcpy(dev_initLabelsQ, *initLabelsQ, countNonDegenIntQ*sizeof(int), cudaMemcpyHostToDevice);
 
  // negative alpha values are not handled explicitly since they are original vertices
  // ******No need to copy alpha values since they are only used to sort edge wise******
  // cudaMemcpy(alphaSortedIndicies, dev_alphaSortedIndicies, countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);

  gpuCalculateInitLabel<<<dimGrid, dimBlock>>>(
      sizeP,  dev_psP2,
      dev_intersectionsP, dev_intersectionsQ, dev_alphaValuesP,
      dev_neighborP,
      countNonDegenIntP, countNonDegenIntQ, dev_initLabelsP, dev_initLabelsQ);

  cudaMemcpy(*intersectionsP, dev_intersectionsP, countNonDegenIntP*2*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(*intersectionsQ, dev_intersectionsQ, countNonDegenIntQ*2*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(*neighborP, dev_neighborP, countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*neighborQ, dev_neighborQ, countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*neighborMapP, dev_neighborMapP, countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*neighborMapQ, dev_neighborMapQ, countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);

  cudaMemcpy(*initLabelsP, dev_initLabelsP, countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*initLabelsQ, dev_initLabelsQ, countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*alphaValuesP, dev_alphaValuesP, countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*alphaValuesQ, dev_alphaValuesQ, countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  
  cudaDeviceSynchronize();

  // int limitP=countNonDegenIntP;
  // int limitQ=countNonDegenIntQ;
  int limitP=10;
  int limitQ=10;

  printf("intersectionP");
  for (int i = 0; i < limitP*2; ++i){
    if(i%2==0) 
      printf("\n%d %d ", i/2, *(*alphaValuesP+(i/2)));
    // printf(" %f ", intersectionsP[i]);
    printf(" %f ", *(*intersectionsP+i));
  }
  printf("\n\nintersectionQ");
  for (int i = 0; i < limitQ*2; ++i){
    if(i%2==0)
      printf("\n%d %d ", i/2, *(*alphaValuesQ+(i/2)));
    printf(" %f ", *(*intersectionsQ+i));
  }
  // printf("\n\nalpha P\n");
  // for (int i = 0; i < limitP; ++i){
  //   printf(" %d>%d ", i, alphaValuesP[i]);
  // }
  // printf("\n\nalpha Q\n");
  // for (int i = 0; i < limitQ; ++i){
  //   printf(" %d>%d ", i, alphaValuesQ[i]);
  // }
  // printf("\n");
  printf("\nneighbor P\n");
  for (int i = 0; i < limitP; ++i){
    printf(" %d-%d ", i, *(*neighborP+i));
  }
  printf("\nneighbor Q\n");
  for (int i = 0; i < limitQ; ++i){
    printf(" %d-%d ", i, *(*neighborQ+i));
  }
  // printf("\n");
  // for (int i = 0; i < limitP; ++i)
  // {
  //   printf(" %d-%d ", i, *(*neighborMapP+i));
  // }
  // printf("\n");
  // for (int i = 0; i < limitQ; ++i)
  // {
  //   printf(" %d-%d ", i, *(*neighborMapQ+i));
  // }
  printf("\nLabel P\n");
  for (int i = 0; i < limitP; ++i){
    printf(" %d>%d ", i, *(*initLabelsP+i));
  }
  printf("\nLable Q\n");
  for (int i = 0; i < limitQ; ++i){
    printf(" %d>%d ", i, *(*initLabelsQ+i));
  }
  printf("\n");


  cudaFree(dev_polyPX);
  cudaFree(dev_polyPY);
  cudaFree(dev_polyQX);
  cudaFree(dev_polyQY);
}