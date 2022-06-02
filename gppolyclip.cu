#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/scan.h>

#include "constants.h"

typedef struct
{
  double x, y;
} point;

// compute twice the signed area of the triange [P,Q,R]
__device__ double A(const point& P, const point& Q, const point& R) {
	return (Q.x-P.x) * (R.y-P.y) - (Q.y-P.y) * (R.x-P.x);
}

// difference of two 2D points
__device__ point sub(const point& a, const point& b) { 
  point r;
  r.x=a.x-b.x;
  r.y=a.y-b.y;
  return r; 
}

// add two 2D points
__device__ point add(const point& a, const point& b) { 
  point r;
  r.x=a.x+b.x;
  r.y=a.y+b.y;
  return r; 
}

// multiply two 2D points
__device__ double mul(const point& a, const point& b) { 
  point r;
  r.x=a.x*b.x;
  r.y=a.y*b.y;
  return (r.x+r.y); 
}

// multiply scalar with 2D points
__device__ point mulScalar(const double c, const point& b) { 
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
__device__ int getNeighborIndex(int id, int *neighborMap, int *neighbor){
  int neighborId=neighborMap[id];
  if(neighborId!=-100){
    return neighbor[neighborId];
  }
  return -1; //no neighbor for this vertex
}

/*
-----------------------------------------------------------------
Function to return intersection  type
Returns the type of the intersection
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
__device__ int getIntersectType(const point& P1, const point& P2, const point& Q1, const point& Q2,  double& alpha, double& beta) {
	double AP1 = A(P1,Q1,Q2);
	double AP2 = A(P2,Q1,Q2);

	if (fabs(AP1-AP2) > EPSILON) {
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
	}
	else
		if (fabs(AP1) < EPSILON) {
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
			if ( (alpha > EPSILON) && (alpha < 1.0-EPSILON) )
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
			if ( (beta > EPSILON) && (beta < 1.0-EPSILON) )
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
Function to count all intersections. 
Return prefix sum arrays.
  *prefix sum of count of all intersection vertices x2 (P and Q)
  *prefix sum of count of all intersection vertices excluding 
   degenerate cases x2 (P and Q)
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
__global__ void deviceCountIntersections(double *polyPX, double *polyPY, double *polyQX,  double *polyQY, int sizeP, int sizeQ, int *psP1, int *psP2, int *psQ1, int *psQ2){
  int id=blockIdx.x+threadIdx.x;
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
    //   printf("id %d count1 %d count2  %d (%f,%f) (%f,%f) %d\n", id, count1, count2, P1.x, P1.y, Q1.x, Q1.y, i);

  }
  count2++; //represent the parent vertex 
  if(id<sizeP){
    psP1[pid]=count1;
    psP2[pid]=count2;
  } else{
    psQ1[pid]=count1;
    psQ2[pid]=count2;
  }
  __syncthreads();
  thrust::exclusive_scan(thrust::device, psP1, psP1 + sizeP+1, psP1);   //sizeP location contains the total size of the count1
  thrust::exclusive_scan(thrust::device, psP2, psP2 + sizeP+1, psP2);
  thrust::exclusive_scan(thrust::device, psQ1, psQ1 + sizeQ+1, psQ1);   //sizeQ location contains the total size of the count1
  thrust::exclusive_scan(thrust::device, psQ2, psQ2 + sizeQ+1, psQ2);
  // printf("id %d count1 %d count2  %d (%f,%f) (%f,%f)\n", id, count1, count2, P1.x, P1.y, P2.x, P2.y);
  // __syncthreads();
  // if(id==0){
  //   for(int ii=0; ii<sizeP; ++ii){
  //     printf("%d *%d ", dev_psP1[ii], dev_psP2[ii]);
  //   }
  //   printf("\nend\n");
  //   for(int ii=0; ii<sizeQ; ++ii){
  //     printf("%d *%d ", dev_psQ1[ii], dev_psQ2[ii]);
  //   }
  //   printf("\nend\n");
  // }
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
__global__ void deviceCalculateIntersections(double *polyPX, double *polyPY, double *polyQX, double *polyQY, int sizeP, int sizeQ, int *psP1, int *psP2, int *psQ1, int *psQ2, double *intersectionsP, double *intersectionsQ, int *neighborP, int *neighborQ, int *neighborMapP, int *neighborMapQ, int sizeNP, int sizeNQ, int *initLabelsP, int *initLabelsQ){
  int id=blockIdx.x+threadIdx.x;
  double alpha;
  double beta;
  point I;
  int count1=0, count2=0, size=sizeQ, indexIntP, indexIntQ;
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
    intersectionsQ[psQ2[pid]*3]=poly1X[pid];       //consider edge for the intersection array
    intersectionsQ[psQ2[pid]*3+1]=poly1Y[pid];
    intersectionsQ[psQ2[pid]*3+2]=-100;     //alpha value define it is a parent, not intersection
    neighborMapQ[psQ2[pid]+count2]=-100;    //default neighbor value. No neighbor
    // printf("id %d loc %d x:%f y:%f\n", id, psQ2[pid], intersectionsQ[psQ2[pid]*3], intersectionsQ[psQ2[pid]*3+1]);
    indexIntQ=getIntersectionStartIndex(pid, psQ1);
  } else {
    intersectionsP[psP2[pid]*3]=poly1X[pid];       //consider edge for the intersection array
    intersectionsP[psP2[pid]*3+1]=poly1Y[pid];
    intersectionsP[psP2[pid]*3+2]=-100;     //alpha value define it is a parent, not intersection
    neighborMapP[psP2[pid]+count2]=-100;    //default neighbor value. No neighbor
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
    if(i && id<sizeP){
      count1++;
      if(i==1 || i==3 || i==5 || i==7)
        count2++;
      switch(i) {
        // case X_INTERSECTION:
        case 1:
          I = add(mulScalar((1.0-alpha), P1), mulScalar(alpha, P2));
          I.x=getValueTolarence(I.x);
          I.y=getValueTolarence(I.y);
          // printf("* %d %d %d %d %f %f\n", (psP2[pid]+count2), indexIntP, count1-1, psP2[pid]+count2, I.x, I.y);
          intersectionsP[(psP2[pid]+count2)*3]=I.x;       //consider edge for the intersection array
          intersectionsP[(psP2[pid]+count2)*3+1]=I.y;
          intersectionsP[(psP2[pid]+count2)*3+2]=alpha;
          neighborMapP[psP2[pid]+count2]=indexIntP+count1-1;
          neighborP[indexIntP+count1-1]=psP2[pid]+count2;                    //neighbor of new vertex
          break;
        // X-overlap
        case 5:
          // printf("** %d %d %d %d\n", (psP2[pid]+count2), indexIntP, count1-1, psP2[pid]+count2);
          intersectionsP[(psP2[pid]+count2)*3]=Q1.x;
          intersectionsP[(psP2[pid]+count2)*3+1]=Q1.y;
          intersectionsP[(psP2[pid]+count2)*3+2]=alpha;
          neighborMapP[psP2[pid]+count2]=indexIntP+count1-1;
          neighborP[indexIntP+count1-1]=psP2[pid]+count2;                    //neighbor of new vertex
          break;
        // case T_INTERSECTION_Q:
        // case T_OVERLAP_Q:
        case 2:
        case 6:
          intersectionsP[psP2[pid]*3+2]=alpha;          //***** error prone. Did not checked in depth
          neighborMapP[psP2[pid]+count2]=indexIntP+count1-1;
          neighborP[indexIntP+count1-1]=psP2[pid]+count2;                    //neighbor of new vertex
        break;
        // case T_INTERSECTION_P:
        // case T_OVERLAP_P:
        case 3:
        case 7:
          // printf("*** %d %d %d %d\n", (psP2[pid]+count2), indexIntP, count1-1, psP2[pid]+count2);
          intersectionsP[(psP2[pid]+count2)*3]=Q1.x;
          intersectionsP[(psP2[pid]+count2)*3+1]=Q1.y;
          intersectionsP[(psP2[pid]+count2)*3+2]=alpha;
          neighborMapP[psP2[pid]+count2]=indexIntP+count1-1;
          neighborP[indexIntP+count1-1]=psP2[pid]+count2;                    //neighbor of new vertex
          break;
        // case V_INTERSECTION:
        // case V_OVERLAP:
        case 4:
        case 8:
          // printf("**** %d %d %d %d\n", (psP2[pid]+count2), indexIntP, count1-1, psP2[pid]+count2);
          intersectionsP[psP2[pid]*3+2]=alpha;          //***** error prone. Did not checked in depth
          neighborMapP[psP2[pid]+count2]=indexIntP+count1-1;
          neighborP[indexIntP+count1-1]=psP2[pid]+count2;                    //neighbor of new vertex
          break;
      } 
    } else if(i && id>=sizeP){
      initLabelsQ[(psQ2[pid]+count2)]=-100;    //make init label to default -100 
      count1++;
        if(i==1 || i==3 || i==5 || i==7)
        // if(i==1 || i==2 || i==5 || i==6)
          count2++;
      switch(i) {
        //
        // X-intersection
        //
        // case X_INTERSECTION:
        case 1:
          I = add(mulScalar((1.0-alpha), P1), mulScalar(alpha, P2));
          I.x=getValueTolarence(I.x);
          I.y=getValueTolarence(I.y);
          // printf("/* %d %d %d %d %f %f\n", (psQ2[pid]+count2), indexIntQ, count1-1, psQ2[pid]+count2, I.x, I.y);
          intersectionsQ[(psQ2[pid]+count2)*3]=I.x;       //consider edge for the intersection array
          intersectionsQ[(psQ2[pid]+count2)*3+1]=I.y;
          intersectionsQ[(psQ2[pid]+count2)*3+2]=alpha;
          neighborMapQ[psQ2[pid]+count2]=indexIntQ+count1-1;
          neighborQ[indexIntQ+count1-1]=psQ2[pid]+count2;                    //neighbor of new vertex
          break;
        // case X_OVERLAP:
        case 5:
          // printf("/** %d %d %d %d\n", (psQ2[pid]+count2), indexIntQ, count1-1, psQ2[pid]+count2);
          intersectionsQ[(psQ2[pid]+count2)*3]=Q1.x;    
          intersectionsQ[(psQ2[pid]+count2)*3+1]=Q1.y;
          intersectionsQ[(psQ2[pid]+count2)*3+2]=beta;
          neighborMapQ[psQ2[pid]+count2]=indexIntQ+count1-1;
          neighborQ[indexIntQ+count1-1]=psQ2[pid]+count2;                    //neighbor of new vertex
          break;
        // case T_INTERSECTION_Q:
        // case T_OVERLAP_Q: 
        // was 2, 6
        case 3:
        case 7:
          // printf("/*** %d %d %d %d\n", (psQ2[pid]+count2), indexIntQ, count1-1, psQ2[pid]+count2);
          intersectionsQ[(psQ2[pid]+count2)*3]=Q1.x;
          intersectionsQ[(psQ2[pid]+count2)*3+1]=Q1.y;
          intersectionsQ[(psQ2[pid]+count2)*3+2]=alpha;   
          neighborMapQ[psQ2[pid]+count2]=indexIntQ+count1-1;
          neighborQ[indexIntQ+count1-1]=psQ2[pid]+count2;                    //neighbor of new vertex
          break;
        // case T_INTERSECTION_P:
        // case T_OVERLAP_P:
        // was 3, 7
        case 2:
        case 6:
          intersectionsQ[psQ2[pid]*3+2]=alpha;          //***** error prone. Did not checked in depth
          neighborMapQ[psQ2[pid]+count2]=indexIntQ+count1-1;
          neighborQ[indexIntQ+count1-1]=psQ2[pid]+count2;                    //neighbor of new vertex
        break;
        // case V_INTERSECTION:
        // case V_OVERLAP:
        case 4:
        case 8:
          // printf("/**** %d %d %d %d\n", (psQ2[pid]+count2), indexIntQ, count1-1, psQ2[pid]+count2);
          intersectionsQ[psQ2[pid]*3+2]=alpha;          //***** error prone. Did not checked in depth
          neighborMapQ[psQ2[pid]+count2]=indexIntQ+count1-1;
          neighborQ[indexIntQ+count1-1]=psQ2[pid]+count2;                    //neighbor of new vertex
          break;
      } 
    }
  }

  // Apply initial label
  __syncthreads();
  // noly worls with the intersectionP part and copy lables to Q later
  if(id>=sizeP) return;
  int start=psP2[pid], end=psP2[pid+1];
  int tmpId, nId, pMNId, pPNId;
  point pM, pP, qM, qP, current;
  int qMType, qPType, tmpIniLabel;

  for (int i = start; i < end; i++){
    initLabelsP[i]=-100;
    if(intersectionsP[i*3+2]!=-100){    //consider intersections only
      current.x=intersectionsP[i*3]; 
      current.y=intersectionsP[i*3+1]; 
      tmpId=getCircularId(i-1, sizeNP);
      // determine local configuration at this intersection vertex
      pM.x=intersectionsP[tmpId*3];                // P-, predecessor of I on P
      pM.y=intersectionsP[tmpId*3+1];                // P-, predecessor of I on P
      if(intersectionsP[tmpId*3+2]!=-100)
        pMNId=getNeighborIndex(tmpId, neighborMapP, neighborQ); //get neighbor id of P_m vertex
      else pMNId=-100;

      tmpId=getCircularId(i+1, sizeNP);
      pP.x=intersectionsP[tmpId*3];                // P+, successor of I on P
      pP.y=intersectionsP[tmpId*3+1];                // P+, successor of I on P
      if(intersectionsP[tmpId*3+2]!=-100)
        pPNId=getNeighborIndex(tmpId, neighborMapP, neighborQ); //get neighbor id of P_p vertex
      else pPNId=-100;

      nId=getNeighborIndex(i, neighborMapP, neighborQ);
      tmpId=getCircularId(nId-1, sizeNQ);
      qM.x=intersectionsQ[tmpId*3];     // Q-, predecessor of I on Q
      qM.y=intersectionsQ[tmpId*3+1];     // Q-, predecessor of I on Q
      qMType=oracle(pMNId, pPNId, tmpId, qM, pM, current, pP);

      tmpId=getCircularId(nId+1, sizeNQ);
      qP.x=intersectionsQ[tmpId*3];     // Q+, successor of I on P
      qP.y=intersectionsQ[tmpId*3+1];     // Q+, successor of I on P
      qPType=oracle(pMNId, pPNId, tmpId, qP, pM, current, pP);

      tmpIniLabel=getInitialLabel(qMType, qPType);
      initLabelsP[i]=tmpIniLabel;
      initLabelsQ[nId]=tmpIniLabel;
      // printf("%d %d (%f, %f) (%f, %f) (%f, %f) (%f, %f)\n", i, nId, pM.x, pM.y, pP.x, pP.y, qM.x, qM.y, qP.x, qP.y);
      // printf(">>> %d %d %d %d\n", i, qMType, qPType, getInitialLabel(qMType, qPType));
    }
  }
}

/* this method is replaced by deviceCalculateIntersections() remove this in future to clean the code
/*
-----------------------------------------------------------------
Function to calculate all intersections
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
/*
__global__ void intersect(double *polyPX, double *polyPY, double *polyQX,  double *polyQY, double *dev_intersectionsP, double *dev_intersectionsQ, int sizeP, int sizeQ){
  int id=threadIdx.x;
  double alpha;
  double beta;
  point I;

  if(id>sizeP*sizeQ) return;

  int count[9] = {0,0,0,0,0,0,0,0,0};
  point P1, P2, Q1, Q2;
  int pid=id/sizeQ;
  int qid=(id+1)%sizeQ;

  // neighbor types = {P1->0, Q1->1, I_P->2, I_Q->3}
  // arrays to save neighbor information
  int *neighborP;
  int *neighborQ;
  neighborP=(int *) malloc(sizeP*sizeQ*sizeof(int));
  neighborQ=(int *) malloc(sizeP*sizeQ*sizeof(int));

  P1.x = polyPX[pid];
  P1.y = polyPY[pid];
  P2.x = polyPX[pid+1];
  P2.y = polyPY[pid+1];
  Q1.x = polyQX[qid];
  Q1.y = polyQY[qid];
  Q2.x = polyQX[qid+1];
  Q2.y = polyQY[qid+1];

  // reset P2 vertex of last edge to first vertex
  if(qid == sizeQ-1){
    Q2.x = polyQX[0];
    Q2.y = polyQY[0];
  }
  if(pid == sizeP-1){
    P2.x = polyPX[0];
    P2.y = polyPY[0];
  }

  // printf("%f %f %f %f %f %f %f %f \n", P1.x, P1.y, P2.x, P2.y, Q1.x, Q1.y, Q2.x, Q2.y);
  //
  // determine intersection or overlap type
  //
  // IntersectionType i = intersect(edgeP, edgeQ, alpha, beta);
  int i = getIntersectType(P1, P2, Q1, Q2, alpha, beta);
  printf("type %d\n", i);

  // shared variable. removed temporarily
  // count[i]++;

  // vertex* P1 = edgeP.one;
  // vertex* Q1 = edgeQ.one;

  switch(i) {
  //
  // X-intersection
  //
  // case X_INTERSECTION:
  case 1:
    // I = (1.0-alpha)*edgeP.one->p + alpha*edgeP.two->p;
    I = add(mulScalar((1.0-alpha), P1), mulScalar(alpha, P2));
    // I_P = new vertex(I,alpha);
    // I_Q = new vertex(I,beta);
    // insertVertex(I_P, edgeP);
    // insertVertex(I_Q, edgeQ);
    // link(I_P, I_Q);
    // printf("innn %f %f\n", I.x, I.y);
    dev_intersectionsP[id*3]=I.x;       //consider edge for the intersection array
    dev_intersectionsP[id*3+1]=I.y;
    dev_intersectionsP[id*3+2]=alpha;
    neighborP[id]=3;                    //consider I_Q for the neighbor 

    dev_intersectionsQ[id*3]=I.x;
    dev_intersectionsQ[id*3+1]=I.y;
    dev_intersectionsQ[id*3+2]=beta;
    neighborQ[id]=2;
    break;
  
  //
  // X-overlap
  //
  // case X_OVERLAP:
  case 5:
    // I_Q = new vertex(P1->p, beta);
    // insertVertex(I_Q, edgeQ);
    // link(P1, I_Q);

    // I_P = new vertex(Q1->p, alpha);
    // insertVertex(I_P, edgeP);
    // link(I_P, Q1);
    dev_intersectionsQ[id*3]=P1.x;    
    dev_intersectionsQ[id*3+1]=P1.y;
    dev_intersectionsQ[id*3+2]=beta;
    neighborQ[id]=0;  

    dev_intersectionsP[id*3]=Q1.x;
    dev_intersectionsP[id*3+1]=Q1.y;
    dev_intersectionsP[id*3+2]=alpha;
    neighborP[id]=1;
    break;

  //
  // T-intersection or T_overlap on Q
  //
  // case T_INTERSECTION_Q:
  // case T_OVERLAP_Q:
  case 2:
  case 6:
    // I_Q = new vertex(P1->p, beta);
    // insertVertex(I_Q, edgeQ);
    // link(P1, I_Q);

    dev_intersectionsQ[id*3]=P1.x;
    dev_intersectionsQ[id*3+1]=P1.y;
    dev_intersectionsQ[id*3+2]=beta;   
    neighborQ[id]=0;
    break;

  //
  // T-intersection or T-overlap on P
  //
  // case T_INTERSECTION_P:
  // case T_OVERLAP_P:
  case 3:
  case 7:
    // I_P = new vertex(Q1->p, alpha);
    // insertVertex(I_P, edgeP);
    // link(I_P, Q1);

    dev_intersectionsP[id*3]=Q1.x;
    dev_intersectionsP[id*3+1]=Q1.y;
    dev_intersectionsP[id*3+2]=alpha;
    neighborP[id]=1;
    break;

  //
  // V-intersection or V-overlap
  //
  // case V_INTERSECTION:
  // case V_OVERLAP:
  case 4:
  case 8:
    // link(P1,Q1);
    neighborP[id]=1;
    neighborQ[id]=0;
    break;
  } 

}*/

/*
-----------------------------------------------------------------
Function to get next id of a given vertex 
Runs in GPU
Called from Device
polyType 0 is polygon P and 1 if polygon Q
-------------------------------------------------------------------
*/
// __device__ int getNextID(int threadID, int polyType, long size){
//   if(polyType == 0) return (size+threadID/size+1)%size;

//   return (size+threadID%size+1)%size;
// }

/*
-----------------------------------------------------------------
Function to label all instersection points
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
/*
__device__ void labelIntersectionPoints(double *polyPX, double *polyPY, double *polyQX,  double *polyQY, double *dev_intersectionsP, double *dev_intersectionsQ, int *neighborP, int *neighborQ, int sizeP, int sizeQ){
  int threadID=threadIdx.x;
  // if the current thread does not have an intersection point to label, exit the function
  if()
    return;

  // pass type of the current point (P/Q/intersection_p/intersection_Q)
  int previousID=getPreviousID(threadID, 2, sizeQ);
  int nextID=getNextID(threadID, 2, sizeQ);

  // intersection point related to current thread
  point I;
  I.x=intersectionsP[id*3];
  I.y=intersectionsP[id*3+1];
  // determine local configuration at this intersection vertex
  point P_m, P_p, Q_m, Q_p
  P_m.x = polyPX[previousID]; //I->prev;                // P-, predecessor of I on P
  P_m.y = polyPY[previousID]; //I->prev;                // P-, predecessor of I on P
  P_p.x = polyPX[nextID]; //I->next;                // P+, successor of I on P
  P_p.y = polyPY[nextID]; //I->next;                // P+, successor of I on P
  
  int neighborPreviousID=getPreviousID(threadID, neighborP[threadID], sizeQ);
  int neighborNextID=getNextID(threadID, neighborP[threadID], sizeQ);
  // neighbor changes based on the neighbor arrays
  // Q_m.x = polyQX[id]; //I->neighbour->prev;     // Q-, predecessor of I on Q
  // Q_m.y = polyQY[id]; //I->neighbour->prev;     // Q-, predecessor of I on Q
  // Q_p.x = polyQX[id+1]; //I->neighbour->next;     // Q+, successor of I on P
  // Q_p.y = polyQY[id+1]; //I->neighbour->next;     // Q+, successor of I on P

  if(neighborP[threadID] == 0){
    Q_m.x = polyPX[neighborPreviousID]; //I->neighbour->prev;     // Q-, predecessor of I on Q
    Q_m.y = polyPY[neighborPreviousID]; //I->neighbour->prev;     // Q-, predecessor of I on Q
    Q_p.x = polyPX[neighborNextID]; //I->neighbour->next;     // Q+, successor of I on P
    Q_p.y = polyPY[neighborNextID]; //I->neighbour->next;     // Q+, successor of I on P
  } else if(neighborP[threadID] == 1){
    Q_m.x = polyQX[neighborPreviousID]; //I->neighbour->prev;     // Q-, predecessor of I on Q
    Q_m.y = polyQY[neighborPreviousID]; //I->neighbour->prev;     // Q-, predecessor of I on Q
    Q_p.x = polyQX[neighborNextID]; //I->neighbour->next;     // Q+, successor of I on P
    Q_p.y = polyQY[neighborNextID]; //I->neighbour->next;     // Q+, successor of I on P
  } else if(neighborP[threadID] == 2){
    Q_m.x = polyPX[neighborPreviousID]; //I->neighbour->prev;     // Q-, predecessor of I on Q
    Q_m.y = polyPY[neighborPreviousID]; //I->neighbour->prev;     // Q-, predecessor of I on Q
    Q_p.x = polyPX[neighborNextID]; //I->neighbour->next;     // Q+, successor of I on P
    Q_p.y = polyPY[neighborNextID]; //I->neighbour->next;     // Q+, successor of I on P
  } else if(neighborP[threadID] == 3){
    Q_m.x = polyQX[neighborPreviousID]; //I->neighbour->prev;     // Q-, predecessor of I on Q
    Q_m.y = polyQY[neighborPreviousID]; //I->neighbour->prev;     // Q-, predecessor of I on Q
    Q_p.x = polyQX[neighborNextID]; //I->neighbour->next;     // Q+, successor of I on P
    Q_p.y = polyQY[neighborNextID]; //I->neighbour->next;     // Q+, successor of I on P
  }

  // check positions of Q- and Q+ relative to (P-, I, P+)
  int Q_m_type = oracle(Q_m, P_m, I, P_p);
  int Q_p_type = oracle(Q_p, P_m, I, P_p);

}
*/
__global__ void hellworld(int a, int *x) {
  printf("HIIIII \n");
  *x=a+10;
}

/*
-----------------------------------------------------------------
Function to get circular id of a given id 
Runs in CPU
Called from Host
-------------------------------------------------------------------
*/
int getCircularIdHost(int id, int maxCount){
  if(maxCount==id) return 0;
  else if(id==-1) return maxCount-1;
  else return id;
}

// count how many intersection points and prefix sums
void countIntersections(double *polyPX, double *polyPY, double *polyQX,  double *polyQY, int sizeP, int sizeQ, int *countNonDegenIntP, int *countNonDegenIntQ, double **intersectionsP, double **intersectionsQ, int **initLabelsP, int **initLabelsQ){
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

  int threadsPerBlock = 256;
  int blocksPerGrid = (2*(sizeP+sizeQ) + threadsPerBlock - 1) / threadsPerBlock;
  dim3 dimBlock(threadsPerBlock, 1, 1), dimGrid(blocksPerGrid, 1, 1); 

  deviceCountIntersections<<<dimGrid, dimBlock>>>(dev_polyPX, dev_polyPY, dev_polyQX, dev_polyQY, sizeP, sizeQ, dev_psP1, dev_psP2, dev_psQ1, dev_psQ2);

  cudaMemcpy(&psP1, dev_psP1, (sizeP+1)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psP2, dev_psP2, (sizeP+1)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psQ1, dev_psQ1, (sizeQ+1)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&psQ2, dev_psQ2, (sizeQ+1)*sizeof(int), cudaMemcpyDeviceToHost);

  // for (int i = 0; i < sizeQ+1; ++i)
  // {
  //   printf(" %d ", psQ2[i]);
  // }
  // printf("--- \n");

  // for (int i = 0; i < sizeQ+1; ++i)
  // {
  //   printf(" %d ", psQ1[i]);
  // }
  // printf("--- \n");
  // Phase 2: Calcualte intersections and save them in the arrays. Make neighbor connections
  // int countIntersections=psP1[sizeP-1];
  // int countNonDegenIntP=psP2[sizeP-1]+(psP1[sizeP-1]-psP1[sizeP-2])+1;
  // int countNonDegenIntQ=psQ2[sizeQ-1]+(psQ1[sizeQ-1]-psQ1[sizeQ-2])+1;

  int countIntersections=psP1[sizeP];
  *countNonDegenIntP=psP2[sizeP];
  *countNonDegenIntQ=psQ2[sizeQ];
  // printf("non degen size %d %d\n", countNonDegenIntP, countNonDegenIntQ);
  // double intersectionsP[*countNonDegenIntP*3], intersectionsQ[*countNonDegenIntQ*3];
  int neighborP[countIntersections], neighborQ[countIntersections], neighborMapP[*countNonDegenIntP], neighborMapQ[*countNonDegenIntQ];
  // int initLabelsP[*countNonDegenIntP], initLabelsQ[*countNonDegenIntQ];
  double *dev_intersectionsP, *dev_intersectionsQ;
  int *dev_neighborP, *dev_neighborQ, *dev_neighborMapP, *dev_neighborMapQ, *dev_initLabelsP, *dev_initLabelsQ;

  *intersectionsP=(double *)malloc(*countNonDegenIntP*3*sizeof(double));
  *intersectionsQ=(double *)malloc(*countNonDegenIntQ*3*sizeof(double));
  *initLabelsP=(int *)malloc(*countNonDegenIntP*sizeof(int));
  *initLabelsQ=(int *)malloc(*countNonDegenIntQ*sizeof(int));

  // Allocate memory in device 
  cudaMalloc((void **) &dev_intersectionsP, *countNonDegenIntP*3*sizeof(double));
  cudaMalloc((void **) &dev_intersectionsQ, *countNonDegenIntQ*3*sizeof(double));
  cudaMalloc((void **) &dev_neighborP, countIntersections*sizeof(int));
  cudaMalloc((void **) &dev_neighborQ, countIntersections*sizeof(int));
  cudaMalloc((void **) &dev_neighborMapP, *countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_neighborMapQ, *countNonDegenIntQ*sizeof(int));
  cudaMalloc((void **) &dev_initLabelsP, *countNonDegenIntP*sizeof(int));
  cudaMalloc((void **) &dev_initLabelsQ, *countNonDegenIntQ*sizeof(int));

  deviceCalculateIntersections<<<dimGrid, dimBlock>>>(dev_polyPX, dev_polyPY, dev_polyQX, dev_polyQY, sizeP, sizeQ, dev_psP1, dev_psP2, dev_psQ1, dev_psQ2, dev_intersectionsP, dev_intersectionsQ, dev_neighborP, dev_neighborQ, dev_neighborMapP, dev_neighborMapQ, *countNonDegenIntP, *countNonDegenIntQ, dev_initLabelsP, dev_initLabelsQ);

  cudaMemcpy(*intersectionsP, dev_intersectionsP, *countNonDegenIntP*3*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(*intersectionsQ, dev_intersectionsQ, *countNonDegenIntQ*3*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(&neighborP, dev_neighborP, countIntersections*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&neighborQ, dev_neighborQ, countIntersections*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&neighborMapP, dev_neighborMapP, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&neighborMapQ, dev_neighborMapQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*initLabelsP, dev_initLabelsP, *countNonDegenIntP*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(*initLabelsQ, dev_initLabelsQ, *countNonDegenIntQ*sizeof(int), cudaMemcpyDeviceToHost);

  cudaDeviceSynchronize();

  // load labeled intersection arrays to PP and QQ
  /*int prevId;
  polygon P[countNonDegenIntP], Q[countNonDegenIntQ];
  vertex tmpVertex;
  for(int i=0; i<countNonDegenIntP; ++i){
    tmpVertex.p.x=intersectionsP[i*3];
    tmpVertex.p.y=intersectionsP[i*3+1];
    tmpVertex.alpha=intersectionsP[i*3+2];              
    tmpVertex.label=initLabelsP[i]; 
    // bool source; 
    if(intersectionsP[i*3+2]!=-100){
      tmpVertex.intersection=true;  
      P[i]->neighbor=Q[neighborMapP[i]];
    }
    else tmpVertex.intersection=false;  
    P[i]=tmpVertex;
    prevId=getCircularId(i-1, countNonDegenIntP);
    P[i]->prev=P[prevId];
    P[prevId]->next=P[i];
  }*/

  // for(int i=0; i<countNonDegenIntQ; ++i){
  //   tmpVertex.next=;
  //   tmpVertex.neighbour=;
  // }
/*
  printf("intersectionP");
  for (int i = 0; i < *countNonDegenIntP*3; ++i)
  {
    if(i%3==0)
      printf("\n%d", i/3);
    // printf(" %f ", intersectionsP[i]);
    printf(" %f ", *(*intersectionsP+i));
  }
  printf("\n\nintersectionQ");
  for (int i = 0; i < *countNonDegenIntQ*3; ++i)
  {
    if(i%3==0)
      printf("\n%d", i/3);
    printf(" %f ", *(*intersectionsQ+i));
  }
  printf("\n");/*
  for (int i = 0; i < countIntersections; ++i)
  {
    printf(" %d ", neighborP[i]);
  }
  printf("\n");
  for (int i = 0; i < countIntersections; ++i)
  {
    printf(" %d ", neighborQ[i]);
  }
  printf("\n");
  for (int i = 0; i < countNonDegenIntP; ++i)
  {
    printf(" %d-%d ", i, neighborMapP[i]);
  }
  printf("\n");
  for (int i = 0; i < countNonDegenIntQ; ++i)
  {
    printf(" %d-%d ", i, neighborMapQ[i]);
  }
  printf("\n");
  for (int i = 0; i < countNonDegenIntP; ++i)
  {
    printf(" %d>%d ", i, initLabelsP[i]);
  }
  printf("\n");
  for (int i = 0; i < countNonDegenIntQ; ++i)
  {
    printf(" %d>%d ", i, initLabelsQ[i]);
  }
  printf("\n");
*/


  cudaFree(dev_polyPX);
  cudaFree(dev_polyPY);
  cudaFree(dev_polyQX);
  cudaFree(dev_polyQY);
}
// old method replaced by countIntersections() 2nd phase
/*
void calculateIntersections(double *polyPX, double *polyPY, double *polyQX,  double *polyQY, int sizeP, int sizeQ){

  // int j;
  // for(j=0; j<sizeP; j++){
  //   printf("-> %f ", polyPX[j]);
  // }
  // printf("\n");

  double *dev_polyPX, *dev_polyPY, *dev_polyQX, *dev_polyQY;
  double *dev_intersectionsP, *dev_intersectionsQ, *dev_neighbors;
  double intersectionsP[sizeP*sizeQ*3], intersectionsQ[sizeP*sizeQ*3];

  // Allocate memory in device 
  cudaMalloc((void **) &dev_polyPX, sizeP*sizeof(double));
  cudaMalloc((void **) &dev_polyPY, sizeP*sizeof(double));
  cudaMalloc((void **) &dev_polyQX, sizeQ*sizeof(double));
  cudaMalloc((void **) &dev_polyQY, sizeQ*sizeof(double));
  cudaMalloc((void **) &dev_intersectionsP, 3*sizeP*sizeQ*sizeof(double));
  cudaMalloc((void **) &dev_intersectionsQ, 3*sizeP*sizeQ*sizeof(double));

  // Copy input vectors from host memory to GPU buffers.
  cudaMemcpy(dev_polyPX, polyPX, sizeP*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_polyPY, polyPY, sizeP*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_polyQX, polyQX, sizeQ*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_polyQY, polyQY, sizeQ*sizeof(double), cudaMemcpyHostToDevice);

  printf("in cpu before\n");

  // need to chage shape of kernel lanch******
  int threadsPerBlock = 256;
  int blocksPerGrid = (2*(sizeP+sizeQ) + threadsPerBlock - 1) / threadsPerBlock;
  dim3 dimBlock(threadsPerBlock, 1, 1), dimGrid(blocksPerGrid, 1, 1); 
  intersect <<<dimGrid, dimBlock>>> (dev_polyPX, dev_polyPY, dev_polyQX, dev_polyQY, dev_intersectionsP, dev_intersectionsQ, sizeP, sizeQ);

  cudaMemcpy(&intersectionsP, dev_intersectionsP, 3*sizeP*sizeQ*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(&intersectionsQ, dev_intersectionsQ, 3*sizeP*sizeQ*sizeof(double), cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();

  printf("in cpu after\n");

  int j;
  for(j=0; sizeP*sizeQ>j; ++j){
    printf(">> %f %f %f\n", intersectionsP[j*3], intersectionsP[j*3+1], intersectionsP[j*3+2]);
  }

  cudaFree(dev_polyPX);
  cudaFree(dev_polyPY);
  cudaFree(dev_polyQX);
  cudaFree(dev_polyQY);
}*/

void testhello() {
  printf("in cpu before\n");
  int x;
  int *dev_x;
  cudaMalloc((void**)&dev_x, sizeof(int));

  hellworld <<<1,2>>> (11, dev_x);

  cudaMemcpy(&x, dev_x, sizeof(int), cudaMemcpyDeviceToHost);
  printf("***== %d\n", x);

  cudaDeviceSynchronize();

  printf("in cpu after\n");
}

// int main(){
//   testhello();
//   // cudaDeviceSynchronize();

//   return 0;
// }