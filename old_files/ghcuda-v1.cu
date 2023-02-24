#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define EPSILON  0.000000001 

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
Function to return intersection  type
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/
// __device__ int getIntersectType(const edge& edgeP, const edge& edgeQ, double& alpha, double& beta) {
__device__ int getIntersectType(const point& P1, const point& P2, const point& Q1, const point& Q2,  double& alpha, double& beta) {
  // const point2D& P1 = edgeP.one->p;
  // const point2D& P2 = edgeP.two->p;
  // const point2D& Q1 = edgeQ.one->p;
  // const point2D& Q2 = edgeQ.two->p;

	double AP1 = A(P1,Q1,Q2);
	double AP2 = A(P2,Q1,Q2);

	if (fabs(AP1-AP2) > EPSILON) {
		// from here: [P1,P2] and [Q1,Q2] are not parallel

		//
		// analyse potential intersection
		//

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

    //
		// distinguish intersection types
    //
    
		if (alpha_in_0_1 && beta_in_0_1) return (1);
			// return (X_INTERSECTION);
			

		if (alpha_is_0 && beta_in_0_1) return (2);
			// return (T_INTERSECTION_Q);

		if (beta_is_0 && alpha_in_0_1) return (3);
			// return (T_INTERSECTION_P);

		if (alpha_is_0 && beta_is_0) return (4);
			// return (V_INTERSECTION);
	}
	else
		if (fabs(AP1) < EPSILON) {
			// from here: [P1,P2] and [Q1,Q2] are collinear

			//
			// analyse potential overlap
			//

			// point2D dP = P2-P1;
			// point2D dQ = Q2-Q1;
			// point2D PQ = Q1-P1;
			
      point dP = sub(P2, P1);
			point dQ = sub(Q2, Q1);
			point PQ = sub(Q1, P1);


			// compute alpha and beta
			// alpha = (PQ*dP) / (dP*dP);
			// beta = -(PQ*dQ) / (dQ*dQ);

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

      //
			// distinguish intersection types
      //

			if (alpha_in_0_1 && beta_in_0_1) return (5);
				// return (X_OVERLAP);

			if (alpha_not_in_0_1 && beta_in_0_1) return (6);
				// return (T_OVERLAP_Q);

			if (beta_not_in_0_1 && alpha_in_0_1) return (7);
				// return (T_OVERLAP_P);

			if (alpha_is_0 && beta_is_0) return (8);
				// return (V_OVERLAP);
		}

  return (0);
	// return (NO_INTERSECTION); 
}

/*
-----------------------------------------------------------------
Function to calculate all intersections
Runs in GPU
Called from Host
-------------------------------------------------------------------
*/
 
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

}

/*
-----------------------------------------------------------------
Function to get previous id of a given vertex 
Runs in GPU
Called from Device
polyType 0 is polygon P and 1 if polygon Q if 2 intersection P if 3 intersection Q
-------------------------------------------------------------------
*/

__device__ int getPreviousID(int threadID, int polyType, long size){
  if(polyType == 0) return (size+threadID/size-1)%size;       //polygon P
  else if(polyType == 1) return (size+threadID%size-1)%size;  //polygon Q
  else if(polyType == 2) return threadID/size;                //intersectionP
  else if(polyType == 3) return threadID%size;                //interscetionQ
}

/*
-----------------------------------------------------------------
Function to get next id of a given vertex 
Runs in GPU
Called from Device
polyType 0 is polygon P and 1 if polygon Q
-------------------------------------------------------------------
*/

__device__ int getNextID(int threadID, int polyType, long size){
  if(polyType == 0) return (size+threadID/size+1)%size;

  return (size+threadID%size+1)%size;
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

int oracle(int threadID, point* Q, point* P1, point* P2, point* P3) {

  // is Q linked to P1 ?
  // if (P1->intersection && (P1->neighbour == Q))
  if(neighbor)
    return(2);

  // is Q linked to P2 ?
  // if (P3->intersection && (P3->neighbour == Q))
  if(neighborP[threadID] == )
    return(3);

  // check relative position of Q with respect to chain (P1,P2,P3)
  double s1 = A(Q, P1, P2);
  double s2 = A(Q, P2, P3);
  double s3 = A(P1, P2, P3);

  if (s3 > 0) { 
    // chain makes a left turn
    if (s1 > 0 && s2 > 0)
      return(0);
    else
      return(1);
  }
  else {
    // chain makes a right turn (or is straight)
    if (s1 < 0 && s2 < 0)
      return(1);
    else
      return(0);
  }
}

/*
-----------------------------------------------------------------
Function to label all instersection points
Runs in GPU
Called from Device
-------------------------------------------------------------------
*/

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

__global__ void hellworld(int a, int *x) {
  printf("HIIIII \n");
  *x=a+10;
}
 
void calculateIntersections(double *polyPX, double *polyPY, double *polyQX,  double *polyQY, int sizeP, int sizeQ){

  // int j;
  // for(j=0; j<sizeP; j++){
  //   printf("-> %f ", polyPX[j]);
  // }
  // printf("\n");

  double *dev_polyPX, *dev_polyPY, *dev_polyQX, *dev_polyQY;
  double *dev_intersectionsP, *dev_intersectionsQ;
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
  intersect <<<1, sizeP*sizeQ>>> (dev_polyPX, dev_polyPY, dev_polyQX, dev_polyQY, dev_intersectionsP, dev_intersectionsQ, sizeP, sizeQ);

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
}

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