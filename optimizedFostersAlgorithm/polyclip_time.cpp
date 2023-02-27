////////////////////////////////////////////////////////////////////////
//
// polyclip.cpp
//
// computes the intersection (or union) of two non-self-intersecting complex
// polygons, with possibly multiple and nested components, even in case of
// degenerate intersections (vertex on edge, overlapping edges, etc.)
//
// compile:
//
// g++ -std=c++11 polyclip_time.cpp -o polyclip
// g++ -O2 -std=c++11 polyclip_time.cpp -o polyclip
//
// usage:
//
// polyclip [-union] input1.poly input2.poly output.poly

// ./polyclip ../examples2/Fig8-P.poly ../examples2/Fig8-Q.poly results/Fig8-R.poly
// ./polyclip ../examples/Fig14-P.poly ../examples/Fig14-Q.poly results/Fig14-R.poly
// ./polyclip ../examples/Fig15-P.poly ../examples/Fig15-Q.poly results/Fig15-R.poly
// ./polyclip ../examples/Fig16-P.poly ../examples/Fig16-Q.poly results/Fig16-R.poly
// ./polyclip ../examples/Fig17-P.poly ../examples/Fig17-Q.poly results/Fig17-R.poly
// ./polyclip ../examples/Fig18-P.poly ../examples/Fig18-Q.poly results/Fig18-R.poly
// ./polyclip ../examples/Fig19-P.poly ../examples/Fig19-Q.poly results/Fig19-R.poly

//real-world dataset experiment 
/*
./polyclip ../examples2/s.txt ../examples2/c.txt results/Fig-large-R.poly
*/
// synthetic dataset experiment
/*
./polyclip ../../datasets/GH-Synthetic/syntheticO_n/lakes_851348.txt ../../datasets/GH-Synthetic/syntheticO_n/lakes-synthetic_1_1.txt results/Fig-large-R.poly
./polyclip ../examples2/s.txt ../../datasets/GH-Synthetic/syntheticO_n/s-synthetic_15_10.txt results/Fig-large-R.poly
./polyclip ../examples2/c.txt ../../datasets/GH-Synthetic/syntheticO_n/c-synthetic_100_0.txt results/Fig-large-R.poly

./polyclip ../../datasets/GH-Synthetic/syntheticAlltoALL/worst-syntheticP_500.txt ../../datasets/GH-Synthetic/syntheticAlltoALL/worst-syntheticQ_500.txt results/Fig-large-R.poly
./polyclip ../../datasets/GH-Synthetic/syntheticAlltoALL/worst-syntheticP_700.txt ../../datasets/GH-Synthetic/syntheticAlltoALL/worst-syntheticQ_700.txt results/Fig-large-R.poly
./polyclip ../../datasets/GH-Synthetic/syntheticAlltoALL/worst-syntheticP_1000.txt ../../datasets/GH-Synthetic/syntheticAlltoALL/worst-syntheticQ_1000.txt results/Fig-large-R.poly
./polyclip ../../datasets/GH-Synthetic/syntheticAlltoALL/worst-syntheticP_1500.txt ../../datasets/GH-Synthetic/syntheticAlltoALL/worst-syntheticQ_1500.txt results/Fig-large-R.poly
./polyclip ../../datasets/GH-Synthetic/syntheticAlltoALL/worst-syntheticP_2000.txt ../../datasets/GH-Synthetic/syntheticAlltoALL/worst-syntheticQ_2000.txt results/Fig-large-R.poly
*/

// ./polyclip ../examples2/Fig8-P.poly ../examples2/Fig8-Q.poly results/Fig8-R.poly

////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <chrono>
#include <algorithm>


using namespace std;
using namespace std::chrono;

#include "point2D.h"
#include "polygon.h"
#include "readSCPolygon.h"

////////////////////////////////////////////////////////////////////////

#define EPSILON  0.000000001            // tolerance

vector<polygon> PP, QQ;                 // two input polygons
vector<vector<int>> PP_Disjoints, QQ_Disjoints;                 // two polygons after cmbr
vector<polygon> RR;                     // output polygon

#include "readShapefile.cpp"

bool UNION = false;                     // global switch for computing union instead of intersection

bool DEBUG = true;   //flag to enable timing
bool DEBUG2 = false;   //flag to enable timing
bool PF = false;       //flag to enable print information
bool PR = false;       //flag to print results in a file

////////////////////////////////////////////////////////////////////////
//
// compute intersection of [P1,P2] and [Q1,Q2]
//
enum IntersectionType {                 // types of intersection (detected in the first phase)
  NO_INTERSECTION,
  X_INTERSECTION,
  T_INTERSECTION_Q,
  T_INTERSECTION_P,
  V_INTERSECTION,
  X_OVERLAP,
  T_OVERLAP_Q,
  T_OVERLAP_P,
  V_OVERLAP
};

/*
-----------------------------------------------------------------
Function to check if there is a overlap between given 2 edges 
Returns 1 if there is a overlap; else 0
-------------------------------------------------------------------
*/
int LSMF(const point2D& P1, const point2D& P2, const point2D& Q1, const point2D& Q2){
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


// get CMBR for PP and QQ
void getCMBR(double *cmbr){
  vector<double> PPMBR;
  vector<double> QQMBR;
  // double *cmbr; //minx, miny, maxx, maxy
  double minX, minY, maxX, maxY;

  PPMBR=getMBR(PP[0]);
  QQMBR=getMBR(QQ[0]);

  if(PF){
    cout<<"MBR_P ["<<PPMBR[0]<<", "<<PPMBR[1]<<", "<<PPMBR[2]<<", "<<PPMBR[3]<<endl;
    cout<<"MBR_Q ["<<QQMBR[0]<<", "<<QQMBR[1]<<", "<<QQMBR[2]<<", "<<QQMBR[3]<<endl;
  }

  // check intersection between MBRs
  if(PPMBR[0]>QQMBR[2] || PPMBR[2]<QQMBR[0]){
    printf("No Overlap between polygons\n");
    exit(0);
  }
  if(PPMBR[1]>QQMBR[3] || PPMBR[3]<QQMBR[1]){
    printf("No Overlap between polygons\n");
    exit(0);
  }

  cmbr[0]=max(PPMBR[0], QQMBR[0]);
  cmbr[1]=max(PPMBR[1], QQMBR[1]);
  cmbr[2]=min(PPMBR[2], QQMBR[2]);
  cmbr[3]=min(PPMBR[3], QQMBR[3]);
  if(PF){
    cout<<"CMBR ["<<cmbr[0]<<", "<<cmbr[1]<<", "<<cmbr[2]<<", "<<cmbr[3]<<endl;
  }
}

void CMF(){
  double *cmbr;
  double minX, minY, maxX, maxY;
  vertex* V1;
  bool isPrevEdgeValid=true;
  vector<int> tmp;
  int vid=0;

  cmbr=(double *)malloc(4*sizeof(double));
  getCMBR(cmbr);

  for (polygon& P : PP){ 
    vid=0;
    for (edge edgeP : P.edges(SOURCE)){
      minX=min(edgeP.one->p.x, edgeP.two->p.x);
      minY=min(edgeP.one->p.y, edgeP.two->p.y);
      maxX=max(edgeP.one->p.x, edgeP.two->p.x);
      maxY=max(edgeP.one->p.y, edgeP.two->p.y);
      if((minX>cmbr[2] || maxX<cmbr[0]) || (minY>cmbr[3] || maxY<cmbr[1])){
        if(!isPrevEdgeValid){
          P.removeVertex(edgeP.one);
          vid--;
        }else{
          tmp.push_back(vid);
        }
        isPrevEdgeValid=false;
      } else {
        if(!isPrevEdgeValid){
          tmp.push_back(vid);
        }
        isPrevEdgeValid=true;
      }
      vid++;
    }
    PP_Disjoints.push_back(tmp);
  }
  isPrevEdgeValid=true;
  for (polygon& Q : QQ){
    vid=0;
    for (edge edgeQ : Q.edges(SOURCE)){
      minX=min(edgeQ.one->p.x, edgeQ.two->p.x);
      minY=min(edgeQ.one->p.y, edgeQ.two->p.y);
      maxX=max(edgeQ.one->p.x, edgeQ.two->p.x);
      maxY=max(edgeQ.one->p.y, edgeQ.two->p.y);
      if((minX>cmbr[2] || maxX<cmbr[0]) || (minY>cmbr[3] || maxY<cmbr[1])){
        if(!isPrevEdgeValid){
          Q.removeVertex(edgeQ.one);
          vid--;
        }else{
          tmp.push_back(vid);
        }
        isPrevEdgeValid=false;
      } else {
        if(!isPrevEdgeValid){
          tmp.push_back(vid);
        }
        isPrevEdgeValid=true;
      }
      vid++;
    }
    QQ_Disjoints.push_back(tmp);
  }
  if(PF){
    int i=0;
    cout<<"AFTER CF ==> New P size: "<<PP[i].size;
    cout<<" New Q size: "<<QQ[i].size<<endl;
  }
}

IntersectionType intersect(const edge& edgeP, const edge& edgeQ, double& alpha, double& beta) {
  const point2D& P1 = edgeP.one->p;
  const point2D& P2 = edgeP.two->p;
  const point2D& Q1 = edgeQ.one->p;
  const point2D& Q2 = edgeQ.two->p;

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
    
		if (alpha_in_0_1 && beta_in_0_1)
			return (X_INTERSECTION);

		if (alpha_is_0 && beta_in_0_1)
			return (T_INTERSECTION_Q);

		if (beta_is_0 && alpha_in_0_1)
			return (T_INTERSECTION_P);

		if (alpha_is_0 && beta_is_0)
			return (V_INTERSECTION);
	}
	else
		if (fabs(AP1) < EPSILON) {
			// from here: [P1,P2] and [Q1,Q2] are collinear

			//
			// analyse potential overlap
			//

			point2D dP = P2-P1;
			point2D dQ = Q2-Q1;
			point2D PQ = Q1-P1;

			// compute alpha and beta
			alpha = (PQ*dP) / (dP*dP);
			beta = -(PQ*dQ) / (dQ*dQ);

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

			if (alpha_in_0_1 && beta_in_0_1)
				return (X_OVERLAP);

			if (alpha_not_in_0_1 && beta_in_0_1)
				return (T_OVERLAP_Q);

			if (beta_not_in_0_1 && alpha_in_0_1)
				return (T_OVERLAP_P);

			if (alpha_is_0 && beta_is_0)
				return (V_OVERLAP);
		}

	return (NO_INTERSECTION);
}
//
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
//
// compute all intersections
//
void computeIntersections2() {
  // cout << "Computing intersections...\n\n";
  
  point2D I;
  double alpha;
  double beta;
  vertex* I_P;
  vertex* I_Q;

  int count[9] = {0,0,0,0,0,0,0,0,0};
  
  //
  // loop over the source edges of P and Q
  //
  for (polygon& P : PP) 
    for (edge edgeP : P.edges(SOURCE))
      for (polygon& Q : QQ) 
        for (edge edgeQ : Q.edges(SOURCE)) {
          if(LSMF(edgeP.one->p , edgeP.two->p, edgeQ.one->p, edgeQ.two->p)){
            //
            // determine intersection or overlap type
            //
            IntersectionType i = intersect(edgeP, edgeQ, alpha, beta);
            // cout << "" << i << " " << edgeP.one->p << " " << edgeP.two->p<< " " << edgeQ.one->p << " " << edgeQ.two->p << " " << alpha << " " << beta << endl;
            count[i]++;
      
            vertex* P1 = edgeP.one;
            vertex* Q1 = edgeQ.one;
      
            switch(i) {
            //
            // X-intersection
            //
            case X_INTERSECTION:
              I = (1.0-alpha)*edgeP.one->p + alpha*edgeP.two->p;
              I_P = new vertex(I,alpha);
              I_Q = new vertex(I,beta);
              insertVertex(I_P, edgeP);
              insertVertex(I_Q, edgeQ);
              link(I_P, I_Q);
              break;
      
            //
            // X-overlap
            //
            case X_OVERLAP:
              I_Q = new vertex(P1->p, beta);
              insertVertex(I_Q, edgeQ);
              link(P1, I_Q);
      
              I_P = new vertex(Q1->p, alpha);
              insertVertex(I_P, edgeP);
              link(I_P, Q1);
              break;
      
            //
            // T-intersection or T_overlap on Q
            //
            case T_INTERSECTION_Q:
            case T_OVERLAP_Q:
              I_Q = new vertex(P1->p, beta);
              insertVertex(I_Q, edgeQ);
              link(P1, I_Q);
              break;
      
            //
            // T-intersection or T-overlap on P
            //
            case T_INTERSECTION_P:
            case T_OVERLAP_P:
              I_P = new vertex(Q1->p, alpha);
              insertVertex(I_P, edgeP);
              link(I_P, Q1);
              break;
      
            //
            // V-intersection or V-overlap
            //
            case V_INTERSECTION:
            case V_OVERLAP:
              link(P1,Q1);
              break;
            }
          }
        }
  if(PF){       
    cout << "\n\n--------------------------------\n"; 
    cout << "... " << count[1] << " non-degenerate and ";
    cout << count[2]+count[3]+count[4]+count[5]+count[6]+count[7]+count[8] << " degenerate intersections found" << endl;
    cout << "       (" << count[2]+count[3] << " T-intersections, ";
    cout << count[4] << " V-intersections,\n        ";
    cout << count[5] << " X-overlaps, ";
    cout << count[6]+count[7] << " T-overlaps, ";
    cout << count[8] << " V-overlaps)" << endl;
    cout << "... " << count[1]+count[5]+count[3]+count[7] << " vertices added to P" << endl;
    cout << "... " << count[1]+count[5]+count[2]+count[6] << " vertices added to Q" << endl;
  }
}
//
////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////
//
// compute all intersections
//
void computeIntersections() {
  if(PF)
    cout << "Computing intersections...\n\n";
  
  point2D I;
  double alpha;
  double beta;
  vertex* I_P;
  vertex* I_Q;

  int count[9] = {0,0,0,0,0,0,0,0,0};
  int pid=-1, peid=0, pcurrent=0, qid=-1, qeid=0, qcurrent=0;
  bool check1=true, check2=true;
  
  //
  // loop over the source edges of P and Q
  //
  for (polygon& P : PP){
    pid++;
    peid=0;
    pcurrent=0;
    for (edge edgeP : P.edges(SOURCE)){
      if(PP_Disjoints[pid].size()>0) if(peid==PP_Disjoints[pid][pcurrent]){
        if(PP_Disjoints[pid].size()-1>pcurrent) pcurrent++;
        check1=false;
      }
      if(check1){
        for (polygon& Q : QQ){ 
          qid++;
          qeid=0;
          qcurrent=0;
          for (edge edgeQ : Q.edges(SOURCE)){
            if(QQ_Disjoints[qid].size()>0) if(qeid==QQ_Disjoints[qid][qcurrent]){
              if(QQ_Disjoints[qid].size()-1>qcurrent) qcurrent++;
              check2=false;
            }
            if(check2){
              if(LSMF(edgeP.one->p , edgeP.two->p, edgeQ.one->p, edgeQ.two->p)){
                //
                // determine intersection or overlap type
                //
                IntersectionType i = intersect(edgeP, edgeQ, alpha, beta);
                count[i]++;
          
                vertex* P1 = edgeP.one;
                vertex* Q1 = edgeQ.one;
          
                switch(i){
                //
                // X-intersection
                //
                case X_INTERSECTION:
                  I = (1.0-alpha)*edgeP.one->p + alpha*edgeP.two->p;
                  I_P = new vertex(I,alpha);
                  I_Q = new vertex(I,beta);
                  insertVertex(I_P, edgeP);
                  insertVertex(I_Q, edgeQ);
                  link(I_P, I_Q);
                  break;
          
                //
                // X-overlap
                //
                case X_OVERLAP:
                  I_Q = new vertex(P1->p, beta);
                  insertVertex(I_Q, edgeQ);
                  link(P1, I_Q);
          
                  I_P = new vertex(Q1->p, alpha);
                  insertVertex(I_P, edgeP);
                  link(I_P, Q1);
                  break;
          
                //
                // T-intersection or T_overlap on Q
                //
                case T_INTERSECTION_Q:
                case T_OVERLAP_Q:
                  I_Q = new vertex(P1->p, beta);
                  insertVertex(I_Q, edgeQ);
                  link(P1, I_Q);
                  break;
          
                //
                // T-intersection or T-overlap on P
                //
                case T_INTERSECTION_P:
                case T_OVERLAP_P:
                  I_P = new vertex(Q1->p, alpha);
                  insertVertex(I_P, edgeP);
                  link(I_P, Q1);
                  break;
          
                //
                // V-intersection or V-overlap
                //
                case V_INTERSECTION:
                case V_OVERLAP:
                  link(P1,Q1);
                  break;
                }
              }
            }
          }
        }
      }
    }
    
  }
  if(PF){ 
    cout << "... " << count[1] << " non-degenerate and ";
    cout << count[2]+count[3]+count[4]+count[5]+count[6]+count[7]+count[8] << " degenerate intersections found" << endl;
    cout << "       (" << count[2]+count[3] << " T-intersections, ";
    cout << count[4] << " V-intersections,\n        ";
    cout << count[5] << " X-overlaps, ";
    cout << count[6]+count[7] << " T-overlaps, ";
    cout << count[8] << " V-overlaps)" << endl;
    cout << "... " << count[1]+count[5]+count[3]+count[7] << " vertices added to P" << endl;
    cout << "... " << count[1]+count[5]+count[2]+count[6] << " vertices added to Q" << endl;
  }
}
//
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
//
// checks the relative position of Q with respect to (P1,P2,P3)
//
enum RelativePositionType {             // types of relative positions of Q w.r.t. (P1,P2,P3)
  LEFT,
  RIGHT,
  IS_P_m,
  IS_P_p
};

RelativePositionType oracle(vertex* Q, vertex* P1, vertex* P2, vertex* P3) {

	// is Q linked to P1 ?
	if ( P1->intersection && (P1->neighbour == Q) )
		return(IS_P_m);

	// is Q linked to P2 ?
	if ( P3->intersection && (P3->neighbour == Q) )
		return(IS_P_p);

	// check relative position of Q with respect to chain (P1,P2,P3)
	double s1 = A( Q->p, P1->p, P2->p);
	double s2 = A( Q->p, P2->p, P3->p);
	double s3 = A(P1->p, P2->p, P3->p);

	if (s3 > 0) { 
		// chain makes a left turn
		if (s1 > 0 && s2 > 0)
			return(LEFT);
		else
			return(RIGHT);
	}
	else {
		// chain makes a right turn (or is straight)
		if (s1 < 0 && s2 < 0)
			return(RIGHT);
		else
			return(LEFT);
	}
}
//
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
//
// label all intersections
//
void labelIntersections() {
  if(PF)
    cout << "\nLabelling intersections...\n\n";
  
  high_resolution_clock::time_point start1, start2, start3, start4, start5, start6, start7, end1, end2, end3, end4, end5, end6, end7;
	//
	// 1) initial classification
  //

  int count[2] = {0,0};

  if(DEBUG2)
    start1 = high_resolution_clock::now();
  // loop over intersection vertices of P
  for (polygon& P : PP) 
    for (vertex* I : P.vertices(INTERSECTION)) {

  		// determine local configuration at this intersection vertex
  		vertex* P_m = I->prev;                // P-, predecessor of I on P
  		vertex* P_p = I->next;                // P+, successor of I on P
  		vertex* Q_m = I->neighbour->prev;     // Q-, predecessor of I on Q
  		vertex* Q_p = I->neighbour->next;     // Q+, successor of I on P
  
  		// check positions of Q- and Q+ relative to (P-, I, P+)
  		RelativePositionType Q_m_type = oracle(Q_m, P_m, I, P_p);
  		RelativePositionType Q_p_type = oracle(Q_p, P_m, I, P_p);
  
  		//
  		// check non-overlapping cases
  		//
  		if ((Q_m_type == LEFT  && Q_p_type == RIGHT) ||
  				(Q_m_type == RIGHT && Q_p_type == LEFT )) {
  			I->label = CROSSING;
        count[0]++;
      }
  
  		if ((Q_m_type == LEFT  && Q_p_type == LEFT ) ||
  				(Q_m_type == RIGHT && Q_p_type == RIGHT)) {
  			I->label = BOUNCING;
        count[1]++;
      }
  		//
  		// check overlapping cases
  		//
  		if ( ( (Q_p_type == IS_P_p) && (Q_m_type == RIGHT) ) ||
  		     ( (Q_m_type == IS_P_p) && (Q_p_type == RIGHT) ) )
  			I->label = LEFT_ON;
  
  		if ( ( (Q_p_type == IS_P_p) && (Q_m_type == LEFT) ) ||
  		     ( (Q_m_type == IS_P_p) && (Q_p_type == LEFT) ) )
  			I->label = RIGHT_ON;
  
  		if ( ( (Q_p_type == IS_P_p) && (Q_m_type == IS_P_m) ) ||
  		     ( (Q_m_type == IS_P_p) && (Q_p_type == IS_P_m) ) )
  			I->label = ON_ON;
  
  		if ( ( (Q_m_type == IS_P_m) && (Q_p_type == RIGHT) ) ||
  		     ( (Q_p_type == IS_P_m) && (Q_m_type == RIGHT) ) )
  			I->label = ON_LEFT;
  
  		if ( ( (Q_m_type == IS_P_m) && (Q_p_type == LEFT) ) ||
  		     ( (Q_p_type == IS_P_m) && (Q_m_type == LEFT) ) )
  			I->label = ON_RIGHT;
  	}
  if(DEBUG2)
    end1 = high_resolution_clock::now();
  if(PF)
    cout << "... " << count[0] << " crossing and " << count[1] << " bouncing intersection vertices" << endl;

	//
	// 2) classify intersection chains
	//

  count[0] = count[1] = 0;
  if(DEBUG2)
    start2 = high_resolution_clock::now();
  // loop over intersection vertices of P
  for (polygon& P : PP) 
    for (vertex* I : P.vertices(INTERSECTION)) {
      
      // start of an intersection chain ?
  		if (I->label == LEFT_ON ||
  				I->label == RIGHT_ON) {
  
  			// remember status of the first chain vertex and vertex itself
  			RelativePositionType x; 
  			if (I->label == LEFT_ON)
  				x = LEFT;
  			else
  				x = RIGHT;
        vertex* X = I;
  
  			// proceed to end of intersection chain and mark all visited vertices as NONE
  			do {
  				I->label = NONE;
  				I = I->next;
  			} while (I->label == ON_ON);
  
  			RelativePositionType y; 
  			if (I->label == ON_LEFT)
  				y = LEFT;
  			else
  				y = RIGHT;
  
        // determine type of intersection chain
        IntersectionLabel chainType;
        if (x != y) { 
        	chainType = DELAYED_CROSSING;
          count[0]++;
        } 
        else {
        	chainType = DELAYED_BOUNCING;
          count[1]++;
        }
        
        // mark both ends of an intersection chain with chainType (i.e., as DELAYED_*)
        X->label = chainType;
        I->label = chainType;
  		}
    }
  if(DEBUG2)
    end2 = high_resolution_clock::now();
  if(PF)
    cout << "... " << count[0] << " delayed crossings and " << count[1] << " delayed bouncings" << endl;

	//
	// 3) copy labels from P to Q
	//
  if(DEBUG2)
    start3 = high_resolution_clock::now();
  // loop over intersection vertices of P
  for (polygon& P : PP) 
    for (vertex* I : P.vertices(INTERSECTION))
      I->neighbour->label = I->label;
  if(DEBUG2)
    end3 = high_resolution_clock::now();
  //
  // 3.5) check for special cases
  //
  if(DEBUG2)
    start4 = high_resolution_clock::now();
  set<polygon*> noIntersection[2];
  set<polygon*> identical[2];

  count[0] = count[1] = 0;

  for (int i=0; i<2; ++i) {
    vector<polygon>* P_or_Q = &PP;      // if i=0, then do it for P w.r.t. Q
    vector<polygon>* Q_or_P = &QQ;

    if (i==1) {                         // if i=1, then do it for Q w.r.t. P
      P_or_Q = &QQ; 
      Q_or_P = &PP;
    }
  
    // loop over all components of P (or Q)
    for (polygon& P : *P_or_Q)
      if (P.noCrossingVertex(UNION)) {
        //
        // P_ has no crossing vertex (but may have bounces or delayed bounces, except for UNION),
        // hence it does not intersect with Q_or_P
        //
        noIntersection[i].insert(&P);   // remember component, and ignore it later in step 4
        
        // is P identical to some component of and Q_or_P?
        if (P.allOnOn()) {
          identical[i].insert(&P);      // -> remember for further processing below
        }
        else {  
          // is P inside Q_or_P?
          bool isInside = false;
          point2D p = P.getNonIntersectionPoint();
          for (polygon& Q : *Q_or_P)
            if ( Q.pointInPoly(p) )
              isInside = !isInside;              
          if (isInside ^ UNION) {
            RR.push_back(P);             // -> add P to the result
            count[0]++;
          }
        }
      }
  }

  // handle components of P that are identical to some component of Q
  for (polygon* P : identical[0]) {
    // is P a hole?
    bool P_isHole = false;  
    for (polygon& P_ : PP)
      if ( ( P_.root != P->root ) && (P_.pointInPoly(P->root->p)) )
        P_isHole = !P_isHole;

    for (polygon* Q : identical[1])
      for (vertex* V : Q->vertices(ALL))
        if (V == P->root->neighbour) {  // found Q that matches P
          // is Q a hole?
          bool Q_isHole = false;  
          for (polygon& Q_ : QQ)
            if ( ( Q_.root != Q->root ) && (Q_.pointInPoly(Q->root->p)) )
              Q_isHole = !Q_isHole;

          // if P and Q are both holes or both are not holes
          if (P_isHole == Q_isHole) {
            RR.push_back(*P);           // -> add P to the result
            count[1]++;
          }          
          goto next_P;
        }
    next_P: ;
  }
  if(DEBUG2)
    end4 = high_resolution_clock::now();
  if(PF)
   cout << "... " << count[0] << " interior and " << count[1] << " identical components added to result\n";

	//
	// 4) set entry/exit flags
	//
  if(DEBUG2)
    start5 = high_resolution_clock::now();
  set<vertex*> split[2];                // split vertex candidates for P and Q
  set<vertex*> crossing[2];             // CROSSING vertex candidates for P and Q

  for (int i=0; i<2; ++i) {
    vector<polygon>* P_or_Q = &PP;      // if i=0, then do it for P w.r.t. Q
    vector<polygon>* Q_or_P = &QQ;

    if (i==1) {                         // if i=1, then do it for Q w.r.t. P
      P_or_Q = &QQ; 
      Q_or_P = &PP;
    }
  
    // loop over all components of P (or Q)
    for (polygon& P : *P_or_Q) {

      // ignore P if it does not intersect with Q_or_P (detected in step 3.5 above)
      if(noIntersection[i].find(&P) != noIntersection[i].end())
        continue;
        
      // start at a non-intersection vertex of P
      vertex* V = P.getNonIntersectionVertex();
  
      // check if it is inside or outside Q (or P)
      // and set ENTRY/EXIT status accordingly          
      EntryExitLabel status = ENTRY;
      for (polygon& Q : *Q_or_P)
        if (Q.pointInPoly(V->p))
          toggle(status);

      //
      // starting at V, loop over those vertices of P, that are either
      // a crossing intersection or marked as ends of an intersection chain
      //
      bool first_chain_vertex = true;     // needed for dealing with crossing chains
      
      for (vertex* I : P.vertices(INTERSECTION, V)) {
        //
        // in the case of normal crossings, we...
        //
        if (I->label == CROSSING) {
          // mark vertex with current ENTRY/EXIT status
    			I->enex = status;
          // toggle status from ENTRY to EXIT or vice versa
    			toggle(status);
        }
        
        //
        // identify split vertex candidates (INTERIOR bouncing vertices)
        //
        if ( (I->label == BOUNCING) && ((status == EXIT) ^ UNION) )
          split[i].insert(I);

        //
        // in the case of a delayed crossing chain, we
        // mark both end points of the chain with the current ENTRY/EXIT status,
        // toggling the status only at the end last chain vertex,
        // and, in case of a delayed EXIT  crossing, the first vertex
        //  or, in case of a delayed ENTRY crossing, the last  vertex,
        // of the chain as CROSSING
        //
        if (I->label == DELAYED_CROSSING) {
          // mark vertex with current ENTRY/EXIT status
          I->enex = status;
  
          if (first_chain_vertex) {       // are we at the first vertex of a delayed crossing chain?
            if ((status == EXIT) ^ UNION)
              I->label = CROSSING;        // mark first vertex as CROSSING
            first_chain_vertex = false;
          } 
          else {                          // here we are at the last vertex of a delayed crossing chain
            if ((status == ENTRY) ^ UNION)
              I->label = CROSSING;        // mark last vertex as CROSSING
            first_chain_vertex = true;
            
            // toggle status from ENTRY to EXIT or vice versa (only for last chain vertex)
            toggle(status);
          }
        }
  
        //
        // in the case of a delayed bouncing chain, we
        // mark both end points of the chain with the current ENTRY/EXIT status
        // toggling the status at both end points of the chain,
        // and, in case of a delayed INTERIOR bouncing, both end points
        // of the chain as CROSSING candidates
        //
        if (I->label == DELAYED_BOUNCING) {
          // mark vertex with current ENTRY/EXIT status
          I->enex = status;
  
          if (first_chain_vertex) {       // are we at the first vertex of a delayed crossing chain?
            if ((status == EXIT) ^ UNION)
              crossing[i].insert(I);      // mark first EXIT vertex as CROSSING candidate
            first_chain_vertex = false;
          } 
          else {                          // here we are at the last vertex of a delayed crossing chain
            if ((status == ENTRY) ^ UNION)
              crossing[i].insert(I);      // mark last ENTRY vertex as CROSSING candidate
            first_chain_vertex = true;
            
          }
          // toggle status from ENTRY to EXIT or vice versa (for first AND last chain vertex)
          toggle(status);
        }
  		}
    }  
  }
  if(DEBUG2)
    end5 = high_resolution_clock::now();
	//
	// 5) handle split vertex pairs
	//
  if(DEBUG2)
    start6 = high_resolution_clock::now();
  count[0] = 0;
  
  // loop over P's split candidates
  for (vertex* I_P : split[0]) {
    vertex* I_Q = I_P->neighbour;
    
    // check if the neighbour on Q is also a split candidate
    if (split[1].find(I_Q) != split[1].end()) {
      count[0]++;
      //
      // split vertex pair
      //
      
      // duplicate vertices
			vertex* V_P = new vertex(I_P->p);
			vertex* V_Q = new vertex(I_Q->p);
      
      // compute areas to compare local orientation
      double sP = A( I_P->prev->p, I_P->p, I_P->next->p);
      double sQ = A( I_Q->prev->p, I_Q->p, I_Q->next->p);

      // link vertices correctly
      if (sP*sQ > 0) {                  // same local orientation
        link(I_P, V_Q);
        link(I_Q, V_P);
      }
      else {                            // different local orientation
        link(V_P, V_Q);
      }

      // add duplicate vertices to P and Q    
      insertVertex(V_P, I_P);
      insertVertex(V_Q, I_Q);

      // mark all four vertices correctly
      if (!UNION) {
        I_P->enex = EXIT;
        V_P->enex = ENTRY;
        I_Q->enex = EXIT;
        V_Q->enex = ENTRY;
      }
      else {
        I_P->enex = ENTRY;
        V_P->enex = EXIT;
        I_Q->enex = ENTRY;
        V_Q->enex = EXIT;
      }

      I_P->label = CROSSING;
      V_P->label = CROSSING;
      I_Q->label = CROSSING;
      V_Q->label = CROSSING;
    }
  }
  if(DEBUG2)
    end6 = high_resolution_clock::now();
  if(PF)
    cout << "... " << count[0] << " bouncing vertex pairs split" << endl;
  
	//
	// 6) handle CROSSING vertex candidates
	//
  if(DEBUG2)
    start7 = high_resolution_clock::now();
  // loop over P's CROSSING candidates
  for (vertex* I_P : crossing[0]) {
    vertex* I_Q = I_P->neighbour;
    
    // check if the neighbour on Q is also a CROSSING candidate
    if (crossing[1].find(I_Q) != crossing[1].end()) {
      //
      // mark CROSSING candidate pair as such
      //
      I_P->label = CROSSING;
      I_Q->label = CROSSING;
    }
  }
  if(DEBUG2){
    end7 = high_resolution_clock::now();

    auto duration1 = duration_cast<microseconds>(end1 - start1);
    auto duration2 = duration_cast<microseconds>(end2 - start2);
    auto duration3 = duration_cast<microseconds>(end3 - start3);
    auto duration4 = duration_cast<microseconds>(end4 - start4);
    auto duration5 = duration_cast<microseconds>(end5 - start5);
    auto duration6 = duration_cast<microseconds>(end6 - start6);
    auto duration7 = duration_cast<microseconds>(end7 - start7);
    cout << " Time: initial classification: " << fixed << duration1.count() << setprecision(10) << " microseconds" << endl;
    cout << " Time: classify intersection chains: " << fixed << duration2.count() << setprecision(10) << " microseconds" << endl;
    cout << " Time: copy labels from P to Q: " << fixed << duration3.count() << setprecision(10) << " microseconds" << endl;
    cout << " Time: check for special cases: " << fixed << duration4.count() << setprecision(10) << " microseconds" << endl;
    cout << " Time: set entry/exit flags: " << fixed << duration5.count() << setprecision(10) << " microseconds" << endl;
    cout << " Time: handle split vertex pairs: " << fixed << duration6.count() << setprecision(10) << " microseconds" << endl;
    cout << " Time: handle CROSSING vertex candidates: " << fixed << duration7.count() << setprecision(10) << " microseconds" << endl;

  }
}
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// create the result polygon
//
void createResult() {
  if(PF)
    cout << "\nCreating result...\n";
  //
  // for all crossing vertices
  //
  // NOTE: all crossing vertices that are visited while contructing a
  //       component of the result polygon are marked as "not intersection",
  //       so that they cannot serve as start vertex of another component
  // 

  for (polygon& P : PP) 
    for (vertex* I : P.vertices(CROSSING_INTERSECTION)) {
      polygon R;                         // result polygon component

      vertex* V = I;                      // start traversal at I
      V->intersection = false;            // mark visited vertices
  
      do {
        EntryExitLabel status = V->enex;
        toggle(status);
    		do {                              // traverse P (or Q) and add vertices to R, until...
          if ((status == EXIT) ^ UNION)
            V = V->next;                  // move forward  from an ENTRY vertex to the next EXIT  vertex
          else
            V = V->prev;                  // move backward from an EXIT  vertex to the next ENTRY vertex
          V->intersection = false;        // mark visited vertices
  
    			// add vertex to result polygon 
          R.newVertex(V->p);
        } while ( !(V->enex == status)    // ... we arrive at a vertex with opposite entry/exit flag, or
                  && (V != I) );          // at the initial vertex I
  
        if (V != I) {
          V = V->neighbour;               // switch from P to Q or vice versa
          V->intersection = false;        // mark visited vertices
        }
      } while (V != I);                   // the result polygon component is complete, 
                                          // if we are back to the initial vertex I
      RR.push_back(R);
    }
}
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// cleans up the result polygon by removing unnecessary collinear vertices
//
void cleanUpResult() {
  if(PF)
   cout << "\nPost-processing...\n\n";
  int count = 0;
  for (polygon& R : RR) {
    while ( (R.root != NULL) && (fabs(A(R.root->prev->p,R.root->p,R.root->next->p)) < EPSILON) ) {
      R.removeVertex(R.root);
      count++;
    }
    if (R.root != NULL)
      for (vertex* V : R.vertices(ALL))
        if (fabs(A(V->prev->p,V->p,V->next->p)) < EPSILON) {
          R.removeVertex(V);
          count++;
        }
  }
  if(PF)
    cout << "... " << count << " vertices removed\n\n";
}
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// prints statistics (#components, #vertices) for a complex polygon
//
void printInfo(vector<polygon>& PP) {
  int sum = 0;
  cout << "has "<< PP.size() << " component";
  if (PP.size() > 1) 
    cout << "s";
  cout << " with ";
  for (int i=0; i<PP.size(); i++) {
    int count = 0;
    for (vertex* V : PP[i].vertices(ALL))
      count++;
    cout << count;
    if (i<PP.size()-1) 
      cout << " + ";
    sum += count;
  }
  if (PP.size() > 1)
    cout << " = " << sum;
  cout << " vertices\n\n";
}
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// load complex polygon P from file "s" and print statistics
//
void loadPolygon(vector<polygon>& PP, string s) {
  ifstream from(s);
  do {
    polygon P;
    from >> P;
    if (!from.fail())
      PP.push_back(P);
  } while (!from.eof());
  from.close();
  if(PF)
    printInfo(PP);
}
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// save polygons to array
//
void gpc_read_polygon(FILE *fp, string polyName){
  int v, size; 
  point2D point;
  polygon P;
  double x, y;
  fscanf(fp, "%d", &size);
          
  for (v=0; v<size; v++){
    fscanf(fp, "%lf,%lf", &x, &y);
    point.x=x;
    point.y=y;
    P.newVertex(point, true);
    // printf(" %f %f, ", p->contour[0].vertex[v].x,
    //                    p->contour[0].vertex[v].y );
  }
  // P.numVertices=*size;
  if(polyName=="PP") PP.push_back(P);
  else QQ.push_back(P);
  if(PF) printf("Polygon with %d vertices reading completed!\n", size);
}

////////////////////////////////////////////////////////////////////////
//
// save complex polygon P to file "s" and print statistics
//
void savePolygon(vector<polygon>& PP, string s) {
  ofstream to(s);
  for (polygon& P : PP)
    to << P;
  to.close();
  if(PF)
    printInfo(PP);
}
//
////////////////////////////////////////////////////////////////////////


int main(int argc, char* argv[])
{
  // /*
  // check input parameters
  if (argc < 4) {
    if(PF)
      cout << "insufficient number of parameters. exe-file input1 path input2 path result path" << endl;
    exit(0);
  }
  
  int argn = 1;
  if (string(argv[1]) == "-union") {
    if(PF)
      cout << "\n!!! computing UNION instead of INTERSECTION !!!\n";
    UNION = true;
    argn++;
  }
// */
  high_resolution_clock::time_point start, end, start1, end1, start2, end2, start3, end3;

  // -------------------------------------------------------
  // read input polygons
  // -------------------------------------------------------
  /*
  if(PF){ 
    cout << "\nP "; loadPolygon(PP,string(argv[argn++]));
    cout << "\nQ "; loadPolygon(QQ,string(argv[argn++]));
  } else {
    loadPolygon(PP,string(argv[argn++]));
    loadPolygon(QQ,string(argv[argn++]));
  }
  */
  // -------------------------------------------------------

  // -------------------------------------------------------
  // use with large polygons with size, point configuration data files
  // -------------------------------------------------------
  FILE *pfile, *qfile;
  pfile=fopen(argv[argn++], "r");
  qfile=fopen(argv[argn++], "r");
  gpc_read_polygon(pfile, "PP");
  gpc_read_polygon(qfile, "QQ");

  if(argc > 4){
    if(string(argv[argn+1]) == "save") PR=true;
    if(string(argv[argn+1]) == "debug") {
      PF=true;
      DEBUG2=true;
    }
  }
  if(argc > 5)
    if(string(argv[argn+2]) == "debug") {
      PF=true;
      DEBUG2=true;
    }
  // -------------------------------------------------------


// original 
  // cout << "\nP "; loadPolygon(PP,string(argv[argn++]));
  // cout <<   "Q "; loadPolygon(QQ,string(argv[argn++]));

	// int PPID=9560, QQID=9560;
  // loadPolygonFromShapeFile(PPTmp, string("/mnt/d//Datasets/sports/sports"), PPID+1);
	// loadPolygonFromShapeFile(QQTmp, string("/mnt/d//Datasets/sports/sports"), QQID+1);
/*
	int PPID=36; //ne_10m_ocean
	// int PPID=0; //ne_10m_ocean
	// int PPID=2742; //ne_10m_ocean
	// int PPID=2741; //ne_10m_ocean
  loadPolygonFromShapeFile2(PPTmp, string("../../datasets/ne_10m_ocean.csv"), PPID+1);
  // int QQID=521; //continents
  // int QQID=1661; //continents
  // int QQID=1048; //continents
  // int QQID=1081; //continents
  // int QQID=1193; //continents
	// loadPolygonFromShapeFile2(QQTmp, string("../../datasets/continents.csv"), QQID+1);
  int QQID=1; //ne_10m_land
  // int QQID=4; //ne_10m_land
	loadPolygonFromShapeFile2(QQTmp, string("../../datasets/ne_10m_land.csv"), QQID+1);
  
  if(PF){
    cout << "PP Polygon size " << PPTmp[PPID].size;
    cout << " QQ Polygon size " << QQTmp[QQID].size << endl;
  }

  PP.push_back(PPTmp[PPID]);
  QQ.push_back(QQTmp[QQID]);
  // QQ.push_back(PPTmp[PPID]);
*/
  start = high_resolution_clock::now();

  // CMF();

  if(DEBUG)
    start1 = high_resolution_clock::now();
  // phase 1
  computeIntersections2();
  if(DEBUG)
    end1 = high_resolution_clock::now(); 

  if(DEBUG)
    start2 = high_resolution_clock::now();
  // phase 2
  labelIntersections();
  if(DEBUG)
    end2 = high_resolution_clock::now(); 

  if(DEBUG)
    start3 = high_resolution_clock::now();
  // phase 3
  createResult();
  if(DEBUG)
    end3 = high_resolution_clock::now(); 
 
  end = high_resolution_clock::now();
  // post-processing
  cleanUpResult();
  
  // write output polygon
  if(PR){
    // cout << "R "; 
    // savePolygon(RR,string(argv[argn]));
    cout << "R "; 
    savePolygon(RR,string(argv[argn]));
    // savePolygon(RR,string("results/ocean1.poly"));
  }else{
    // savePolygon(RR,string(argv[argn]));
    // savePolygon(RR,string("results/ocean1.poly"));
  }


  if(DEBUG && PF){
    // Calculating total time taken by the program.
    auto duration = duration_cast<microseconds>(end - start);
    auto duration1 = duration_cast<microseconds>(end1 - start1);
    auto duration2 = duration_cast<microseconds>(end2 - start2);
    auto duration3 = duration_cast<microseconds>(end3 - start3);

    cout << "All time in microseconds\nTime: Total : " << fixed
    << duration.count() << setprecision(10) << endl;
    cout << "Time: computeIntersections(): " << fixed
    << duration1.count() << setprecision(10) << endl;
    cout << "Time: labelIntersections(): " << fixed
    << duration2.count() << setprecision(10) << endl;
    cout << "Time: createResult(): " << fixed
    << duration3.count() << setprecision(10) << endl;
    cout << "\nintersection calc exec %: " << fixed << (duration1.count()*100.0/duration.count()) << setprecision(10) << endl;
    cout << "labeling calc exec %: " << (duration2.count()*100.0/duration.count()) << endl;
    cout << "result create calc exec %: " << (duration3.count()*100.0/duration.count()) << endl;    
  }
  else{
    // Calculating total time taken by the program.
    auto duration = duration_cast<microseconds>(end - start);
    auto duration1 = duration_cast<microseconds>(end1 - start1);
    auto duration2 = duration_cast<microseconds>(end2 - start2);
    auto duration3 = duration_cast<microseconds>(end3 - start3);

    cout << setprecision(10) << duration1.count() << ", " <<
     duration2.count() << ", " <<
     duration3.count() << ", " << duration.count();
  }
}
