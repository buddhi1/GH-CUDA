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
// g++ -std=c++11 polyclip.cpp -o polyclip.exe
//
// usage:
//
// polyclip [-union] input1.poly input2.poly output.poly
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <set>

using namespace std;

#include "constants.h"
#include "point2D.h"
#include "polygon.h"
#include "readSCPolygon.h"

////////////////////////////////////////////////////////////////////////

#define EPSILON  0.000000001            // tolerance

vector<polygon> PP, QQ; 
vertex **PPVertexPointers, **QQVertexPointers;                // two input polygons
vector<polygon> RR;                     // output polygon

bool UNION = false;                     // global switch for computing union instead of intersection

////////////////////////////////////////////////////////////////////////
//
// compute intersection of [P1,P2] and [Q1,Q2]
//
enum IntersectionType {                 // types of intersection (detected in the first phase)
  NO_INTERSECTION, //0
  X_INTERSECTION,  //1
  T_INTERSECTION_Q, //2
  T_INTERSECTION_P, //3
  V_INTERSECTION, //4
  X_OVERLAP,      //5
  T_OVERLAP_Q,    //6
  T_OVERLAP_P,    //7
  V_OVERLAP       //8
};

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
void computeIntersections() {
  if(DEBUG_INFO_PRINT) cout << "Computing intersections...\n\n";
  
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
          //
    			// determine intersection or overlap type
          //
          IntersectionType i = intersect(edgeP, edgeQ, alpha, beta);
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
  if(DEBUG_INFO_PRINT){  
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
  if(DEBUG_INFO_PRINT) cout << "\nLabelling intersections...\n\n";
  
	//
	// 1) initial classification
  //

  int count[2] = {0,0};
  /*
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

  cout << "... " << count[0] << " crossing and " << count[1] << " bouncing intersection vertices" << endl;
  */
	//
	// 2) classify intersection chains
	//

  count[0] = count[1] = 0;

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
  				I->neighbour->label = NONE; //////////***********new
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
        X->neighbour->label = chainType;  //*********new 
        I->neighbour->label = chainType;  //*********new
        // cout << "lbl2 "<<I->p.x<<","<<I->p.y<<"-"<<chainType<<" "<<X->p.x<<","<<X->p.y<<"-"<<chainType<<endl;
  		}
    }

  if(DEBUG_INFO_PRINT) cout << "... " << count[0] << " delayed crossings and " << count[1] << " delayed bouncings" << endl;

  // step 3 is commented and done in new steps in stesp 2 nested loop and GPU
	//*********************************************************
	// 3) copy labels from P to Q
	//
  // if(DEBUG_INFO_PRINT) cout <<"Label step 3 start"<<endl;

  // loop over intersection vertices of P
  // for (polygon& P : PP) 
  //   for (vertex* I : P.vertices(INTERSECTION)){
  //     I->neighbour->label = I->label;
  //   }

  // if(DEBUG_INFO_PRINT) cout <<"Label step 3 completed"<<endl;
  //*******************************************************************

  // 
  // 3.5) check for special cases
  //
  
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
  
  if(DEBUG_INFO_PRINT) cout << "... " << count[0] << " interior and " << count[1] << " identical components added to result\n";

	//
	// 4) set entry/exit flags
	//

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

	//
	// 5) handle split vertex pairs
	//

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

  if(DEBUG_INFO_PRINT) cout << "... " << count[0] << " bouncing vertex pairs split" << endl;
  
	//
	// 6) handle CROSSING vertex candidates
	//

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
}
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// create the result polygon
//
void createResult() {
  if(DEBUG_INFO_PRINT) cout << "\nCreating result...\n";
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
  if(DEBUG_INFO_PRINT) cout << "\nPost-processing...\n\n";
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
  if(DEBUG_INFO_PRINT) cout << "... " << count << " vertices removed\n\n";
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
  if(DEBUG_INFO_PRINT) printInfo(PP);
}
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// save complex polygon P to file "s" and print statistics
//
void savePolygon(vector<polygon>& PP, string s) {
  ofstream to(s);
  for (polygon& P : PP)
    to << P;
  to.close();
  if(DEBUG_INFO_PRINT) printInfo(PP);
}
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// save polygons to array
//
void gpc_read_polygon(FILE *fp, double **px, double **py, int *size, string polyName){
  int v; 
  point2D point;
  polygon P;
  fscanf(fp, "%d", size);

  MALLOC(*px, *size*sizeof(double), "vertex x values", double);
  MALLOC(*py, *size*sizeof(double), "vertex y values", double);
          
  for (v=0; v<*size; v++){
    fscanf(fp, "%lf,%lf", (*px+v), (*py+v));
    point.x=*(*px+v);
    point.y=*(*py+v);
    P.newVertex(point, true);
    // printf(" %f %f, ", p->contour[0].vertex[v].x,
    //                    p->contour[0].vertex[v].y );
  }
  // P.numVertices=*size;
  if(polyName=="PP") PP.push_back(P);
  else QQ.push_back(P);
  if(DEBUG_INFO_PRINT) printf("Polygon %s with %d vertices reading completed!\n", polyName.c_str(), *size);
}

//
////////////////////////////////////////////////////////////////////////

// int main(int argc, char* argv[])
// {
//   // check input parameters
//   if (argc < 4) {
//     cout << "insufficient number of parameters" << endl;
//     exit(0);
//   }
  
//   int argn = 1;
//   if (string(argv[1]) == "-union") {
//     cout << "\n!!! computing UNION instead of INTERSECTION !!!\n";
//     UNION = true;
//     argn++;
//   }

//   // read input polygons
//   cout << "\nP "; loadPolygon(PP,string(argv[argn++]));
//   cout <<   "Q "; loadPolygon(QQ,string(argv[argn++]));
  
//   // phase 1
//   computeIntersections();  

//   // phase 2
//   labelIntersections();

//   // phase 3
//   createResult();
  
//   // post-processing
//   cleanUpResult();
  
//   // write output polygon
//   cout << "R "; savePolygon(RR,string(argv[argn]));
// }
