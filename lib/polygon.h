////////////////////////////////////////////////////////////////////////
//
// polygon.h
//
// contains everything that the polyclip algorithm needs for
// handling polygons
//
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>

using namespace std;

////////////////////////////////////////////////////////////////////////
//
// a simple class to handle polygon vertices
//
////////////////////////////////////////////////////////////////////////

enum IntersectionLabel {      // for the classification of intersection vertices in the second phase
  NONE,
  CROSSING,
  BOUNCING,
  LEFT_ON,
  RIGHT_ON,
  ON_ON,
  ON_LEFT,
  ON_RIGHT,
  DELAYED_CROSSING,
  DELAYED_BOUNCING
};

enum EntryExitLabel {         // for marking intersection vertices as "entry" or "exit"
  EXIT,
  ENTRY,
  NEITHER
};

class vertex
{
public:
  point2D p;                  // coordinates of the vertex
  vertex* prev;               // pointer to previous vertex
  vertex* next;               // pointer to next vertex
  vertex* neighbour;          // pointer to neighbouring vertex for intersection vertices
  bool source;                // to mark source vertices of the polygon
  bool intersection;          // to mark intersection vertices
  double alpha;               // to describe relative edge position of an intersection vertex
  IntersectionLabel label;    // type of intersection vertex
  EntryExitLabel enex;        // entry/exit "flag"

  //
  // default-constructor for generating a vertex with coordinates (0,0)
  //
	vertex() : 
    p(), prev(NULL), next(NULL), neighbour(NULL),
    intersection(false), source(false),
    alpha(-1.0), label(NONE), enex(NEITHER)
  {};

  //
  // constructor for generating a vertex with explicit coordinates
  //
  vertex(double x, double y) :
    p(x,y), prev(NULL), next(NULL), neighbour(NULL),
    intersection(false), source(false),
    alpha(-1.0), label(NONE), enex(NEITHER)
  {};

  //
  // constructor for generating a vertex from a 2D point;
  // if alpha is given (and then in [0,1], thus > -1.0),
  // then the vertex will be marked as an intersection vertex
  //
  vertex(point2D& q, double alpha = -1.0) : 
    p(q), prev(NULL), next(NULL), neighbour(NULL),
    source(false),
    alpha(alpha), label(NONE), enex(NEITHER)
  {
    intersection = (alpha > -1.0);
  };
};

//
// link vertex P to vertex Q, using the neighbour pointer
// and mark both vertices as intersection vertices
//
void link(vertex* P, vertex* Q) {
  P->neighbour = Q;
  Q->neighbour = P;
  P->intersection = true;
  Q->intersection = true;
}

//
// insert a vertex between two existing vertices
// if alpha is given (and then in [0,1], thus > -1.0),
// then the vertex will be inserted, based on alpha
//
void insertVertex (vertex* V, vertex* curr, double alpha = -1.0) {
  if (alpha > -1.0)
		do {
			curr = curr->next;
		} while (!curr->source && curr->alpha < alpha);
  else
    curr = curr->next;
    
	curr->prev->next = V;
	V->prev = curr->prev;
	V->next = curr;
	curr->prev = V;
}

////////////////////////////////////////////////////////////////////////
//
// a simple class to handle polygon edges
//
////////////////////////////////////////////////////////////////////////

class edge 
{
public:
  vertex* one;
  vertex* two;

  //
  // constructor for generating an edge from vertex P to vertex Q
  //
  edge (vertex* P, vertex* Q) : one(P), two(Q) { };
};

//
// insert an intersection vertex between the endpoints of edge E
//
void insertVertex (vertex* V, edge& E) {
  insertVertex(V, E.one, V->alpha);
}

////////////////////////////////////////////////////////////////////////
//
// iterators to loop over vertices and edges, with the option to
// restrict to source vertices or edges and to intersection vertices
//
////////////////////////////////////////////////////////////////////////

enum IteratorType {
  SOURCE,
  INTERSECTION,
  CROSSING_INTERSECTION,
  ALL
};

class vertexIterator 
{ 
private:
  class iterator {
  public:
    iterator(vertex* root, IteratorType IterType) :
      root(root), V(NULL), iterType(IterType) 
    { 
      if (root == NULL)
        return;
  
      if (nextVertex() == NULL)         // no (source/intersection) vertex found
        root = V = NULL;                // -> mark iterator as "end"
    }
  
    const iterator& operator++() { 
      nextVertex(); 
      return *this; 
    }

    vertex* operator*() {
      return V;
    }
  
    bool operator!=(const iterator& other) const {
      return (root != other.root) || (V != other.V);
    }

  private:
    vertex* root;
    vertex* V;
    IteratorType iterType;

    //
    // find the next vertex
    // if iterType is ALL, then it is just the next vertex
    // if iterType is SOURCE, then it is the next source vertex
    // if iterType is INTERSECTION, then it is the next intersection vertex
    // if iterType is CROSSING_INTERSECTION, then it is the next intersection vertex with CROSSING label
    //
    vertex* nextVertex() {
      bool nextFound = false;
  
      if (V == NULL) {                  // find first (source/intersection) vertex
        V = root;
        switch(iterType) {
        case ALL:
          nextFound = true;
          break;
        case SOURCE:
          if (V->source)
            nextFound = true;
          break;
        case INTERSECTION:
          if (V->intersection)
            nextFound = true;
          break;
        case CROSSING_INTERSECTION:
          if (V->intersection && (V->label == CROSSING))
            nextFound = true;
          break;
        }
      }
          
      while (!nextFound) {              // find next (source/intersection) vertex
        switch(iterType) {
        case ALL:
          V = V->next;
          break;
        case SOURCE:
          do {
            V = V->next;
          } while (!V->source && V != root);
          break;
        case INTERSECTION:
          do {
            V = V->next;
          } while (!V->intersection && V != root);
          break;
        case CROSSING_INTERSECTION:
          do {
            V = V->next;
          } while ( ( !V->intersection || (V->label != CROSSING) ) && V != root);
          break;
        }
  
        if (V == root) {                // back at the root vertex?
          root = V = NULL;              // -> mark iterator as "end"
          return(V);
        }
        
        switch(iterType) {
        case ALL:
          nextFound = true;
          break;
        case SOURCE:
          if (V->source)
            nextFound = true;
          break;
        case INTERSECTION:
          if (V->intersection)
            nextFound = true;
          break;
        case CROSSING_INTERSECTION:
          if (V->intersection && (V->label == CROSSING))
            nextFound = true;
          break;
        }
      }      
      return(V);
    }
  };

public:
  vertexIterator() : root(NULL) {};
  
  iterator begin() { return iterator(root, iterType); }
  iterator end()   { return iterator(NULL, iterType); }

  vertex* root;
  IteratorType iterType;
};


class edgeIterator 
{ 
private:
  class iterator {
  public:
    iterator(vertex* root, IteratorType IterType) :
      root(root), one(NULL), two(NULL), iterType(IterType) 
    { 
      if (root == NULL)
        return;
  
      if (nextEdge() == NULL)           // no source edge found
        root = one = two = NULL;        // -> mark iterator as "end"
    }
  
    const iterator& operator++() { nextEdge(); return *this; }

    edge operator*() {
      return edge(one,two);
    }
  
    bool operator!=(const iterator& other) const {
      return (root != other.root) || (one != other.one) || (two != other.two);
    }

  private:
    vertex* root;
    vertex* one;
    vertex* two;
    IteratorType iterType;

    //
    // find the next vertex, starting at curr
    // if iterType is ALL, then it is just the next vertex
    // if iterType is SOURCE, then it is the next source vertex
    //
    vertex* nextVertex(vertex* curr) {
      if (curr == NULL)
        return(NULL);
        
      switch(iterType) {
      case ALL:
        curr = curr->next;
        break;
        
      case SOURCE:
        do {
          curr = curr->next;
        } while (!curr->source);
        break;
      }
      
      return(curr);
    }
    
    //
    // find the next edge
    //
    vertex* nextEdge() {
      if (root == NULL)                 // empty polygon?
        return (NULL);
     
      if (one == NULL) {                // find one (source) vertex
        one = root;                     // note: root is always a (source) vertex
        two = nextVertex(one);
        if (two == one)                 // just one (source) vertex
          return(NULL);                 // -> no (source) edges
        return(one);
      }
      
      if (two == root) {                // back at the root vertex?
        root = one = two = NULL;        // -> mark iterator as "end"
        return(NULL);
      }
  
      one = two;
      two = nextVertex(one);
  
      return (one);
    }    
  };

public:
  edgeIterator() : root(root) {};
  
  iterator begin() { return iterator(root, iterType); }
  iterator end()   { return iterator(NULL, iterType); }

  vertex* root;
  IteratorType iterType;
};


////////////////////////////////////////////////////////////////////////
//
// a simple polygon class with essential functionality
//
////////////////////////////////////////////////////////////////////////

class polygon
{
public:
  //
  // default-constructor for generating an empty polygon
  //
  polygon() : root(NULL) {}
  int size=0;
  //
  // create a new vertex and append it to the polygon
  //
  void newVertex(point2D& v, bool source = false) {
    vertex* V = new vertex(v);
    V->source = source;
    
    if (root == NULL) {
      // add vertex as the very first vertex of the polygon
      V->next = V;
      V->prev = V;
      root = V;
      size++;
    } 
    else {
      // add vertex at the end 
      V->prev = root->prev;
      V->next = root;
      root->prev->next = V;
      root->prev = V;
      size++;
    } 
  }

  //
  // removes the vertex V from the polygon
  //
  void removeVertex (vertex* V) {
    if (root == V) {
      root = V->next;
      if (root->next == root)
        root = NULL;
    }
    V->prev->next = V->next;
    V->next->prev = V->prev;
    delete(V);
  }
  
  //
  // check, if the point R lies inside the polygon or not
  // cf. Algorithm 6 in Hormann & Agathos [2001]
  //
  bool pointInPoly(point2D& R) {
  	int w = 0;
    for (edge E : edges(ALL)) {
			point2D& P0 = E.one->p;
			point2D& P1 = E.two->p;

			if ( (P0.y < R.y) != (P1.y < R.y) )
				if (P0.x >= R.x) {
					if (P1.x > R.x)
						w = w + 2 * (P1.y > P0.y) - 1;
					else
						if ( (A(P0,P1,R) > 0) == (P1.y > P0.y) )
							w = w + 2 * (P1.y > P0.y) - 1;
				}
				else
					if (P1.x > R.x)
						if ( (A(P0,P1,R) > 0) == (P1.y > P0.y) )
							w = w + 2 * (P1.y > P0.y) - 1;
		}
  
  	return ( (w % 2) != 0 );
  }

  //
  // check, if all vertices have the ON_ON label
  //
  bool allOnOn() {
    for (vertex* V : vertices(ALL))
      if (V->label != ON_ON)
        return(false);
    return(true);
  }

  //
  // check, if the polygon does not contain any crossing intersection vertex
  // or crossing intersection chain or (if we want to compute the union instead
  // of the intersection) a bouncing vertex or a bouncing intersection chain
  //
  bool noCrossingVertex(bool union_case = false) {
    for (vertex* V : vertices(ALL)) 
      if (V->intersection) {
        if ( (V->label == CROSSING) || (V->label == DELAYED_CROSSING) )
          return(false);

        if (union_case && ( (V->label == BOUNCING) || (V->label == DELAYED_BOUNCING) ) )
          return(false);
      }
    return(true);
  }

  //
  // return a non-intersection point
  //
  point2D getNonIntersectionPoint() {
    for (vertex* V : vertices(ALL))
      if (!V->intersection)
        return(V->p);

    // no non-intersection vertex found -> find suitable edge midpoint
    for (vertex* V : vertices(ALL))
      // make sure that edge from V to V->next is not collinear with other polygon
      if ( (V->next->neighbour != V->neighbour->prev) && (V->next->neighbour != V->neighbour->next) )
        // return edge midpoint
        return(0.5*(V->p + V->next->p));
  }
  
  //
  // return and insert a non-intersection vertex
  //
  vertex* getNonIntersectionVertex() {
    for (vertex* V : vertices(ALL))
      if (!V->intersection)
        return(V);

    // no non-intersection vertex found -> generate and return temporary vertex
    for (vertex* V : vertices(ALL))
      // make sure that edge from V to V->next is not collinear with other polygon
      if ( (V->next->neighbour != V->neighbour->prev) && (V->next->neighbour != V->neighbour->next) ) {
        // add edge midpoint as temporary vertex
        point2D p = 0.5*(V->p + V->next->p);
        vertex* T = new vertex(p);
        insertVertex(T, V);
        return(T);
      }
    return(NULL);
  }

  //
  // return iterator to loop over (source/intersection) vertices,
  // starting at first (by default, it starts at root)
  //
  vertexIterator& vertices(IteratorType iterType, vertex* first = NULL) { 
    vertexIter.iterType = iterType;

    if (first == NULL)
      vertexIter.root = root;
    else
      vertexIter.root = first;      

    return vertexIter; 
  }
    
  //
  // return iterator to loop over (source) edges
  //
  edgeIterator& edges(IteratorType iterType) { 
    edgeIter.iterType = iterType;
    edgeIter.root = root;
    return edgeIter; 
  }

  vertex* root;               // root vertex of the polygon
  // int numVertices;            // number of vertices in the polygon

protected:
  vertexIterator vertexIter;  // vertex iterator
  edgeIterator edgeIter;      // edge iterator
};

//
// read polygon from a stream
//
inline istream& operator>>(istream& s, polygon& P) { 
  point2D v;
  char c;
  bool readOn = true;
  int count=0;
  
  do {
    s >> v >> c;
  	if (s.fail())
  		readOn = false;
  	else {
      count++;
      P.newVertex(v, true);

      if (c == ';')
        readOn = false;
  	}
  } while (readOn && !s.eof());
  // P.numVertices=count;

  return (s); 
}

//
// write polygon to a stream
//
inline ostream& operator<<(ostream& s, polygon& P) { 
  if (P.root == NULL)
    return (s);
  for (vertex* V : P.vertices(ALL)) {
    s << V->p;
    if (V != P.root->prev)
      s << "," << endl;
  }
  return (s << ";\n"); 
}

//
// toggle status from ENTRY to EXIT or vice versa
//
void toggle(EntryExitLabel& status) {
  if (status == ENTRY) {
    status = EXIT;
    return;
  }
  if (status == EXIT) {
    status = ENTRY;
    return;
  }
}
