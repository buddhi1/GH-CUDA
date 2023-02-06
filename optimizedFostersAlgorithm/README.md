This folder contains the implementation of [1]. We downloaded this from the author publication and optimized it using the LSMF filter to fairly compare the sequential run time against our parallel code run time.

We use this optimized version as the baseline to compare speedup for the test cases.  

[1] Foster, Erich L., Kai Hormann, and Romeo Traian Popa. "Clipping simple polygons with degenerate intersections." Computers & Graphics: X 2 (2019): 100007.


1) Introduction
---------------

This is a C++ implementation of our extension of the
Greiner-Hormann clipping algorithm, which computes the
intersection (or union) of two non-self-intersecting complex
polygons, with possibly multiple and nested components, even in
case of degenerate intersections (vertex on edge, overlapping
edges, etc.).


2) Compiling
------------

The code can be compiled by executing

  g++ -std=c++11 polyclip.cpp -o polyclip.exe

or
  gcc -g polyclip.cpp -lstdc++ -std=c++11 -o polyclip.exe


3) Executing
------------

To compute the intersection (or union) of two polygon, execute

  polyclip [-union] input1.poly input2.poly output.poly

where "-union" is the optional switch for computing the union
instead of the intersection, "input?.poly" are the names of the
two input polygon files and "output.poly" is the name of the
result polygon file.


4) File Format
--------------

The "*.poly" files must have the following structure. Each line
contains two numbers (int or double), the x and the y coordinates
of a vertex, followed by a "," or a ";", where the "," is used to
separate the vertices of a polygon component and ";" marks the end
of the component. For example, the following 7 lines:

0 0,
1 0,
0 1;
-0.5 -0.5,
1.5 -0.5,
1.5 1.5,
-0.5 1.5;

describe a polygon with 2 components, a right triangle inside a
square. All vertices in one file must be different from each
other.


5) Admitted Input
-----------------

The following features are allowed in the input polygons:

 - the vertex order in each component can be CW or CCW
 - components can be nested (AKA holes)
 - the two input polygons are allowed to have degenerate
   intersections (vertex on edge, overlapping edges, etc.)
   with each other

The following features are not allowed in the input polygons:

 - the polygons should not self-intersect (including degenerate
   self-intersections like vertex on vertex, vertex on edge),
   although the result will be correct as long as the self-
   intersection does not lie on the other polygon


6) Robustness
-------------

The implementation is based on C++ floating point numbers with
double precision and therefore not robust. The EPSILON parameter
(set to 0.000000001) is used as a tolerance for equality checks,
and two numbers are considered equal if their difference is less
than EPSILON.
