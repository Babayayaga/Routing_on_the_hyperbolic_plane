# Routing on the hyperbolic plane
This repository contains the code for the bachelor thesis _Routing on the Hyperbolic Plane_ at the University of Stuttgart. The implementation is written in C++ and utilizes the constrained Delaunay triangulation class of CGAL (version 5.6.1). The visualization is built using Qt 5, and the images below were captured from the applet.

## Abstract

 The Euclidean shortest path problem is a problem in computational geometry,
 which requires finding a shortest path between two points in the Euclidean
 plane avoiding a set of polygonal obstacles. This thesis explores the analog of
 this problem on a surface with different geometry, the hyperbolic plane. To
 remain within the Euclidean framework, we introduce the Beltrami-Klein and
 Poincaré disk models- two prominent map projections of the hyperbolic plane
 into the Euclidean plane. Using these models, we adapt and implement the
 canonical approach that solves the Euclidean shortest path problem for the
 hyperbolic plane. This approach involves triangulating the polygonal domain
 to derive the visibility graph from the set of obstacles. A query between two
 points can then be answered by applying shortest-path algorithms, such as A*
 or Dijkstra, on this graph. Furthermore, we evaluate an approximation scheme
 that runs these shortest-path algorithms on a subgraph of the triangulation.
 Due to the absence of polygonal land maps in the hyperbolic plane, an
 algorithm to generate them has been developed and implemented.

## The applet
In addition to the implemented algorithms, the applet provides the following functions
- freely moving around within the Poincaré disk and Beltrami-Klein models to explore polygonal domains
- manually drawing of polygonal domains and storing/loading the created/generated files
- displaying various visibility/routing properties interactively
- exporting the viewport as .svg or .png

<img src="https://github.com/user-attachments/assets/c74927c3-6a4e-46ee-ab8c-5ae0fb7d053e" width="230">
<img src="https://github.com/user-attachments/assets/f24cd5ae-f9bd-41f7-9df6-ac22e8955da4" width="230">

<img src = "https://github.com/user-attachments/assets/456ead60-7d66-4148-a694-6db457693da8" width = "230">
<img src = "https://github.com/user-attachments/assets/9683c3b4-6416-4fec-9bdd-af2eacfbb33d" width = "230">
<img src = "https://github.com/user-attachments/assets/6094d5b6-7f74-4c01-a8f8-cbaa4f3ee476" width = "230">
<img src = "https://github.com/user-attachments/assets/d6063c26-3552-42f6-a2b9-5fae6663017f" width = "230">

<img src = "https://github.com/user-attachments/assets/6446ea86-00d4-4683-8a27-18fa84d3aee0" width = "230">
<img src = "https://github.com/user-attachments/assets/47258cd5-2872-4e44-a02e-f9660f744cce" width = "230">
<img src = "https://github.com/user-attachments/assets/d43c9477-93fb-43f4-b00a-8784c86ee8db" width = "230">


