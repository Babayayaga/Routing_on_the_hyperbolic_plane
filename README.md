# Routing on the hyperbolic plane
This repository contains the code for the thesis `Routing on the hyperbolic plane`. It solves and visualizes the analoge of the [Euclidean shortest 
path problem](https://en.wikipedia.org/wiki/Euclidean_shortest_path) on the hyperbolic plane, calculation and visualization are carried out in the 
[Beltrami-Klein model](https://en.wikipedia.org/wiki/Beltrami%E2%80%93Klein_model) and the [Poincaré disk model](https://en.wikipedia.org/wiki/Poincar%C3%A9_disk_model). The problem 
is solved by first computing a [constrained Delaunay triangulation (CGAL)](https://doc.cgal.org/latest/Triangulation_2/index.html#Section_2D_Triangulations_Constrained_Delaunay)
as navigation mesh. Then, we use a [triangular expansion algorithm](https://www.researchgate.net/figure/Triangular-expansion-algorithm-recursion-entering-triangle-through-edge-e_fig1_260873070) 
to construct a visibility graph, now A* or Dijkstra are applied. For large problems, we restrict the visibiliy graph to be
a subgraph of the triangulation, the results are non-optimal. A non-optimal path between start- and endpoint can be improved by computing point-to-point
visibility tests. These visibility tests use the line_walk method, provided by the [Trianulation_2 class](https://pages.github.com/) of CGAL.
## Generating a polygonal domain
