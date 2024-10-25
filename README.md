# Routing on the hyperbolic plane
This repository contains the code for the thesis `Routing on the hyperbolic plane`. It solves and visualizes the analoge of the [Euclidean shortest 
path problem](https://en.wikipedia.org/wiki/Euclidean_shortest_path) on the hyperbolic plane, calculation and visualization are carried out in the 
[Beltrami-Klein model](https://en.wikipedia.org/wiki/Beltrami%E2%80%93Klein_model) and the [Poincaré disk model](https://en.wikipedia.org/wiki/Poincar%C3%A9_disk_model). A routing problem in a hyperbolic polygonal domain,
is solved by first computing a [constrained Delaunay triangulation (CGAL)](https://doc.cgal.org/latest/Triangulation_2/index.html#Section_2D_Triangulations_Constrained_Delaunay)
as navigation mesh. Then, we use a [triangular expansion algorithm](https://www.researchgate.net/figure/Triangular-expansion-algorithm-recursion-entering-triangle-through-edge-e_fig1_260873070) 
to construct a visibility graph, now A* or Dijkstra are applied. For large problems, we can choose to restrict the visibiliy graph to be
a subgraph of the triangulation, the results are non-optimal. A non-optimal path between start- and endpoint can be improved by computing point-to-point
visibility tests. These visibility tests use the line_walk method, provided by the [Trianulation_2 class](https://pages.github.com/) of CGAL.
## Traits for the constrained Delaunay triangulation
To be able to use CGALs Triangulation classes in the Beltrami-Klein and Poincaré disk model, we specified with the Traits classes `BeltramiKleinTraits` `PoincareDiskTraits`, how they behave. Both classes inherit by default from the [`Exact_predicates_inexact_constructions_kernel`](https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html). The `BeltramiKleinTraits` overwrite the `Side_of_oriented_circle_2` predicate. The `PoincareDiskTraits` overwrite the `Orientation_2` predicate and the `Construct_segment_2` construction. In the Beltrami-Klein model segments between two points are Euclidean segments. So the required `Orientation_2` predicate is the Euclidean one, thus exact. Whereas in the Poincaré disk model, a segment between two points is the circular arc contained in the circle perpendicular to the unit circle. We can only get the orientation of three points, by constructing this circular arc through two of them, and checking whether the third points is in the circles iterior. Since this involves a construction, the `Exact_predicates_exact_constructions_kernel`
is needed for exact results. Alltogether, computing the orientation of three points is faster in the Beltrami-Klein model, thus all algorithms in this project run faster in it (see thesis for detailed comparision). 
## Generating a polygonal domain
Due to the lack of "real" routing problems on the hyperbolic plane, we implemented the random generation of polygons (see `RandomDomainGenerator`). First, we sample a number of random points within the disk $D$ with Euclidean radius $0 < R_e < 1$. A point is sampled by sampling its angle uniformly, while the radial coordinate density is given by $\rho(r)=\sinh(r) / (\cosh(R)-1)$, where $R$ is the hyperbolic radius of $D$. On the point set we compute a Delaunay triangulation and choose the `in_domain` value ([`Delaunay_mesh_face_base_2`](https://doc.cgal.org/latest/Mesh_2/classCGAL_1_1Delaunay__mesh__face__base__2.html)) for its triangles randomly. We apply some "smoothing" steps, among others the following one: Every triangle not `in_domain` with more than one neighboring `in_domain` triangle, is set to be `in_domain`. And vice versa.
From this we infer the constrained edges of the triangulation and insert them.

<img src="https://github.com/user-attachments/assets/c74927c3-6a4e-46ee-ab8c-5ae0fb7d053e" width="300">
<img src="https://github.com/user-attachments/assets/f24cd5ae-f9bd-41f7-9df6-ac22e8955da4" width="300">

(Routing scenarios with ~20k points, left Beltrami-Klein model, right Poincaré disk)

