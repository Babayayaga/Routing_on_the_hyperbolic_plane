# Routing on the hyperbolic plane
This repository contains the code for the thesis `Routing on the hyperbolic plane`. It solves and visualizes the analoge of the [Euclidean shortest 
path problem](https://en.wikipedia.org/wiki/Euclidean_shortest_path) on the hyperbolic plane, calculation and visualization are carried out in the 
[Beltrami-Klein model](https://en.wikipedia.org/wiki/Beltrami%E2%80%93Klein_model) and the [Poincaré disk model](https://en.wikipedia.org/wiki/Poincar%C3%A9_disk_model). Here, a routing problem in a hyperbolic polygonal domain
is solved by first computing a [constrained Delaunay triangulation (CGAL)](https://doc.cgal.org/latest/Triangulation_2/index.html#Section_2D_Triangulations_Constrained_Delaunay)
as navigation mesh. Then, we use a [triangular expansion algorithm](https://www.researchgate.net/figure/Triangular-expansion-algorithm-recursion-entering-triangle-through-edge-e_fig1_260873070) 
to construct a visibility graph, now A* or Dijkstra are applied. For large problems, we can choose to restrict the visibiliy graph to be
a subgraph of the triangulation, the results are non-optimal. A non-optimal path between start- and endpoint can be improved by computing point-to-point
visibility tests with points on the path. These visibility tests use the `line_walk` method, provided by the [Triangulation_2 class](https://pages.github.com/) of CGAL.
## Traits for the constrained Delaunay triangulation
To be able to use CGALs Triangulation classes in the Beltrami-Klein and Poincaré disk model, we specified with the Traits classes `BeltramiKleinTraits` `PoincareDiskTraits` how they behave. Both classes inherit by default from the [`Exact_predicates_inexact_constructions_kernel`](https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html). The `BeltramiKleinTraits` class overwrites the `Side_of_oriented_circle_2` predicate. The `PoincareDiskTraits` class overwrites the `Orientation_2` predicate and the `Construct_segment_2` construction. In the Beltrami-Klein model, segments between two points are Euclidean segments. So the required `Orientation_2` predicate is the Euclidean one, thus guaranteed to be exact. Whereas, in the Poincaré disk model, a segment between two points is the circular arc that is contained in the unique circle through both points, perpendicular to the unit circle. Calculating the orientation of three points involves constructing this circular arc through two of them, and checking whether the third point is in the circles iterior. Since this involves a construction, the `Exact_predicates_exact_constructions_kernel`
is needed for exact results. Alltogether, computing the orientation of three points is faster in the Beltrami-Klein model, thus all algorithms in this project run faster in it (see thesis for detailed comparision). 
## Generating a polygonal domain
Due to the lack of "real" routing problems on the hyperbolic plane, we implemented random generation of polygons (see `RandomDomainGenerator`). First, we sample a number of random points within the disk $D$ with Euclidean radius $0 < R_e < 1$. A point is sampled by sampling its angle uniformly, while the radial coordinate density is given by $\rho(r)=\sinh(r) / (\cosh(R)-1)$, where $R$ is the hyperbolic radius of $D$. On the point set we compute a Delaunay triangulation and choose the `in_domain` value ([`Delaunay_mesh_face_base_2`](https://doc.cgal.org/latest/Mesh_2/classCGAL_1_1Delaunay__mesh__face__base__2.html)) for its triangles randomly. We apply some "smoothing" steps, among others the following one: Every triangle not `in_domain` with more than one neighboring `in_domain` triangle, is set to be `in_domain`. And vice versa.
In the end, we infer the constrained edges of the triangulation and insert them as constraints.

<img src="https://github.com/user-attachments/assets/c74927c3-6a4e-46ee-ab8c-5ae0fb7d053e" width="300">
<img src="https://github.com/user-attachments/assets/f24cd5ae-f9bd-41f7-9df6-ac22e8955da4" width="300">

(Routing scenarios with ~20k points, left Beltrami-Klein model, right Poincaré disk)

## Visualization
For displaying, we use Qts `QMainWindow`, it contains a `QGraphicsScene` and a `QGraphicsView`. The two custom `QGraphicsItems`, namly `TriangulationGraphicsItem` and `RoutingScenarioGraphicsItem` are added to the scene, they are painted when their `paint()` method is called by the scene. The `RoutingScenarioGraphicsItem` is parent of the `TriangulationGraphicsItem`. All routing specific variables are stored inside the chief class of this project, `RoutingScenario`. The `RoutingScenarioGraphicsItem` has acess to them and visualizes them, if enabled. The triangulation is painted by the `TriangulationGraphicsItem`, constrained edges can be drawn seperatly, as well as only the edges triangulating the free space. Whenever `TriangulationGraphicsItem` or `RoutingScenarioGraphicsItem` want to paint an edge between two nodes, they pass their intern `QPainter` to `HyperbolicPainter`. From a given start and end point, `HyperbolicPainter` can draw and construct Euclidean lines or circular arcs connecting both points, by calling the `construct_segment_2_object` from our traits-classes. We note that the `drawArc()` method of the `QPainter` takes the spanning angle as an integer, this means that for arcs contained in large circles, the arc doesn't seem to connect the two points, because it's so inaccurate. Therefore, we can set the `drawing approximation` parameter, to specify the minimum radius of an arc, that is drawn as a Euclidean segment. 

<img src = "https://github.com/user-attachments/assets/456ead60-7d66-4148-a694-6db457693da8" width = "230">
<img src = "https://github.com/user-attachments/assets/9683c3b4-6416-4fec-9bdd-af2eacfbb33d" width = "230">
<img src = "https://github.com/user-attachments/assets/6094d5b6-7f74-4c01-a8f8-cbaa4f3ee476" width = "230">
<img src = "https://github.com/user-attachments/assets/d6063c26-3552-42f6-a2b9-5fae6663017f" width = "230">

We can infer closed polygonal chains, that may only intersect in their vertices, from the constrained edges of the triangulation. This decomposition into polygonal chains is used to store and load in a `RoutingScenario`, to disregard small obstacles and to draw the iterior of an obstacle with a `QPainterPath`. The fill rule of `QPainterPath` is set to `Qt::OddEvenFill`.

<img src = "https://github.com/user-attachments/assets/6446ea86-00d4-4683-8a27-18fa84d3aee0" width = "230">
<img src = "https://github.com/user-attachments/assets/47258cd5-2872-4e44-a02e-f9660f744cce" width = "230">
<img src = "https://github.com/user-attachments/assets/d43c9477-93fb-43f4-b00a-8784c86ee8db" width = "230">
<img src = "https://github.com/user-attachments/assets/c4f24ee7-48fc-410c-b7a9-2fadf095de96" width = "230">

(Beltrami-Klein model left, Poincaré disk model right: Decomposition of the constrained edges into poylgonal chains. Starting both times with a triangulation on 1mio. points, aka the level of detail. On the left, minimum chain size is set to 10 and ~300k points resulted to be on obstacles. On the right, minimum chain size is set to 1000 and in total ~100k points on obstacles.)
