#ifndef ROUTINGSCENARIO_H
#define ROUTINGSCENARIO_H

/*
 * Solves a routing scenario, using an underlying constrained Delaunay triangulation.
 * For optimal paths, we compute a visibility graph on the vertices of polygons (obstacles).
 * The visibility graph is computed via an all-to-all visibility check via a
 * method similiar to the line_walk algorithm (see CGAL Triangulation_2), or
 * with a much faster triangle expansion algorithm. Then A* or Dijkstra.
 * For non-optimal results, we use the triangulation as a visibility graph and optimize the resulting path.
 *
 * A mapping between the Triangulation_2 Vertex_handles and indices is stored, we're using the indices
 * for describing the visibility graph.
 *
 * Obstacles of the routing scenario are stored as closed polygonal chains, that may only
 * intersect in its vertices. The saving is required for painting the interior of obstacles
 * (see QPainterPath and RoutingScenarioGraphicsItem).
 * A routing scenario can be saved using the polygonal chain decomposition.
 *
 * The algorithms only use the constrained edges of the triangulation and the in_domain values of
 * the triangle faces.
 */

#include <iostream>
#include <queue>
#include <CGAL/Timer.h>
#include "Triangulation_line_face_iterator.h"

namespace CGAL::Qt {
    template<typename T>
    class Routing_scenario {
        typedef typename T::Geom_traits Traits;
        typedef typename T::Point_2 Point_2;
        typedef typename T::Edge Edge;
        typedef typename T::Vertex_handle Vertex_handle;
        typedef typename T::Face_handle Face_handle;
        typedef typename T::Finite_vertices_iterator Finite_vertices_iterator;
        typedef typename T::Finite_faces_iterator Finite_faces_iterator;
        typedef typename T::Vertex_circulator Vertex_circulator;
        typedef typename T::Edge_circulator Edge_circulator;
        typedef typename T::Face_circulator Face_circulator;
        typedef typename T::Line_face_circulator Line_face_circulator;
        typedef typename T::Constrained_edges_iterator Constrained_edges_iterator;
        typedef Beltrami_klein_traits<Exact_predicates_inexact_constructions_kernel> Beltrami_klein_traits;
        typedef Triangulation_line_face_iterator<T> Line_face_iterator;

    public:
        Routing_scenario();

        //set new point to start point
        void set_start_point(const Point_2 p) {
            start_point = p;
            start_node_handle = insert_point(p);
        }

        void set_destination_point(const Point_2 p) {
            destination_point = p;
            destination_node_handle = insert_point(p);
        }

        //set already existing vertex to start
        void set_point_to_start(const Vertex_handle vp) {
            start_point = vp->point();
            start_node_handle = vp;
            path.clear();
            dijkstra_result.clear();
            defined_dijkstra = false;
            defined_path = false;
        }

        void set_point_to_destination(const Vertex_handle vp) {
            path.clear();
            defined_path = false;
            destination_point = vp->point();
            destination_node_handle = vp;
        }

        int number_of_vertices() {
            return t->number_of_vertices();
        }

        int number_of_faces() {
            return t->number_of_faces();
        }

        int edges_visibility_graph() {
            return adjacencies.size() / 2;
        }

        Vertex_handle insert_point(Point_2 p) {
            defined_domain = false;
            clear_graphs();
            Vertex_handle vh = t->insert(p);
            if (!vertex_index_map.contains(vh)) {
                vertex_index_map.insert({vh, number_of_vertices() - 1});
                index_vertex_map.push_back(vh);
            }
            return vh;
        }

        template<typename InputIterator>
        void insert_points(InputIterator first, InputIterator last) {
            defined_domain = false;
            clear_graphs();
            std::vector<Point_2> points(first, last);
            t->insert(points.begin(), points.end());
            vertex_index_map.clear();
            index_vertex_map.clear();
            int counter = 0;
            for (Finite_vertices_iterator fi = t->finite_vertices_begin(); fi != t->finite_vertices_end(); ++fi) {
                vertex_index_map.insert({fi, counter});
                index_vertex_map.push_back(fi);
                ++counter;
            }
        }

        template<typename InputIterator>
        void insert_obstacle(InputIterator first, InputIterator last, const bool closed) {
            defined_domain = false;
            clear_decomposition();
            clear_graphs();
            std::vector<Point_2> points(first, last);
            std::vector<Vertex_handle> obstacle;
            for (int i = 0; i < points.size(); i++) {
                Vertex_handle vh = t->insert(points[i]);
                obstacle.push_back(vh);
                if (!vertex_index_map.contains(vh)) {
                    vertex_index_map.insert({vh, number_of_vertices() - 1});
                    index_vertex_map.push_back(vh);
                }
            }
            obstacles.insert(obstacle);
            t->insert_constraint(points.begin(), points.end(), closed);
        }

        void insert_edge_constraint(const Vertex_handle va, const Vertex_handle vb) {
            defined_domain = false;
            clear_decomposition();
            clear_graphs();
            if (!vertex_index_map.contains(va)) {
                vertex_index_map.insert({va, number_of_vertices() - 1});
                index_vertex_map.push_back(va);
            }
            if (!vertex_index_map.contains(vb)) {
                vertex_index_map.insert({vb, number_of_vertices() - 1});
                index_vertex_map.push_back(vb);
            }
            t->insert_constraint(va, vb);
        }

        std::list<std::vector<Vertex_handle>> get_obstacles() {
            if (!defined_decomposition) {
                discover_polygon_decomposition();
            }
            std::list<std::vector<Vertex_handle> > set;
            for (std::vector<Vertex_handle> obstacle: obstacles) {
                std::vector<Vertex_handle> vertices;
                for (Vertex_handle vh: obstacle) {
                    vertices.push_back(vh);
                }
                set.push_back(vertices);
            }
            return set;
        }

        //remove vertices, requiered to be in triangulation
        void remove_vertices(std::set<Vertex_handle> set) {
            clear_routing_scenario();
            for (Vertex_handle vh: set) {
                t->remove(vh);
            }
            int counter = 0;
            for (Finite_vertices_iterator fvi = t->finite_vertices_begin(); fvi != t->finite_vertices_end(); ++fvi) {
                vertex_index_map.insert({fvi, counter});
                index_vertex_map.push_back(fvi);
                ++counter;
            }
        }

        void remove_constrained_edge(Face_handle fh, const int index) {
            clear_decomposition();
            t->remove_constrained_edge(fh, index);
        }

        void clear() {
            clear_routing_scenario();
            t->clear();
        }

        //computes path to destination_node from dijkstra_result
        //returns true if nodes are reachable
        bool get_path_from_dijkstra() {
            const int i = vertex_index_map[start_node_handle];
            path.clear();
            if (dijkstra_result[vertex_index_map[destination_node_handle]].first == -1) {
                defined_path = false;
                return false;
            }
            int node = vertex_index_map[destination_node_handle];
            do {
                path.push_back(node);
                node = dijkstra_result[node].first;
            } while (node != i);
            path.push_back(i);
            defined_path = true;
            return true;
        }

        Point_2 get_point_on_path(const int path_position) {
            return index_vertex_map[path[path_position]]->point();
        }

        Vertex_handle get_handle_on_path(const int path_position) {
            return index_vertex_map[path[path_position]];
        }

        std::vector<int> get_indices_path() {
            return path;
        }

        std::vector<int> get_point_path() {
            std::vector<Point_2> point_path;
            for (const int index: path) {
                point_path.push_back(index_vertex_map[path[index]]->point());
            }
            return point_path;
        }

        double get_path_length() {
            double sum = 0;
            for (int i = 0; i < path.size() - 1; ++i) {
                sum += hyperbolic_distance(index_vertex_map[path[i]]->point(),
                                           index_vertex_map[path[i + 1]]->point());
            }
            return sum;
        }

        void print_out_visibility_graph() {
            std::cout << "visibility graph: " << std::endl;
            //each vertex is part of the visibility graph
            for (Finite_vertices_iterator vp = t->finite_vertices_begin();
                 vp != t->finite_vertices_end(); vp++) {
                int index = vertex_index_map[vp];
                for (int i = offsets[index]; i < offsets[index + 1]; i++) {
                    std::cout << "vertex_handle with id: " << index
                            << " has neighbor: " << adjacencies[i] << std::endl;
                }
            }
            std::cout << "visibility graph end." << std::endl;
        }

        void print_out_offsets() const {
            std::cout << "offset vector: " << std::endl;
            for (const int i: offsets) {
                std::cout << i << std::endl;
            }
            std::cout << "offset vector end." << std::endl;
        }

        void print_out_adjacencies() {
            std::cout << "adjacency vector: " << std::endl;
            for (int i: adjacencies) {
                std::cout << index_vertex_map[i]->point() << std::endl;
            }
            std::cout << "adjacency vector end." << std::endl;
        }

        void print_out_pred() {
            std::cout << "predecessor vector: " << std::endl;
            for (int i = 0; i < number_of_vertices(); i++) {
                std::cout << "index: " << i << " , pred: " << dijkstra_result[i].first << std::endl;
            }
            std::cout << "predecessor vector end." << std::endl;
        }

        void print_out_vertex_index_mapping() {
            std::cout << "vertex->index mapping: " << std::endl;
            for (std::pair<Vertex_handle, int> pair: vertex_index_map) {
                std::cout << "par: (" << pair.first->point() << " , " << pair.second << ")" << std::endl;
            }
            std::cout << "vertex->index mapping end." << std::endl;
        }

        void print_out_index_vertex_mapping() {
            std::cout << "index->vertex mapping: " << std::endl;
            int counter = 0;
            for (Vertex_handle vh: index_vertex_map) {
                std::cout << "id: " << counter << " , " << vh->point() << ")" << std::endl;
                counter++;
            }
            std::cout << "index->vertex mapping end." << std::endl;
        }

        void test_suite() {
            std::cout << std::endl;
            std::cout << "TEST SUITE START" << std::endl;
            clear_graphs();

            Timer timer;
            timer.start();
            build_visibility_graph_naive();
            timer.stop();
            const double time_1 = timer.time();
            timer.reset();

            const std::vector<int> offsets_copy = offsets;
            const std::vector<int> adjacencies_copy = adjacencies;
            clear_graphs();

            timer.start();
            build_visibility_graph();
            timer.stop();
            const double time_2 = timer.time();
            timer.reset();

            std::cout << "1. Algorithm: Time: " << time_1 << std::endl;
            std::cout << "2. Algorithm: Time: " << time_2 << std::endl;
            std::cout << "1. Algorithm: Size of adjacency list = " << adjacencies_copy.size() << std::endl;
            std::cout << "2. Algorithm: Size of adjacency list = " << adjacencies.size() << std::endl;

            for (int i = 0; i < offsets.size(); ++i) {
                if (offsets_copy[i] != offsets[i]) {
                    std::cout << "Unequal results in offset vector! At index: " << i << std::endl;
                    break;
                }
            }

            //test shortest path
            //although if visibility graphs are equal, result are equal
            /*if (active_node_handle == nullptr) {
                set_point_active(index_vertex_map[0]);
            }
            if (destination_node_handle == nullptr) {
                set_point_to_destination(index_vertex_map[index_vertex_map.size() - 1]);
            }
            clear_graphs();
            build_visibility_graph_naive();
            compute_distances();
            const std::pair<double, int> results_1 = a_star();
            const std::vector<int> path_copy = path;
            clear_graphs();
            build_visibility_graph();
            compute_distances();
            const std::pair<double, int> results_2 = a_star();

            //compare reachability and length
            //-1 length means unreachable
            std::cout << "1. Algorithm: Path length = " << results_1.first << std::endl;
            std::cout << "2. Algorithm: Path length = " << results_2.first << std::endl;
            std::cout << "1. Algorithm: Path size = " << path_copy.size() << std::endl;
            std::cout << "2. Algorithm: Path size = " << path.size() << std::endl;

            for (int i = 0; i < path.size(); i++) {
                if (path[i] != path_copy[i]) {
                    std::cout << "Paths unequal!" << std::endl;
                    break;
                }
            }*/
            clear_graphs();
            std::cout << "TEST SUITE END" << std::endl;
            std::cout << std::endl;
        }

        Point_2 beltrami_klein_to_poincare(const Point_2 p) {
            const double abs = p.x() * p.x() + p.y() * p.y();
            double d = 1 + sqrt(1 - abs);
            return Point_2(p.x() / d, p.y() / d);
        }

        Point_2 poincare_to_beltrami_klein(const Point_2 p) {
            double d = 1 + p.x() * p.x() + p.y() * p.y();
            return Point_2(2 * p.x() / d, 2 * p.y() / d);
        }

        //debugging methods
        void compute_visible_vertices_from_start_node() {
            visibles_start_node.clear();
            for (Vertex_handle vh: what_points_can_p_see(start_node_handle)) {
                visibles_start_node.insert(vh);
            }
        }

        bool compute_intersected_faces();

        std::set<Vertex_handle> what_points_can_p_see(Vertex_handle vp);

        void dijkstra();

        //return true if destination node reachable from start node
        bool a_star();

        void build_visibility_graph_naive();

        bool can_p_see_q(Vertex_handle vp, Vertex_handle vq);

        void build_visibility_graph();

        void use_triangulation_as_visibility_graph();

        void path_optimization();

        //sets in_domain for every triangle in triangulation
        void discover_components();

        //underlying constrained delaunay triangulation
        T *t;

        std::set<std::vector<Vertex_handle> > obstacles;

        Point_2 start_point, destination_point;
        Vertex_handle start_node_handle;
        Vertex_handle destination_node_handle;

        std::unordered_map<Vertex_handle, int> vertex_index_map;
        std::vector<Vertex_handle> index_vertex_map;
        std::vector<int> offsets;
        std::vector<int> adjacencies;
        std::vector<double> distances;
        bool defined_visibility_graph;
        bool defined_domain;
        bool defined_decomposition;

        //at position i, predecessor and cost of i are stored
        std::vector<std::pair<int, double> > dijkstra_result;
        //get queries from dijkstra_result if it's computed
        bool defined_dijkstra;
        bool defined_path;

        //for debugging, crossed_faces, faces intersected by direct connection
        //between start and destination node
        std::set<Face_handle> crossed_faces;
        //for debugging, nodes visibile from start_node
        std::set<Vertex_handle> visibles_start_node;

    private:
        void clear_graphs() {
            adjacencies.clear();
            dijkstra_result.clear();
            distances.clear();
            offsets.clear();
            path.clear();
            defined_visibility_graph = false;
            defined_dijkstra = false;
            defined_path = false;
        }

        void clear_routing_scenario() {
            start_node_handle = nullptr;
            destination_node_handle = nullptr;
            vertex_index_map.clear();
            index_vertex_map.clear();
            crossed_faces.clear();
            visibles_start_node.clear();
            clear_graphs();
            defined_domain = false;
            clear_decomposition();
        }

        void clear_decomposition() {
            defined_decomposition = false;
            obstacles.clear();
        }

        double hyperbolic_distance(const Point_2 p_, const Point_2 q_) {
            Point_2 p = p_;
            Point_2 q = q_;

            //if we are in Beltrami Klein model, first translate p_ and q_ to Poincare disk model
            if (typeid(Beltrami_klein_traits) == typeid(Traits)) {
                p = beltrami_klein_to_poincare(p);
                q = beltrami_klein_to_poincare(q);
            }

            const double px = to_double(p.x());
            const double py = to_double(p.y());
            const double qx = to_double(q.x());
            const double qy = to_double(q.y());
            const double a = 1 - px * qx - py * qy;
            const double b = py * qx - px * qy;
            const double d = a * a + b * b;
            double real = (qx - px) * a + (qy - py) * b;
            double img = (qy - py) * a - (qx - px) * b;
            real = real / d;
            img = img / d;

            double r = std::sqrt(to_double(real * real + img * img));

            if (typeid(Beltrami_klein_traits) == typeid(Traits)) {
                r = 2 * r / (1 + r * r);
            }

            return std::log((1 + r) / (1 - r));
        }

        void discover_polygon_decomposition();

        void compute_distances();

        //computes visible vertices from c that lie in the angle LCR
        //L - point of l, R - point of r, C - point of c
        void visibles_from_angle(Edge e, Vertex_handle l, Vertex_handle r, Vertex_handle c, int *adjacency_counter);

        std::vector<int> path;
    };

    template<typename T>
    Routing_scenario<T>::Routing_scenario() : defined_domain(false), defined_visibility_graph(false),
                                              start_node_handle(nullptr), destination_node_handle(nullptr),
                                              defined_dijkstra(false), defined_decomposition(false),
                                              defined_path(false) {
        t = new T();
    }

    template<typename T>
    bool Routing_scenario<T>::compute_intersected_faces() {
        crossed_faces.clear();
        if (!defined_domain) {
            discover_components();
        }
        if (t->dimension() >= 2) {
            Line_face_circulator lfc = t->line_walk(start_point, destination_point, start_node_handle->face()),
                    done(lfc);
            if (lfc != nullptr) {
                if (lfc->has_vertex(destination_node_handle)) {
                    Face_handle fh = lfc->neighbor(lfc->cw(lfc->index(start_node_handle)));
                    if (lfc->is_in_domain() && fh->is_in_domain()) {
                        return false;
                    }
                    return true;
                }
                do {
                    if (lfc->is_in_domain()) {
                        return false;
                    }
                    if (!t->is_infinite(lfc)) {
                        crossed_faces.insert(lfc);
                    }
                    ++lfc;
                } while (!lfc->has_vertex(destination_node_handle));
                //for the last face
                if (!lfc->is_in_domain() && !t->is_infinite(lfc)) {
                    crossed_faces.insert(lfc);
                    return true;
                }
                return false;
            }
        }
        return true;
    }

    //returns all vertices, visible from p
    template<typename T>
    std::set<typename T::Vertex_handle> Routing_scenario<T>::what_points_can_p_see(Vertex_handle vp) {
        std::set<Vertex_handle> set;
        if (!defined_domain) {
            discover_components();
        }
        for (Finite_vertices_iterator vq = t->finite_vertices_begin();
             vq != t->finite_vertices_end(); ++vq) {
            if (can_p_see_q(vp, vq)) {
                set.insert(vq);
            }
        }
        return set;
    }

    template<typename T>
    bool Routing_scenario<T>::can_p_see_q(const Vertex_handle vp, const Vertex_handle vq/*, int* steps*/) {
        if (vp == vq) {
            return false;
        }
        //Line_face_circulator lfc = t->line_walk(p, q, vp->face()), done(lfc);
        Line_face_iterator lfc = Line_face_iterator(vp, t, vq->point());
        if (lfc != nullptr) {
            if (lfc->has_vertex(vq)) {
                //++*steps;
                //get neighbor triangle (share same edge)
                Face_handle fh = lfc->neighbor(lfc->cw(lfc->index(vp)));
                if (lfc->is_in_domain() && fh->is_in_domain()) {
                    return false;
                }
                return true;
            }
            do {
                //++*steps;
                if (lfc->is_in_domain()) {
                    return false;
                }
                ++lfc;
            } while (/*!lfc->has_vertex(vq)*/ lfc->vertex(lfc.i) != vq);
            if (lfc->is_in_domain()) {
                return false;
            }
        }
        return true;
    }

    //iterate through all pairs of vertices and check visibility of each pair
    template<typename T>
    void Routing_scenario<T>::build_visibility_graph_naive() {
        clear_graphs();
        if (!defined_domain) {
            discover_components();
        }
        int adjacency_counter = 0;
        offsets.push_back(0);
        for (int i = 0; i < number_of_vertices(); ++i) {
            for (int j = 0; j < number_of_vertices(); ++j) {
                if (i != j) {
                    if (i < j) {
                        if (can_p_see_q(index_vertex_map[i], index_vertex_map[j]/*, &steps*/)) {
                            adjacencies.push_back(j);
                            ++adjacency_counter;
                        }
                    } else {
                        for (int a = offsets[j]; a < offsets[j + 1]; ++a) {
                            if (adjacencies[a] == i) {
                                adjacencies.push_back(j);
                                ++adjacency_counter;
                                break;
                            }
                            if (adjacencies[a] > i) {
                                break;
                            }
                        }
                    }
                }
            }
            offsets.push_back(adjacency_counter);
        }
        compute_distances();
        defined_visibility_graph = true;
    }

    template<typename T>
    void Routing_scenario<T>::build_visibility_graph() {
        clear_graphs();
        if (!defined_domain) {
            discover_components();
        }
        int adjacency_counter = 0;
        offsets.push_back(0);

        //std::cout << "start building visibility graph" << std::endl;

        for (int i = 0; i < number_of_vertices(); ++i) {
            //std::cout << "i: " << i << std::endl;
            Vertex_handle center = index_vertex_map[i];
            //check visibility for all incident vertices
            Edge_circulator ec = t->incident_edges(center), done(ec);
            do {
                if (t->is_infinite(*ec)) {
                    continue;
                }
                Face_handle fh = ec->first;
                if (fh->is_in_domain() && t->mirror_edge(*ec).first->is_in_domain()) {
                    continue;
                }
                Vertex_handle nh = fh->vertex(fh->ccw(ec->second));
                ++adjacency_counter;
                adjacencies.push_back(vertex_index_map[nh]);
            } while (++ec != done);
            //for each incident face, we find the visible vertices enclosed by the rays starting at center
            Face_circulator fc = t->incident_faces(center), first(fc);
            do {
                //std::cout << "new face" << std::endl;
                if (t->is_infinite(fc) || fc->is_in_domain()) {
                    continue;
                }
                int index = fc->index(center);
                //Point_2 left = fc->vertex(fc->cw(index))->point();
                //Point_2 right = fc->vertex(fc->ccw(index))->point();
                Vertex_handle left = fc->vertex(fc->cw(index));
                Vertex_handle right = fc->vertex(fc->ccw(index));
                visibles_from_angle(Edge(fc, index), left, right, center/*->point()*/, &adjacency_counter);
            } while (++fc != first);
            offsets.push_back(adjacency_counter);
        }
        defined_visibility_graph = true;
        compute_distances();
    }

    template<typename T>
    void Routing_scenario<T>::visibles_from_angle(Edge e, Vertex_handle l, Vertex_handle r, Vertex_handle c,
                                                  int *adjacency_counter) {
        Edge edge = t->mirror_edge(e);
        Face_handle fh = edge.first;
        if (t->is_infinite(fh) || fh->is_in_domain()) {
            //std::cout << "recursion end" << std::endl;
            return;
        }

        /*std::cout << "left: " << l->point() << std::endl;
        std::cout << "right: " << r->point() << std::endl;
        std::cout << "adjacency counter: " << *adjacency_counter << std::endl;*/
        int index = edge.second;
        Vertex_handle vh = fh->vertex(index);
        Point_2 p = vh->point();

        //std::cout << "begin orientation" << std::endl;
        Orientation o_left = t->orientation(c->point(), l->point(), p);
        //std::cout << "finished orientation" << std::endl;
        if (o_left == LEFT_TURN) {
            //std::cout << "LEFT_TURN " << l->point() << std::endl;
            visibles_from_angle(Edge(fh, fh->ccw(index)), l, r, c, adjacency_counter);
            //std::cout << "l: " << l->point() << std::endl;
            return;
        }
        if (o_left == COLLINEAR) {
            //std::cout << "COLLINEAR " << l->point()  << std::endl;
            ++*adjacency_counter;
            adjacencies.push_back(vertex_index_map[vh]);
            visibles_from_angle(Edge(fh, fh->ccw(index)), l, r, c, adjacency_counter);
            return;
        }

        Orientation o_right = t->orientation(c->point(), r->point(), p);
        if (o_right == RIGHT_TURN) {
            //std::cout << "RIGHT_TURN " << l->point() << std::endl;
            visibles_from_angle(Edge(fh, fh->cw(index)), l, r, c, adjacency_counter);
            return;
        }
        if (o_right == COLLINEAR) {
            //std::cout << "COLLINEAR " << l->point() << std::endl;
            ++*adjacency_counter;
            adjacencies.push_back(vertex_index_map[vh]);
            visibles_from_angle(Edge(fh, fh->cw(index)), l, r, c, adjacency_counter);
            return;
        }

        ++*adjacency_counter;
        adjacencies.push_back(vertex_index_map[vh]);
        visibles_from_angle(Edge(fh, fh->cw(index)), l, vh, c, adjacency_counter);
        visibles_from_angle(Edge(fh, fh->ccw(index)), vh, r, c, adjacency_counter);
    }

    //only checking visibility for all incident vertices
    template<typename T>
    void Routing_scenario<T>::use_triangulation_as_visibility_graph() {
        clear_graphs();
        if (!defined_domain) {
            discover_components();
        }
        int adjacency_counter = 0;
        offsets.push_back(0);
        for (int i = 0; i < number_of_vertices(); ++i) {
            Vertex_handle vh = index_vertex_map[i];
            Edge_circulator ec = t->incident_edges(vh), done(ec);
            do {
                if (t->is_infinite(*ec)) {
                    continue;
                }
                Face_handle fh = ec->first;
                if (fh->is_in_domain() && t->mirror_edge(*ec).first->is_in_domain()) {
                    continue;
                }
                Vertex_handle nh = fh->vertex(fh->ccw(ec->second));
                ++adjacency_counter;
                adjacencies.push_back(vertex_index_map[nh]);
            } while (++ec != done);
            offsets.push_back(adjacency_counter);
        }
        compute_distances();
        defined_visibility_graph = true;
    }

    //to be called after use_triangulation_as_visibility_graph
    //optimizes path by checking visibilities of pairs of points
    template<typename T>
    void Routing_scenario<T>::path_optimization() {
        if (path.size() < 3) {
            return;
        }
        std::vector<int> shortend_path;
        shortend_path.push_back(path[0]);
        int index = 0;
        while (index != path.size() - 1) {
            int check = path.size() - 1;
            while (!can_p_see_q(index_vertex_map[path[index]], index_vertex_map[path[check]])) {
                --check;
            }
            shortend_path.push_back(path[check]);
            index = check;
        }
        shortend_path.push_back(path[path.size() - 1]);
        path = shortend_path;
    }

    //uses lazy-deletion, (see https://stackoverflow.com/questions/9209323/easiest-way-of-using-min-priority-queue-with-key-update-in-c)
    //note that in our scenario, visibility graph is not always connected
    template<typename T>
    void Routing_scenario<T>::dijkstra() {
        auto cmp = [](const std::pair<int, double> &a, const std::pair<int, double> &b) -> bool {
            return a.second > b.second;
        };
        std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >,
            decltype(cmp)> queue(cmp);
        std::vector<double> costs(number_of_vertices());
        std::vector<int> predecessors(number_of_vertices());
        std::vector<bool> closed(number_of_vertices());
        const int start_index = vertex_index_map[start_node_handle];
        std::fill(closed.begin(), closed.end(), false);
        std::fill(predecessors.begin(), predecessors.end(), -1);
        std::fill(costs.begin(), costs.end(), DBL_MAX);
        predecessors[start_index] = start_index;
        costs[start_index] = 0;
        queue.push(std::make_pair(start_index, 0));
        while (!queue.empty()) {
            const std::pair<int, double> node = queue.top();
            const int index = node.first;
            queue.pop();
            if (!closed[index]) {
                closed[index] = true;
                //std::cout << "index: " << index << std::endl;
                for (int i = offsets[index]; i < offsets[index + 1]; i++) {
                    int ni = adjacencies[i];
                    //std::cout << "neighbor index: " << ni << std::endl;
                    if (!closed[ni]) {
                        double new_d = node.second + distances[i];
                        //std::cout << "costs[index]: " << node.second << std::endl;
                        //std::cout << "distances[i]: " << distances[i] << std::endl;
                        if (new_d < costs[ni]) {
                            costs[ni] = new_d;
                            //std::cout << "set costs to be: " << costs[ni] << std::endl;
                            predecessors[ni] = index;
                            //std::cout << "set predecessor of: " << ni << " to be: " << predecessors[ni] << std::endl;
                            queue.push(std::make_pair(ni, new_d));
                            //std::cout << "queue top: " << queue.top() << std::endl;
                        }
                    }
                }
            }
        }
        defined_dijkstra = true;
        std::vector<std::pair<int, double> > target(number_of_vertices());
        for (int i = 0; i < number_of_vertices(); i++) {
            target[i] = std::make_pair(predecessors[i], costs[i]);
        }
        dijkstra_result = target;
    }

    template<typename T>
    bool Routing_scenario<T>::a_star() {
        auto cmp = [](const std::pair<int, double> a, const std::pair<int, double> b) -> bool {
            return a.second > b.second;
        };
        std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >,
            decltype(cmp)> queue(cmp);
        std::vector<double> costs(number_of_vertices());
        std::vector<double> heuristic(number_of_vertices());
        std::vector<int> predecessors(number_of_vertices());
        std::vector<bool> closed(number_of_vertices());
        const int start_index = vertex_index_map[start_node_handle];
        const int destination_index = vertex_index_map[destination_node_handle];
        std::fill(closed.begin(), closed.end(), false);
        std::fill(predecessors.begin(), predecessors.end(), -1);
        std::fill(costs.begin(), costs.end(), DBL_MAX);
        //heuristic is segment length between i and destination_node
        for (int i = 0; i < number_of_vertices(); i++) {
            if (i != destination_index) {
                heuristic[i] = hyperbolic_distance(index_vertex_map[i]->point(), destination_point);
            } else {
                heuristic[i] = 0;
            }
        }
        predecessors[start_index] = start_index;
        costs[start_index] = 0;
        queue.push(std::make_pair(start_index, 0));
        int relaxed_vertices_counter = 0;
        while (!queue.empty()) {
            std::pair<int, double> node = queue.top();
            int index = node.first;
            if (index == destination_index) {
                break;
            }
            queue.pop();
            if (!closed[index]) {
                relaxed_vertices_counter++;
                closed[index] = true;
                for (int i = offsets[index]; i < offsets[index + 1]; i++) {
                    int ni = adjacencies[i];
                    if (!closed[ni]) {
                        double new_d = costs[index] + distances[i];
                        if (new_d < costs[ni]) {
                            costs[ni] = new_d;
                            predecessors[ni] = index;
                            queue.push(std::make_pair(ni, new_d + heuristic[ni]));
                        }
                    }
                }
            }
        }

        //std::cout << "visited vertices: " < relaxed_vertices_counter << std::endl;

        //if destination is unreachable
        if (predecessors[destination_index] == -1) {
            defined_path = false;
            return false;
        }
        //if reachable, then set the path
        int node = destination_index;
        path.clear();
        do {
            path.push_back(node);
            node = predecessors[node];
        } while (node != start_index);
        path.push_back(start_index);
        defined_path = true;
        return true;
    }

    template<typename T>
    void Routing_scenario<T>::discover_polygon_decomposition() {
        obstacles.clear();

        std::vector<Vertex_handle> visited_points;
        std::set<Edge> edges;

        for (Constrained_edges_iterator cei = t->constrained_edges_begin(); cei != t->constrained_edges_end(); ++cei) {
            edges.insert(*cei);
        }

        while (!edges.empty()) {
            Edge edge = *edges.begin();
            Vertex_handle start_vh = edge.first->vertex(edge.first->cw(edge.second));
            Vertex_handle current_vh = start_vh;
            do {
                visited_points.push_back(current_vh);
                Edge_circulator ec = t->incident_edges(current_vh), done(ec);
                do {
                    Edge mirrored_edge = t->mirror_edge(*ec);
                    if (edges.contains(*ec) || edges.contains(mirrored_edge)) {
                        const Face_handle fh = ec->first;
                        int index = fh->index(current_vh);
                        current_vh = fh->vertex((3 - (index + ec->second)) % 3);
                        edges.erase(*ec);
                        edges.erase(mirrored_edge);
                        break;
                    }
                } while (++ec != done);
            } while (current_vh != start_vh);
            obstacles.insert(visited_points);
            visited_points.clear();
        }
        defined_decomposition = true;
    }

    template<typename T>
    void Routing_scenario<T>::discover_components() {
        if (t->dimension() != 2)
            return;

        Face_handle start = t->infinite_face();
        std::list<Face_handle> queue;
        queue.push_back(start);
        std::set<Face_handle> traversed;
        traversed.insert(start);
        start->set_in_domain(false);

        while (!queue.empty()) {
            Face_handle fh = queue.front();
            queue.pop_front();

            for (int i = 0; i < 3; i++) {
                Face_handle fi = fh->neighbor(i);
                if (!traversed.contains(fi)) {
                    queue.push_back(fi);
                    if (t->is_constrained(Edge(fh, i))) {
                        fi->set_in_domain(!fh->is_in_domain());
                    } else {
                        fi->set_in_domain(fh->is_in_domain());
                    }
                    traversed.insert(fi);
                }
            }
        }
        defined_domain = true;
    }

    template<typename T>
    void Routing_scenario<T>::compute_distances() {
        for (int i = 0; i < number_of_vertices(); i++) {
            Vertex_handle vh = index_vertex_map[i];
            for (int j = offsets[i]; j < offsets[i + 1]; j++) {
                int ni = adjacencies[j];
                Vertex_handle nh = index_vertex_map[ni];
                if (i < ni) {
                    distances.push_back(hyperbolic_distance(vh->point(), nh->point()));
                } else {
                    for (int a = offsets[ni]; a < offsets[ni + 1]; a++) {
                        if (adjacencies[a] == i) {
                            distances.push_back(distances[a]);
                        }
                    }
                }
            }
        }
    }
}


#endif //ROUTINGSCENARIO_H
