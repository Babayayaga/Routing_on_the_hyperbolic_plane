#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include "include/internal/BeltramiKleinTraits.h"
#include "include/internal/PoincareDiskTraits.h"
#include "include/internal/RandomDomainGenerator.h"
#include "include/internal/Qt/Input/PointInput.h"
#include "include/internal/RoutingScenario.h"


int main() {
    typedef CGAL::Beltrami_klein_traits<CGAL::Exact_predicates_inexact_constructions_kernel> Beltrami_klein_traits;
    typedef CGAL::Poincare_disk_traits<CGAL::Exact_predicates_inexact_constructions_kernel> Poincare_disk_traits;
    //set either to Poincare_disk_traits or Beltrami_klein_traits
    typedef Beltrami_klein_traits K;
    typedef CGAL::Delaunay_mesh_face_base_2<K> Face_base;
    typedef CGAL::Triangulation_vertex_base_2<K> Vertex_base;
    typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base> TDS;
    typedef CGAL::No_constraint_intersection_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> T;
    typedef CGAL::Qt::Routing_scenario<T> RoutingScenario;
    typedef K::Point_2 Point_2;
    auto routing_scenario = RoutingScenario();
    auto random_generator = CGAL::Qt::Random_domain_generator<T>(&routing_scenario);

    int number_of_points;
    int min;
    double radius;
    double threshold;
    bool remove_points_in_free_space;
    bool end;
    bool new_domain;
    do {
        std::cout << "-----Random generation-----" << std::endl;
        std::cout << "Number of points: ";
        std::cin >> number_of_points;
        std::cout << "Size smallest polygonal chain: ";
        std::cin >> min;
        std::cout << "Radius: ";
        std::cin >> radius;
        std::cout << "Threshold: ";
        std::cin >> threshold;
        std::cout << "Disregard points with no incident constraints (0/1)? ";
        std::cin >> remove_points_in_free_space;
        std::cout << "Generating..." << std::endl;
        random_generator.generate_random_domain(number_of_points, radius, threshold, 5, 5, 5, true,
                                                remove_points_in_free_space);
        std::cout << "Result: Number of points: " << routing_scenario.number_of_vertices() << ". Number of faces: " <<
                routing_scenario.number_of_faces() << std::endl;

        do {
            int start_id;
            int end_id;
            int algorithm;
            bool a_star_or_dijkstra;
            bool optimize_path;
            std::cout << "-----Routing-----" << std::endl;
            std::cout << "Select id of start node: ";
            std::cin >> start_id;
            std::cout << "Select id of end node: ";
            std::cin >> end_id;
            routing_scenario.set_point_to_start(routing_scenario.index_vertex_map[start_id]);
            routing_scenario.set_point_to_destination(routing_scenario.index_vertex_map[end_id]);

            std::cout << "Choose algorithm: " << std::endl;
            std::cout <<
                    "(0) All-to-all visibility check (1) Triagular expansion algorithm (2) Routing on triangulation: ";
            std::cin >> algorithm;
            if (algorithm == 0) {
                routing_scenario.build_visibility_graph_naive();
            } else if (algorithm == 1) {
                routing_scenario.build_visibility_graph();
            } else {
                routing_scenario.use_triangulation_as_visibility_graph();
            }
            std::cout << "Computing visibility graph..." << std::endl;
            std::cout << "Result: Number of edges: " << routing_scenario.edges_visibility_graph() << std::endl;
            std::cout << "Compute A* or Dijkstra (0/1)? ";
            std::cin >> a_star_or_dijkstra;
            if (a_star_or_dijkstra) {
                routing_scenario.dijkstra();
                std::cout << "Computing Dijkstra..." << std::endl;
                routing_scenario.get_path_from_dijkstra();
            } else {
                routing_scenario.a_star();
                std::cout << "Computing A*..." << std::endl;
            }
            std::cout << "Result: Path length: " << routing_scenario.get_path_length()
                    << ". Number of nodes: " << routing_scenario.get_indices_path().size() << std::endl;

            std::cout << "Optimize path (0/1)?  ";
            std::cin >> optimize_path;
            if (optimize_path) {
                routing_scenario.path_optimization();
                std::cout << "Optimize path..." << std::endl;
            }

            std::cout << "Continue with current polygonal domain (0/1)? ";
            std::cin >> new_domain;
        } while (!new_domain);
            std::cout << "Again (0/1)? ";
            std::cin >> end;
    } while (!end);
    exit(0);
}