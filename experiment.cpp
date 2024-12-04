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
        std::cout << "Generating..." << std::endl;
        random_generator.generate_random_domain(number_of_points, radius, threshold, 5, 5, 5, true, false, 0);
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
            CGAL::Timer timer;
            if (algorithm == 0) {
                std::cout << "Computing visibility graph..." << std::endl;
                timer.start();
                routing_scenario.build_visibility_graph_naive();
                timer.stop();
            } else if (algorithm == 1) {
                std::cout << "Computing visibility graph..." << std::endl;
                timer.start();
                routing_scenario.build_visibility_graph();
                timer.stop();
            } else {
                int additional_points = 0;
                std::cout << "Insert uniformly distributed random points? N: ";
                std::cin >> additional_points;
                if (additional_points > 0) {
                    std::cout << "Inserting additional points..." << std::endl;
                    random_generator.insert_uniformly_distributed_points(additional_points, radius);
                    routing_scenario.discover_components();
                }
                std::cout << "Computing visibility graph..." << std::endl;
                timer.start();
                routing_scenario.use_triangulation_as_visibility_graph();
                timer.stop();
            }
            std::cout << "Computing visibility graph took: " << timer.time() << " seconds. Number of edges: "
                    << routing_scenario.edges_visibility_graph() << " Number of nodes: "
                    << routing_scenario.number_of_vertices() << std::endl;
            timer.reset();
            std::cout << "Compute A* or Dijkstra (0/1)? ";
            std::cin >> a_star_or_dijkstra;
            if (a_star_or_dijkstra) {
                std::cout << "Computing Dijkstra..." << std::endl;
                routing_scenario.dijkstra();
                routing_scenario.get_path_from_dijkstra();
                std::cout << "Average path length Dijkstra: " << routing_scenario.average_path_length_dijkstra() <<
                        std::endl;
            } else {
                std::cout << "Computing A*..." << std::endl;
                routing_scenario.a_star();
            }

            std::cout << "Optimize path (0/1)?  ";
            std::cin >> optimize_path;
            if (optimize_path) {
                routing_scenario.path_optimization();
            }

            std::cout << "Continue with current polygonal domain (0/1)? ";
            std::cin >> new_domain;
        } while (new_domain);
        std::cout << "Again (0/1)? ";
        std::cin >> end;
    } while (end);


    /*int error = 0;
    for(int i = 0; i < 1000; ++i) {
        random_generator.generate_random_domain(20*1000, 5, 0.1, 5, 5, 5, true, false, 0);
        if(!routing_scenario.test_suite()) {
            error = 1;
            break;
        }
    }
    std::cout << "finished with " << error << std::endl;
    exit(0);*/
}
