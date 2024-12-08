#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include "include/internal/BeltramiKleinTraits.h"
#include "include/internal/PoincareDiskTraits.h"
#include "include/internal/RandomDomainGenerator.h"
#include "include/internal/Qt/Input/PointInput.h"
#include "include/internal/RoutingScenario.h"


int main() {
    typedef CGAL::Beltrami_klein_traits<CGAL::Exact_predicates_inexact_constructions_kernel> Beltrami_klein_traits;
    typedef CGAL::Poincare_disk_traits<> Poincare_disk_traits;
    //set either to Poincare_disk_traits or Beltrami_klein_traits
    typedef Beltrami_klein_traits B;
    typedef Poincare_disk_traits P;

    typedef CGAL::Delaunay_mesh_face_base_2<P> Face_base_P;
    typedef CGAL::Triangulation_vertex_base_2<P> Vertex_base_P;
    typedef CGAL::Triangulation_data_structure_2<Vertex_base_P, Face_base_P> TDS_P;

    typedef CGAL::Delaunay_mesh_face_base_2<B> Face_base_B;
    typedef CGAL::Triangulation_vertex_base_2<B> Vertex_base_B;
    typedef CGAL::Triangulation_data_structure_2<Vertex_base_B, Face_base_B> TDS_B;
    typedef CGAL::No_constraint_intersection_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<P, TDS_P, Itag> T_P;
    typedef CGAL::Constrained_Delaunay_triangulation_2<B, TDS_B, Itag> T_B;
    typedef CGAL::Qt::Routing_scenario<T_P> RoutingScenario_P;
    typedef CGAL::Qt::Routing_scenario<T_B> RoutingScenario_B;
    typedef P::Point_2 Point_2;
    auto routing_scenario_p = RoutingScenario_P();
    auto random_generator_p = CGAL::Qt::Random_domain_generator<T_P>(&routing_scenario_p);

    auto routing_scenario_b = RoutingScenario_B();
    auto random_generator_b = CGAL::Qt::Random_domain_generator<T_B>(&routing_scenario_b);


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
        random_generator_p.generate_random_domain(number_of_points, radius, threshold, 5, 5, 5, true, false, 0);
        std::cout << "Result: Number of points: " << routing_scenario_p.number_of_vertices() << ". Number of faces: " <<
                routing_scenario_p.number_of_faces() << std::endl;

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
            routing_scenario_p.set_point_to_start(routing_scenario_p.index_vertex_map[start_id]);
            routing_scenario_p.set_point_to_destination(routing_scenario_p.index_vertex_map[end_id]);

            std::cout << "Choose algorithm: " << std::endl;
            std::cout <<
                    "(0) All-to-all visibility check (1) Triagular expansion algorithm (2) Routing on triangulation: ";
            std::cin >> algorithm;
            CGAL::Timer timer;
            if (algorithm == 0) {
                std::cout << "Computing visibility graph..." << std::endl;
                timer.start();
                routing_scenario_p.build_visibility_graph_naive();
                timer.stop();
            } else if (algorithm == 1) {
                std::cout << "Computing visibility graph..." << std::endl;
                timer.start();
                routing_scenario_p.build_visibility_graph();
                timer.stop();
            } else {
                int additional_points = 0;
                std::cout << "Insert uniformly distributed random points? N: ";
                std::cin >> additional_points;
                if (additional_points > 0) {
                    std::cout << "Inserting additional points..." << std::endl;
                    random_generator_p.insert_uniformly_distributed_points(additional_points, radius);
                    routing_scenario_p.discover_components();
                }
                std::cout << "Computing visibility graph..." << std::endl;
                timer.start();
                routing_scenario_p.use_triangulation_as_visibility_graph();
                timer.stop();
            }
            std::cout << "Computing visibility graph took: " << timer.time() << " seconds. Number of edges: "
                    << routing_scenario_p.edges_visibility_graph() << " Number of nodes: "
                    << routing_scenario_p.number_of_vertices() << std::endl;
            timer.reset();
            std::cout << "Compute A* or Dijkstra (0/1)? ";
            std::cin >> a_star_or_dijkstra;
            if (a_star_or_dijkstra) {
                std::cout << "Computing Dijkstra..." << std::endl;
                routing_scenario_p.dijkstra();
                routing_scenario_p.get_path_from_dijkstra();
                std::cout << "Average path length Dijkstra: " << routing_scenario_p.average_path_length_dijkstra() <<
                        std::endl;
            } else {
                std::cout << "Computing A*..." << std::endl;
                routing_scenario_p.a_star();
            }

            std::cout << "Optimize path (0/1)?  ";
            std::cin >> optimize_path;
            if (optimize_path) {
                routing_scenario_p.path_optimization();
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
    exit(0);
*/

    //polygonal domain generation test
    //poincare
    //std::cout << "Poincare" << std::endl;
    /*std::vector<int> trials = {100, 100, 100, 100, 100, 100, 100, 100, 100, 50, 1};
    int number_of_points = 0;
    CGAL::Timer timer;
    std::list<int> tests = {5000, 10000, 50000, 100000, 250000, 500000, 1000000,
        5000000, 10000000, 20000000, 25000000};
    int test_counter = 0;*/
    /*for (int test : tests) {
        std::cout << "test: " << test << std::endl;
        const double radius = std::acosh(test / (2 * CGAL_PI * 100) + 1);
        std::cout << "radius: " << radius << std::endl;
        timer.start();
        for (int i = 0; i < trials[test_counter]; ++i) {
            random_generator_p.generate_random_domain(test, radius, 0.55, 5, 5, 50, true, false, 0);
            number_of_points += routing_scenario_p.number_of_vertices();
        }
        timer.stop();
        const double time = timer.time();
        const double average_time = time / trials[test_counter];
        timer.reset();
        std::cout << "average time: " << average_time << std::endl;
        std::cout << "average number_of_points: " << number_of_points / trials[test_counter] << std::endl;
        number_of_points = 0;
        std::cout << std::endl;
        ++test_counter;
    }*/
    //beltrami
    /*test_counter = 0;
    std::cout << "BELTRAMI" << std::endl;
    for (int test : tests) {
        std::cout << "test: " << test << std::endl;
        const double radius = std::acosh(test / (2 * CGAL_PI * 100) + 1);
        std::cout << "radius: " << radius << std::endl;
        timer.start();
        for (int i = 0; i < trials[test_counter]; ++i) {
            if(i%10==0) {
                std::cout << "fffff" << std::endl;
            }
            random_generator_b.generate_random_domain(test, radius, 0.55, 5, 5, 50, true, false, 0);
            number_of_points += routing_scenario_b.number_of_vertices();
        }
        timer.stop();
        const double time = timer.time();
        const double average_time = time / trials[test_counter];
        timer.reset();
        std::cout << "average time: " << average_time << std::endl;
        std::cout << "average number_of_points: " << number_of_points / trials[test_counter] << std::endl;
        number_of_points = 0;
        std::cout << std::endl;
        ++test_counter;
    }*/

    //exact_predicates_exact_constructions_test
    /*const int trials = 100;
    std::list<int> tests = {5000, 10000, 50000, 100000, 250000, 500000};
    std::list <double> radii = {2, 3, 4, 5, 6, 7, 8, 9};
    for(int test : tests) {
        random_generator.generate_random_domain(test, radius, 0.55, 5, 5, 50, true, false, 0);
    }*/

    //

}
