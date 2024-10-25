#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include "include/internal/BeltramiKleinTraits.h"
#include "include/internal/PoincareDiskTraits.h"
#include "include/internal/RandomDomainGenerator.h"
#include "include/internal/Qt/Input/PointInput.h"
#include "include/internal/RoutingScenario.h"


int main() {
    typedef CGAL::Beltrami_klein_traits<CGAL::Exact_predicates_inexact_constructions_kernel> K;
    typedef CGAL::Delaunay_mesh_face_base_2<K> Face_base;
    typedef CGAL::Triangulation_vertex_base_2<K> Vertex_base;
    typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base> TDS;
    typedef CGAL::No_constraint_intersection_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> T;
    typedef CGAL::Qt::Routing_scenario<T> RoutingScenario;
    typedef K::Point_2 Point_2;
    RoutingScenario routing_scenario = RoutingScenario();
    auto rg = CGAL::Qt::Random_domain_generator<T>(&routing_scenario);
    rg.generate_random_domain(20000000, 0.999, 0.01, 5, 5, 5, true, false);
    std::cout << "Number of vertices: " << routing_scenario.number_of_vertices() << std::endl;

    routing_scenario.set_start_point(CGAL::ORIGIN);
    routing_scenario.set_destination_point(Point_2(0, 0.99));

    CGAL::Timer timer;
    timer.start();
    routing_scenario.use_triangulation_as_visibility_graph();
    timer.stop();
    std::cout << "Building visibility graph took: " << timer.time() << " seconds." << std::endl;
    timer.reset();
    //routing_scenario.use_triangulation_as_visibility_graph();
    //routing_scenario.build_visibility_graph_naive();
    timer.start();
    routing_scenario.a_star();
    timer.stop();
    std::cout << "A* took: " << timer.time() << " seconds." << std::endl;
    timer.reset();

    std::cout << "FINISH" << std::endl;
}
