#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include "include/internal/BeltramiKleinTraits.h"
#include "include/internal/PoincareDiskTraits.h"
#include "include/internal/RandomDomainGenerator.h"
#include "include/internal/Qt/Input/PointInput.h"
#include "include/internal/RoutingScenario.h"

typedef CGAL::Beltrami_klein_traits<CGAL::Exact_predicates_inexact_constructions_kernel> Beltrami_klein_traits;
typedef CGAL::Poincare_disk_traits<CGAL::Exact_predicates_inexact_constructions_kernel> Poincare_disk_traits;

//set either to Poincare_disk_traits or Beltrami_klein_traits
typedef Poincare_disk_traits K;

typedef CGAL::Delaunay_mesh_face_base_2<K> Face_base;
typedef CGAL::Triangulation_vertex_base_2<K> Vertex_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base> TDS;
typedef CGAL::No_constraint_intersection_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> T;
typedef CGAL::Qt::Routing_scenario<T> RoutingScenario;

typedef K::Point_2 Point_2;
typedef T::Vertex_handle Vertex_handle;

auto routing_scenario = RoutingScenario();
auto random_generator = CGAL::Qt::Random_domain_generator<T>(&routing_scenario);

void benchmark_random_generation() {
    const std::vector<double> r_h_list = {5, 6, 7, 8, 9, 10, 11, 12};

    const std::vector<int> trials = {100, 100, 100, 100, 50, 25, 10, 10};

    CGAL::Timer timer;
    int r_h_pointer = 0;
    for (const double r_h: r_h_list) {
        std::cout << "------------r_h: " << r_h << std::endl;
        int n = (int) (100 * 2 * CGAL_PI * (std::cosh(r_h) - 1));
        std::cout << "n: " << n << std::endl;
        double sum_time = 0;
        int n_sum = 0;

        for (int i = 0; i < trials[r_h_pointer]; ++i) {
            timer.start();
            random_generator.generate_random_domain(n, r_h, 0.55, 5, 5, 10, true, false, 0);
            timer.stop();
            n_sum += routing_scenario.number_of_vertices();
            sum_time += timer.time();
            timer.reset();
        }

        std::cout << "average time: " << sum_time / trials[r_h_pointer] << std::endl;
        std::cout << "average n: " << n_sum / trials[r_h_pointer] << std::endl;
        ++r_h_pointer;
    }
}

void load_polygonal_domain(std::string file_name) {
    std::ifstream ifs("../resources/polygonal_domains/" + file_name);
    if (ifs.fail()) {
        std::cout << "Failed to open file." << std::endl;
        return;
    }
    routing_scenario.clear();

    std::vector<Point_2> obstacle;
    int length;
    std::string model;
    ifs >> model;
    bool b;
    int projection = -1;
    if (typeid(K) == typeid(Beltrami_klein_traits) && model == "P") {
        std::cout << "Translating points from Poincare disk model to Beltrami-Klein model" << std::endl;
        projection = 1;
    }
    if (typeid(K) == typeid(Poincare_disk_traits) && model == "B") {
        std::cout << "Translating points from Beltrami-Klein model to Poincare disk model" << std::endl;
        projection = 0;
    }

    while (ifs >> length) {
        for (int i = 0; i < length; ++i) {
            double x;
            ifs >> x;
            double y;
            ifs >> y;

            if (projection == -1) {
                obstacle.push_back(Point_2(x, y));
            } else if (projection == 0) {
                obstacle.push_back(routing_scenario.beltrami_klein_to_poincare(Point_2(x, y)));
            } else {
                obstacle.push_back(routing_scenario.poincare_to_beltrami_klein(Point_2(x, y)));
            }
        }
        routing_scenario.insert_obstacle(obstacle.begin(), obstacle.end(), true);
        obstacle.clear();
    }
    routing_scenario.discover_components();

    std::cout << "Loading successfully" << std::endl;
    std::cout << "number of vertices: " << routing_scenario.number_of_vertices() << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
}

void save_polygonal_domain() {
    std::string file_name;
    std::cout << "Enter file name: ";
    std::cin >> file_name;
    if (!file_name.empty()) {
        std::ofstream ofs("../resources/polygonal_domains/" + file_name);
        if (ofs.fail()) {
            std::cout << "Failed to find directory." << std::endl;
            return;
        }
        if (typeid(K) == typeid(Beltrami_klein_traits)) {
            ofs << "B" << std::endl;
        } else {
            ofs << "P" << std::endl;
        }

        CGAL::Timer timer;
        timer.start();

        std::list<std::vector<Vertex_handle> > obstacles = routing_scenario.get_obstacles();
        for (std::vector<Vertex_handle> obstacle: obstacles) {
            const int length = obstacle.size();
            ofs << length << std::endl;
            for (const Vertex_handle vh: obstacle) {
                ofs << std::fixed << std::setprecision(14) << vh->point().x() << " " << vh->point().y() << std::endl;
            }
        }
        std::cout << "saving successfully" << std::endl;
    }
}

void generate_polygonal_domain() {
    //domain gen. parameters
    int n;
    int m;
    double r_h;
    double t;
    double pd;

    bool save;
    std::cout << "----random generation----" << std::endl;
    std::cout << "insert random gen. parameters:" << std::endl;

    std::cout << "r_h: ";
    std::cin >> r_h;
    std::cout << "pd: ";
    std::cin >> pd;
    if (pd != 0) {
        n = (int) (pd * 2 * CGAL_PI * (std::cosh(r_h) - 1));
        std::cout << "n: " << n << std::endl;
    } else {
        std::cout << "n: ";
        std::cin >> n;
    }
    std::cout << "m: ";
    std::cin >> m;
    std::cout << "t: ";
    std::cin >> t;

    random_generator.generate_random_domain(n, r_h, t, 5, 5, m, true, false, 0);
    std::cout << "number of vertices: " << routing_scenario.number_of_vertices() << std::endl;

    std::cout << "Save (0/1) ?";
    std::cin >> save;
    if (save) {
        save_polygonal_domain();
    }
}

void benchmark_routing_on_triangulation() {
    bool end;
    do {
        int trials;
        std::cout << "trials: ";
        std::cin >> trials;

        CGAL::Random random;
        std::vector<std::pair<Vertex_handle, Vertex_handle> > queries;
        std::vector<double> path_lengths;
        //std::vector<std::vector<int> > shortest_paths;
        for (int i = 0; i < trials; ++i) {
            const int start_index = random.get_int(0, routing_scenario.number_of_vertices() - 1);
            int dest_index = start_index;
            while (dest_index == start_index) {
                dest_index = random.get_int(0, routing_scenario.number_of_vertices() - 1);
            }
            Vertex_handle start = routing_scenario.index_vertex_map[start_index];
            Vertex_handle dest = routing_scenario.index_vertex_map[dest_index];
            queries.push_back(std::make_pair(start, dest));
        }

        std::cout << std::endl;
        std::cout << "--routing on visibility graph--" << std::endl;
        std::cout << std::endl;
        CGAL::Timer timer;
        timer.start();
        routing_scenario.build_visibility_graph();
        timer.stop();
        const double build_v_time = timer.time();
        std::cout << "building visbility graph took: " << build_v_time << std::endl;
        std::cout << "number of edges visbility graph: " << routing_scenario.edges_visibility_graph() << std::endl;
        timer.reset();
        double a_star_sum_time = 0;
        double dijkstra_sum_time = 0;
        double a_star_sum_path_length = 0, dijkstra_sum_length = 0;
        double reachable_counter = 0;
        for (std::pair<Vertex_handle, Vertex_handle> query: queries) {
            routing_scenario.set_point_to_start(query.first);
            routing_scenario.set_point_to_destination(query.second);

            timer.start();
            const bool reachable = routing_scenario.a_star();
            timer.stop();

            if (reachable) {
                //shortest_paths.push_back(routing_scenario.get_indices_path());
                double path_length = routing_scenario.get_path_length(routing_scenario.get_indices_path());
                a_star_sum_time += timer.time();
                a_star_sum_path_length += path_length;
                ++reachable_counter;
                path_lengths.push_back(path_length);
            } else {
                //shortest_paths.push_back({});
                path_lengths.push_back(DBL_MAX);
            }
            timer.reset();

            timer.start();
            routing_scenario.dijkstra();
            timer.stop();

            dijkstra_sum_time += timer.time();
            dijkstra_sum_length += routing_scenario.average_path_length_dijkstra();

            timer.reset();
        }
        //std::cout << "reachable paths: " << reachable_counter << std::endl;
        std::cout << std::endl;
        std::cout << "average A* time: " << a_star_sum_time / reachable_counter << std::endl;
        //std::cout << "average A* path length: " << a_star_sum_path_length / reachable_counter << std::endl;
        //std::cout << "one-to-one total time: " << build_v_time + a_star_sum_time / reachable_counter << std::endl;
        std::cout << std::endl;
        std::cout << "average Dijkstra time: " << dijkstra_sum_time / trials << std::endl;
        //std::cout << "average Dijkstra path length: " << dijkstra_sum_length / trials << std::endl;
        //std::cout << "one-to-all total time: " << build_v_time + dijkstra_sum_time / trials << std::endl;
        std::cout << std::endl;

        bool end1;
        do {
            //routing on triangulation
            bool b1;
            std::cout << "Insert extra points (0/1) ? ";
            std::cin >> b1;
            if (b1) {
                int amount;
                //bool blue_noise;
                int candidates = 0;
                double radius;
                const int old_number_of_vertices = routing_scenario.number_of_vertices();
                std::cout << "number of extra points: ";
                std::cin >> amount;
                /*std::cout << "blue noise sampling (0/1) ?";
                std::cin >> blue_noise;
                if(blue_noise) {
                    std::cout << "number of candidates: ";
                    std::cin >> candidates;
                }*/
                std::cout << "sample radius: ";
                std::cin >> radius;

                timer.start();
                random_generator.insert_uniformly_distributed_points(amount, radius, false, 0);
                timer.stop();
                std::cout << "inserting extra points took: " << timer.time() << std::endl;
                std::cout << "points inserted into free space: " << routing_scenario.number_of_vertices() -
                        old_number_of_vertices << std::endl;
                timer.reset();
            }

            bool b2;
            int opti = -1;
            std::cout << "Use path optimization (0/1) ? ";
            std::cin >> b2;
            if (b2) {
                std::cout << "Use greedy optimization   (0)" << std::endl;
                std::cout << "Use full optimization     (1)" << std::endl;
                std::cout << "Which optimizazion: ";
                std::cin >> opti;
            }

            std::cout << std::endl;
            std::cout << "--routing on subgraph--" << std::endl;
            std::cout << std::endl;
            routing_scenario.discover_components();
            timer.start();
            routing_scenario.use_triangulation_as_visibility_graph();
            timer.stop();
            const double build_s_time = timer.time();
            std::cout << "building subgraph took: " << build_s_time << std::endl;
            std::cout << "number of edges subgraph: " << routing_scenario.edges_visibility_graph() << std::endl;
            timer.reset();
            a_star_sum_time = 0;
            dijkstra_sum_time = 0;
            double opti_time = 0;
            double a_star_approx_sum_path_length = 0, dijkstra_approx_sum_path_length = 0;
            std::vector<double> approx_path_lengths;
            //std::vector<std::vector<int> > approx_paths;
            reachable_counter = 0;
            for (std::pair<Vertex_handle, Vertex_handle> query: queries) {
                routing_scenario.set_point_to_start(query.first);
                routing_scenario.set_point_to_destination(query.second);

                timer.start();
                const bool reachable = routing_scenario.a_star();
                timer.stop();

                std::vector<int> path = routing_scenario.get_indices_path();

                if (reachable) {
                    ++reachable_counter;
                    if (opti != -1) {
                        CGAL::Timer opti_timer;
                        opti_timer.start();
                        if (opti == 0) {
                            path = routing_scenario.greedy_optimization(path);
                        } else if (opti == 1) {
                            path = routing_scenario.path_optimization(path);
                        }
                        opti_timer.stop();
                        opti_time += opti_timer.time();
                    }
                    //approx_paths.push_back(routing_scenario.get_indices_path());
                    double path_length = routing_scenario.get_path_length(path);
                    a_star_sum_time += timer.time();
                    a_star_approx_sum_path_length += path_length;
                    approx_path_lengths.push_back(path_length);
                } else {
                    //approx_paths.push_back({});
                    approx_path_lengths.push_back(DBL_MAX);
                }
                timer.reset();

                timer.start();
                routing_scenario.dijkstra();
                timer.stop();

                dijkstra_sum_time += timer.time();
                dijkstra_approx_sum_path_length += routing_scenario.average_path_length_dijkstra();

                timer.reset();
            }

            double min = DBL_MAX, max = 0;
            for (int i = 0; i < trials; ++i) {
                if (path_lengths[i] != DBL_MAX) {
                    const double ratio = approx_path_lengths[i] / path_lengths[i];
                    if (ratio > max) {
                        max = ratio;
                    }
                    if (ratio < min) {
                        min = ratio;
                    }
                }
            }
            if (opti != -1) {
                std::cout << "average path optimization took: " << opti_time / reachable_counter << std::endl;
            }
            std::cout << std::endl;
            //std::cout << "reachable paths: " << reachable_counter << std::endl;
            std::cout << "average A* time: " << a_star_sum_time / reachable_counter << std::endl;
            std::cout << "average Dijkstra time: " << dijkstra_sum_time / trials << std::endl;
            //std::cout << "average Dijkstra path length: " << dijkstra_approx_sum_path_length / trials << std::endl;
            //std::cout << "-> Dijkstra quality of approx. paths is: " << dijkstra_approx_sum_path_length / dijkstra_sum_length << std::endl;
            //std::cout << "average A* path length: " << a_star_approx_sum_path_length / reachable_counter << std::endl;
            std::cout << "-> A* quality of approx. paths is: " << a_star_approx_sum_path_length / a_star_sum_path_length
                    << std::endl;
            std::cout << "min ratio: " << min << std::endl;
            std::cout << "max ratio: " << max << std::endl;
            //std::cout << "one-to-one total time: " << build_s_time + a_star_sum_time / reachable_counter << std::endl;
            std::cout << std::endl;
            //std::cout << "one-to-all total time: " << build_s_time + dijkstra_sum_time / trials << std::endl;
            std::cout << std::endl;

            if (b1) {
                routing_scenario.remove_all_unconstrained_points();
                std::cout << "removed unconstrained points" << std::endl;
            }

            std::cout << "Again (0/1)? ";
            std::cin >> end1;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
        } while (end1);
        std::cout << "Again (0/1)? ";
        std::cin >> end;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    } while (end);
}

void big_benchmark_routing_on_triangulation() {
    const std::vector<std::string> domains = {
        /*"10_100_10_55", "11_100_10_55", "12_100_10_55",
        "10_100_10_60", "11_100_10_60", "12_100_10_60",
        "10_100_10_62", "11_100_10_62", "12_100_10_62"*/
        "10_100_10_65", "11_100_10_65"/*, "12_100_10_65"*/
    };

    const std::vector<int> trials = {
        /*1000, 500, 250,
        1000, 500, 250,*/
        1000, 500/*, 250*/
    };

    for (int i = 0; i < domains.size(); ++i) {
        load_polygonal_domain(domains[i]);
        std::vector<std::pair<Vertex_handle, Vertex_handle> > queries;
        std::vector<double> path_lengths;
        //loading query file
        std::ifstream ifs("../resources/polygonal_domains/query_" + domains[i]);
        if (ifs.fail()) {
            std::cout << "Failed to open file." << std::endl;
            return;
        }
        for (int m = 0; m < trials[i]; ++m) {
            int start;
            ifs >> start;
            int dest;
            ifs >> dest;
            queries.push_back(std::make_pair(routing_scenario.index_vertex_map[start],
                                             routing_scenario.index_vertex_map[dest]));
        }
        std::cout << std::endl;
        std::cout << "--routing on visibility graph--" << std::endl;
        std::cout << std::endl;
        CGAL::Timer timer;
        timer.start();
        routing_scenario.build_visibility_graph();
        timer.stop();
        const double build_v_time = timer.time();
        std::cout << "building visbility graph took: " << build_v_time << std::endl;
        std::cout << "number of edges visbility graph: " << routing_scenario.edges_visibility_graph() << std::endl;
        timer.reset();
        double a_star_sum_time = 0;
        double dijkstra_sum_time = 0;
        double a_star_sum_path_length = 0, dijkstra_sum_length = 0;
        double reachable_counter = 0;
        for (std::pair<Vertex_handle, Vertex_handle> query: queries) {
            routing_scenario.set_point_to_start(query.first);
            routing_scenario.set_point_to_destination(query.second);

            timer.start();
            const bool reachable = routing_scenario.a_star();
            timer.stop();

            if (reachable) {
                //shortest_paths.push_back(routing_scenario.get_indices_path());
                double path_length = routing_scenario.get_path_length(routing_scenario.get_indices_path());
                a_star_sum_time += timer.time();
                a_star_sum_path_length += path_length;
                ++reachable_counter;
                path_lengths.push_back(path_length);
            } else {
                //shortest_paths.push_back({});
                path_lengths.push_back(DBL_MAX);
            }
            timer.reset();

            timer.start();
            routing_scenario.dijkstra();
            timer.stop();

            dijkstra_sum_time += timer.time();
            dijkstra_sum_length += routing_scenario.average_path_length_dijkstra();

            timer.reset();
        }
        std::cout << std::endl;
        std::cout << "average A* time: " << a_star_sum_time / reachable_counter << std::endl;
        std::cout << std::endl;
        std::cout << "average Dijkstra time: " << dijkstra_sum_time / trials[i] << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "--routing on subgraph--" << std::endl;
        std::cout << std::endl;
        timer.start();
        routing_scenario.use_triangulation_as_visibility_graph();
        timer.stop();
        const double build_s_time = timer.time();
        std::cout << "building subgraph took: " << build_s_time << std::endl;
        std::cout << "number of edges subgraph: " << routing_scenario.edges_visibility_graph() << std::endl;
        timer.reset();

        a_star_sum_time = 0;
        dijkstra_sum_time = 0;
        double opti_time1 = 0, opti_time2 = 0;
        double a_star_approx_sum_path_length = 0, dijkstra_approx_sum_path_length = 0;
        double approx_path_length1 = 0, approx_path_length2 = 0;
        std::vector<double> approx_path_lengths;
        std::vector<double> optimized_path_lengths1;
        std::vector<double> optimized_path_lengths2;
        //std::vector<std::vector<int> > approx_paths;
        reachable_counter = 0;
        for (std::pair<Vertex_handle, Vertex_handle> query: queries) {
            routing_scenario.set_point_to_start(query.first);
            routing_scenario.set_point_to_destination(query.second);

            timer.start();
            const bool reachable = routing_scenario.a_star();
            timer.stop();

            if (reachable) {
                ++reachable_counter;
                //approx_paths.push_back(routing_scenario.get_indices_path());
                std::vector<int> path = routing_scenario.get_indices_path();

                double path_length = routing_scenario.get_path_length(path);
                a_star_sum_time += timer.time();
                a_star_approx_sum_path_length += path_length;
                approx_path_lengths.push_back(path_length);

                timer.reset();

                timer.start();
                std::vector<int> optimized_path1 = routing_scenario.greedy_optimization(path);
                timer.stop();
                opti_time1 += timer.time();
                timer.reset();

                timer.start();
                std::vector<int> optimized_path2 = routing_scenario.path_optimization(path);
                timer.stop();
                opti_time2 += timer.time();
                timer.reset();

                double path_length1 = routing_scenario.get_path_length(optimized_path1);
                approx_path_length1 += path_length1;
                optimized_path_lengths1.push_back(path_length1);

                double path_length2 = routing_scenario.get_path_length(optimized_path2);
                approx_path_length2 += path_length2;
                optimized_path_lengths2.push_back(path_length2);
            } else {
                //approx_paths.push_back({});
                approx_path_lengths.push_back(DBL_MAX);
                optimized_path_lengths1.push_back(DBL_MAX);
                optimized_path_lengths2.push_back(DBL_MAX);
            }
            timer.reset();

            timer.start();
            routing_scenario.dijkstra();
            timer.stop();

            dijkstra_sum_time += timer.time();
            dijkstra_approx_sum_path_length += routing_scenario.average_path_length_dijkstra();

            timer.reset();
        }

        double min = DBL_MAX, max = 0;
        for (int j = 0; j < trials[i]; ++j) {
            if (path_lengths[j] != DBL_MAX) {
                const double ratio = approx_path_lengths[j] / path_lengths[j];
                if (ratio > max) {
                    max = ratio;
                }
                if (ratio < min) {
                    min = ratio;
                }
            }
        }

        double min1 = DBL_MAX, max1 = 0;
        for (int j = 0; j < trials[i]; ++j) {
            if (path_lengths[j] != DBL_MAX) {
                const double ratio = optimized_path_lengths1[j] / path_lengths[j];
                if (ratio > max1) {
                    max1 = ratio;
                }
                if (ratio < min1) {
                    min1 = ratio;
                }
            }
        }

        double min2 = DBL_MAX, max2 = 0;
        for (int j = 0; j < trials[i]; ++j) {
            if (path_lengths[j] != DBL_MAX) {
                const double ratio = optimized_path_lengths2[j] / path_lengths[j];
                if (ratio > max2) {
                    max2 = ratio;
                }
                if (ratio < min2) {
                    min2 = ratio;
                }
            }
        }
        std::cout << std::endl;
        std::cout << "A* time: " << a_star_sum_time / reachable_counter << std::endl;
        std::cout << "Dijkstra time: " << dijkstra_sum_time / trials[i] << std::endl;
        std::cout << std::endl;
        std::cout << "avg: " << a_star_approx_sum_path_length / a_star_sum_path_length
                << std::endl;
        std::cout << "min: " << min << std::endl;
        std::cout << "max: " << max << std::endl;
        std::cout << std::endl;
        std::cout << "PATH OPTIMIZATION" << std::endl;
        std::cout << "GREEDY" << std::endl;
        std::cout << "avg: " << approx_path_length1 / a_star_sum_path_length << std::endl;
        std::cout << "min: " << min1 << std::endl;
        std::cout << "max: " << max1 << std::endl;
        std::cout << "time: " << opti_time1 / reachable_counter << std::endl;
        std::cout << std::endl;
        std::cout << "FULL" << std::endl;
        std::cout << "avg: " << approx_path_length2 / a_star_sum_path_length << std::endl;
        std::cout << "min: " << min2 << std::endl;
        std::cout << "max: " << max2 << std::endl;
        std::cout << "time: " << opti_time2 / reachable_counter << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }
}

void benchmark_extra_points() {
    const std::vector<std::string> domains = {
        /*"10_100_10_55", "11_100_10_55", "12_100_10_55",
        "10_100_10_60", "11_100_10_60", "12_100_10_60",
        "10_100_10_62", "11_100_10_62", "12_100_10_62"*/
        "10_100_10_65", "11_100_10_65" /*, "12_100_10_65"*/
    };

    const std::vector<int> trials = {
        /*1000, 500, 250,
        1000, 500, 250,*/
        /*1000, 500,*/ 250
    };

    const std::vector<double> radius = {
        /*10, 11, 12,
        10, 11, 12,*/
        /*10, 11,*/ 12
    };

    const std::vector<std::vector<double> > extra_points_factor = {
        /*{1.0 / 64, 1.0 / 32, 1.0 / 16, 1.0 / 8}, {1.0 / 64, 1.0 / 32, 1.0 / 16, 1.0 / 8},
        {1.0 / 64, 1.0 / 32, 1.0 / 16, 1.0 / 8},
        {1.0 / 64, 1.0 / 32, 1.0 / 16, 1.0 / 8}, {1.0 / 64, 1.0 / 32, 1.0 / 16, 1.0 / 8},
        {1.0 / 64, 1.0 / 32, 1.0 / 16, 1.0 / 8},*/
        /*{1.0 / 64, 1.0 / 32, 1.0 / 16, 1.0 / 8}, {1.0 / 64, 1.0 / 32, 1.0 / 16, 1.0 / 8},*/
        {1.0 / 64, 1.0 / 32, 1.0 / 16, 1.0 / 8}
    };

    CGAL::Timer timer;

    for (int i = 0; i < domains.size(); ++i) {
        std::cout << "DOMAIN NUMBER: " << i << std::endl;
        load_polygonal_domain(domains[i]);
        std::vector<std::pair<Vertex_handle, Vertex_handle> > queries;
        std::vector<double> path_lengths;
        int reachable_counter = 0;
        double a_star_sum_path_length = 0;
        //loading query file
        std::ifstream ifs("../resources/polygonal_domains/query_" + domains[i]);
        if (ifs.fail()) {
            std::cout << "Failed to open file." << std::endl;
            return;
        }
        for (int m = 0; m < trials[i]; ++m) {
            int start;
            ifs >> start;
            int dest;
            ifs >> dest;
            queries.push_back(std::make_pair(routing_scenario.index_vertex_map[start],
                                             routing_scenario.index_vertex_map[dest]));
        }

        routing_scenario.build_visibility_graph();
        for (std::pair<Vertex_handle, Vertex_handle> query: queries) {
            routing_scenario.set_point_to_start(query.first);
            routing_scenario.set_point_to_destination(query.second);

            const bool reachable = routing_scenario.a_star();

            if (reachable) {
                double path_length = routing_scenario.get_path_length(routing_scenario.get_indices_path());
                ++reachable_counter;
                path_lengths.push_back(path_length);
                a_star_sum_path_length += path_length;
            } else {
                path_lengths.push_back(DBL_MAX);
            }
        }

        for (double factor: extra_points_factor[i]) {
            std::cout << "-------------------------factor: " << factor << std::endl;
            timer.start();
            const int amount = (int) 2.0 * CGAL_PI * (std::cosh(radius[i]) - 1) * factor * 100.0;
            random_generator.insert_uniformly_distributed_points(amount, radius[i], false, 0);
            timer.stop();
            std::cout << "amount: " << amount << std::endl;
            std::cout << "inserting took: " << timer.time() << std::endl;
            timer.reset();
            std::cout << std::endl;
            routing_scenario.discover_components();
            timer.start();
            routing_scenario.use_triangulation_as_visibility_graph();
            timer.stop();
            const double build_s_time = timer.time();
            std::cout << "building subgraph took: " << build_s_time << std::endl;
            std::cout << "number of edges subgraph: " << routing_scenario.edges_visibility_graph() << std::endl;
            std::cout << "average vertex degree: " << routing_scenario.average_vertex_degree_visibility_graph() <<
                    std::endl;
            timer.reset();
            double a_star_sum_time = 0;
            double opti_time = 0;
            double a_star_approx_sum_path_length = 0;
            std::vector<double> approx_path_lengths;
            //std::vector<std::vector<int> > approx_paths;
            reachable_counter = 0;
            for (std::pair<Vertex_handle, Vertex_handle> query: queries) {
                routing_scenario.set_point_to_start(query.first);
                routing_scenario.set_point_to_destination(query.second);

                timer.start();
                const bool reachable = routing_scenario.a_star();
                timer.stop();
                a_star_sum_time += timer.time();

                std::vector<int> path = routing_scenario.get_indices_path();

                if (reachable) {
                    ++reachable_counter;
                    CGAL::Timer opti_timer;
                    opti_timer.start();
                    path = routing_scenario.path_optimization(path);
                    opti_timer.stop();
                    opti_time += opti_timer.time();
                    double path_length = routing_scenario.get_path_length(path);
                    a_star_approx_sum_path_length += path_length;
                    approx_path_lengths.push_back(path_length);
                } else {
                    approx_path_lengths.push_back(DBL_MAX);
                }
                timer.reset();
            }

            double min = DBL_MAX, max = 0;
            for (int j = 0; j < trials[i]; ++j) {
                if (path_lengths[j] != DBL_MAX) {
                    const double ratio = approx_path_lengths[j] / path_lengths[j];
                    if (ratio > max) {
                        max = ratio;
                    }
                    if (ratio < min) {
                        min = ratio;
                    }
                }
            }

            std::cout << std::endl;
            std::cout << "average A* time: " << a_star_sum_time / reachable_counter << std::endl;
            std::cout << "full optimization took: " << opti_time / reachable_counter << std::endl;
            std::cout << "avg: " << a_star_approx_sum_path_length / a_star_sum_path_length
                    << std::endl;
            std::cout << "min: " << min << std::endl;
            std::cout << "max: " << max << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            routing_scenario.remove_all_unconstrained_points();
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }
}

void benchmark_tea() {
    int trials;
    std::cout << "trials: ";
    std::cin >> trials;

    CGAL::Timer sum;
    sum.start();
    //double number_o_tests = 0;
    //double sum_times_orientation_tests = 0;
    //double min_o_time = DBL_MAX, max_o_time = 0;
    double distance_comp_time = 0;
    double min = DBL_MAX, max = 0;
    for (int i = 0; i < trials; ++i) {
        CGAL::Timer timer;
        timer.start();
        const std::pair<double, double> pair = routing_scenario.build_visibility_graph();
        timer.stop();

        //sum_times_orientation_tests +=pair.first;
        //number_o_tests += pair.second;

        distance_comp_time += pair.second;

        /*if(min_o_time > pair.first) {
            min_o_time = pair.first;
        }

        if(max_o_time < pair.first) {
            max_o_time = pair.first;
        }*/

        if (timer.time() < min) {
            min = timer.time();
        }
        if (timer.time() > max) {
            max = timer.time();
        }
    }
    sum.stop();
    std::cout << "edges visibility graph: " << routing_scenario.edges_visibility_graph() << std::endl;
    std::cout << "average vertex degree: " << routing_scenario.average_vertex_degree_visibility_graph() << std::endl;
    std::cout << "average time: " << sum.time() / trials << std::endl;
    std::cout << "max. time: " << max << std::endl;
    std::cout << "min. time: " << min << std::endl;
    std::cout << "average distance comp time: " << distance_comp_time / trials << std::endl;
    //std::cout << "average number of o-tests: " << number_o_tests / trials << std::endl;
    //std::cout << "average time o-test: " << sum_times_orientation_tests / trials << std::endl;
    //std::cout << "max. time o-test: " << max_o_time << std::endl;
    //std::cout << "min. time o-test: " << min_o_time << std::endl;

    bool b1;
    std::cout << "Again (0/1) ? ";
    std::cin >> b1;
    if (b1) {
        benchmark_tea();
    }
}

void generate_query_list() {
    int size;
    std::cout << "query list size: ";
    std::cin >> size;

    CGAL::Random random;
    std::vector<std::pair<int, int> > queries;
    for (int i = 0; i < size; ++i) {
        const int start_index = random.get_int(0, routing_scenario.number_of_vertices() - 1);
        int dest_index = start_index;
        while (dest_index == start_index) {
            dest_index = random.get_int(0, routing_scenario.number_of_vertices() - 1);
        }
        queries.push_back(std::make_pair(start_index, dest_index));
    }

    std::string file_name;
    std::cout << "save under file name: ";
    std::cin >> file_name;

    if (!file_name.empty()) {
        std::ofstream ofs("../resources/polygonal_domains/" + file_name);
        if (ofs.fail()) {
            std::cout << "Failed to find directory." << std::endl;
            return;
        }

        std::list<std::vector<Vertex_handle> > obstacles = routing_scenario.get_obstacles();
        for (std::pair<int, int> query: queries) {
            ofs << query.first << " " << query.second << std::endl;
        }
        std::cout << "saving successfully" << std::endl;
    }
}

void benchmark_on_domain() {
    int action;
    std::cout << "Benchmark TEA                         (0)" << std::endl;
    std::cout << "Benchmark routing on triangulation    (1)" << std::endl;
    std::cout << "Generate point-to-point query list    (2)" << std::endl;
    std::cout << "Leave                                 (3)" << std::endl;
    std::cout << "Select action: ";
    std::cin >> action;
    if (action == 0) {
        benchmark_tea();
    } else if (action == 1) {
        benchmark_routing_on_triangulation();
    } else if (action == 2) {
        generate_query_list();
    }
}

void big_benchmark_tea() {
    const std::list<double> t_list = {
        0.55, 0.60 , 0.62
    };
    const std::list<double> r_h_list = {10.5, 11.25, 11.5, 11.75};
    const std::vector<int> trial_list = {50, 25, 25, 25};

    CGAL::Timer timer;

    for (const double t: t_list) {
        int r_h_counter = 0;
        for (const double r_h: r_h_list) {
            const int n = (int) (100 * 2 * CGAL_PI * (std::cosh(r_h) - 1));

            double avg_N = 0;
            double avg_edges_V = 0;
            double min_time = DBL_MAX, max_time = 0, avg_time = 0;
            double avg_avg_vertex_degree_V = 0, min_avg_degree = DBL_MAX, max_avg_degree = 0;
            for (int i = 0; i < trial_list[r_h_counter]; ++i) {

                if(i%5==0) {
                    std::cout << "zzz" << std::endl;
                }

                random_generator.generate_random_domain(n, r_h, t, 5, 5, 10, true, false, 0);

                avg_N += routing_scenario.number_of_vertices();

                timer.start();
                routing_scenario.build_visibility_graph();
                timer.stop();

                const double time = timer.time();
                avg_time += time;
                if (time > max_time) {
                    max_time = time;
                }
                if (time < min_time) {
                    min_time = time;
                }
                timer.reset();

                const double avg_degree = routing_scenario.average_vertex_degree_visibility_graph();
                avg_avg_vertex_degree_V += avg_degree;
                if (avg_degree < min_avg_degree) {
                    min_avg_degree = avg_degree;
                }
                if (avg_degree > max_avg_degree) {
                    max_avg_degree = avg_degree;
                }

                avg_edges_V += routing_scenario.edges_visibility_graph();
            }
            std::cout << "new row" << std::endl;
            std::cout << "t: " << t << " r_h: " << r_h << std::endl;
            std::cout << "avg_N: " << avg_N / trial_list[r_h_counter] << " avg_edges: "
                    << avg_edges_V / trial_list[r_h_counter] << " avg_degree: "
                    << avg_avg_vertex_degree_V / trial_list[r_h_counter] << std::endl;
            std::cout << "min_avg_degree: " << min_avg_degree << " max_avg_degree: " << max_avg_degree << std::endl;
            std::cout << "min_time: " << min_time << " max_time: "
                    << max_time << " avg_time: "
                    << avg_time / trial_list[r_h_counter] << std::endl;
            ++r_h_counter;
        }
    }
}

int main() {
    bool again;
    int action;
    do {
        std::cout << "Benchmark polygonal domain gen.   (0)" << std::endl;
        std::cout << "Load polygonal domain             (1)" << std::endl;
        std::cout << "Generate polygonal domain         (2)" << std::endl;
        std::cout << "Benchmark TEA                     (3)" << std::endl;
        std::cout << "Benchmark Routing on triang.      (4)" << std::endl;
        std::cout << "Benchmark extra points            (5)" << std::endl;
        std::cout << "Select action: ";
        std::cin >> action;
        if (action == 0) {
            benchmark_random_generation();
        } else if (action == 1) {
            std::cout << std::endl;
            std::string file_name;
            std::cout << "Enter file name: ";
            std::cin >> file_name;
            load_polygonal_domain(file_name);
            std::cout << std::endl;
            benchmark_on_domain();
        } else if (action == 2) {
            std::cout << std::endl;
            generate_polygonal_domain();
            std::cout << std::endl;
            benchmark_on_domain();
        } else if (action == 3) {
            std::cout << std::endl;
            big_benchmark_tea();
            std::cout << std::endl;
        } else if (action == 4) {
            std::cout << std::endl;
            big_benchmark_routing_on_triangulation();
            std::cout << std::endl;
        } else {
            std::cout << std::endl;
            benchmark_extra_points();
            std::cout << std::endl;
        }

        std::cout << "Again (0/1)? ";
        std::cin >> again;
        std::cout << std::endl;
    } while (again);
}
