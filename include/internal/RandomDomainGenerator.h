#ifndef RANDOMWORLDGENERATOR_H
#define RANDOMWORLDGENERATOR_H

/*
 * generates random polygonal domains in the Poincare disk model, or
 * the Beltrami-Klein model. First we sample random points in the Poincare disk
 * model, therefore a Euclidean radius R in the Poincare disk is choosen, all sampled points lie within a disk
 * centered at 0 with radius R.
 * The angle of a random point is distributed uniformly and the radius is distributed with
 * sinh(x)/(cosh(R)-1), where x is the sampled hyperbolic radius of the point. To sample a
 * hyperbolic radius we take Euclidean er and compute the corresponding hyperbolic radius. With this
 * procedure and simplifying we arive at probability_density(er, area), where area is the area of the
 * hyperbolic circle centered at 0 and Euclidean radius R.
 * To get samples from probability_density(er, area), rejection and inverse-sampling can be used.
 */

#include <CGAL/point_generators_2.h>
#include "internal/RoutingScenario.h"
#include <random>

namespace CGAL::Qt {
    template<typename T>
    class Random_domain_generator {
        typedef typename T::Point_2 Point_2;
        typedef typename T::Geom_traits Geom_traits;
        typedef typename T::Finite_faces_iterator Finite_faces_iterator;
        typedef typename T::Finite_edges_iterator Finite_edges_iterator;
        typedef typename T::Finite_vertices_iterator Finite_vertices_iterator;
        typedef typename T::Constrained_edges_iterator Constrained_edges_iterator;
        typedef typename T::Face_circulator Face_circulator;
        typedef typename T::Edge_circulator Edge_circulator;
        typedef typename T::Face_handle Face_handle;
        typedef typename T::Vertex_handle Vertex_handle;
        typedef typename T::Vertex_circulator Vertex_circulator;
        typedef typename T::Edge Edge;

    public:
        explicit Random_domain_generator(Routing_scenario<T> *r);

        static double probability_density(const double er, const double area) {
            /*const double area = std::cosh(std::log((1 + radius) / (1 - radius))) - 1;*/
            //return 2 * er / (1 - er * er) / area;
            //return 2*er; //Euclidean
            return std::sinh(er) / area;
        }

        static double inverse_antiderivative(const double x, const double area) {
            //return std::sqrt(1 - std::exp(-x * area));
            return std::acosh(x * area + 1);
        }

        static double antiderivative(const double x, const double area) {
            //return -std::log(1 - x * x) / area;
            return (std::cosh(x) - 1) / area;
        }

        void generate_random_domain(const int number_of_points, const double radius,
                                    const double threshold, const int erosions, const int dilations, const int min,
                                    const bool erosion_before_dilation, const bool blue_noise, const int candidates) {
            //std::cout << std::endl;
            //std::cout << "-----Domain generation-----" << std::endl;
            Timer sum;
            sum.start();
            r->clear();

            T *t = r->t;

            Timer timer;
            timer.start();
            double sampling_time, insertion_time;
            if(blue_noise) {
                this->  blue_noise(candidates, number_of_points, radius);
                timer.stop();
                sampling_time = timer.time();
                insertion_time = timer.time();
                //std::cout << "Blue noise toook: " << timer.time() << " seconds." << std::endl;
                timer.reset();
            } else {
                std::vector<Point_2> points = inverse_sampling(number_of_points, radius);
                timer.stop();
                sampling_time = timer.time();
                //std::cout << "Sampling points took: " << timer.time() << " seconds." << std::endl;
                timer.reset();
                timer.start();
                r->insert_points(points.begin(), points.end());
                timer.stop();
                insertion_time = timer.time();
                //std::cout << "Inserting points took: " << timer.time() << " seconds." << std::endl;
                timer.reset();
            }

            timer.start();
            set_faces_in_domain_method1(t, threshold);
            timer.stop();
            const double boolean_time = timer.time();
            //std::cout << "Initial in_domain assignment took: " << timer.time() << " seconds." << std::endl;
            timer.reset();

            timer.start();
            make_smoother(t, erosion_before_dilation, erosions, dilations);
            timer.stop();
            const double smoothing_time = timer.time();
            //std::cout << "Erosion and dilation took: " << timer.time() << " seconds." << std::endl;
            timer.reset();

            timer.start();
            discover_edges(t);
            timer.stop();
            const double edges_time = timer.time();
            //std::cout << "Discovering and inserting constrained edges took: " << timer.time() << " seconds." <<
                    //std::endl;
            timer.reset();

            double remove_obstacle_time = 0;
            if (min > 3) {
                timer.start();
                remove_small_obstacles(t, min);
                timer.stop();
                remove_obstacle_time = timer.time();
                //std::cout << "Removing obstacles with less than min points took: " << timer.time() << " seconds." <<
                  //      std::endl;
                timer.reset();
            }

            timer.start();
            r->remove_all_unconstrained_points();
            timer.stop();
            const double remove_points_time = timer.time();
            //std::cout << "Removing points with no incident constrained edges took: " << timer.time() << " seconds." <<
                    //std::endl;
            timer.reset();

            timer.start();
            r->discover_components();
            timer.stop();
            const double discover_time = timer.time();
            //std::cout << "Discovering components took: " << timer.time() << " seconds." << std::endl;
            timer.reset();

            sum.stop();
            const double sum_time = sum.time();
            /*std::cout << "sampling points: " << sampling_time << " %: " << sampling_time / sum_time  << std::endl;
            std::cout << "inserting points: " << insertion_time << " %: " << insertion_time / sum_time  << std::endl;
            std::cout << "boolean assignment: " << boolean_time << " %: " << boolean_time / sum_time  << std::endl;
            std::cout << "smoothing took: " << smoothing_time << " %: " << smoothing_time / sum_time  << std::endl;
            std::cout << "discover edges: " << edges_time << " %: " << edges_time / sum_time  << std::endl;
            std::cout << "remove small obstacles: " << remove_obstacle_time << " %: " << remove_obstacle_time / sum_time  << std::endl;
            std::cout << "remove points: " << remove_points_time << " %: " << remove_points_time / sum_time  << std::endl;
            std::cout << "discover components: " << discover_time << " %: " << discover_time / sum_time  << std::endl;
            std::cout << "-----Generating domain took: " << sum.time() << " seconds.-----" << std::endl;
            std::cout << std::endl;*/
        }

        void insert_uniformly_distributed_points(const int n, const double radius, const bool b, const int k) {
            /*Timer timer;
            timer.start();*/
            if(!b) {
                std::vector<Point_2> points = inverse_sampling(n, radius);
                r->insert_points(points.begin(), points.end());
                //r->remove_unconstrained_points_in_obstacle_interior();
            } else {
                blue_noise(k, n, radius);
                r->remove_unconstrained_points_in_obstacle_interior();
            }
            /*timer.stop();
            std::cout << "Inserting additional points took: " << timer.time() << " seconds." << std::endl;
            timer.reset();*/
        }

        //remove polygonal chains with less than min vertices
        void remove_small_obstacles(T *t, const int min) {
            std::list<std::vector<Vertex_handle> > obstacles = r->get_obstacles();
            std::set<Edge> set;
            for (std::vector<Vertex_handle> obstacle: obstacles) {
                if (obstacle.size() < min) {
                    //remove edges of obstacle
                    for (int i = 0; i < obstacle.size(); ++i) {
                        Face_circulator fc = t->incident_faces(obstacle[i]), done(fc);
                        do {
                            int index = fc->index(obstacle[i]);
                            if (fc->vertex(fc->ccw(index)) == obstacle[(i + 1) % obstacle.size()]) {
                                r->remove_constrained_edge(fc, fc->cw(index));
                                break;
                            }
                        } while (++fc != done);
                    }
                }
            }
        }

        //give every vertex a value from [-1,1], then decide for each triangle the in_domain value
        void set_faces_in_domain_method1(T *t, const double threshold) {
            //std::map<Vertex_handle, float> map;
            std::vector<double> value;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> value_dist(0, 1);
            //Random random;

            for (int i = 0; i < r->number_of_vertices(); ++i) {
                value.push_back(value_dist(gen));
            }

            //smoothing values
            /*for (int i = 0; i < r->number_of_vertices(); ++i) {
                double sum = 0;
                double neighbors = 0;
                Vertex_handle vh = r->index_vertex_map[i];
                Vertex_circulator vc = r->t->incident_vertices(vh), done(vc);
                do {
                    if (!r->t->is_infinite(vc)) {
                        const double neighbor_value = value[r->vertex_index_map[vc]];
                        sum += neighbor_value;
                        ++neighbors;
                    }
                } while (++vc != done);
                value[i] = sum / neighbors;
            }*/

            for (Finite_faces_iterator fi = t->finite_faces_begin(); fi != t->finite_faces_end(); ++fi) {
                const double x = value[r->vertex_index_map[fi->vertex(0)]];
                const double y = value[r->vertex_index_map[fi->vertex(1)]];
                const double z = value[r->vertex_index_map[fi->vertex(2)]];
                if ((x + y + z) / 3 /*value_dist(gen)*/ > threshold) {
                    fi->set_in_domain(true);
                } else {
                    fi->set_in_domain(false);
                }
            }
        }

        /*void set_faces_in_domain_method2(T *t, const double threshold) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> value_dist(0, 1);

            for (Finite_vertices_iterator fvi = t->finite_vertices_begin(); fvi != t->finite_vertices_end(); ++fvi) {
                Face_circulator fc = t->incident_faces(fvi), first(fc);
                do {
                    if (t->is_infinite(fc)) {
                        continue;
                    }
                    int index = fc->index(fvi);
                    Vertex_handle left = fc->vertex(fc->cw(index));
                    Vertex_handle right = fc->vertex(fc->ccw(index));
                    visit_k_nearest(Edge(fc, index), left, right, fvi);
                } while (++fc != first);
            }
        }*/

        void make_smoother(T *t, const bool b, const int erosions, const int dilations) {
            if (b) {
                for (int i = 0; i < dilations; i++) {
                    dilation(2, t);
                }
                for (int i = 0; i < erosions; i++) {
                    erosion(2, t);
                }
            } else {
                for (int i = 0; i < erosions; i++) {
                    erosion(2, t);
                }
                for (int i = 0; i < dilations; i++) {
                    dilation(2, t);
                }
            }
        }

        //set in_domain=true for triangles with two in_domain=true neighbors
        static void dilation(const int level, T *t) {
            std::set<Face_handle> changed;
            for (Finite_faces_iterator fi = t->finite_faces_begin(); fi != t->finite_faces_end(); ++fi) {
                if (!fi->is_in_domain()) {
                    int count = 0;
                    if (fi->neighbor(0)->is_in_domain()) {
                        ++count;
                    }
                    if (fi->neighbor(1)->is_in_domain()) {
                        ++count;
                    }
                    if (fi->neighbor(2)->is_in_domain()) {
                        ++count;
                    }
                    if (count >= level) {
                        changed.insert(fi);
                    }
                }
            }
            for (const Face_handle fh: changed) {
                fh->set_in_domain(true);
            }
        }

        //set in_domain=false for triangles with two in_domain=false neighbors
        static void erosion(const int level, T *t) {
            std::set<Face_handle> changed;
            for (Finite_faces_iterator fi = t->finite_faces_begin(); fi != t->finite_faces_end(); ++fi) {
                if (fi->is_in_domain()) {
                    int count = 0;
                    if (!fi->neighbor(0)->is_in_domain()) {
                        count++;
                    }
                    if (!fi->neighbor(1)->is_in_domain()) {
                        count++;
                    }
                    if (!fi->neighbor(2)->is_in_domain()) {
                        count++;
                    }
                    if (count >= level) {
                        changed.insert(fi);
                    }
                }
            }
            for (const Face_handle fh: changed) {
                fh->set_in_domain(false);
            }
        }

        void blue_noise(const int k, const int n, const double radius) {
            std::vector<Point_2> first_list = inverse_sampling(3, radius);
            r->insert_points(first_list.begin(), first_list.end());
            for(int i = 3; i < n; ++i) {
                //std::cout << i << "-th iteration" << std::endl;
                double min_max_distance = 0;
                int best_index;
                std::vector<Point_2> candidates = inverse_sampling(k, radius);
                for(int j = 0; j < k; ++j) {
                    Point_2 candidate = candidates[j];
                    Face_handle face = r->t->locate(candidate);
                    Vertex_handle v0 = face->vertex(0);
                    Vertex_handle v1 = face->vertex(1);
                    Vertex_handle v2 = face->vertex(2);
                    double d0, d1, d2;
                    if (!r->t->is_infinite(v0)) {
                        d0 = r->hyperbolic_distance(v0->point(), candidate);
                    } else {
                        d0 = DBL_MAX;
                    }
                    if (!r->t->is_infinite(v1)) {
                        d1 = r->hyperbolic_distance(v1->point(), candidate);
                    } else {
                        d1 = DBL_MAX;
                    }
                    if (!r->t->is_infinite(v2)) {
                        d2 = r->hyperbolic_distance(v2->point(), candidate);
                    } else {
                        d2 = DBL_MAX;
                    }
                    double min;
                    if(d1 > d0) {
                        if(d2 > d0) {
                            min = d0;
                        }
                        else {
                            min = d2;
                        }
                    } else {
                        if(d2 > d1) {
                            min = d1;
                        }
                        else {
                            min = d2;
                        }
                    }
                    if(min > min_max_distance) {
                        min_max_distance = min;
                        best_index = j;
                    }
                }
                r->insert_point(candidates[best_index]);
            }
        }

        std::vector<Point_2> inverse_sampling(const int number_of_points, const double radius) {
            std::vector<Point_2> samples;
            std::random_device rd;
            std::mt19937 gen(rd());
            const double area = std::cosh(radius) - 1;
            const double x_max = antiderivative(radius, area);
            //std::cout << "area " << area << std::endl;
            //std::cout << "x_max " << x_max << std::endl;
            std::uniform_real_distribution<> x_distr(0, x_max);
            std::uniform_real_distribution<> angle_distr(0, 2 * CGAL_PI);
            for (int i = 0; i < number_of_points; i++) {
                double r = inverse_antiderivative(x_distr(gen), area);
                const double angle = angle_distr(gen);

                //if Beltrami Klein model is used
                if (typeid(Geom_traits) == typeid(Beltrami_klein_traits<>)) {
                    r = std::tanh(r);
                } else {
                    r = std::tanh(r / 2);
                }

                double x_ = std::cos(angle) * r;
                double y_ = std::sin(angle) * r;
                //std::cout << "sampled point: x " << x_ << ", y: " << y_ << std::endl;
                samples.push_back(Point_2(x_, y_));
            }
            return samples;
        }

        std::vector<Point_2> rejection_sampling(const int number_of_points, const double radius) {
            int rejected_counter = 0; //for testing how often point is rejected
            int number_of_tries = 0;
            std::vector<Point_2> samples;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> x_distr(0, radius);
            const double area = std::cosh(radius) - 1;
            const double y_max = probability_density(radius, area);
            std::uniform_real_distribution<> y_distr(0, y_max);
            std::uniform_real_distribution<> angle_distr(0, 2 * CGAL_PI);
            while (samples.size() < number_of_points) {
                double x = x_distr(gen);
                const double y = y_distr(gen);
                if (y <= probability_density(x, area)) {
                    const double angle = angle_distr(gen);

                    //from hyperbolic radius to Euclidean radius in Poincaré disk model
                    //if Beltrami Klein model is used
                    if (typeid(Geom_traits) == typeid(Beltrami_klein_traits<>)) {
                        x = std::tanh(x);
                    } else {
                        x = std::tanh(x / 2);
                    }

                    double x_ = std::cos(angle) * x;
                    double y_ = std::sin(angle) * x;

                    //std::cout << "sampled point: x " << x_ << ", y: " << y_ << std::endl;
                    samples.push_back(Point_2(x_, y_));
                } else {
                    ++rejected_counter;
                }
                ++number_of_tries;
            }
            std::cout << "rejected " << rejected_counter << " points out of " << number_of_tries << " attempts." <<
                    std::endl;
            return samples;
        }

        //opposite of discover_components() algorithm
        //discovers constrained edges
        void discover_edges(T *t) {
            for (Finite_edges_iterator fe = t->finite_edges_begin(); fe != t->finite_edges_end(); ++fe) {
                const Face_handle fh = fe->first;
                const int i = fe->second;
                if (fh->neighbor(i)->is_in_domain() != fh->is_in_domain()) {
                    r->insert_edge_constraint(fh->vertex(fh->cw(i)), fh->vertex(fh->ccw(i)));
                }
            }
        }

    private:
        Routing_scenario<T> *r;
    };

    template<typename T>
    Random_domain_generator<
        T>::Random_domain_generator(Routing_scenario<T> *r) : r(r) {
    }
}

#endif //RANDOMWORLDGENERATOR_H
