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
        typedef typename T::Edge Edge;

    public:
        explicit Random_domain_generator(Routing_scenario<T> *r);

        static double probability_density(const double er, const double area) {
            return 2 * er / (1 - er * er) / area;
            //return 2*er; //Euclidean
        }

        static double inverse_antiderivative(const double x, const double area) {
            return std::sqrt(1 - std::exp(-x * area));
        }

        static double antiderivative(const double x, const double area) {
            return -std::log(1 - x * x) / area;
        }

        void generate_random_domain(const int number_of_points, const double radius,
                                    const double threshold, const int erosions, const int dilations, const int min,
                                    const bool b) {
            std::cout << std::endl;
            std::cout << "-----Domain generation-----" << std::endl;
            Timer sum;
            sum.start();

            T *t = r->t;

            Timer timer;
            timer.start();
            std::vector<Point_2> points = inverse_sampling(number_of_points, radius);
            timer.stop();
            std::cout << "Sampling points took: " << timer.time() << " seconds." << std::endl;
            timer.reset();

            timer.start();
            t->insert(points.begin(), points.end());
            timer.stop();
            std::cout << "Inserting points took: " << timer.time() << " seconds." << std::endl;
            timer.reset();

            timer.start();
            set_faces_in_domain(t, threshold);
            timer.stop();
            std::cout << "Initial in_domain assignment took: " << timer.time() << " seconds." << std::endl;
            timer.reset();

            timer.start();
            make_smoother(t, b, erosions, dilations);
            timer.stop();
            std::cout << "Erosion and dilation took: " << timer.time() << " seconds." << std::endl;
            timer.reset();

            timer.start();
            discover_edges(t);
            timer.stop();
            std::cout << "Discovering and inserting constrained edges took: " << timer.time() << " seconds." <<
                    std::endl;
            timer.reset();

            if (min > 3) {
                timer.start();
                remove_small_obstacles(t, min);
                //since we are removing constrained edges, this step is needed
                //when removing constrained edge we don't know the new in_domain value of both incident faces
                //can also be neglected, but then remove_unconstrained_points_in_obstacle_iterior() gives slightly wrong results
                r->discover_components();
                timer.stop();
                std::cout << "Removing obstacles with less than min points took: " << timer.time() << " seconds." <<
                        std::endl;
                timer.reset();
            }

            timer.start();
            r->remove_all_unconstrained_points();
            timer.stop();
            std::cout << "Removing points with no incident constrained edges took: " << timer.time() << " seconds." <<
                    std::endl;
            timer.reset();

            //we could also remove triangles, where all edges are constrained edges
            //and triangles where two sides are constrained edges

            timer.start();
            r->discover_components();
            timer.stop();
            std::cout << "Discovering components took: " << timer.time() << " seconds." << std::endl;
            timer.reset();

            sum.stop();
            std::cout << "-----Generating domain took: " << sum.time() << " seconds.-----" << std::endl;
            std::cout << std::endl;
        }

        void insert_uniformly_distributed_points(const int n, const double radius) {
            /*Timer timer;
            timer.start();*/
            std::vector<Point_2> points = inverse_sampling(n, radius);
            r->insert_points(points.begin(), points.end());
            r->remove_unconstrained_points_in_obstacle_interior();
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
        void set_faces_in_domain(T *t, const double threshold) {
            std::map<Vertex_handle, float> map;
            Random random;
            for (Finite_vertices_iterator fi = t->finite_vertices_begin(); fi != t->finite_vertices_end(); ++fi) {
                map.insert({fi, random.get_double(-1, 1)});
            }

            for (Finite_faces_iterator fi = t->finite_faces_begin(); fi != t->finite_faces_end(); ++fi) {
                const double x = map[fi->vertex(0)];
                const double y = map[fi->vertex(1)];
                const double z = map[fi->vertex(2)];
                if ((x + y + z) / 3 > threshold) {
                    fi->set_in_domain(true);
                } else {
                    fi->set_in_domain(false);
                }
            }
        }

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

        std::vector<Point_2> inverse_sampling(const int number_of_points, const double radius) {
            std::vector<Point_2> samples;
            std::random_device rd;
            std::mt19937 gen(rd());
            const double area = std::cosh(std::log((1 + radius) / (1 - radius))) - 1;
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
                    r = 2 * r / (1 + r * r);
                }

                double x_ = std::cos(angle) * r;
                double y_ = std::sin(angle) * r;
                //std::cout << "sampled point: x " << x_ << ", y: " << y_ << std::endl;
                samples.push_back(Point_2(x_, y_));
            }
            return samples;
        }

        std::vector<Point_2> rejection_sampling(const int number_of_points, const double radius) {
            std::vector<Point_2> samples;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> x_distr(0, radius);
            //one could simplify the area computation
            const double area = std::cosh(std::log((1 + radius) / (1 - radius))) - 1;
            const double y_max = probability_density(radius, area);
            std::uniform_real_distribution<> y_distr(0, y_max);
            std::uniform_real_distribution<> angle_distr(0, 2 * CGAL_PI);
            //std::cout << "y_max: " << y_max << std::endl;
            while (samples.size() < number_of_points) {
                double x = /*radius * std::sqrt(*/x_distr(gen)/*)*/;
                const double y = y_distr(gen);
                //std::cout << "x: " << x << " y: " << y << std::endl;
                //std::cout << "std::log((x + 1) / (1 - x)) = " << std::log((x + 1) / (1 - x)) << std::endl;
                if (y <= probability_density(x, area)) {
                    const double angle = angle_distr(gen);

                    //if Beltrami Klein model is used
                    if (typeid(Geom_traits) == typeid(Beltrami_klein_traits<>)) {
                        x = 2 * x / (1 + x * x);
                    }

                    double x_ = std::cos(angle) * x;
                    double y_ = std::sin(angle) * x;
                    //std::cout << "sampled point: x " << x_ << ", y: " << y_ << std::endl;
                    samples.push_back(Point_2(x_, y_));
                }
            }
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
