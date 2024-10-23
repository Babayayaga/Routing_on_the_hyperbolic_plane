#ifndef ROUTINGSCENARIOGRAPHICSITEM_H
#define ROUTINGSCENARIOGRAPHICSITEM_H

/*
 * displays a Routing_scenario.
 * the interior of obstacles of the Routing_scenario are painted via a QPainterPath
 * the decomposition into closed polygonal chains can be displayed.
 * the RoutingScenarioGraphicsItem has a child, the TriangulationGraphicsItem,
 * to display the (constrained) edges of the ConstrainedDelaunayTriangulation_2,
 * that underlies the Routing_scenario.
 *
 * A transformation map (t_map) is stored, it maps a Vertex_handle to a point.
 * Like with our TriangulationGraphicsItem, the QPainterPath and the decomposition edges aren't recomputed in
 * every paint(), instead they get only recomputed if the underlying triangulation changes (changed()),
 * or when a transformation is applied. So painting is faster.
 */

#include "internal/Qt/HyperbolicPainter.h"
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>
#include "internal/RoutingScenario.h"
#include "internal/Qt/GraphicItems/TriangulationGraphicsItem.h"

#include <QGraphicsScene>
#include <QPainter>

namespace CGAL::Qt {
    template<typename T>
    class Routing_scenario_graphics_item : public QGraphicsItem {
        typedef typename CGAL::Qt::TriangulationGraphicsItem<T> TriangulationGraphicsItem;
        typedef typename Routing_scenario<T> Routing_scenario;
        typedef typename T::Point_2 Point_2;
        typedef typename T::Geom_traits Geom_traits;
        typedef typename Geom_traits::Triangle_2 Triangle_2;
        typedef typename T::Face_handle Face_handle;
        typedef typename T::Vertex_handle Vertex_handle;
        typedef typename T::Finite_vertices_iterator Finite_vertices_iterator;
        typedef typename Hyperbolic_painter<Geom_traits>::QHyperbolic_segment QHyperbolic_segment;

    public:
        explicit Routing_scenario_graphics_item(Routing_scenario *r);

        void changed();

        void repaint();

        void transformation();

        QRectF boundingRect() const override;

        void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) override;

        void set_show_decomposition(const bool b) {
            show_decomposition = b;
            if (b) {
                construct_decomposition();
                for(int i = 0; i < decomposition.size(); ++i) {
                    decomposition_colors.push_back(QColor(QRandomGenerator::global()->bounded(255),
                                              QRandomGenerator::global()->bounded(255),
                                              QRandomGenerator::global()->bounded(255)));
                }
            }
            repaint();
        }

        void set_show_obstacle_interior(const bool b) {
            show_obstacle_interior = b;
            if (b) {
                construct_obstacle_path();
            }
            repaint();
        }

        void construct_obstacle_path() {
            std::list<std::vector<Point_2>> obstacles;
            for (std::vector<Vertex_handle> obstacle: r->get_obstacles()) {
                std::vector<Point_2> points;
                for (Vertex_handle vh: obstacle) {
                    if (transformed) {
                        points.push_back(t_map[vh]);
                    } else {
                        points.push_back(vh->point());
                    }
                }
                obstacles.push_back(points);
            }
            obstacle_path = hyperbolic_painter.construct_path(obstacles, approximation_radius);
        }

        void construct_decomposition() {
            decomposition.clear();
            for (std::vector<Vertex_handle> obstacle: r->get_obstacles()) {
                std::vector<QHyperbolic_segment> obstacle_segments;
                for (int i = 0; i < obstacle.size(); ++i) {
                    QHyperbolic_segment segment;
                    if (transformed) {
                        segment = hyperbolic_painter.construct_segment(
                            t_map[obstacle[i]], t_map[obstacle[(i + 1) % obstacle.size()]], approximation_radius);
                        obstacle_segments.push_back(segment);
                    } else {
                        segment = hyperbolic_painter.construct_segment(
                            obstacle[i]->point(), obstacle[(i + 1) % obstacle.size()]->point(), approximation_radius);
                    }
                    obstacle_segments.push_back(segment);
                }
                decomposition.push_back(obstacle_segments);
            }
        }

        void set_obstacle_brush(const QBrush &brush) {
            obstacle_brush = brush;
            repaint();
        }

        void set_show_point_to_point_visibility(const bool b) {
            show_point_to_point_visibility = b;
            repaint();
        }

        void set_show_point_to_all_visibility(const bool b) {
            show_point_to_all_visibility = b;
            repaint();
        }

        void set_show_visibility_graph(const bool b) {
            show_visibility_graph = b;
            repaint();
        }

        void set_show_dijkstra_tree(const bool b) {
            show_dijkstra_tree = b;
            repaint();
        }

        void set_show_path(const bool b) {
            show_path = b;
            repaint();
        }

        void set_approximation_radius(const double r) {
            approximation_radius = r;
        }

        void set_show_origin(const bool b) {
            show_origin = b;
            repaint();
        }

        void reset_transformation() {
            transformed = false;
            triangulation_graphics_item->reset_transformation();
            t_map.clear();
            location = ORIGIN;
            if (show_obstacle_interior) {
                construct_obstacle_path();
            }
            if(show_decomposition) {
                construct_decomposition();
            }
            repaint();
        }

        //p gets send to origin
        Point_2 translation_to_origin(const Point_2 p, const Point_2 q) {
            Point_2 p_ = p;
            Point_2 q_ = q;

            //check if we are in beltrami klein model
            if (typeid(Geom_traits) == typeid(Beltrami_klein_traits<>)) {
                p_ = r->beltrami_klein_to_poincare(p_);
                q_ = r->beltrami_klein_to_poincare(q_);
            }

            const double px = to_double(p_.x());
            const double py = to_double(p_.y());
            const double qx = to_double(q_.x());
            const double qy = to_double(q_.y());

            const double a = 1 - px * qx - py * qy;
            const double b = py * qx - px * qy;

            const double d = a * a + b * b;
            double real = (qx - px) * a + (qy - py) * b;
            double img = (qy - py) * a - (qx - px) * b;
            real = real / d;
            img = img / d;

            Point_2 result = Point_2(real, img);

            //go back to beltrami klein model
            if (typeid(Geom_traits) == typeid(Beltrami_klein_traits<>)) {
                result = r->poincare_to_beltrami_klein(result);
            }

            return result;
        }

        Point_2 translation_to_origin_inverse(const Point_2 p, const Point_2 q) {
            Point_2 p_ = p;
            Point_2 q_ = q;

            if (typeid(Geom_traits) == typeid(Beltrami_klein_traits<>)) {
                p_ = r->beltrami_klein_to_poincare(p_);
                q_ = r->beltrami_klein_to_poincare(q_);
            }

            const double px = to_double(p_.x());
            const double py = to_double(p_.y());
            const double qx = to_double(q_.x());
            const double qy = to_double(q_.y());

            const double a = 1 + px * qx + py * qy;
            const double b = px * qy - py * qx;

            const double d = a * a + b * b;
            double real = (qx + px) * a + (qy + py) * b;
            double img = (qy + py) * a - (qx + px) * b;
            real = real / d;
            img = img / d;

            Point_2 result = Point_2(real, img);

            if (typeid(Geom_traits) == typeid(Beltrami_klein_traits<>)) {
                result = r->poincare_to_beltrami_klein(result);
            }

            return result;
        }

        Point_2 rotate_point_around_origin(const Point_2 p, const double angle) {
            const double x = p.x();
            const double y = p.y();
            const double c = cos(angle);
            const double s = sin(angle);
            const double x_ = x * c - y * s;
            const double y_ = y * c + x * s;
            return Point_2(x_, y_);
        }

        void set_focus(const Point_2 q) {
            t_map.clear();
            for (Finite_vertices_iterator fi = r->t->finite_vertices_begin(); fi != r->t->finite_vertices_end(); ++fi) {
                t_map.insert({fi, translation_to_origin(q, fi->point())});
            }
            location = translation_to_origin_inverse(q, ORIGIN);
            transformed = true;
        }

        void composite_focusing(const Point_2 q) {
            if (transformed) {
                for (Finite_vertices_iterator fi = r->t->finite_vertices_begin(); fi != r->t->finite_vertices_end(); ++
                     fi) {
                    Point_2 p = t_map[fi];
                    t_map.erase(fi);
                    t_map.insert({fi, translation_to_origin(q, p)});
                }
            } else {
                for (Finite_vertices_iterator fi = r->t->finite_vertices_begin(); fi != r->t->finite_vertices_end(); ++
                     fi) {
                    t_map.insert({fi, translation_to_origin(q, fi->point())});
                }
            }
            transformed = true;
            triangulation_graphics_item->set_tmap(&t_map);
        }

        void composite_rotation(const double angle) {
            if (transformed) {
                for (Finite_vertices_iterator fi = r->t->finite_vertices_begin(); fi != r->t->finite_vertices_end(); ++
                     fi) {
                    Point_2 p = t_map[fi];
                    t_map.erase(fi);
                    t_map.insert({fi, rotate_point_around_origin(p, angle)});
                }
            } else {
                for (Finite_vertices_iterator fi = r->t->finite_vertices_begin(); fi != r->t->finite_vertices_end(); ++
                     fi) {
                    t_map.insert({fi, rotate_point_around_origin(fi->point(), angle)});
                }
            }
            transformed = true;
            triangulation_graphics_item->set_tmap(&t_map);
        }

        void clear();

        Routing_scenario *r;
        TriangulationGraphicsItem *triangulation_graphics_item;
        Hyperbolic_painter<Geom_traits> hyperbolic_painter;
        QPen epen, vpen1, vpen2, vpen3;
        QBrush obstacle_brush;
        Converter<Geom_traits> convert;
        bool show_point_to_point_visibility, show_visibility_graph, show_point_to_all_visibility;
        bool show_dijkstra_tree, show_path;
        bool show_obstacle_interior;
        bool show_origin;
        bool show_decomposition;
        QPainterPath obstacle_path;
        std::vector<std::vector<QHyperbolic_segment>> decomposition;
        std::vector<QColor> decomposition_colors;
        double approximation_radius;

        //holds current transformation
        std::map<Vertex_handle, Point_2> t_map;
        bool transformed;
        Point_2 location;
    };

    template<typename T>
    Routing_scenario_graphics_item<T>::Routing_scenario_graphics_item(Routing_scenario *r)
        : r(r), hyperbolic_painter(0), show_point_to_all_visibility(false)
          , show_point_to_point_visibility(false), show_visibility_graph(false), show_dijkstra_tree(false),
          show_path(false), obstacle_brush(QColor(140, 140, 140, 100)), show_obstacle_interior(false),
          show_origin(false)
          , approximation_radius(10), transformed(false), location(ORIGIN), show_decomposition(false) {
        triangulation_graphics_item = new TriangulationGraphicsItem(r->t);
        triangulation_graphics_item->setFlag(QGraphicsItem::ItemUsesExtendedStyleOption);
        triangulation_graphics_item->setParentItem(this);
        triangulation_graphics_item->setFlag(ItemStacksBehindParent);

        epen.setWidthF(0.000);
        epen.setBrush(::Qt::darkGreen);

        vpen1.setStyle(::Qt::SolidLine);
        vpen1.setWidth(8);
        vpen1.setBrush(::Qt::magenta);
        vpen1.setCapStyle(::Qt::RoundCap);
        vpen1.setJoinStyle(::Qt::RoundJoin);
        vpen2 = vpen1;
        vpen2.setColor(::Qt::red);
        vpen3 = vpen1;
        vpen3.setColor(::Qt::darkGreen);
    }

    template<typename T>
    void Routing_scenario_graphics_item<T>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                                                  QWidget *widget) {
        painter->save();

        hyperbolic_painter = Hyperbolic_painter<Geom_traits>(painter);

        if (show_obstacle_interior) {
            painter->setPen(::Qt::NoPen);
            painter->setBrush(obstacle_brush);
            painter->drawPath(obstacle_path);
            painter->setPen(::Qt::NoBrush);
        }

        if (show_visibility_graph) {
            painter->setPen(epen);
            for (int i = 0; i < r->number_of_vertices(); i++) {
                Vertex_handle vh = r->index_vertex_map[i];
                for (int j = r->offsets[i]; j < r->offsets[i + 1]; j++) {
                    int neighbor_index = r->adjacencies[j];
                    if (i < neighbor_index) {
                        Vertex_handle neighbor = r->index_vertex_map[neighbor_index];
                        QHyperbolic_segment segment;

                        if (!transformed) {
                            segment = hyperbolic_painter.construct_segment(
                                neighbor->point(), vh->point(), approximation_radius);
                        } else {
                            segment = hyperbolic_painter.construct_segment(
                                t_map[neighbor], t_map[vh], approximation_radius);
                        }
                        hyperbolic_painter.draw_hyperbolic_segment(segment);
                    }
                }
            }
        }

        if (show_point_to_point_visibility && r->start_node_handle != nullptr && r->destination_node_handle !=
            nullptr) {
            painter->setBrush(QBrush(::Qt::lightGray));
            painter->setPen(::Qt::NoPen);
            for (Face_handle face: r->crossed_faces) {
                QPainterPath triangle_path;

                if (!transformed) {
                    Triangle_2 triangle = r->t->triangle(face);
                    triangle_path = hyperbolic_painter.construct_polygon_path({
                        triangle.vertex(0), triangle.vertex(1), triangle.vertex(2), triangle.vertex(0)
                    });
                } else {
                    triangle_path = hyperbolic_painter.construct_polygon_path({
                        t_map[face->vertex(0)], t_map[face->vertex(1)], t_map[face->vertex(2)], t_map[face->vertex(0)]
                    });
                }
                painter->drawPath(triangle_path);
            }
            painter->setBrush(::Qt::NoBrush);
            epen.setStyle(::Qt::DashLine);
            painter->setPen(epen);
            QHyperbolic_segment segment;

            if (!transformed) {
                segment = hyperbolic_painter.construct_segment(
                    r->start_point, r->destination_point, 100);
            } else {
                segment = hyperbolic_painter.construct_segment(
                    t_map[r->start_node_handle], t_map[r->destination_node_handle], 100);
            }
            hyperbolic_painter.draw_hyperbolic_segment(segment);
            epen.setStyle(::Qt::SolidLine);
        }

        //show all points visible from start node
        if (show_point_to_all_visibility) {
            painter->setPen(epen);
            for (Vertex_handle vh: r->visibles_start_node) {
                QHyperbolic_segment segment;
                if (!transformed) {
                    segment = hyperbolic_painter.construct_segment(
                        r->start_point, vh->point(), approximation_radius);
                } else {
                    segment = hyperbolic_painter.construct_segment(
                        t_map[r->start_node_handle], t_map[vh], approximation_radius);
                }
                hyperbolic_painter.draw_hyperbolic_segment(segment);
            }
        }

        if (show_dijkstra_tree) {
            painter->setPen(epen);
            const int start_index = r->vertex_index_map[r->start_node_handle];
            for (int i = 0; i < r->dijkstra_result.size(); i++) {
                int p_index = r->dijkstra_result[i].first;
                if (i != start_index && p_index != -1) {
                    QHyperbolic_segment segment;
                    if (!transformed) {
                        segment = hyperbolic_painter.construct_segment(r->index_vertex_map[i]->point(),
                                                                       r->index_vertex_map[p_index]->point(),
                                                                       approximation_radius);
                    } else {
                        segment = hyperbolic_painter.construct_segment(t_map[r->index_vertex_map[i]],
                                                                       t_map[r->index_vertex_map[p_index]],
                                                                       approximation_radius);
                    }
                    hyperbolic_painter.draw_hyperbolic_segment(segment);
                }
            }
        }

        if (show_path) {
            painter->setPen(epen);
            for (int i = 0; i < r->get_indices_path().size() - 1; i++) {
                QHyperbolic_segment segment;

                if (!transformed) {
                    segment = hyperbolic_painter.construct_segment(
                        r->get_point_on_path(i),
                        r->get_point_on_path(i + 1), approximation_radius);
                } else {
                    segment = hyperbolic_painter.construct_segment(
                        t_map[r->get_handle_on_path(i)],
                        t_map[r->get_handle_on_path(i + 1)], approximation_radius);
                }
                hyperbolic_painter.draw_hyperbolic_segment(segment);
            }
        }

        if (show_decomposition) {
            for (int i = 0; i < decomposition.size(); ++i) {
                epen.setColor(decomposition_colors[i]);
                for (QHyperbolic_segment segment : decomposition[i]) {
                    painter->setPen(epen);
                    hyperbolic_painter.draw_hyperbolic_segment(segment);
                }
            }
            epen.setBrush(::Qt::darkGreen);
        }

        painter->restore();

        //if defined, draw start/destination_point
        QTransform matrix = painter->worldTransform();
        painter->resetTransform();
        if (r->start_node_handle != nullptr) {
            painter->setPen(vpen1);
            if (!transformed) {
                painter->drawPoint(matrix.map(convert(r->start_point)));
            } else {
                painter->drawPoint(matrix.map(convert(t_map[r->start_node_handle])));
            }
        }
        if (r->destination_node_handle != nullptr) {
            painter->setPen(vpen2);
            if (!transformed) {
                painter->drawPoint(matrix.map(convert(r->destination_point)));
            } else {
                painter->drawPoint(matrix.map(convert(t_map[r->destination_node_handle])));
            }
        }
        if (show_origin) {
            painter->setPen(vpen3);
            painter->drawPoint(matrix.map(convert(ORIGIN)));
        }
    }

    template<typename T>
    void Routing_scenario_graphics_item<T>::clear() {
        decomposition.clear();
        decomposition_colors.clear();
        reset_transformation();
        triangulation_graphics_item->clear();
        show_visibility_graph = false;
        show_point_to_all_visibility = false;
        show_point_to_point_visibility = false;
        repaint();
    }

    template<typename T>
    QRectF Routing_scenario_graphics_item<T>::boundingRect() const {
        return {-1, -1, 2, 2};
    }

    //only is called, when the underlying triangulation changes
    template<typename T>
    void Routing_scenario_graphics_item<T>::changed() {
        //std::cout << "changed" << std::endl;
        decomposition.clear();
        decomposition_colors.clear();
        if (transformed) {
            for (Finite_vertices_iterator fvi = r->t->finite_vertices_begin(); fvi != r->t->finite_vertices_end(); ++fvi) {
                Point_2 p = translation_to_origin(location, fvi->point());
                t_map.insert({fvi, p});
            }
            triangulation_graphics_item->set_tmap(&t_map);
        }
        if (show_obstacle_interior) {
            construct_obstacle_path();
        }
        if (show_decomposition) {
            construct_decomposition();
            for(int i = 0; i < decomposition.size(); ++i) {
                decomposition_colors.push_back(QColor(QRandomGenerator::global()->bounded(255),
                                          QRandomGenerator::global()->bounded(255),
                                          QRandomGenerator::global()->bounded(255)));
            }
        }
        triangulation_graphics_item->changed();
        repaint();
    }

    template<typename T>
    void Routing_scenario_graphics_item<T>::transformation() {
        //std::cout << "transformation" << std::endl;
        transformed = true;
        triangulation_graphics_item->transformation();
        decomposition.clear();
        if (show_obstacle_interior) {
            construct_obstacle_path();
        }
        if (show_decomposition) {
            construct_decomposition();
        }
        repaint();
    }

    template<typename T>
    void Routing_scenario_graphics_item<T>::repaint() {
        triangulation_graphics_item->set_approximation_radius(approximation_radius);
        triangulation_graphics_item->update();
        update();
    }
}
#endif
