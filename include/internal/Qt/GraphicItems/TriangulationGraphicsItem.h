#ifndef CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H
#define CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H

/*
 * displays a ConstrainedDelaunayTriangulation_2
 */

#include "internal/Qt/HyperbolicPainter.h"
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>

namespace CGAL::Qt {
    template<typename T>
    class TriangulationGraphicsItem : public QGraphicsItem {
        typedef typename T::Geom_traits Geom_traits;
        typedef typename T::Point_2 Point_2;
        typedef typename T::Face_handle Face_handle;
        typedef typename T::Vertex_handle Vertex_handle;
        typedef typename Geom_traits::Triangle_2 Triangle_2;
        typedef typename Hyperbolic_painter<Geom_traits>::QHyperbolic_segment QHyperbolic_segment;

    public:
        explicit TriangulationGraphicsItem(T *t_);

        void changed();

        void repaint();

        void transformation();

        QRectF boundingRect() const override;

        void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) override;

        void set_vertices_pen(const QPen &pen) {
            vertices_pen = pen;
        }

        void set_edges_pen(const QPen &pen) {
            edges_pen = pen;
        }

        void set_constraints_pen(const QPen &pen) {
            constraints_pen = pen;
        }

        void set_visible_nodes(const bool b) {
            show_vertices = b;
            repaint();
        }

        void set_visible_edges(const bool b) {
            show_edges = b;
            if (show_edges) {
                construct_edges();
            }
            repaint();
        }

        void set_show_triangulation_between_obstacles(const bool b) {
            show_triangulation_between_obstacles = b;
            if (show_triangulation_between_obstacles) {
                construct_edges_between_obstacles();
            }
            repaint();
        }

        void set_visible_constraints(const bool b) {
            show_constraints = b;
            if (show_constraints) {
                construct_constrained_edges();
            }
            repaint();
        }

        void set_approximation_radius(const double r) {
            approximation_bound = r;
        }

        void set_tmap(std::map<Vertex_handle, Point_2> *map) {
            t_map = map;
            transformed = true;
            repaint();
        }

        void reset_transformation() {
            transformed = false;
            edges.clear();
            constrained_edges.clear();
            if (show_constraints) {
                construct_constrained_edges();
            }
            if (show_edges) {
                construct_edges();
            }
            if(show_triangulation_between_obstacles) {
                construct_edges_between_obstacles();
            }
            repaint();
        }

        void clear() {
            edges.clear();
            constrained_edges.clear();
            transformed = false;
        }

    protected:
        void paint_vertices(QPainter *painter);

        void construct_edges();

        void construct_edges_between_obstacles();

        void construct_constrained_edges();

        T *t;
        Hyperbolic_painter<Geom_traits> hyperbolic_painter;
        QPen vertices_pen;
        QPen edges_pen;
        QPen constraints_pen;
        bool show_edges, show_triangulation_between_obstacles;
        bool show_vertices;
        bool show_constraints;
        Converter<Geom_traits> convert;
        double approximation_bound;
        bool transformed;
        std::map<Vertex_handle, Point_2> *t_map;
        std::list<QHyperbolic_segment> edges;
        std::list<QHyperbolic_segment> constrained_edges;
    };


    template<typename T>
    TriangulationGraphicsItem<T>::TriangulationGraphicsItem(T *t_)
        : t(t_), hyperbolic_painter(0), show_edges(false), show_triangulation_between_obstacles(false),
          show_vertices(false), show_constraints(false),
          approximation_bound(10), transformed(false) {
        set_vertices_pen(QPen(::Qt::red, 2, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
        set_edges_pen(QPen(::Qt::black, 0, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
        set_constraints_pen(QPen(::Qt::darkMagenta, 0.01, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
    }

    template<typename T>
    QRectF TriangulationGraphicsItem<T>::boundingRect() const {
        return {-1, -1, 2, 2};
    }

    template<typename T>
    void
    TriangulationGraphicsItem<T>::paint_vertices(QPainter *painter) {
        Converter<Geom_traits> convert;
        painter->setPen(vertices_pen);
        QTransform matrix = painter->worldTransform();
        painter->resetTransform();
        for (typename T::Finite_vertices_iterator it = t->finite_vertices_begin();
             it != t->finite_vertices_end(); ++it) {
            QPointF point;
            if (transformed) {
                point = matrix.map(convert((*t_map)[it]));
                painter->drawPoint(point);
            } else {
                point = matrix.map(convert(it->point()));
                painter->drawPoint(point);
            }
        }
    }

    template<typename T>
    void TriangulationGraphicsItem<T>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                                             QWidget *widget) {
        painter->setPen(edges_pen);
        hyperbolic_painter = Hyperbolic_painter<Geom_traits>(painter);

        if (show_edges || show_triangulation_between_obstacles) {
            for (QHyperbolic_segment edge: edges) {
                hyperbolic_painter.draw_hyperbolic_segment(edge);
            }
        }
        if (show_constraints) {
            painter->setPen(constraints_pen);
            for (QHyperbolic_segment edge: constrained_edges) {
                hyperbolic_painter.draw_hyperbolic_segment(edge);
            }
        }

        if (show_vertices) {
            paint_vertices(painter);
        }
    }

    template<typename T>
    void TriangulationGraphicsItem<T>::construct_edges() {
        edges.clear();
        if (transformed) {
            for (typename T::Finite_edges_iterator eit = t->finite_edges_begin();
                 eit != t->finite_edges_end(); ++eit) {
                Vertex_handle vh = eit->first->vertex(eit->first->cw(eit->second));
                Vertex_handle vi = eit->first->vertex(eit->first->ccw(eit->second));
                QHyperbolic_segment segment = hyperbolic_painter.construct_segment(
                    (*t_map)[vh], (*t_map)[vi], approximation_bound);
                edges.push_back(segment);
            }
        } else {
            for (typename T::Finite_edges_iterator eit = t->finite_edges_begin();
                 eit != t->finite_edges_end(); ++eit) {
                Point_2 p = eit->first->vertex(eit->first->cw(eit->second))->point();
                Point_2 q = eit->first->vertex(eit->first->ccw(eit->second))->point();
                QHyperbolic_segment segment = hyperbolic_painter.construct_segment(p, q, approximation_bound);
                edges.push_back(segment);
            }
        }
    }

    //works only if all triangles got their in_domain values
    template<typename T>
    void TriangulationGraphicsItem<T>::construct_edges_between_obstacles() {
        edges.clear();
        if (transformed) {
            for (typename T::Finite_edges_iterator eit = t->finite_edges_begin();
                 eit != t->finite_edges_end(); ++eit) {
                Face_handle fh = eit->first;
                Face_handle mirrored = t->mirror_edge(*eit).first;
                if(!(fh->is_in_domain() && mirrored->is_in_domain())) {
                    Vertex_handle vh = eit->first->vertex(eit->first->cw(eit->second));
                    Vertex_handle vi = eit->first->vertex(eit->first->ccw(eit->second));
                    QHyperbolic_segment segment = hyperbolic_painter.construct_segment(
                        (*t_map)[vh], (*t_map)[vi], approximation_bound);
                    edges.push_back(segment);
                }
            }
        } else {
            for (typename T::Finite_edges_iterator eit = t->finite_edges_begin();
                 eit != t->finite_edges_end(); ++eit) {
                Face_handle fh = eit->first;
                Face_handle mirrored = t->mirror_edge(*eit).first;
                if(!(fh->is_in_domain() && mirrored->is_in_domain())) {
                    Point_2 p = eit->first->vertex(eit->first->cw(eit->second))->point();
                    Point_2 q = eit->first->vertex(eit->first->ccw(eit->second))->point();
                    QHyperbolic_segment segment = hyperbolic_painter.construct_segment(p, q, approximation_bound);
                    edges.push_back(segment);
                }
            }
        }
    }

    template<typename T>
    void TriangulationGraphicsItem<T>::construct_constrained_edges() {
        constrained_edges.clear();
        if (transformed) {
            for (typename T::Constrained_edges_iterator cei = t->constrained_edges_begin();
                 cei != t->constrained_edges_end(); ++cei) {
                Vertex_handle vh = cei->first->vertex(cei->first->cw(cei->second));
                Vertex_handle vi = cei->first->vertex(cei->first->ccw(cei->second));
                QHyperbolic_segment segment = hyperbolic_painter.construct_segment((*t_map)[vh], (*t_map)[vi],
                    approximation_bound);
                constrained_edges.push_back(segment);
            }
        } else {
            for (typename T::Constrained_edges_iterator cei = t->constrained_edges_begin();
                 cei != t->constrained_edges_end(); ++cei) {
                Point_2 p = cei->first->vertex(cei->first->cw(cei->second))->point();
                Point_2 q = cei->first->vertex(cei->first->ccw(cei->second))->point();
                QHyperbolic_segment segment = hyperbolic_painter.construct_segment(p, q, approximation_bound);
                constrained_edges.push_back(segment);
            }
        }
    }

    template<typename T>
    void TriangulationGraphicsItem<T>::changed() {
        edges.clear();
        constrained_edges.clear();
        if (t->number_of_vertices() > 0) {
            if (show_edges) {
                construct_edges();
            }
            if (show_triangulation_between_obstacles) {
                construct_edges_between_obstacles();
            }
            if (show_constraints) {
                construct_constrained_edges();
            }
        }
        repaint();
    }

    template<typename T>
    void TriangulationGraphicsItem<T>::transformation() {
        transformed = true;
        edges.clear();
        constrained_edges.clear();
        if (show_edges) {
            construct_edges();
        } else if (show_triangulation_between_obstacles) {
            construct_edges_between_obstacles();
        }
        if (show_constraints) {
            construct_constrained_edges();
        }
        repaint();
    }

    template<typename T>
    void TriangulationGraphicsItem<T>::repaint() {
        update();
    }
}
#endif
