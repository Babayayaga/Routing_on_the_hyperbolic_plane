#ifndef CGAL_HYPERBOLIC_PAINTER_OSTREAM_H
#define CGAL_HYPERBOLIC_PAINTER_OSTREAM_H

/*
 * used to construct and draw hyperbolic segments and the QPainterPath for
 * the obstacle interiors, with Qt
 */

#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/Converter.h>
#include "internal/PoincareDiskTraits.h"

namespace CGAL::Qt {
    template<typename Traits>
    class Hyperbolic_painter {
        typedef Poincare_disk_traits<> Gt;
        typedef Gt::Hyperbolic_segment_2 Hyperbolic_segment_2;
        typedef Gt::Circular_arc_2 Circular_arc_2;
        typedef Gt::Euclidean_segment_2 Euclidean_segment_2;
        typedef Gt::Point_2 Point_2;
        typedef Gt::Circle_2 Circle_2;
        typedef Gt::Triangle_2 Triangle_2;
        typedef Gt::Construct_segment_2 Construct_segment_2;
        typedef Gt::FT FT;

        Gt gt;
        constexpr static double APPROXIMATION_BOUND = 100;
        constexpr static double COEFF = 180 * 16 / CGAL_PI;
        QPainter *qp{};
        Converter<Gt> convert;

    public:
        struct QArc {
            Circle_2 circle;
            int angle;
            int span;
        };

        typedef std::variant<QArc, Euclidean_segment_2> QHyperbolic_segment;

        explicit Hyperbolic_painter(QPainter *p, QRectF rect = QRectF()) : qp(p), convert(rect) {
        }

        static std::pair<double, double> getAnglesFromArc(const Circular_arc_2 *arc) {
            const Circle_2 &circ = arc->supporting_circle();
            const Point_2 &center = circ.center();
            const Point_2 &source = arc->source();
            const Point_2 &target = arc->target();
            double asource = std::atan2(-to_double(source.y() - center.y()),
                                        to_double(source.x() - center.x()));
            double atarget = std::atan2(-to_double(target.y() - center.y()),
                                        to_double(target.x() - center.x()));
            std::swap(asource, atarget);
            double aspan = atarget - asource;
            if (aspan < 0) {
                aspan += 2 * CGAL_PI;
            }
            return std::make_pair(int(asource * COEFF), int(aspan * COEFF));
        }

        QHyperbolic_segment construct_segment(const Point_2 p, const Point_2 q,
                                              const double approximation_bound = APPROXIMATION_BOUND) {
            //if we are in beltrami klein model always draw euclidean segment
            if (typeid(Traits) == typeid(Beltrami_klein_traits<>)) {
                return Euclidean_segment_2(p, q);
            }
            Hyperbolic_segment_2 hyperbolic_segment = gt.construct_segment_2_object()(p, q);
            if (const Euclidean_segment_2 *seg = boost::get<Euclidean_segment_2>(&hyperbolic_segment)) {
                return Euclidean_segment_2(seg->source(), seg->target());
            }
            const Circular_arc_2 &arc = boost::get<const Circular_arc_2 &>(hyperbolic_segment);
            if (arc.squared_radius() > approximation_bound) {
                return Euclidean_segment_2(arc.source(), arc.target());
            }
            const std::pair<double, double> angles = getAnglesFromArc(&arc);
            return QArc(arc.supporting_circle(), angles.first, angles.second);
        }

        void draw_hyperbolic_segment(QHyperbolic_segment segment) {
            if (std::holds_alternative<QArc>(segment)) {
                QArc arc = std::get<QArc>(segment);
                qp->drawArc(convert(arc.circle.bbox()), arc.angle, arc.span);
            } else {
                Euclidean_segment_2 seg = std::get<Euclidean_segment_2>(segment);
                const Point_2 &source = seg.source();
                const Point_2 &target = seg.target();
                const QPointF src(to_double(source.x()), to_double(source.y()));
                const QPointF tgt(to_double(target.x()), to_double(target.y()));
                qp->drawLine(src, tgt);
            }
        }

        QPainterPath construct_polygon_path(std::vector<Point_2> obstacle,
                                            const double approximation_bound = APPROXIMATION_BOUND) const {
            QPainterPath obstaclePath;
            const Construct_segment_2 construct_segment_2 = gt.construct_segment_2_object();

            Point_2 oldPoint = obstacle.front();
            obstaclePath.moveTo(QPointF(to_double(oldPoint.x()), to_double(oldPoint.y())));
            for (auto it = ++obstacle.begin();
                 it != obstacle.end(); ++it) {
                Point_2 current = *it;

                if (typeid(Traits) == typeid(Beltrami_klein_traits<>)) {
                    QPointF qtarget(to_double(current.x()), to_double(current.y()));
                    obstaclePath.lineTo(qtarget);
                    oldPoint = current;
                    continue;
                }

                Hyperbolic_segment_2 hyperbolicSegment = construct_segment_2(oldPoint, current);
                if (const Circular_arc_2 *arc = boost::get<Circular_arc_2>(&hyperbolicSegment)) {
                    if (arc->squared_radius() <= approximation_bound) {
                        std::pair<double, double> angles = getAnglesFromArc(arc);
                        if (current == arc->target()) {
                            obstaclePath.arcTo(convert(arc->supporting_circle().bbox()),
                                               angles.first / 16.0 + angles.second / 16.0, -angles.second / 16.0);
                        } else if (current == arc->source()) {
                            obstaclePath.arcTo(convert(arc->supporting_circle().bbox()), angles.first / 16.0,
                                               angles.second / 16.0);
                        }
                    } else {
                        QPointF qtarget(to_double(current.x()), to_double(current.y()));
                        obstaclePath.lineTo(qtarget);
                    }
                } else if (boost::get<Euclidean_segment_2>(&hyperbolicSegment)) {
                    QPointF qtarget(to_double(current.x()), to_double(current.y()));
                    obstaclePath.lineTo(qtarget);
                }
                oldPoint = current;
            }
            return obstaclePath;
        }

        QPainterPath construct_path(std::list<std::vector<Point_2> > obstacles,
                                    const double approximation_bound = APPROXIMATION_BOUND) const {
            QPainterPath obstaclePath;
            for (std::list<std::vector<Point_2>>::iterator obstacle_iterator = obstacles.begin();
                 obstacle_iterator != obstacles.end(); ++obstacle_iterator) {
                std::vector<Point_2> obstacle = *obstacle_iterator;
                obstacle.push_back(obstacle.front());
                obstaclePath.addPath(construct_polygon_path(obstacle, approximation_bound));
            }
            return obstaclePath;
        }
    };
}

#endif
