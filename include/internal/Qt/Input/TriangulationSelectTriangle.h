#ifndef TRIANGULATIONSELECTTRIANGLE_H
#define TRIANGULATIONSELECTTRIANGLE_H

#include "internal/Qt/Input/QObjectInput.h"
#include <CGAL/Qt/Converter.h>

#include <QGraphicsSceneMouseEvent>
#include <QEvent>

#include "internal/RoutingScenario.h"
#include "internal/Qt/GraphicItems/HyperbolicTriangleGraphicsItem.h"

namespace CGAL::Qt {
    template<typename T>
    class TriangulationSelectTriangle : public QObjectInput {
    public:
        typedef typename T::Geom_traits K;
        typedef typename K::Triangle_2 Triangle_2;
        typedef typename K::Circle_2 Circle_2;
        typedef typename T::Face_handle Face_handle;
        typedef typename T::Vertex_handle Vertex_handle;
        typedef typename T::Point_2 Point_2;

        TriangulationSelectTriangle(QGraphicsScene *s, Routing_scenario<T> *r, Routing_scenario_graphics_item<T> *rgi,
                                    QObject *parent);

        Hyperbolic_triangle_graphics_item<T> *active_triangle;

    protected:
        void mouseMoveEvent(QGraphicsSceneMouseEvent *event);

        void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);

        void keyPressedEvent(QKeyEvent *event);

        bool eventFilter(QObject *obj, QEvent *event) override;

        Routing_scenario<T> *routing_scenario;
        Routing_scenario_graphics_item<T> *routing_scenario_graphics_item;
        Converter<K> convert;
        QGraphicsScene *scene_;
        Point_2 p;
    };

    template<typename T>
    TriangulationSelectTriangle<T>::
    TriangulationSelectTriangle(QGraphicsScene *s, Routing_scenario<T> *r, Routing_scenario_graphics_item<T> *rgi,
                                QObject *parent) : QObjectInput(parent), routing_scenario(r),
                                                   routing_scenario_graphics_item(rgi),
                                                   scene_(s), p(ORIGIN) {
        active_triangle = new Hyperbolic_triangle_graphics_item<T>();
        active_triangle->setZValue(10);
        active_triangle->setVisible(true);
    }

    template<typename T>
    void
    TriangulationSelectTriangle<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent *event) {
        p = convert(event->scenePos());
        if (std::sqrt(p.x() * p.x() + p.y() * p.y() < 1) && routing_scenario->number_of_vertices() > 2) {
            if (routing_scenario_graphics_item->transformed) {
                p = routing_scenario_graphics_item->translation_to_origin_inverse(routing_scenario_graphics_item->location, p);
            }
            Face_handle face = routing_scenario->t->locate(p);
            if (!routing_scenario->t->is_infinite(face)) {
                Point_2 v0 = face->vertex(0)->point();
                Point_2 v1 = face->vertex(1)->point();
                Point_2 v2 = face->vertex(2)->point();
                Vertex_handle nearest = face->vertex(0);
                if (CGAL::squared_distance(p, v0) > CGAL::squared_distance(p, v1)) {
                    nearest = face->vertex(1);
                }
                if (CGAL::squared_distance(p, nearest->point()) > CGAL::squared_distance(p, v2)) {
                    nearest = face->vertex(2);
                }
                if (event->button() == ::Qt::LeftButton) {
                    active_triangle->setVisible(false);
                    Q_EMIT(changedActiveNode(CGAL::make_object(nearest)));
                } else if (event->button() == ::Qt::RightButton) {
                    active_triangle->setVisible(false);
                    Q_EMIT(changedDestinationNode(CGAL::make_object(nearest)));
                }
            }
        }
    }

    template<typename T>
    void TriangulationSelectTriangle<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event) {
        p = convert(event->scenePos());
    }

    template<typename T>
    void TriangulationSelectTriangle<T>::keyPressedEvent(QKeyEvent *event) {
        if (event->key() == ::Qt::Key_Shift
            && std::sqrt(p.x() * p.x() + p.y() * p.y() < 1)
            && routing_scenario->number_of_vertices() > 2) {
            if (routing_scenario_graphics_item->transformed) {
                p = routing_scenario_graphics_item->translation_to_origin_inverse(routing_scenario_graphics_item->location, p);
            }
            Face_handle f = routing_scenario->t->locate(p);
            if (!routing_scenario->t->is_infinite(f)) {
                active_triangle->setTriangle(
                    Triangle_2(routing_scenario_graphics_item->translation_to_origin(
                                   routing_scenario_graphics_item->location, f->vertex(0)->point()),
                               routing_scenario_graphics_item->translation_to_origin(
                                   routing_scenario_graphics_item->location, f->vertex(1)->point()),
                               routing_scenario_graphics_item->translation_to_origin(
                                   routing_scenario_graphics_item->location, f->vertex(2)->point())));
                active_triangle->setVisible(true);
                return;
            }
        }
        active_triangle->setVisible(false);
    }

    template<typename T>
    bool
    TriangulationSelectTriangle<T>::eventFilter(QObject *obj, QEvent *event) {
        if (event->type() == QEvent::GraphicsSceneMouseRelease) {
            auto *mouseEvent = dynamic_cast<QGraphicsSceneMouseEvent *>(event);
            mouseReleaseEvent(mouseEvent);
            //return true;
        }
        if (event->type() == QEvent::KeyPress) {
            auto *keyEvent = dynamic_cast<QKeyEvent *>(event);
            keyPressedEvent(keyEvent);
            //return true;
        }
        if (event->type() == QEvent::GraphicsSceneMouseMove) {
            auto *mouseEvent = dynamic_cast<QGraphicsSceneMouseEvent *>(event);
            mouseMoveEvent(mouseEvent);
            //return false;
        }
        return QObject::eventFilter(obj, event);
    }
}

#endif
