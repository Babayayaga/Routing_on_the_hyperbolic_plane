#ifndef POLYLINEINPUT_H
#define POLYLINEINPUT_H
#include <CGAL/export/Qt.h>

#include <QPolygonF>
#include <QPointF>
#include <QKeyEvent>

#include "internal/Qt/Input/QObjectInput.h"
#include <CGAL/Qt/Converter.h>
#include <QGraphicsLineItem>
#include <QGraphicsItemGroup>

#include "../GraphicItems/HyperbolicSegmentGraphicsItem.h"

#include <boost/variant.hpp>

namespace CGAL::Qt {
    template<typename T>
    class Polyline_input : public QObjectInput {
    public:
        typedef typename T::Geom_traits K;
        typedef typename T::Point_2 Point_2;
        typedef typename K::Circle_2 Circle_2;
        typedef typename K::Construct_segment_2 Construct_segment_2;
        typedef typename K::Segment_2 Segment_2;

        Polyline_input(QObject *parent, QGraphicsScene *s, Routing_scenario_graphics_item<T>* rgi);

    protected:
        bool eventFilter(QObject *obj, QEvent *event) override;

        bool mousePressEvent(const QGraphicsSceneMouseEvent *event);

        void mouseMoveEvent(const QGraphicsSceneMouseEvent *event);

        void rubberbands(Point_2 p);

        void generate_polygon();

        std::vector<typename K::Point_2> points;
        Converter<K> convert;
        QGraphicsItemGroup *path_item;
        Hyperbolic_segment_graphics_item<T> *b, *e;
        QPointF sp;
        QGraphicsScene *scene_;
        Routing_scenario_graphics_item<T>* routing_scenario_graphics_item;
    };

    template<typename T>
    Polyline_input<T>::Polyline_input(QObject *parent, QGraphicsScene *s, Routing_scenario_graphics_item<T>* rgi)
        : QObjectInput(parent), path_item(nullptr),
          b(nullptr), e(nullptr), scene_(s), routing_scenario_graphics_item(rgi) {
    }

    template<typename T>
    bool Polyline_input<T>::mousePressEvent(const QGraphicsSceneMouseEvent *event) {
        Point_2 p = convert(event->scenePos());
        if (p.x() * p.x() + p.y() * p.y() < 1) {
            if (event->modifiers()) {
                return false;
            }
            if (event->button() != ::Qt::RightButton
                && event->button() != ::Qt::LeftButton) {
                return false;
            }
            points.push_back(p);
            if (event->button() == ::Qt::RightButton) {
                generate_polygon();
                points.clear();
                if (b) {
                    scene_->removeItem(b);
                    delete b;
                    b = nullptr;
                }
                if (path_item) {
                    scene_->removeItem(path_item);
                    delete path_item;
                    path_item = nullptr;
                }
                if (e) {
                    scene_->removeItem(e);
                    delete e;
                    e = nullptr;
                }
                return true;
            }
            if (event->button() == ::Qt::LeftButton) {
                if (!path_item) {
                    path_item = new QGraphicsItemGroup();
                    path_item->setZValue(1);
                    scene_->addItem(path_item);
                }
                if (e) {
                    path_item->addToGroup(e->clone());
                }
                return true;
            }
        }
        return false;
    }

    template<typename T>
    void Polyline_input<T>::mouseMoveEvent(const QGraphicsSceneMouseEvent *event) {
        Point_2 p = convert(event->scenePos());
        if (p.x() * p.x() + p.y() * p.y() < 1) {
            rubberbands(p);
        }
    }

    template<typename T>
    void Polyline_input<T>::rubberbands(const Point_2 p) {
        if (points.empty()) {
            return;
        }
        if (!b) {
            b = new Hyperbolic_segment_graphics_item<T>();
            scene_->addItem(b);
            b->setZValue(1);
        }
        if (!e) {
            e = new Hyperbolic_segment_graphics_item<T>();
            scene_->addItem(e);
            e->setZValue(1);
        }
        b->setPoints(p, points.front());
        e->setPoints(p, points.back());
    }

    template<typename T>
    bool Polyline_input<T>::eventFilter(QObject *obj, QEvent *event) {
        if (event->type() == QEvent::GraphicsSceneMousePress) {
            auto *mouseEvent = dynamic_cast<QGraphicsSceneMouseEvent *>(event);
            if (!mousePressEvent(mouseEvent)) {
                return QObject::eventFilter(obj, event);
            }
        } else if (event->type() == QEvent::GraphicsSceneMouseMove) {
            auto *mouseEvent = dynamic_cast<QGraphicsSceneMouseEvent *>(event);
            mouseMoveEvent(mouseEvent);
            return QObject::eventFilter(obj, event);
        }
        return QObject::eventFilter(obj, event);
    }

    template<typename T>
    void Polyline_input<T>::generate_polygon() {
        if (points.size() > 2) {
            if(routing_scenario_graphics_item->location != ORIGIN) {
                std::vector<Point_2> t_points;
                for(Point_2 p : points) {
                    t_points.push_back(routing_scenario_graphics_item->translation_to_origin_inverse(routing_scenario_graphics_item->location, p));
                }
                Q_EMIT(generate(CGAL::make_object(t_points)));
            } else {
                Q_EMIT(generate(CGAL::make_object(points)));
            }
        }
    }
}
#endif // POLYLINEINPUT_H
