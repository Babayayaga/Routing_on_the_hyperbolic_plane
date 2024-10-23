#ifndef CGAL_QT_POINT_INPUT_AND_TRIANGULATION_CONFLICT_ZONE
#define CGAL_QT_POINT_INPUT_AND_TRIANGULATION_CONFLICT_ZONE

#include "internal/Qt/Input/QObjectInput.h"
#include <CGAL/Qt/Converter.h>

#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <QObject>
#include <CGAL/Object.h>

#include "internal/RoutingScenario.h"
#include "internal/Qt/GraphicItems/RoutingScenarioGraphicsItem.h"

namespace CGAL::Qt {
    template<typename T>
    class PointInput : public QObjectInput {
    public:
        typedef typename T::Geom_traits K;
        typedef typename K::Circle_2 Circle_2;
        typedef typename T::Face_handle Face_handle;
        typedef typename T::Point_2 Point_2;

        PointInput(QGraphicsScene *s, Routing_scenario_graphics_item<T> *rgi, QObject *parent);

    protected:
        void mousePressEvent(const QGraphicsSceneMouseEvent *event);

        void mouseReleaseEvent(const QGraphicsSceneMouseEvent *event);

        bool eventFilter(QObject *obj, QEvent *event) override;

        Routing_scenario_graphics_item<T> *routing_scenario_graphics_item;
        Converter<K> convert;
        QGraphicsScene *scene_;
        Point_2 p;
    };

    template<typename T>
    PointInput<T>::
    PointInput(QGraphicsScene *s, Routing_scenario_graphics_item<T> *rgi, QObject *parent)
        : QObjectInput(parent), routing_scenario_graphics_item(rgi), scene_(s) {
    }

    template<typename T>
    void PointInput<T>::mousePressEvent(const QGraphicsSceneMouseEvent *event) {
        p = convert(event->scenePos());
    }

    template<typename T> void
    PointInput<T>::mouseReleaseEvent(const QGraphicsSceneMouseEvent *event) {
        //float rounded_x = std::round(CGAL::to_double(p.x()) * 1000.0) / 1000.0;
        //float rounded_y = std::round(CGAL::to_double(p.y()) * 1000.0) / 1000.0;
        //p = Point_2(rounded_x, rounded_y);

        if(p.x() * p.x() + p.y() * p.y() < 1) {
            if (routing_scenario_graphics_item->location != ORIGIN) {
                p = routing_scenario_graphics_item->translation_to_origin_inverse(routing_scenario_graphics_item->location, p);
            }
            if (event->button() == ::Qt::LeftButton) {
                Q_EMIT(changedActiveNode(CGAL::make_object(p)));
            }
            else if(event->button() == ::Qt::RightButton) {
                if (event->modifiers() & ::Qt::ShiftModifier) {
                    Q_EMIT (generate(CGAL::make_object(p)));
                } else {
                    Q_EMIT(changedDestinationNode(CGAL::make_object(p)));
                }
            }
        }
    }

    template<typename T> bool
    PointInput<T>::eventFilter(QObject *obj, QEvent *event) {
        if (event->type() == QEvent::GraphicsSceneMousePress) {
            auto *mouseEvent = dynamic_cast<QGraphicsSceneMouseEvent *>(event);
            mousePressEvent(mouseEvent);
            return true;
        }
        if (event->type() == QEvent::GraphicsSceneMouseRelease) {
            auto *mouseEvent = dynamic_cast<QGraphicsSceneMouseEvent *>(event);
            mouseReleaseEvent(mouseEvent);
            return true;
        }
        return QObject::eventFilter(obj, event);
    }
} // namespace CGAL::Qt


#endif
