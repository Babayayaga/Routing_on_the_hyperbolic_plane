#ifndef WASDINPUT_H
#define WASDINPUT_H

#include "internal/Qt/Input/QObjectInput.h"
#include "internal/Qt/GraphicItems/RoutingScenarioGraphicsItem.h"
#include <QEvent>

namespace CGAL::Qt {
    template<typename T>
    class WASDInput : public QObjectInput {
    public:
        typedef typename T::Geom_traits K;
        typedef typename K::Circle_2 Circle_2;
        typedef typename T::Face_handle Face_handle;
        typedef typename T::Point_2 Point_2;

        WASDInput(Routing_scenario_graphics_item<T> *r, QObject *parent);

        void set_walking_speed(const double speed) {
            walking_speed = speed;
        }

        void set_rotation_speed(const double speed) {
            rotation_speed = speed;
        }

        double walking_speed = 0.05;
        double rotation_speed = 1.0 / 20;

        void keyPressedEvent(const QKeyEvent *event);

        bool eventFilter(QObject *obj, QEvent *event) override;

        Routing_scenario_graphics_item<T> *routing_scenario_graphics_item;
    };

    template<typename T>
    WASDInput<T>::
    WASDInput(Routing_scenario_graphics_item<T> *rgi, QObject *parent)
        : QObjectInput(parent), routing_scenario_graphics_item(rgi) {
    }

    template<typename T>
    void WASDInput<T>::keyPressedEvent(const QKeyEvent *event) {
        if (event->key() == ::Qt::Key_W) {
            routing_scenario_graphics_item->composite_focusing(Point_2(0, walking_speed));
        }
        if (event->key() == ::Qt::Key_A) {
            routing_scenario_graphics_item->composite_rotation(-CGAL_PI * rotation_speed);
        }
        if (event->key() == ::Qt::Key_S) {
            routing_scenario_graphics_item->composite_focusing(Point_2(0, -walking_speed));
        }
        if (event->key() == ::Qt::Key_D) {
            routing_scenario_graphics_item->composite_rotation(CGAL_PI * rotation_speed);
        }
        Q_EMIT(transformed());
    }

    template<typename T>
    bool
    WASDInput<T>::eventFilter(QObject *obj, QEvent *event) {
        if (event->type() == QEvent::KeyPress) {
            auto *keyEvent = dynamic_cast<QKeyEvent *>(event);
            keyPressedEvent(keyEvent);
            return true;
        }
        return QObject::eventFilter(obj, event);
    }
}
#endif
