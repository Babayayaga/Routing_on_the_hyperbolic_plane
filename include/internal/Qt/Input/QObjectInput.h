#ifndef GRAPHICSVIEWINPUT_H
#define GRAPHICSVIEWINPUT_H

#include <CGAL/export/Qt.h>
#include <CGAL/auto_link/Qt.h>

#include <QObject>

#include "../../PoincareDiskTraits.h"
#include <CGAL/Point_2.h>

namespace CGAL::Qt {

    class QObjectInput : public QObject {
        Q_OBJECT

    public:
        explicit QObjectInput(QObject *parent)
            : QObject(parent) {
        }

    Q_SIGNALS:
        void generate(CGAL::Object o);
        void modelChanged();
        void changedActiveNode(CGAL::Object o);
        void changedDestinationNode(CGAL::Object o);
        void transformed();

    public Q_SLOTS:
        virtual void processInput(CGAL::Object /*o*/) {
        }
    };
} // namespace CGAL::Qt
#endif //GRAPHICSVIEWINPUT_H
