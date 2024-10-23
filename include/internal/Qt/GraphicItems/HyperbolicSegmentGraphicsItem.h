#ifndef QGRAPHICSHYPERBOLICSEGMENTITEM_H
#define QGRAPHICSHYPERBOLICSEGMENTITEM_H

#include <QGraphicsItem>
#include <CGAL/Qt/Converter.h>
#include "../../PoincareDiskTraits.h"
#include "internal/Qt/HyperbolicPainter.h"

namespace CGAL::Qt {
    template<typename T>
    class Hyperbolic_segment_graphics_item : public QGraphicsItem {
        typedef typename T::Geom_traits Geom_traits;
        typedef typename T::Point_2 Point_2;
        typedef typename Hyperbolic_painter<Geom_traits>::QHyperbolic_segment QHyperbolic_segment;

    public:
        Point_2 p1;
        Point_2 p2;
        Hyperbolic_painter<Geom_traits> hyperbolic_painter;

        Hyperbolic_segment_graphics_item();

        void setPoints(Point_2 p1, Point_2 p2);

        Hyperbolic_segment_graphics_item *clone();

        QRectF boundingRect() const override;

        void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) override;
    };

    template<typename T>
    Hyperbolic_segment_graphics_item<T>::Hyperbolic_segment_graphics_item() : hyperbolic_painter(0) {
    }

    template<typename T>
    QRectF Hyperbolic_segment_graphics_item<T>::boundingRect() const {
        return {-1, -1, 2, 2};
    }

    template<typename T>
    void Hyperbolic_segment_graphics_item<T>::setPoints(Point_2 p1, Point_2 p2) {
        this->p1 = p1;
        this->p2 = p2;
        update();
    }

    template<typename T>
    Hyperbolic_segment_graphics_item<T> *Hyperbolic_segment_graphics_item<T>::clone() {
        auto *clone = new Hyperbolic_segment_graphics_item();
        clone->setPoints(p1, p2);
        return clone;
    }

    template<typename T>
    void Hyperbolic_segment_graphics_item<T>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                                                    QWidget *widget) {
        hyperbolic_painter = Hyperbolic_painter<Geom_traits>(painter);
        painter->setPen(QPen(::Qt::red, 0.01, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));

        QHyperbolic_segment segment = hyperbolic_painter.construct_segment(p1, p2, 50);

        hyperbolic_painter.draw_hyperbolic_segment(segment);
    }
}

#endif
