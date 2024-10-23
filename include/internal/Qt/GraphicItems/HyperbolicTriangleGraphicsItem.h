#ifndef HYPERBOLICTRIANGLEGRAPHICSITEM_H
#define HYPERBOLICTRIANGLEGRAPHICSITEM_H
#include <QGraphicsItem>
#include <CGAL/Qt/Converter.h>
#include "PoincareDiskTraits.h"
#include "internal/Qt/HyperbolicPainter.h"

namespace CGAL::Qt {
    template<typename T>
    class Hyperbolic_triangle_graphics_item : public QGraphicsItem {
        typedef typename T::Geom_traits Geom_traits;
        typedef typename Geom_traits::Triangle_2 Triangle_2;
        typedef typename T::Point_2 Point_2;

    public:
        Hyperbolic_painter<Geom_traits> hyperbolic_painter;
        Triangle_2 triangle;

        Hyperbolic_triangle_graphics_item();

        void setTriangle(Triangle_2 triangle_2);

        QRectF boundingRect() const override;

        void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) override;
    };

    template<typename T>
    Hyperbolic_triangle_graphics_item<T>::Hyperbolic_triangle_graphics_item() : hyperbolic_painter(0) {}

    template<typename T>
    QRectF Hyperbolic_triangle_graphics_item<T>::boundingRect() const {
        return {-1, -1, 2, 2};
    }

    template<typename T>
    void Hyperbolic_triangle_graphics_item<T>::setTriangle(Triangle_2 triangle_2) {
        this->triangle = triangle_2;
        update();
    }

    template<typename T>
    void Hyperbolic_triangle_graphics_item<T>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                                                    QWidget *widget) {
        hyperbolic_painter = Hyperbolic_painter<Geom_traits>(painter);
        painter->setBrush(QBrush(QColor(100, 00,00,70)));
        painter->setPen(::Qt::NoPen);
        const QPainterPath triangle_path = hyperbolic_painter.construct_polygon_path({
                triangle.vertex(0), triangle.vertex(1), triangle.vertex(2), triangle.vertex(0)
            });
        painter->drawPath(triangle_path);
    }
}
#endif //HYPERBOLICTRIANGLEGRAPHICSITEM_H
