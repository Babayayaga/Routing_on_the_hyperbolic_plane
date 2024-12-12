#ifndef BELTRAMIKLEINTRAITS_H
#define BELTRAMIKLEINTRAITS_H

/*
 * Traits for a Constrained_Delaunay_Triangulation_2 in the BeltramiKlein model
 * we inherit from the predefined Exact_predicates_inexact_constructions_kernel and overwrite
 * the Side_of_oriented_circle_2 predicate, since we're using a Delaunay triangulation and circles in the
 * Beltrami-Klein model aren't Euclidean circles.
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {

    template<typename Traits>
    class Side_of_oriented_circle_2 {
        typedef typename Traits::Point_2 Point_2;
        typedef typename Traits::Euclidean_side_of_oriented_circle_2 Euclidean_side_of_oriented_circle_2;

    public:
        explicit Side_of_oriented_circle_2(const Traits &gt = Traits()) : _gt(gt) {}

        Oriented_side operator()(const Point_2 p, const Point_2 q, const Point_2 r, const Point_2 query) {
            //translate p,q,r,s to Poincare disk model,
            //in this model circles = Euclidian circles
            const double abs_p = p.x() * p.x() + p.y() * p.y();
            double d_p = 1 + sqrt(1 - abs_p);
            const Point_2 p_ = Point_2(p.x() / d_p, p.y() / d_p);

            const double abs_q = q.x() * q.x() + q.y() * q.y();
            double d_q = 1 + sqrt(1 - abs_q);
            const Point_2 q_ = Point_2(q.x() / d_q, q.y() / d_q);

            const double abs_r = r.x() * r.x() + r.y() * r.y();
            double d_r = 1 + sqrt(1 - abs_r);
            const Point_2 r_ = Point_2(r.x() / d_r, r.y() / d_r);

            const double abs_query = query.x() * query.x() + query.y() * query.y();
            double d_query = 1 + sqrt(1 - abs_query);
            const Point_2 query_ = Point_2(query.x() / d_query, query.y() / d_query);

            return _gt.euclidean_side_of_oriented_circle_2_object()(p_, q_, r_, query_);
        }

    private:
        const Traits &_gt;
    };

    template<typename Kernel = Exact_predicates_inexact_constructions_kernel>
    class Beltrami_klein_traits : public Kernel {
        typedef Beltrami_klein_traits Self;
        typedef Kernel Base;

    public:
        typedef typename Kernel::Side_of_oriented_circle_2 Euclidean_side_of_oriented_circle_2;
        typedef CGAL::Side_of_oriented_circle_2<Self> Side_of_oriented_circle_2;

        explicit Beltrami_klein_traits(const Base &kernel = Base()) : Base(kernel) {}

        Euclidean_side_of_oriented_circle_2 euclidean_side_of_oriented_circle_2_object() const {
            return Base::side_of_oriented_circle_2_object();
        }

        Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const {
            return Side_of_oriented_circle_2(*this);
        }
    };
}

#endif
