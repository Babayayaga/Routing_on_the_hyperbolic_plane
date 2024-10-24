#ifndef MYTRAITS_H
#define MYTRAITS_H

/*
 * Traits for a Constrained_Delaunay_Triangulation_2 in the Poincare disk model
 * we inherit from the predefined Exact_predicates_inexact_constructions_kernel and overwrite
 * the Orientation_2 predicate and the Construct_segment_2 construction.
 * To determine the orientation of three points, we need to know the circle perpendicular to the unit disk,
 * this circle needs to be constructed.
 * So for our Orientation_2 predicate to be exact, we would need to inherit from Exact_predicates_exact_constructions_kernel
 * That's possible but will make the program slower. For exact purposes we implemented the BeltramiKleinTraits.
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Hyperbolic_triangulation_2/internal/Exact_complex.h>
#include <CGAL/Origin.h>
#include <boost/variant.hpp>

namespace CGAL {
    template<typename Traits>
    class Circular_arc_2 {
        typedef typename Traits::FT FT;
        typedef Exact_complex<FT> Cplx;
        typedef typename Traits::Point_2 Point_2;
        typedef typename Traits::Circle_2 Circle;
        typedef typename Traits::Orientation_2 Orientation_2;

        Circle _c;
        Point_2 _s, _t;
        const Traits &_gt;

    public:
        Circular_arc_2(const Circle &c,
                       const Point_2 &source, const Point_2 &target,
                       const Traits &gt = Traits())
            : _c(c), _s(source), _t(target), _gt(gt) {
        }

        Circle supporting_circle() const { return _c; }
        Point_2 source() const { return _s; }
        Point_2 target() const { return _t; }
        FT squared_radius() const { return _c.squared_radius(); }
        Point_2 center() const { return _c.center(); }
        Bbox_2 bbox(void) const { return _gt.construct_bbox_2_object()(*this); }
    };

    template<typename Traits>
    class Orientation_2 {
        typedef typename Traits::FT FT;
        typedef typename Traits::Point_2 Point_2;
        typedef typename Traits::Hyperbolic_segment_2 Hyperbolic_segment_2;
        typedef typename Traits::Euclidean_segment_2 Euclidean_segment_2;
        typedef typename Traits::Circle_2 Circle_2;
        typedef typename Traits::Circular_arc_2 Circular_arc_2;
        typedef typename Traits::Euclidean_Orientation_2 Euclidean_Orientation_2;
        typedef typename Traits::Weighted_point_2 Weighted_point_2;
        typedef Exact_complex<FT> Cplx;

    public:
        explicit Orientation_2(const Traits &gt = Traits()) : _gt(gt) {
        }

        Orientation operator()(const Point_2 p, const Point_2 q, const Point_2 query) {
            if (p != q && q != query && query != p) {
                Circle_2 poincare(Point_2(FT(0), FT(0)), FT(1));
                Point_2 O(FT(0), FT(0));
                Orientation bb = _gt.euclidean_orientation_2_object()(p, q, O);
                if (bb == COLLINEAR) {
                    return _gt.euclidean_orientation_2_object()(p, q, query);
                }

                Cplx p1(p.x(), p.y());
                Cplx q1(q.x(), q.y());
                Cplx O1(0, 0);
                Cplx inv;
                if (p1 == O1) {
                    inv = q1.invert_in_unit_circle();
                } else {
                    inv = p1.invert_in_unit_circle();
                }
                Point_2 ip(inv.real(), inv.imag());
                Circle_2 _c = Circle_2(p, q, ip);

                Orientation eo = _gt.euclidean_orientation_2_object()(p, _c.center(), q);
                Bounded_side bs = _c.bounded_side(query);

                if (bs == ON_BOUNDED_SIDE) {
                    if (eo == LEFT_TURN) {
                        return RIGHT_TURN;
                    }
                    return LEFT_TURN;
                }
                if (bs == ON_UNBOUNDED_SIDE) {
                    if (eo == LEFT_TURN) {
                        return LEFT_TURN;
                    }
                    return RIGHT_TURN;
                }
                return COLLINEAR;
            }
            return COLLINEAR;
        }

    private:
        const Traits &_gt;
    };

    template<typename Traits>
    class Construct_hyperbolic_segment_2 {
        typedef typename Traits::FT FT;
        typedef typename Traits::Point_2 Point_2;
        typedef typename Traits::Weighted_point_2 Weighted_point_2;
        typedef typename Traits::Euclidean_segment_2 Euclidean_segment_2;
        typedef typename Traits::Hyperbolic_segment_2 Hyperbolic_segment_2;
        typedef typename Traits::Circle_2 Circle_2;
        typedef typename Traits::Circular_arc_2 Circular_arc_2;
        typedef typename Traits::Construct_weighted_circumcenter_2 Construct_weighted_circumcenter_2;

    public:
        explicit Construct_hyperbolic_segment_2(const Traits &gt = Traits()) : _gt(gt) {
        }

        typedef Hyperbolic_segment_2 result_type;

        Hyperbolic_segment_2 operator()(const Point_2 &p,
                                        const Point_2 &q) const {
            Origin o;
            if (_gt.collinear_2_object()(p, q, Point_2(o)))
                return Euclidean_segment_2(p, q);
            Weighted_point_2 wp(p);
            Weighted_point_2 wq(q);
            Weighted_point_2 wo(Point_2(o), FT(1));
            Point_2 center = _gt.construct_weighted_circumcenter_2_object()(wp, wo, wq);
            FT sq_radius = _gt.compute_squared_distance_2_object()(p, center);
            Circle_2 circle = _gt.construct_circle_2_object()(center, sq_radius);
            if (_gt.euclidean_orientation_2_object()(p, q, center) == LEFT_TURN)
                return Circular_arc_2(circle, p, q);
            return Circular_arc_2(circle, q, p);
        }

    private:
        const Traits &_gt;
    };

    template<typename Kernel = Exact_predicates_inexact_constructions_kernel>
    class Poincare_disk_traits : public Kernel {
        typedef Poincare_disk_traits Self;
        typedef Kernel Base;

    public:
        typedef typename Kernel::FT FT;
        typedef typename Kernel::Point_2 Point_2;
        typedef typename Kernel::Circle_2 Circle_2;
        typedef Circular_arc_2<Self> Circular_arc_2;
        typedef typename Kernel::Segment_2 Euclidean_segment_2;
        typedef boost::variant<Circular_arc_2, Euclidean_segment_2> Hyperbolic_segment_2;
        typedef typename Kernel::Triangle_2 Triangle_2;
        typedef Hyperbolic_segment_2 Segment_2;

        typedef Construct_hyperbolic_segment_2<Self> Construct_segment_2;
        typedef typename Kernel::Orientation_2 Euclidean_Orientation_2;
        typedef Orientation_2<Self> Orientation_2;

        explicit Poincare_disk_traits(const Base &kernel = Base()) : Base(kernel) {
        }

        Euclidean_Orientation_2 euclidean_orientation_2_object() const {
            return Base::orientation_2_object();
        }

        Orientation_2 orientation_2_object() const {
            return Orientation_2(*this);
        }

        Construct_segment_2
        construct_segment_2_object() const {
            return Construct_segment_2(*this);
        }
    };
}
#endif //MYTRAITS_H
