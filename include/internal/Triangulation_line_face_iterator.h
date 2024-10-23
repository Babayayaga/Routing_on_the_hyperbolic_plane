#ifndef TRIANGULATION_LINE_FACE_ITERATOR_H
#define TRIANGULATION_LINE_FACE_ITERATOR_H

/*
 * similiar to the CGAL Triangulation_line_face_iterator,
 * but fitted to our purposes, we want to use the line_walk method
 * and already know at which vertex we're starting at.
 */

#include <CGAL/Triangulation_utils_2.h>

namespace CGAL {
    template<class T>
    class Triangulation_line_face_iterator : public Triangulation_cw_ccw_2 {
        typedef Triangulation_line_face_iterator Line_face_iterator;
        typedef typename T::Geom_traits Gt;
        typedef typename T::Triangulation_data_structure Tds;
        typedef typename Tds::Vertex Vertex;
        typedef typename Tds::Face Face;
        typedef typename Tds::Edge Edge;
        typedef typename Tds::Vertex_handle Vertex_handle;
        typedef typename Tds::Face_handle Face_handle;
        typedef typename Tds::Face_circulator Face_circulator;
        typedef typename T::Point_2 Point_2;

        enum State {
            undefined = -1,
            vertex_vertex,
            vertex_edge,
            edge_vertex,
            edge_edge
        };

    public:
        Line_face_iterator&   operator++() ;
        void increment();
        Face*                   operator->() {return &*pos;}
        Face&                   operator*() { return *pos;}
        Face_handle             handle() {return pos;}
        bool  operator==(const Face_handle& fh) const { return fh == pos; }
        bool  operator!=(const Face_handle& fh) const { return fh != pos; }
        Face_handle pos;
        const T *t;
        State s;
        int i;
        Point_2 p, q;
        Triangulation_line_face_iterator(Vertex_handle v, const T *tr, const Point_2 &dir);
    };

    //setup
    template<class T>
    Triangulation_line_face_iterator<T>::
    Triangulation_line_face_iterator(Vertex_handle v, const T *tr, const Point_2 &dir)
        : pos(), t(tr), s(undefined) {
        p = v->point();
        q = dir;

        Face_circulator fc = t->incident_faces(v);
        Face_circulator done(fc);
        int ic = fc->index(v);
        Vertex_handle vt = fc->vertex(cw(ic));
        while (t->is_infinite(vt) ||
               t->orientation(p, q, vt->point()) != LEFT_TURN) {
            ++fc;
            ic = fc->index(v);
            vt = fc->vertex(cw(ic));
            if (fc == done) {
                t = nullptr;
                return;
            }
        }

        Vertex_handle vr = fc->vertex(ccw(ic));
        Orientation pqr = RIGHT_TURN;
        while ((!t->is_infinite(vr)) &&
               (pqr = t->orientation(p, q, vr->point())) == LEFT_TURN) {
            --fc;
            ic = fc->index(v);
            vr = fc->vertex(ccw(ic));
        }

        // vr can be infinite or finite.
        // If finite [pqr] is COLLINEAR or RIGHT_TURN
        // reset vt and conclude.  vt is still finite and [pqt] still LEFT_TURN
        ic = fc->index(v);
        vt = fc->vertex(cw(ic));
        /*CGAL_assertion(t->orientation(p, q, vt->point()) ==
                       LEFT_TURN);*/
        if (t->is_infinite(vr)) {
            --fc;
            ic = fc->index(v);
            vr = fc->vertex(ccw(ic));
            pqr = t->orientation(p, q, vr->point());
            switch (pqr) {
                case RIGHT_TURN:
                case COLLINEAR:
                    ++fc;
                    ic = fc->index(t->infinite_vertex());
                    pos = fc;
                    s = vertex_vertex;
                    i = ic;
                    break;
                case LEFT_TURN:
                    t = nullptr;
                    break;
            }
        } else if (pqr == COLLINEAR) {
            pos = fc;
            s = vertex_vertex;
            i = ccw(ic);
        } else {
            // pqr==RIGHT_TURN
            pos = fc;
            s = vertex_edge;
            i = ic;
        }
    }

    template <class T>
    Triangulation_line_face_iterator<T>&Triangulation_line_face_iterator<T>::operator++() {
        CGAL_precondition(pos != Face_handle());
        increment();
        return *this;
    }

    template <class T> void
    Triangulation_line_face_iterator<T>::increment() {
        CGAL_precondition(pos != Face_handle());
        if(s == vertex_vertex || s == edge_vertex) {
            Orientation o;
            do{
                Face_handle n = pos->neighbor(cw(i));
                i = n->index(pos);
                pos = n;
                if (pos->vertex(i) == t->infinite_vertex()){
                    o = COLLINEAR;
                    i = cw(i);
                    break;
                }
                o = t->orientation(p, q, pos->vertex(i)->point());
                i = cw(i);
            }while(o == LEFT_TURN);

            if(o == COLLINEAR) {
                s = vertex_vertex;
                i = ccw(i);
            }
            else {
                s = vertex_edge;
            }
        }
        else {
            Face_handle n = pos->neighbor(i);
            int ni = n->index(pos);
            pos = n ;
            Orientation o = t->is_infinite(pos->vertex(ni)) ?
              COLLINEAR :
              t->orientation(p,q,pos->vertex(ni)->point());

            switch(o){
                case LEFT_TURN:
                    s = edge_edge;
                i = ccw(ni);
                break;
                case RIGHT_TURN:
                    s = edge_edge;
                i = cw(ni);
                break;
                default:
                    s = edge_vertex;
                i = ni;
            }
        }
    }

    template < class Triangulation >bool operator==(typename Triangulation::Triangulation_data_structure::Face_handle fh,
           Triangulation_line_face_iterator<Triangulation> fc) {
        return fc==fh;
    }

    template < class Triangulation > bool operator!=(typename Triangulation::Triangulation_data_structure::Face_handle fh,
               Triangulation_line_face_iterator<Triangulation> fc) {
        return fc!=fh;
    }
}

#endif
