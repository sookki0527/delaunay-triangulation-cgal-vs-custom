#ifndef H_TRIANGLE
#define H_TRIANGLE

#include "numeric.h"
#include "vector2.h"


namespace krs {

template<typename T>
struct Triangle
{
	using Type = T;
	using VertexType = Vector2<Type>;
	
	Triangle() = default;
	Triangle(const Triangle&) = default;
	Triangle(Triangle&&) = default;
	Triangle(const VertexType &v1, const VertexType &v2, const VertexType &v3);

	bool circumCircleContains(const VertexType &v) const;

	Triangle &operator=(const Triangle&) = default;
	Triangle &operator=(Triangle&&) = default;
	bool operator ==(const Triangle &t) const;

	template<typename U>
	friend std::ostream &operator <<(std::ostream &str, const Triangle<U> &t);

	const VertexType *a;
	const VertexType *b;
	const VertexType *c;
	bool isBad = false;

	static_assert(std::is_floating_point<Triangle<T>::Type>::value,
		"Type must be floating-point");
};
template<typename T>
    Triangle<T>::Triangle(const VertexType& v1, const VertexType& v2, const VertexType& v3) :
      a(&v1), b(&v2), c(&v3)
    {}
template<typename T>
bool almost_equal(const Triangle<T> &t1, const Triangle<T> &t2)
{
	return	(almost_equal(*t1.a , *t2.a) || almost_equal(*t1.a , *t2.b) || almost_equal(*t1.a , *t2.c)) &&
			(almost_equal(*t1.b , *t2.a) || almost_equal(*t1.b , *t2.b) || almost_equal(*t1.b , *t2.c)) &&
			(almost_equal(*t1.c , *t2.a) || almost_equal(*t1.c , *t2.b) || almost_equal(*t1.c , *t2.c));
}


template<typename T> 
bool Triangle<T>::circumCircleContains(const Vector2<T>& v) const {
    const T ab = a->norm2();
    const T cd = b->norm2();
    const T ef = c->norm2();

    const T ax = a->x;
    const T ay = a->y;
    const T bx = b->x;
    const T by = b->y;
    const T cx = c->x;
    const T cy = c->y;

    const T circum_x = (ab * (cy - by) + cd * (ay - cy) + ef * (by - ay)) /
                       (ax * (cy - by) + bx * (ay - cy) + cx * (by - ay));
    const T circum_y = (ab * (cx - bx) + cd * (ax - cx) + ef * (bx - ax)) /
                       (ay * (cx - bx) + by * (ax - cx) + cy * (bx - ax));

    const Vector2<T> circum(circum_x / 2, circum_y / 2);

    const T circum_radius = a->dist2(circum);
    const T dist = v.dist2(circum);

    return dist <= circum_radius;
}


}

#endif