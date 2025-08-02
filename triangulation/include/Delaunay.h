#include <vector>
#ifndef H_DELAUNAY
#define H_DELAUNAY
#include "vector2.h"
#include "triangle.h"
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define _USE_MATH_DEFINES
#include <cmath>
#include <optional>
#include <limits>
constexpr double INF = 1e18; 
namespace krs{

template <typename T>
class delaunay {
    using VertexType = Vector2<T>;

    using Type = T;
    using TriangleType = Triangle<T>;
    using VertexWithAngle = std::pair<double, VertexType>;  // (angle, (x, y))
    //std::vector<Edge<double>> _edges;
    std::vector<Vector2<T>> _vertices;
    std::vector<std::pair<std::pair<T, T>, std::pair<T, T> >> _edges;
    std::vector<std::vector<std::pair<T, T>>> _triangles;
    std::vector<std::vector<std::pair<T, T>>> res;
    
public:

    delaunay() = default;
    delaunay(const delaunay&) = default;
    delaunay(delaunay&&) = default;


    std::optional<std::pair<std::pair<T, T>, std::pair<T, T>>>
    potentials(std::vector<Vector2<T>>& vertices1,
               std::vector<Vector2<T>>& vertices2, bool flipped);
    const std::vector<std::pair<std::pair<T, T>, std::pair<T, T>>> getEdges();
    const std::vector<std::vector<std::pair<T, T>>> getTriangles();
    void createEdges(std::vector<VertexType>& vertices);
    void slicingVector(std::vector<VertexType>& vertices);
    bool isTriangle(std::vector<std::pair<std::pair<T, T>, std::pair<T, T> >> edges, std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> p3);
    bool onsegment(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> p3) const;
    const int orientation(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> p3) const;
    bool intersect(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> q1, std::pair<T, T> q2) const;
    bool delaunay_intersect(std::pair<T, T> np1, std::pair<T, T> np2) const;

    delaunay& operator=(delaunay&&) = default;

};
    template<typename T>
    bool almost_equal2(std::pair<T, T>& p1, std::pair<T, T>& p2)
    {
        return almost_equal(p1.first, p2.first) && almost_equal(p1.second, p2.second);
    }
     template<typename T>
    bool containsonePoint(std::pair<std::pair<T, T>, std::pair<T, T> > s, std::pair<T, T> p1, std::pair<T, T> p2)
    {

        return (almost_equal2(s.first, p1) || almost_equal2(s.second, p2)) || ((almost_equal2(s.first, p2) || almost_equal2(s.second, p1)));

    }
    template<typename T>
    bool containstwoPoints(std::pair<std::pair<T, T>, std::pair<T, T> > s, std::pair<T, T> p1, std::pair<T, T> p2)
      {
          return (almost_equal2(s.first, p1) || almost_equal2(s.second, p2)) && ((almost_equal2(s.first, p2) || almost_equal2(s.second, p1)));
      }
    template<typename T>
    bool containstwoPoints(std::pair<std::pair<T, T>, std::pair<T, T>> s,
                        std::pair<T, T> p1,
                        Vector2<T> p2)
    {
        std::pair<T, T> pp2 = { p2.x, p2.y };
        return containstwoPoints(s, p1, pp2);  // 기존 함수로 위임
    }

    template<typename T>
    bool containstwoPoints2(std::pair<std::pair<T, T>, std::pair<T, T> > s, std::pair<T, T> p1, std::pair<T, T> p2)
    {
        return (almost_equal2(s.first, p1)&& almost_equal2(s.second, p2)) || ((almost_equal2(s.first, p2) && almost_equal2(s.second, p1)));
    }

    template<typename T>
   const std::vector<std::pair<std::pair<T, T>, std::pair<T, T> >>
       delaunay<T>::getEdges()
   {
       return _edges;
   }

   template<typename T>
    const std::vector<std::vector<std::pair<T, T>>> delaunay<T>:: getTriangles()
    {
        return _triangles;
    }


   template<typename T>
   bool delaunay<T>::onsegment(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> p3) const
   {
       if (p2.first <= std::max(p1.first, p3.first) && p2.first >= std::min(p1.first, p3.first)
           && p2.second <= std::max(p1.second, p3.second) && p2.second >= std::min(p1.second, p3.second) ){
           return true;
       }
       return false;
   }

   template<typename T>
   const int delaunay<T>::orientation(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> p3) const {

       T val = (p2.second - p1.second) * (p3.first - p2.first) -
           (p2.first - p1.first) * (p3.second - p2.second);
       if (val == 0) {
           return 0;
       }
       return (val > 0) ? 1 : 2;
   }


   template<typename T>
   bool delaunay<T>::intersect(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> q1, std::pair<T, T> q2) const
   {
       T o1 = orientation(p1, q1, p2);
       T o2 = orientation(p1, q1, q2);
       T o3 = orientation(p2, q2, p1);
       T o4 = orientation(p2, q2, q1);

       if (o1 != o2 && o3 != o4)
           return true;

       if (o1 == 0 && onsegment(p1, p2, q1)) return true;
       if (o2 == 0 && onsegment(p1, q2, q1)) return true;
       if (o3 == 0 && onsegment(p2, p1, q2)) return true;
       if (o4 == 0 && onsegment(p2, q1, q2)) return true;

       return false; 

   }

     template<typename T>
   bool delaunay<T>::delaunay_intersect(std::pair<T, T> np1, std::pair<T, T> np2) const {

       int c = 0;
       for (auto& e : _edges) {
           if (intersect(np1, e.first, np2, e.second) && !containsonePoint(e, np1, np2)){ 
               return false;
          }

       }return true;
   }
template <typename T>
void delaunay<T>::createEdges(std::vector<Vector2<T>>& vertices) {

    
    size_t s1 = vertices.size() - 1;
    
    if (s1 == 1) {

        T firstx = vertices[0].x;
        T firsty = vertices[0].y;

   
        std::pair<T, T> pair1 = { firstx, firsty };
        std::pair<T, T> pair2 = { vertices[s1].x , vertices[s1].y };
       
        _edges.push_back({ pair1, pair2 });


    }
    else if (s1 == 2) {


        T firstx = vertices[0].x;
        T firsty = vertices[0].y;

      
        const std::pair<double, double> pair1 = { firstx, firsty };
        const std::pair<double, double> pair2 = { vertices[s1-1].x , vertices[s1-1].y };
        const std::pair<double, double> pair3 = { vertices[s1].x , vertices[s1].y };
        
        _edges.push_back({ pair1, pair2 });
        _edges.push_back({ pair2, pair3 });
        _edges.push_back({ pair3, pair1 });

    };

}
template<typename T>
double euclidean_distance(std::pair<T, T> a, std::pair<T, T> b) {
    T dx = a.first - b.first;
    T dy = a.second- b.second;
    return std::sqrt(static_cast<double>(dx * dx + dy * dy));
}

   template<typename T>
   void delaunay<T>::slicingVector(std::vector<Vector2<T>>& vertices)
   {
       if (vertices.size() <= 3) {
            createEdges(vertices);
            if (vertices.size() == 3) {
                auto a = vertices[0];
                auto b = vertices[1];
                auto c = vertices[2];

                T cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
                if (cross < 0) std::swap(b, c); 

                _triangles.push_back({
                    {a.x, a.y},
                    {b.x, b.y},
                    {c.x, c.y}
                });
           } else if (vertices.size() == 2) {
                _triangles.push_back({
                    {vertices[0].x, vertices[0].y},
                    {vertices[1].x, vertices[1].y}
                });
            }
            return;
        }
  
        size_t mid = vertices.size() / 2;
        std::vector<Vector2<T>> left(vertices.begin(), vertices.begin() + mid);
        std::vector<Vector2<T>> right(vertices.begin() + mid, vertices.end());
        std::cout<< vertices.size() << std::endl;
        
        slicingVector(left);
        slicingVector(right);

        std::optional<std::pair<std::pair<T, T>, std::pair<T, T>>> new_LR = potentials(left, right, false);
        std::optional<std::pair<std::pair<T, T>, std::pair<T, T>>> new_LR2 = potentials(left, right, true);
    
   }



    template <typename T> std::optional<std::pair<std::pair<T, T>, std::pair<T, T>>>
    delaunay<T>::potentials(std::vector<Vector2<T>>& vertices1, std::vector<Vector2<T>>& vertices2, bool flipped) {

        if(!flipped){
            std::sort(vertices1.begin(), vertices1.end(), [](Vector2<T> a, Vector2<T> b) { return a.y > b.y; });
            std::sort(vertices2.begin(), vertices2.end(), [](Vector2<T> a, Vector2<T> b) { return a.y > b.y; });
        }
        else{
            std::sort(vertices1.begin(), vertices1.end(), [](Vector2<T> a, Vector2<T> b) { return a.y < b.y; });
            std::sort(vertices2.begin(), vertices2.end(), [](Vector2<T> a, Vector2<T> b) { return a.y < b.y; });
        }
            
            const VertexType p1(vertices1[0].x, vertices1[0].y); 
            const VertexType p2(vertices2[0].x, vertices2[0].y);
            const std::pair<T, T> np1 = { vertices1[0].x, vertices1[0].y };
            const std::pair<T, T> np2 = { vertices2[0].x, vertices2[0].y };

            std::optional<std::pair<std::pair<T, T>, std::pair<T, T>>> new_LR = std::make_pair(np1, np2);
            
            //bottom most line 
            std::pair<std::pair<T, T>, std::pair<T, T> > enp = { np1, np2 };
            if (delaunay_intersect(np1, np2)) {
                _edges.push_back({ np1, np2 });
                new_LR = {np1, np2};
            }
            std::optional<std::pair<std::pair<T, T>, std::pair<T, T>>> prev_LR;
            while(new_LR.has_value()){

                auto [np1, np2] = new_LR.value();
                

                ///******************  RR edge  *************************/

                T start = 0;
                
                std::vector<VertexWithAngle> angle_sorted_vertices_r;

                // Vertices in Right Side
                for (auto& k : vertices2)
                {
                    if (k == p2) continue;  
                    
                    // base LR-edge direction vector (p2 → p1)
                    T vx1 = p1.x - p2.x;
                    T vy1 = p1.y - p2.y;

                    // vector 2: p2 → k direction
                    T vx2 = k.x - p2.x;
                    T vy2 = k.y - p2.y;

                    // dot product & norm
                    double dot = vx1 * vx2 + vy1 * vy2;
                    double ma = std::sqrt(vx1 * vx1 + vy1 * vy1);
                    double mb = std::sqrt(vx2 * vx2 + vy2 * vy2);

                    double cos_theta = std::clamp(dot / (ma * mb), -1.0, 1.0);
                    double angle = std::acos(cos_theta) * 180.0 / M_PI;

                    const VertexType k1(k.x, k.y);
                    angle_sorted_vertices_r.emplace_back(angle, k1);
                

                }
                std::sort(angle_sorted_vertices_r.begin(), angle_sorted_vertices_r.end(),
                        [](const VertexWithAngle& a, const VertexWithAngle& b) {
                            return a.first < b.first;
                        });
    

                std::optional<VertexType> next_lr_candidate_r;
                for(auto & candidate: angle_sorted_vertices_r){
                
                    const TriangleType t(p1, p2, candidate.second);
                    start++;
                    
                        bool has_inner_point = false;
                        for (auto it = angle_sorted_vertices_r.begin() + start; it != angle_sorted_vertices_r.end(); it++) {

                                if (t.circumCircleContains(it->second)) {
                                    has_inner_point = true;
                                    std::pair<double, double> c = {candidate.second.x, candidate.second.y};
                                
                                    _edges.erase(std::remove_if(begin(_edges), end(_edges), [np2, c]
                                    (std::pair<std::pair<T, T>, std::pair<T, T> >& j) {
                                            return containstwoPoints(j, np2, c); }), end(_edges));
                                    break;
                                }
                        }

                        if (!has_inner_point) {
                            next_lr_candidate_r = candidate.second;
                            break;
                        }
                    
                } 
                
                
                ///******************  LL edge  *************************/
                start = 0;
        
                std::vector<VertexWithAngle> angle_sorted_vertices_l;

                for (auto& k : vertices1)
                {
                    if (k == p1) continue;  
                    
                    // base LR-edge direction vector (p1 → p2)
                    T vx1 = p2.x - p1.x;
                    T vy1 = p2.y - p1.y;

                    // vector 2: p1 → k direction
                    T vx2 = k.x - p1.x;
                    T vy2 = k.y - p1.y;

                    // dot product & norm
                    double dot = vx1 * vx2 + vy1 * vy2;
                    double ma = std::sqrt(vx1 * vx1 + vy1 * vy1);
                    double mb = std::sqrt(vx2 * vx2 + vy2 * vy2);

                    double cos_theta = std::clamp(dot / (ma * mb), -1.0, 1.0);
                    double angle = std::acos(cos_theta) * 180.0 / M_PI;
                    const VertexType k2(k.x, k.y);
                angle_sorted_vertices_l.emplace_back(angle, k2);
                
                }

                std::sort(angle_sorted_vertices_l.begin(), angle_sorted_vertices_l.end(),
                        [](const VertexWithAngle& a, const VertexWithAngle& b) {
                            return a.first < b.first;
                        });

                        
                std::optional<VertexType> next_lr_candidate_l;
                for(auto & candidate: angle_sorted_vertices_l){
                    start++;
                    
                    const TriangleType t(p1, p2, candidate.second);

                    
                        bool has_inner_point = false;
                        for (auto it = angle_sorted_vertices_l.begin() + start; it != angle_sorted_vertices_l.end(); it++) {

                                if (t.circumCircleContains(it->second)) {
                                    has_inner_point = true;
                                    std::pair<double, double> c = {candidate.second.x, candidate.second.y};
                                    _edges.erase(std::remove_if(begin(_edges), end(_edges), [np2, c]
                                    (std::pair<std::pair<T, T>, std::pair<T, T> >& j) {
                                            return containstwoPoints(j, np2, c); }), end(_edges));
                                    break;
                                
                                }
                        }

                        if (!has_inner_point) {
                            next_lr_candidate_l = candidate.second;
                            break;
                        }
                    
                } 

           
            if (next_lr_candidate_l.has_value() && next_lr_candidate_r.has_value()) {

                    const VertexType& c1 = next_lr_candidate_l.value();  
                    const VertexType& c2 = next_lr_candidate_r.value();  
                
                    const TriangleType t1(p1, p2, c1);

                    if (t1.circumCircleContains(c2)) {
        
                        std::pair<T, T> next_p1 = { c1.x, c1.y };
                        if (delaunay_intersect(np2, next_p1)) {
                            _edges.push_back({ np2, next_p1 });
                            
                            _triangles.push_back({np1, np2, next_p1});
                            new_LR = {np2, next_p1};
                        }
                    } else {
                    
                        std::pair<T, T> next_p2 = { c2.x, c2.y };
                        if (delaunay_intersect(np1, next_p2)) {
                            _edges.push_back({ np1, next_p2 });
                            _triangles.push_back({np1, np2, next_p2});
                            new_LR = {np1, next_p2};
                        }
                    }
                }

                else if (next_lr_candidate_r.has_value() && !next_lr_candidate_l.has_value()) {
                
                    const VertexType c2 = next_lr_candidate_r.value();
                    const std::pair<T, T> next_p2 = { c2.x, c2.y };

                    if (delaunay_intersect(np1, next_p2)) {
                        _edges.push_back({ np1, next_p2 });
                        _triangles.push_back({np1, np2, next_p2});
                        new_LR = {np1, next_p2};
                    }

                } else if (next_lr_candidate_l.has_value() && !next_lr_candidate_r.has_value()) {
                
                    const VertexType c1 = next_lr_candidate_l.value();
                    const std::pair<T, T> next_p1 = { c1.x, c1.y };

                    if (delaunay_intersect(np2, next_p1)) {
                        _edges.push_back({ np2, next_p1});
                        _triangles.push_back({np1, np2, next_p1});
                        new_LR = {np2, next_p1};
                    }
                } else {
                    new_LR = std::nullopt;  
                }

                
                if (new_LR.has_value()) {
                    if (prev_LR.has_value() && prev_LR.value() == new_LR.value()) {
                        std::cout << "LR edge unchanged. Breaking loop.\n";
                        break;
                    }
                    auto [np1, np2] = new_LR.value();
                    std::cout << "(" << np1.first << ", " << np1.second << ") , ("
                            << np2.first << ", " << np2.second << ")" << std::endl;
                    prev_LR = new_LR;
                }

            
            }
            return prev_LR;

        
        };
}

#endif