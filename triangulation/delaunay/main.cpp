
#include "Delaunay.h"
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <array>
#include <random>

#include <SFML/Graphics.hpp>

#include "vector2.h"
#include "triangle.h"

int main() {


    std::default_random_engine gen(std::random_device{}());
    std::vector<krs::Vector2<double>> points;

    std::uniform_real_distribution<> dist_x(-300, 300);  
    std::uniform_real_distribution<> dist_y(-300, 300);

    double center_x = 400;
    double center_y = 300;

    for (int i = 0; i < 48; ++i) {
        double x = center_x + dist_x(gen);
        double y = center_y + dist_y(gen);
        points.push_back(krs::Vector2<double>{x, y});
    }

    std::sort(points.begin(), points.end(), [](krs::Vector2<double> a, krs::Vector2<double> b) {
    if (a.x == b.x) return a.y < b.y;
            return a.x < b.x;
        });

    for (auto& a : points) {
        std::cout << "[" << a.x << ", " << a.y << "] ";
    }


    krs::delaunay<double> delaunay;


    delaunay.slicingVector(points);
    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> edges = delaunay.getEdges();
    const std::vector<std::vector<std::pair<double, double>>> triangles = delaunay.getTriangles();


    sf::RenderWindow window(sf::VideoMode(800, 600), "Delaunay triangulation");
    window.setFramerateLimit(1);

    // Transform each points of each vector as a rectangle
    for (const auto p : points) {
        sf::RectangleShape s{ sf::Vector2f(4, 4) };
        s.setPosition(static_cast<double>(p.x), static_cast<double>(p.y));
        window.draw(s);
    }

    

    for (const auto& tri : triangles) {
        if (tri.size() == 2) {
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(tri[0].first, tri[0].second)),
                sf::Vertex(sf::Vector2f(tri[1].first, tri[1].second))
            };
            window.draw(line, 2, sf::Lines);
        }
        else if (tri.size() == 3) {
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(tri[0].first, tri[0].second)),
                sf::Vertex(sf::Vector2f(tri[1].first, tri[1].second)),
                sf::Vertex(sf::Vector2f(tri[1].first, tri[1].second)),
                sf::Vertex(sf::Vector2f(tri[2].first, tri[2].second)),
                sf::Vertex(sf::Vector2f(tri[2].first, tri[2].second)),
                sf::Vertex(sf::Vector2f(tri[0].first, tri[0].second))
            };
            window.draw(line, 6, sf::Lines);
    }
}   


    window.display();

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
    }
    return 0;
}