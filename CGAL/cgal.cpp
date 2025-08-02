#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <GL/freeglut.h>
#include <fstream>
#include "threadsafe_queue.h"
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_2<K>  Triangulation;
typedef Triangulation::Edge_iterator  Edge_iterator;
typedef Triangulation::Point          Point;
int width = 500;
int height = 500;
int frameCount = 0;

int N = 4;
// Global variables for ball state
float ballX = 0.0f, ballY = 0.0f;
float ballVX = 1.5f, ballVY = 0.7f;
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct TriangulationTask {
    std::string input_file;
    std::vector<Point> points;
    Triangulation T;
};

std::vector<TriangulationTask> results;
void triangulate(TriangulationTask& task) {
    task.T.insert(task.points.begin(), task.points.end());
}

void threading(){
    thread_safe_queue<TriangulationTask> task_queue;
    std::mutex results_mutex;

    auto worker = [&]() {
        while (true) {
            TriangulationTask task;
            if (!task_queue.try_pop(task)) break;
            triangulate(task);
            std::lock_guard<std::mutex> lock(results_mutex);
            results.push_back(std::move(task)); 
        }
    };

    for(int i = 0; i < N; i++){
        TriangulationTask task;
        task.input_file = "data/voronoi" + std::to_string(i+1) + ".cin";

        std::ifstream in(task.input_file);
        std::istream_iterator<Point> begin(in), end;
        task.points = std::vector<Point>(begin, end);
        task_queue.push(task);
    }

    std::vector<std::thread> threads;
    for (int i = 0; i < 4; ++i)
        threads.emplace_back(worker);
  
    for (auto& t : threads)
        t.join();
    
}

void camera(){
  
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-1000.0, 1000.0, -1000.0, 1000.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

}

void draw_triangles(const Triangulation &T){
    
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);

    for (auto e = T.finite_edges_begin(); e != T.finite_edges_end(); ++e) {
        auto seg = T.segment(*e);
        glVertex2d(seg.source().x(), seg.source().y());
        glVertex2d(seg.target().x(), seg.target().y());
    }

    glEnd();
}

void idle() {
    glutPostRedisplay(); 
}

void Ball(float x, float y, float radius = 20.0f) {
    glColor3f(1.0, 0.0, 0.0); 
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y); 
    for (int i = 0; i <= 100; ++i) {
        float angle = 2.0f * M_PI * i / 100;
        glVertex2f(x + cos(angle) * radius, y + sin(angle) * radius);
    }
    glEnd();
}

void saveFrame(int frameNumber) {
    int w = width, h = height;
    std::vector<unsigned char> pixels(w * h * 3);
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

    std::ostringstream filename;
    filename << "frames/frame_" << std::setw(4) << std::setfill('0') << frameNumber << ".ppm";

    std::ofstream out(filename.str(), std::ios::binary);
    out << "P6\n" << w << " " << h << "\n255\n";
    out.write(reinterpret_cast<char*>(pixels.data()), pixels.size());
    out.close();
}

void display(){
    glClear(GL_COLOR_BUFFER_BIT);


    ballX += ballVX;
    ballY += ballVY;
    float radius = 20.0f;
    if (ballX > 1000.0f - radius || ballX < -1000.0f + radius) ballVX *= -1.0f;
    if (ballY > 1000.0f - radius || ballY < -1000.0f + radius) ballVY *= -1.0f;



    glViewport(0, 0, width, height); 
    camera();
    Ball(ballX, ballY);

 
    int active_idx = -1;
    if (ballX < 0 && ballY > 0) active_idx = 0; 
    else if (ballX > 0 && ballY > 0) active_idx = 1; 
    else if (ballX < 0 && ballY < 0) active_idx = 2; 
    else if (ballX > 0 && ballY < 0) active_idx = 3; 

   
    if (active_idx == 0) {
        glViewport(0, height/2, width/2, height/2);
        camera();
        draw_triangles(results[0].T);
    }
    if (active_idx == 1) {
        glViewport(width/2, height/2, width/2, height/2);
        camera();
        draw_triangles(results[1].T);
    }
    if (active_idx == 2) {
        glViewport(0, 0, width/2, height/2);
        camera();
        draw_triangles(results[2].T);
    }
    if (active_idx == 3) {
        glViewport(width/2, 0, width/2, height/2);
        camera();
        draw_triangles(results[3].T);
    }
    //saveFrame(frameCount++);
    glutSwapBuffers();
}


void reshape(int w, int h)
{
    width = w;
    height = h;
    glViewport(0, 0, w, h);
}



int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(500, 500);
    glutCreateWindow("CGAL Delaunay Triangulation");
   std::cout << "ballVX = " << ballVX << ", ballVY = " << ballVY << std::endl;

    glClearColor(0.0, 0.0, 0.0, 0.0);

    threading();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
	  glutIdleFunc(idle);
    glutMainLoop();
}