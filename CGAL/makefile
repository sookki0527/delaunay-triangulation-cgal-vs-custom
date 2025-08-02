LIBS = -L/mingw64/lib -lfreeglut -lglew32 -lglu32 -lopengl32 \
       -lboost_thread-mt -lmpfr -lgmp -lpthread
INCLUDES = -I/mingw64/include

cgal: cgal.o
	$(CXX) -o $@ $^ $(LIBS)

cgal.o: cgal.cpp
	$(CXX) -c -O3 -Wall -std=c++17 $(INCLUDES) -DUSEGLEW -o $@ $<
