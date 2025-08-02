CC=g++
SOURCEDIR = delaunay
INCLUDEDIR = include
SOURCES = $(wildcard $(SOURCEDIR)/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
SFMLDIR = C:/Users/sookk/SFML-2.6.2
CFLAGS = -c -Wall -std=c++17 -I$(INCLUDEDIR)  -I$(SFMLDIR)/include
LDFLAGS= -L$(SFMLDIR)/lib
LIBFLAGS=-lsfml-graphics -lsfml-window -lsfml-system
EXECUTABLE=main
BINDIR=/usr/bin

all: $(SOURCES) $(EXECUTABLE)
        
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ $(LIBFLAGS)

clean:
	rm *.o $(EXECUTABLE)

install:
	#install -s $(EXECUTABLE) $(BINDIR)
	sudo cp -u $(EXECUTABLE) $(BINDIR)
uninstall:
	sudo rm $(BINDIR)/$(EXECUTABLE)