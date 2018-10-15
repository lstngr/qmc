CXX = g++
CXXFLAGS = -g -O3
OBJECTS = qmc.o mcstep.o ssedata.o

all: qmc

qmc: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

qmc.o: qmc.cpp
	$(CXX) $(CXXFLAGS) -c qmc.cpp -o $@

mcstep.o: mcstep.cpp mcstep.hpp
	$(CXX) $(CXXFLAGS) -c mcstep.cpp -o $@

ssedata.o: ssedata.cpp ssedata.hpp
	$(CXX) $(CXXFLAGS) -c ssedata.cpp -o $@

clean:
	-rm -fr *.o *.gch qmc html/ latex/
