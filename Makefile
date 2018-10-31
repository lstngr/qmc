CXX = g++
CXXFLAGS = -g -O0
OBJECTS = qmc.o mcstep.o ssedata.o lattice.o observable.o omanager.o

all: qmc

qmc: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

qmc.o: qmc.cpp
	$(CXX) $(CXXFLAGS) -c qmc.cpp -o $@

mcstep.o: mcstep.cpp mcstep.hpp
	$(CXX) $(CXXFLAGS) -c mcstep.cpp -o $@

ssedata.o: ssedata.cpp ssedata.hpp
	$(CXX) $(CXXFLAGS) -c ssedata.cpp -o $@

lattice.o: lattice.cpp lattice.hpp
	$(CXX) $(CXXFLAGS) -c lattice.cpp -o $@

observable.o: observable.cpp observable.hpp
	$(CXX) $(CXXFLAGS) -c observable.cpp -o $@

omanager.o: omanager.cpp omanager.hpp
	$(CXX) $(CXXFLAGS) -c omanager.cpp -o $@

clean:
	-rm -fr *.o *.gch qmc html/ latex/
