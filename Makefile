CXX = g++
CXXFLAGS = -g -O0
OBJECTS = qmc.o mcstep.o ssedata.o lattice.o observable.o iomanager.o debugconf.hpp

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

iomanager.o: iomanager.cpp iomanager.hpp
	$(CXX) $(CXXFLAGS) -c iomanager.cpp -o $@

clean:
	-rm -f *.o *.gch qmc

cleandocs:
	-rm -fr html/ latex/
