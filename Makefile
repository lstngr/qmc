all: qmc

qmc: qmc.o qmcstep.o
	g++ -g qmc.o qmcstep.o -o qmc

qmc.o: qmc.cpp
	g++ -g -c qmc.cpp

qmcstep.o: qmcstep.cpp qmcstep.hpp
	g++ -g -c qmcstep.cpp qmcstep.hpp

clean:
	rm *.o *.gch qmc 
