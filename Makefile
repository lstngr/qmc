all: qmc

qmc: qmc.o mcstep.o
	g++ -g qmc.o mcstep.o -o qmc

qmc.o: qmc.cpp
	g++ -g -c qmc.cpp

qmcstep.o: mcstep.cpp mcstep.hpp
	g++ -g -c mcstep.cpp mcstep.hpp

clean:
	rm *.o *.gch qmc 
