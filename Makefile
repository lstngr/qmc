all: qmc

qmc: qmc.o qmcstep.o
	g++ -g qmc.o mcstep.o -o qmc

qmc.o: qmc.cpp
	g++ -g -c qmc.cpp

mcstep.o: mcstep.cpp mcstep.hpp
	g++ -g -c mcstep.cpp mcstep.hpp

ssedata.o: ssedata.cpp ssedata.hpp
	g++ -g -c ssedata.cpp ssedata.hpp

clean:
	rm -f *.o *.gch qmc 
