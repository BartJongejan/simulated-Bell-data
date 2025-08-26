EXECUTABLE = belldata

SRC = belldata.cpp 

belldata: $(SRC) 
	g++ -pedantic -Wall -O3 -DNDEBUG -o $(EXECUTABLE) $(SRC) -lm 


all: belldata

clean:
	rm $(EXECUTABLE)
	rm *.o
