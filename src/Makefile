CXX			= 	g++
CXXFLAGS 	=   -I/usr/local/include -I/opt/local/include -O3 -fopenmp
CXXFLAGS 	+= 	-Wall -pedantic -std=c++11
LDFLAGS 	= 	-L/opt/local/lib   -lgmp -lgmpxx -lcrypto -lm -lmpfr -lflint

SOURCES 	= 	$(wildcard *.cpp)
HEADERS 	= 	$(wildcard *.h)
OBJECTS 	= 	$(SOURCES:.cpp=.o)

all: $(OBJECTS) pack_into_a_lib 

$(OBJECTS): Makefile $(HEADERS)


pack_into_a_lib:
	ar rs yashe.a  *.o > /dev/null 2>&1


.PHONY: clean
clean:
	rm -f $(OBJECTS) yashe.a
