CXX			= 	g++
CXXFLAGS 	=   -I/usr/local/include -I/opt/local/include -O3 -fopenmp
CXXFLAGS 	+= 	-Wall -pedantic -std=c++11
LDFLAGS 	= 	-L/opt/local/lib   -lgmp -lgmpxx -lcrypto -lm -lmpfr -lflint

SOURCES 	= 	$(wildcard *.cpp)
HEADERS 	= 	$(wildcard *.h)
OBJECTS 	= 	$(SOURCES:.cpp=.o)

all: $(OBJECTS) 

$(OBJECTS): Makefile $(HEADERS)

.PHONY: clean
clean:
	rm -f $(OBJECTS)
