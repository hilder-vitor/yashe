#   make [all]  - makes everything.
#   make TARGET - makes the given target.
#   make clean  - removes all files generated by make.

# Please tweak the following variable definitions as needed by your
# project, except GOOGLE_HEADERS, which you can use in your own targets
# but shouldn't modify.

# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GOOGLE_DIR = googletest

SRC_DIR = ../src

INCLUDES = -I/home/vitor/mylibs/gmp/include -I/home/vitor/mylibs/mpfr/include  -I/home/vitor/mylibs/flint/include -I$(SRC_DIR)
LDFLAGS 	= 	-L/opt/local/lib   -L/home/vitor/mylibs/flint/lib/ -L/home/vitor/mylibs/gmp/lib/ -L/home/vitor/mylibs/mpfr/lib/ -lgmp -lgmpxx -lcrypto -lm -lmpfr -lflint -pthread


YASHE_OBJ=$(SRC_DIR)/Plaintext.o $(SRC_DIR)/RealNumberPlaintext.o $(SRC_DIR)/Ciphertext.o $(SRC_DIR)/RealNumberCiphertext.o $(SRC_DIR)/Sampler.o $(SRC_DIR)/Entropy.o $(SRC_DIR)/Yashe.o $(SRC_DIR)/timing.o

CPPFLAGS += -isystem $(GOOGLE_DIR)/include
CXXFLAGS += -g -Wall -Wextra -pthread -std=c++11

GOOGLE_HEADERS = $(GOOGLE_DIR)/include/gtest/*.h \
                $(GOOGLE_DIR)/include/gtest/internal/*.h

GOOGLE_TEST_SRC = $(GOOGLE_DIR)/src/*.cc $(GOOGLE_DIR)/src/*.h $(GOOGLE_HEADERS)

OBJ_GENERATED = yashe_basic_tests.o yashe_vectorial_tests.o yashe_matrix_tests.o yashe_symmetric_matrix_tests.o

all: $(GOOGLE_DIR)/google_test_lib $(OBJ_GENERATED)

# compile the google lib in order to include a single header in our tests' sources codes
$(GOOGLE_DIR)/google_test_lib: $(GOOGLE_TEST_SRC)
	$(CXX) $(CPPFLAGS) -I$(GOOGLE_DIR) $(CXXFLAGS) -c \
            $(GOOGLE_DIR)/src/gtest-all.cc -o $(GOOGLE_DIR)/google_test_lib

yashe_basic_tests.o: Yashe.o $(GOOGLE_DIR)/google_test_lib yashe_basic_tests.cpp $(GOOGLE_HEADERS)
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) yashe_basic_tests.cpp $(GOOGLE_DIR)/google_test_lib $(YASHE_OBJ) $(LDFLAGS) -o yashe_basic_tests.o

yashe_vectorial_tests.o: Yashe.o $(GOOGLE_DIR)/google_test_lib yashe_vectorial_tests.cpp $(GOOGLE_HEADERS)
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) $(SRC_DIR)/vectorutils.h yashe_vectorial_tests.cpp $(GOOGLE_DIR)/google_test_lib $(YASHE_OBJ) $(LDFLAGS) -o yashe_vectorial_tests.o

yashe_matrix_tests.o: Yashe.o $(GOOGLE_DIR)/google_test_lib yashe_matrix_tests.cpp $(GOOGLE_HEADERS)
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) $(SRC_DIR)/vectorutils.h yashe_matrix_tests.cpp $(GOOGLE_DIR)/google_test_lib $(YASHE_OBJ) $(LDFLAGS) -o yashe_matrix_tests.o

yashe_symmetric_matrix_tests.o: Yashe.o $(GOOGLE_DIR)/google_test_lib yashe_symmetric_matrix_tests.cpp $(GOOGLE_HEADERS)
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CXXFLAGS) $(SRC_DIR)/vectorutils.h yashe_symmetric_matrix_tests.cpp $(GOOGLE_DIR)/google_test_lib $(YASHE_OBJ) $(LDFLAGS) -o yashe_symmetric_matrix_tests.o


Yashe.o: ../src/Makefile
	cd ../src/ && make

clean:
	rm -rf *.o $(OBJ_GENERATED)