CXX = c++ -g -w
LDFLAGS = -lopenbabel -lcblas
CFLAGS = -c
#-framework Accelerate

SONAME=-soname

ifeq ($(shell uname -s), Darwin)
    SONAME=-install_name
    LDFLAGS += -L/opt/local/lib/ -lclapack
    CFLAGS += -I/opt/local/include/openbabel-2.0/ -ffast-math -lm
else
	LDFLAGS += -llapack
endif


all: example

example: example.o
	$(CXX) -o example *.o $(LDFLAGS) 

example.o: example.cpp
	$(CXX) $(CFLAGS) example.cpp

clean:
	rm -rf example.o example