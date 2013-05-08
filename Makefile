CXX = c++ -g -w
LDFLAGS = -lopenbabel -lcblas
CFLAGS = -c
#-framework Accelerate

SONAME=-soname

ifeq ($(shell uname -s), Darwin)
    SONAME=-install_name
    LDFLAGS += -L/opt/local/lib/ -lclapack
    CFLAGS += -I/opt/local/include/openbabel-2.0/
else
	LDFLAGS += -llapack
endif


all: example

example: example.o
	$(CXX) $(LDFLAGS) example.o -o example

example.o: example.cpp
	$(CXX) $(CFLAGS) example.cpp

clean:
	rm -rf example.o example