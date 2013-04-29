CC = c++
CFLAGS = -c -I/opt/local/include/openbabel-2.0/
LDFLAGS = -L/opt/local/lib/ -lopenbabel -lcblas

all: example

example: example.o
	$(CC) $(LDFLAGS) example.o -o example

example.o: example.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) example.cpp

clean:
	rm -rf example.o example