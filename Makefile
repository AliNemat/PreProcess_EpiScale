CC=g++ -std=c++11

CFLAGS=-c -Wall

all: program

program: main.o coord.o algorithm.o
		$(CC) main.o coord.o algorithm.o

main.o: main.cpp
		$(CC) $(CFLAGS) main.cpp

coord.o: coord.cpp
		$(CC) $(CFLAGS) coord.cpp

algorithm.o: algorithm.cpp
		$(CC) $(CFLAGS) algorithm.cpp

clean:
		rm -rf  *.o ./a.out
