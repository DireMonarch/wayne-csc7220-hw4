CC      =       /opt/ohpc/pub/mpi/openmpi3-gnu7/3.1.0/bin/mpicc
CCLINK  =       /opt/ohpc/pub/mpi/openmpi3-gnu7/3.1.0/bin/mpicc
SHELL   =       /bin/sh

all: main

main: jim-haslett-csc7220-hw4.c
	$(CC) -lm -o jim-haslett-csc7220-hw4 jim-haslett-csc7220-hw4.c 

clean:
	/bin/rm -f jim-haslett-csc7220-hw4