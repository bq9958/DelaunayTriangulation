#
CC=gcc
COPTS= -O2 -Wall
OUTPUT = mesh

ifeq ($(MODE),DEBUG)
  COPTS = -g -O0 -Wall
  OUTPUT = mesh_debug
endif


mesh:	mesh.o main_mesh.o eigen.o lplib3.o libmesh6.o  
	$(CC) $(COPTS) -o $(OUTPUT) mesh.o  main_mesh.o eigen.o libmesh6.o lplib3.o -lpthread -lm

lplib3.o :	lplib3.c lplib3.h 
	$(CC) -c $(COPTS)  -I. lplib3.c

libmesh6.o :	libmesh6.c libmesh6.h 
	$(CC) -c $(COPTS)  -I. libmesh6.c  

mesh.o :	mesh.c mesh.h libmesh6.h 
	$(CC) -c $(COPTS)  -I. mesh.c

eigen.o :	eigen.c mesh.h 
	$(CC) -c $(COPTS)  -I. eigen.c

main_mesh.o :	main_mesh.c mesh.h
	$(CC) -c $(COPTS)  -I. main_mesh.c  

clean :
	-rm $(OUTPUT) mesh.o main_mesh.o eigen.o lplib3.o libmesh6.o 
