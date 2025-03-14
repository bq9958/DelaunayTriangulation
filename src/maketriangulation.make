#
CC=gcc
#COPTS= -O2 -Wall
COPTS = -g -Wall

triangulation:	mesh.o triangulation.o main_triangulation.o eigen.o lplib3.o libmesh6.o  
	$(CC) $(COPTS) -o triangulation mesh.o triangulation.o main_triangulation.o eigen.o libmesh6.o lplib3.o -lpthread -lm

lplib3.o :	lplib3.c lplib3.h 
	$(CC) -c $(COPTS)  -I. lplib3.c

libmesh6.o :	libmesh6.c libmesh6.h 
	$(CC) -c $(COPTS)  -I. libmesh6.c  

mesh.o :	mesh.c mesh.h libmesh6.h 
	$(CC) -c $(COPTS)  -I. mesh.c

triangulation.o :	triangulation.c triangulation.h mesh.h
	$(CC) -c $(COPTS)  -I. triangulation.c

eigen.o :	eigen.c mesh.h 
	$(CC) -c $(COPTS)  -I. eigen.c

main_triangulation.o :	main_triangulation.c triangulation.h
	$(CC) -c $(COPTS)  -I. main_triangulation.c  

clean :
	-rm triangulation mesh.o main_triangulation.o eigen.o lplib3.o libmesh6.o triangulation.o