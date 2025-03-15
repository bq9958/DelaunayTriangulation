#
CC=gcc
COPTS= -O2 -Wall
BUILD_DIR=../build


$(shell mkdir -p $(BUILD_DIR))

mesh: $(BUILD_DIR)/mesh.o $(BUILD_DIR)/main_mesh.o $(BUILD_DIR)/eigen.o $(BUILD_DIR)/lplib3.o $(BUILD_DIR)/libmesh6.o  
	$(CC) $(COPTS) -o $(BUILD_DIR)/mesh $(BUILD_DIR)/mesh.o $(BUILD_DIR)/main_mesh.o $(BUILD_DIR)/eigen.o $(BUILD_DIR)/libmesh6.o $(BUILD_DIR)/lplib3.o -lpthread -lm

$(BUILD_DIR)/lplib3.o : lplib3.c lplib3.h 
	$(CC) -c $(COPTS) -I. lplib3.c -o $(BUILD_DIR)/lplib3.o

$(BUILD_DIR)/libmesh6.o : libmesh6.c libmesh6.h 
	$(CC) -c $(COPTS) -I. libmesh6.c -o $(BUILD_DIR)/libmesh6.o  

$(BUILD_DIR)/mesh.o : mesh.c mesh.h libmesh6.h 
	$(CC) -c $(COPTS) -I. mesh.c -o $(BUILD_DIR)/mesh.o

$(BUILD_DIR)/eigen.o : eigen.c mesh.h 
	$(CC) -c $(COPTS) -I. eigen.c -o $(BUILD_DIR)/eigen.o

$(BUILD_DIR)/main_mesh.o : main_mesh.c mesh.h
	$(CC) -c $(COPTS) -I. main_mesh.c -o $(BUILD_DIR)/main_mesh.o  

clean :
	-rm -rf $(BUILD_DIR)