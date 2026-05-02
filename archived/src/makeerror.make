#
CC=gcc
LOG_LEVEL ?= 2

#COPTS= -O2 -Wall
COPTS = -g -Wall -DLOG_LEVEL=$(LOG_LEVEL)
BUILD_DIR=../build

$(shell mkdir -p $(BUILD_DIR))

adaptation: $(BUILD_DIR)/mesh.o $(BUILD_DIR)/triangulation.o \
               $(BUILD_DIR)/main_error.o \
			   $(BUILD_DIR)/adaptation.o \
               $(BUILD_DIR)/eigen.o \
			   $(BUILD_DIR)/comprImage.o \
               $(BUILD_DIR)/lplib3.o \
               $(BUILD_DIR)/libmesh6.o \
               $(BUILD_DIR)/dynamicArray.o
	$(CC) $(COPTS) -o $(BUILD_DIR)/error \
	    $(BUILD_DIR)/mesh.o \
	    $(BUILD_DIR)/triangulation.o \
		$(BUILD_DIR)/adaptation.o \
	    $(BUILD_DIR)/main_error.o \
	    $(BUILD_DIR)/eigen.o \
		$(BUILD_DIR)/comprImage.o \
	    $(BUILD_DIR)/libmesh6.o \
	    $(BUILD_DIR)/lplib3.o \
	    $(BUILD_DIR)/dynamicArray.o -lpthread -lm

$(BUILD_DIR)/lplib3.o : lplib3.c lplib3.h 
	$(CC) -c $(COPTS) -I. lplib3.c -o $(BUILD_DIR)/lplib3.o

$(BUILD_DIR)/libmesh6.o : libmesh6.c libmesh6.h 
	$(CC) -c $(COPTS) -I. libmesh6.c -o $(BUILD_DIR)/libmesh6.o  

$(BUILD_DIR)/mesh.o : mesh.c mesh.h libmesh6.h 
	$(CC) -c $(COPTS) -I. mesh.c -o $(BUILD_DIR)/mesh.o

$(BUILD_DIR)/dynamicArray.o : dynamicArray.c dynamicArray.h mesh.h 
	$(CC) -c $(COPTS) -I. dynamicArray.c -o $(BUILD_DIR)/dynamicArray.o

$(BUILD_DIR)/triangulation.o : triangulation.c triangulation.h mesh.h
	$(CC) -c $(COPTS) -I. triangulation.c -o $(BUILD_DIR)/triangulation.o

$(BUILD_DIR)/eigen.o : eigen.c eigen.h mesh.h
	$(CC) -c $(COPTS) -I. eigen.c -o $(BUILD_DIR)/eigen.o

$(BUILD_DIR)/adaptation.o : adaptation.c adaptation.h
	$(CC) -c $(COPTS) -I. adaptation.c -o $(BUILD_DIR)/adaptation.o

$(BUILD_DIR)/comprImage.o : comprImage.c comprImage.h triangulation.h mesh.h debug.h
	$(CC) -c $(COPTS) -I. comprImage.c -o $(BUILD_DIR)/comprImage.o

$(BUILD_DIR)/main_error.o : main_error.c adaptation.h
	$(CC) -c $(COPTS) -I. main_error.c -o $(BUILD_DIR)/main_error.o  

clean :
	-rm -rf $(BUILD_DIR)
