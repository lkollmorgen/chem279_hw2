# This creates executables for tests of the gaussian orbital integration class 

# GNU C++ Compiler
CPP         = g++      

# compiler flags -O3 optimized version -c create a library object
CPPFLAGS    = -O0 -o

# The directory where our gaussian orbital main library include files are stored.
INC_DIR = ../include

# The directory where lib object files are stored
LIB_DIR     = ../lib

# Executable programs
BIN_DIR     = ../bin

# The next line contains a list of executable object files 
# created by this Makefile.
GO_EXECUTABLES = pt1 pt2

pt1: $(INC_DIR)/gaus_overlap.h
	$(CPP) $(CPPFLAGS) pt1 pt1.cpp \
	-I$(INC_DIR) -L$(LIB_DIR) -lgoclass
	mv pt1 $(BIN_DIR)

run_pt1:
	cd $(BIN_DIR); ./pt1

pt2: $(INC_DIR)/gaus_overlap.h
	$(CPP) $(CPPFLAGS) pt2 pt2.cpp \
	-I$(INC_DIR) -L$(LIB_DIR) -lgoclass
	mv pt2 $(BIN_DIR)

run_pt2:
	cd $(BIN_DIR); ./pt2

all:
	cd ../src; make all
	make pt1
	make run_pt1
	make pt2
	make run_pt2
	
cleanall:
	cd $(BIN_DIR); rm $(GO_EXECUTABLES) 

