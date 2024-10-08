# GNU C++ Compiler
CPP         = g++      # C++ Compuler

CPPFLAGS    = -O3 -c    

# The directory where the include files needed to create the library objects are
INC_DIR = ../include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../lib

# # armadillo library
# ARMA_INC	= $(CONDA_PREFIX)/include
# ARMA_LIB	= $(CONDA_PREFIX)/lib -larmadillo

# linkable object files with class info
GO_OBJS = go_class.o

# Lennard Jones library containing Lennard Jones Class
GO_LIB = libgoclass.a

go_class.o:  $(INC_DIR)/gaus_overlap.h
	$(CPP) $(CPPFLAGS) go_class.cpp \
	-I$(INC_DIR)

# The following target creates a static library 
libgoclass.a: $(GO_OBJS)
	ar ru $(GO_LIB) $(GO_OBJS) 
	mv $(GO_LIB) $(LIB_DIR)
	rm $(GO_OBJS)

all: $(GO_LIB)

cleanall:
	cd $(LIB_DIR); rm $(GO_LIB)

