.PHONY: all clean

CXX = g++

CXXFLAGS = -std=c++17 -Wall -O2 -fPIC

PROJECT_DIR = .
ROOT_LIB_DIR = /opt/root/root-install/lib
ROOT_INCLUDE_DIR = /opt/root/root-install/include

ROOT_LIBS = -L $(ROOT_LIB_DIR) -lGeom -lRGL -lGed -lTreePlayer -lCore -lHist -lGraf -lGraf3d -lMathCore -lGpad -lTree -lRint -lRIO -lPostscript -lMatrix -lPhysics -lMinuit -lMinuit2 -lGui -lASImage -lASImageGui -pthread -lm -ldl -rdynamic -lstdc++

LIBS = -fopenmp

SRCS = $(shell find $(PROJECT_DIR)/src/ -name '*.cc') $(PROJECT_DIR)/Mainfile.cc 
OBJS = $(patsubst %.cc,%.o,$(SRCS))

GREEN := $(shell tput -Txterm setaf 2)
BLUE  := $(shell tput -Txterm setaf 6)
RESET := $(shell tput -Txterm sgr0)

all: LEPAM
	
clean: 
	@rm -rf *.a LEPAM $(OBJS)
		
LEPAM: Mainfile.o libLEPAM.a
	$(info $(GREEN)Linking Mainfile.o with libLEPAM.a and ROOT libs$(RESET))
	@${CXX} ${CXXFLAGS} Mainfile.o -o LEPAM -fopenmp -L. -lLEPAM $(ROOT_LIBS)
	@rm -rf $(OBJS)
	$(info $(BLUE)Built the target LEPAM$(RESET))

%.o: %.cc
	$(info $(GREEN)Generating the object $@ $(RESET))
	@$(CXX) $(CXXFLAGS) -c $< -o $@ -I $(PROJECT_DIR)/include -I $(ROOT_INCLUDE_DIR) $(LIBS)

libLEPAM.a: $(OBJS)
	$(info $(GREEN)Generating the library libLEPAM.a$(RESET))
	@ar rcs libLEPAM.a $^