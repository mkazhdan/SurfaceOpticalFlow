SURFACE_OPTICAL_FLOW_TARGET=SurfaceOpticalFlow

SURFACE_OPTICAL_FLOW_SOURCE=SurfaceOpticalFlow/SurfaceOpticalFlow.cpp


CFLAGS += -fpermissive -fopenmp -Wno-deprecated -Wno-unused-result -Wno-format -msse2 -std=c++11
LFLAGS += -lgomp -lz -lpng -ltiff -ljpeg -lpthread -lGL -lGLU -lfftw3 -lfftw3f -lglut

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math
LFLAGS_RELEASE = -O3 

SRC = ./
BIN = Bin/
BIN_O = ./
INCLUDE = /usr/include/ -I.

CC=gcc
CXX=g++
MD=mkdir

SURFACE_OPTICAL_FLOW_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(SURFACE_OPTICAL_FLOW_SOURCE))))

all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: $(BIN)
all: $(BIN)$(SURFACE_OPTICAL_FLOW_TARGET)

debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: $(BIN)
debug: $(BIN)$(SURFACE_OPTICAL_FLOW_TARGET)

clean:
	rm -f $(BIN)$(SURFACE_OPTICAL_FLOW_TARGET)
	rm -f $(SURFACE_OPTICAL_FLOW_OBJECTS)

$(BIN):
	$(MD) -p $(BIN)

$(BIN)$(SURFACE_OPTICAL_FLOW_TARGET): $(SURFACE_OPTICAL_FLOW_OBJECTS)
	$(CXX) -o $@ $(SURFACE_OPTICAL_FLOW_OBJECTS) $(LFLAGS)

$(BIN_O)%.o: $(SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<


