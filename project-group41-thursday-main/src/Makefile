# this makefile is intended for g++ on Linux

CC = g++
# CFLAGS = -c -Wall -Wpedantic
CFLAGS = -c -O3
LDFLAGS = -pthread
SOURCES = experiment.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp  mono.cpp  stereo.cpp  rf_frontend.cpp  rds.cpp   
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = experiment

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm $(OBJECTS)
