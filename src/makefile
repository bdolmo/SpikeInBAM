CXX = g++
CC = gcc
CXXFLAGS = -Wall -g
CFLAGS = -Wall -g
HTSLIB_DIR = htslib

# INFO for STATIC LINKING:
# I had lots of issues when trying to link curl library (libcurl) from htslib. 
# I ended up building htslib without curl support: ./configure --disable-libcurl

INCLUDES = -I$(HTSLIB_DIR)

LIBS = -L$(HTSLIB_DIR) -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu -lhts -lz -lbz2 -llzma -lpthread

ifeq ($(STATIC),1)
    LIBS += -static -ldeflate
endif

CPP_SRCS = SpikeInBAM.cpp BamReader.cpp BamRecord.cpp BamWriter.cpp RefFasta.cpp
C_SRCS = 
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
C_OBJS = $(C_SRCS:.c=.o)
OBJS = $(CPP_OBJS) $(C_OBJS)
DEPS = $(OBJS:.o=.d)
TARGET = spikeinbam

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

-include $(DEPS)

.PHONY: clean run

clean:
	rm -f $(OBJS) $(TARGET) $(DEPS)

run: $(TARGET)
	./$(TARGET)
