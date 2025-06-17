# Variables
CXX = g++
# Optimization flags
CFLAGS = -O3 -march=native -mtune=native
# Debugging flags
DBGFLAGS = -ggdb3 -m64 -fPIC -fexceptions -DNDEBUG -DIL_STD -ldl -fsanitize=undefined -fsanitize=address
TARGET = gen
SRCS = algo.cpp load.cpp RandomNumbers.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = load.h RandomNumbers.h json.hpp

# Default target
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS)

%.o: %.cpp $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@

# Debug target
debug: CFLAGS += $(DBGFLAGS)  # Append flags for debugging
debug: $(TARGET)             # Use the same target for debugging

clean:
	rm -f $(OBJS) $(TARGET)

usage:
	@echo "Usage of the executable:"
	@echo "./gen instance_type=(char*)cityName seed=(int)seed"
	@echo "Example:"
	@echo "./gen instance_type=newyork seed=1"

# Phony targets
.PHONY: all clean debug usage
