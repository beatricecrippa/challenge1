CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++20
SRC_DIR = src
INCLUDE_DIR = include
TARGET = main
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,%.o,$(SRCS))
INCLUDES = -I$(INCLUDE_DIR)
%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
clean:
	$(RM) *.o

distclean: clean
	$(RM) *~