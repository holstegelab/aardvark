CFLAGS=-lm -lz
LDFLAGS=$(CFLAGS)
CXX=g++
CXX_FLAGS=-std=c++17 -Wall -O2 -g

ODIR = build
SDIR = src


OUT = aardvark

SRCFILES := $(wildcard src/*.cpp)

OBJS = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRCFILES))

DEPS = $(OBJS:%.o=%.d)

.PHONY: all clean

all: $(ODIR)/$(OUT)

clean:
	$(RM) $(ODIR)/$(OUT) $(OBJS) $(DEPS)

$(ODIR)/$(OUT): $(OBJS)
	$(CXX) $(CXX_FLAGS) -o $@ $^ $(CFLAGS)

-include $(DEPS)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) $(CXX_FLAGS) -MMD -MP -c -o $@ $< $(LDFLAGS)