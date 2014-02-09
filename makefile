STANDARDFLAG = -std=c++11
CXX = clang++
CXXFLAGS = -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-float-equal $(STANDARDFLAG)
LD  = clang++
DEBUG = -DENABLE_RANGE_CHECK -O3 -march=native
RELEASE = -DNDEBUG -O3 -march=native
TARGET = blatt1.exe
SOURCES := $(wildcard *.cc)        
MYOBJS := $(patsubst %.cc, %.o, $(SOURCES))
KUSKOBJS = cputime.o createVTKFile.o
OBJS = $(MYOBJS) $(KUSKOBJS)

RELEASE_DEBUG_FLAGS = $(RELEASE)

all: $(TARGET)

$(TARGET): .depend $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@
	
depend: .depend

.depend: cmd = $(CXX) -MM -MF depend $(var); cat depend >> .depend;
.depend: $(SOURCES)
	@echo "Generating dependencies..."
	@$(foreach var, $(SOURCES), $(cmd))
	@rm -f depend

-include .depend
    
%.o: %.cc
	$(CXX) $(CXXFLAGS) $(RELEASE_DEBUG_FLAGS) -c $< -o $@
 
%: %.cc
	$(CXX) $(CXXFLAGS) -o $@ $<
 
clean:
	rm -f .depend $(MYOBJS)
