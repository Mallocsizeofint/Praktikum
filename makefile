STANDARDFLAG = -std=c++11
CXX = clang++
CXXFLAGS = -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-float-equal $(STANDARDFLAG)
LD  = clang++
DEBUG = -DENABLE_RANGE_CHECK -O3
RELEASE = -DNDEBUG -O3
TARGET = blatt1.exe
SOURCES := $(wildcard *.cc)        
MYOBJS := $(patsubst %.cc, %.o, $(SOURCES))
KUSKOBJS = cputime.o createVTKFile.o
OBJS = $(MYOBJS) $(KUSKOBJS)

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
	$(CXX) $(CXXFLAGS) $(DEBUG) -c $< -o $@
 
%: %.cc
	$(CXX) $(CXXFLAGS) -o $@ $<
 
clean:
	rm -f .depend *.o
