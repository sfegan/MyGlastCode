SOURCE = $(wildcard *.cpp)
OBJECTS = $(SOURCE:.cpp=.o)

SWIG      = swig -cpperraswarn -python -c++ $(DEFINES) 
PYCONFIG  = python2.7-config

CFLAGS += -I.

all: libVERITAS.a bin _RandomNumbers.so

libVERITAS.a: $(OBJECTS)
	$(AR) r $@ $^

RANDOMNUMBERSHPP=SimpleRNG.hpp RandomNumbers_TNG.hpp

_RandomNumbers.so: RandomNumbers_wrap.o $(RANDOMNUMBERSHPP:.hpp=.o)
	$(CXX)  `$(PYCONFIG) --ldflags` $(LDFLAGS) -shared -o $@ $^ `$(PYCONFIG) --libs`

RandomNumbers_wrap.o: RandomNumbers_wrap.cxx
	$(CXX) -I. `$(PYCONFIG) --includes` -c $<

RandomNumbers_wrap.cxx: RandomNumbers.i RandomNumbers.hpp $(RANDOMNUMBERSHPP)
	$(SWIG) $<

.PHONY: clean bin

bin: 
	$(MAKE) -C bin

clean:
	$(RM) $(OBJECTS) *.a \
		_RandomNumbers.so RandomNumbers.py* \
		RandomNumbers_wrap.cxx RandomNumbers_wrap.o \
		*~

%.o: %.cxx
	$(CXX) $(CFLAGS) -D -c $<

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<
