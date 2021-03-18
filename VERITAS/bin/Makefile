SOURCE = $(wildcard *.cpp)
OBJECTS = $(SOURCE:.cpp=.o)

all: octaveio

octaveio: octaveio.o
	$(CXX) -L/opt/local/lib -L.. -o $@ $^ -lVERITAS -lhdf5

clean:
	$(RM) $(OBJECTS) *.a *~

%.o: %.cpp
	$(CXX) -I.. $(CFLAGS) -c $<
