SOURCE = $(wildcard *.cxx)
OBJECTS = $(SOURCE:.cxx=.o)

all: libtip.a

libtip.a: $(OBJECTS)
	$(AR) r $@ $^

.PHONY: clean

clean:
	$(RM) $(OBJECTS) *.a *~

%.o: %.cxx
	$(CXX) $(CFLAGS) -I/opt/local/include -I.. --include cstring -DBUILD_WITHOUT_ROOT -c $<
