SOURCEF = $(wildcard *.f)
SOURCEC = $(wildcard *.c)
OBJECTS = $(SOURCEF:.f=.o) $(SOURCEC:.c=.o)

all: libMinuit.a # libDynMinuit.so

libMinuit.a: $(OBJECTS)
	$(AR) r $@ $^

libDynMinuit.so: $(OBJECTS:.o=.shared_o)
	$(LD) -shared -o $@ $^

.PHONY: clean

clean:
	$(RM) $(OBJECTS) $(OBJECTS:.o=.shared_o) *.a *.so *~

%.o: %.cxx
	$(CXX) $(CFLAGS) -I.. --include cstring -DBUILD_WITHOUT_ROOT -c $<

%.shared_o: %.cxx
	$(CXX) $(CFLAGS) -fPIC -o $@ -I.. -c $<

%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.shared_o: %.f
	$(FC) $(FFLAGS) -fPIC -o $@ -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<

%.shared_o: %.c
	$(CC) $(CFLAGS) -fPIC -o $@ -c $<
