#OPT     = -g
OPT     = -g -O3 #-mssse3 -mfpmath=sse,387 -ftree-vectorize -march=core2
CFLAGS  = $(OPT) -Wall -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -I.
FFLAGS  = $(OPT)
INCDIRS = -IVERITAS -I/opt/local/include -I${ST_INC_DIR}
LDFLAGS = -L/opt/local/lib -Ltip -LVERITAS -LMinuit -L. -L${ST_LIB_DIR}
LIBS    = -lMyGlastCode -ltip -lVERITAS -lMinuit -lcfitsio -lgfortran \
	  -lxerces-c -lhdf5
LIBSST  = -lirfInterface -lirfLoader -llatResponse -ltestResponse \
          -ldc1aResponse -lirfUtil -lfacilities -lst_facilities -lst_stream \
          -lastro -lCLHEP -lf2c
CC      = gcc-mp-4.8
CXX     = g++-mp-4.8
FC      = gfortran

export CFLAGS FFLAGS CC CXX FC

SOURCES = Catalog.cpp FT1.cpp FT2.cpp MyMinuit.cpp PSF.cpp Analysis.cpp \
	Magic7Dispatcher.cpp SourceLibrary.cpp RegularSpline.cpp IRF.cpp \
	GTI.cpp FITS.cpp PLIntegratedEA.cpp FT2ROI.cpp FT2Exp.cpp \
	FT1ROIFilter.cpp BlockPartition.cpp FunctionXY.cpp
OBJECTS = $(SOURCES:.cpp=.o)
LIBDEPS = VERITAS/libVeritas.a tip/libtip.a Minuit/libMinuit.a libMyGlastCode.a
PROGS = calc_psf calc_event_orient iminuit orbit m7sort plot_spline_model \
	lomblike dump_psf avg_psf flare_probability bayes_block

all: $(PROGS)

libMyGlastCode.a: $(OBJECTS)
	ar r $@ $^

test_corr: test_corr.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) 

bayes_block: bayes_block.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) -lgsl -lgslcblas

flare_probability: flare_probability.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) $(LIBSST) -lgsl -lgslcblas

calc_psf: calc_psf.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

dump_psf: dump_psf.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

avg_psf: avg_psf.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

calc_event_orient: calc_event_orient.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

iminuit: iminuit.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

orbit: orbit.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

m7sort: m7sort.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

plot_spline_model: plot_spline_model.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

lomblike: lomblike.o | $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS) $(LIBSST)

VERITAS/libVeritas.a: VERITAS
tip/libtip.a: tip
Minuit/libMinuit.a: Minuit

SUBDIRS = VERITAS tip Minuit

.PHONY: distclean clean $(SUBDIRS) $(addsuffix -clean,$(SUBDIRS))

$(SUBDIRS):
	$(MAKE) -C $@

$(addsuffix -clean,$(SUBDIRS)):
	$(MAKE) -C $(@:-clean=) clean

clean: 
	$(RM) *.o $(PROGS) libMyGlastCode.a *~

distclean: clean $(addsuffix -clean,$(SUBDIRS))

%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCDIRS) -c $<
