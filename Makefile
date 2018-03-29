CXX       = g++

DEFINES   = -DDISABLE_RANDOM_SEED

# флаги компиляции
CXXFLAGS  =
CXXFLAGS += -std=c++11 -O2
CXXFLAGS += -g
CXXFLAGS += -Wall -Wextra
CXXFLAGS += $(DEFINES)
CXXFLAGS += $(shell pkg-config --cflags blitz)

# флаги сборки (библиотеки)
LDFLAGS   = 
LDFLAGS  += -L./ -llapack -lopenblasp -lgfortran -lpthread -ldcmt
LDFLAGS  += $(shell pkg-config --libs blitz)

SOURCES   = main.cc
BINARY    = autoreg

VISUAL    = visual
VISUAL_SOURCES = visual.cc
VISUAL_LDFLAGS = $(shell pkg-config --libs freeglut) -lGL

INITIALIZER = initializer_mt
INITIALIZER_SOURCES = initializer_mt.cc

$(BINARY): $(SOURCES) *.hh Makefile
	$(CXX) $(CXXFLAGS) $(SOURCES) $(LDFLAGS) -o $(BINARY)

$(VISUAL): Makefile *.hh $(VISUAL_SOURCES)
	$(CXX) $(CXXFLAGS) $(VISUAL_SOURCES) $(VISUAL_LDFLAGS) -o $(VISUAL)

$(INITIALIZER): Makefile  parallel_mt.hh dc.h $(INITIALIZER_SOURCES)
	$(CXX) $(CXXFLAGS) $(INITIALIZER_SOURCES) -L./ -ldcmt -o $(INITIALIZER)

run: ../tests autoreg.model
run: $(BINARY)
	(cd ../tests; $(PWD)/$(BINARY))

debug: ../tests autoreg.model
debug: $(BINARY)
	(cd ../tests; gdb $(PWD)/$(BINARY))

../tests:
	mkdir -p ../tests

../tests/autoreg.model: autoreg.model
	cp ../input/autoreg.model ../tests

clean:
	rm -f $(BINARY) $(VISUAL)
