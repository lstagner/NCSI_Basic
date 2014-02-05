MODELS = guiding_center.o axisymmetric_tokamak.o
INTEGRATORS = runge-kutta.o noncanonical_symplectic.o 

DRIVER_DEPENDS = input_parser.o $(INTEGRATORS) $(MODELS) 

CXXFLAGS = -g -Wall -Wextra -std=c++0x
BOOSTFLAGS = -L /usr/local/lib -lboost_program_options

all: driver

driver : driver.o $(DRIVER_DEPENDS)
	$(CXX) $(CXXFLAGS) -o $@  $< $(DRIVER_DEPENDS) $(BOOSTFLAGS)

driver.o : driver.cc input_parser.h guiding_center.h
	$(CXX) $(CXXFLAGS) -I $(MODELDIR) -I $(INTEGRATORDIR) -c -o $@ $<

input_parser.o : input_parser.cc input_parser.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Static pattern rule for integrators
$(INTEGRATORS) : %.o: %.cc %.h integrator.h guiding_center.h 
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
