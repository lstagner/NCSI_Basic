MODELS = guiding_center.o 
INTEGRATORS = runge-kutta4.o noncanonical_symplectic.o 
FIELDS = axisymmetric_tokamak.o

DRIVER_DEPENDS = input_parser.o $(INTEGRATORS) $(MODELS) $(FIELDS)

#CXXFLAGS = -g -Wall -Wextra -std=c++0x 
CXXFLAGS = -O3
BOOST_FLAGS=-L$(BOOST_LIBRARY_DIR) -lboost_program_options

all: driver

driver : driver.o $(DRIVER_DEPENDS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOST_FLAGS)

input_parser.o : input_parser.cc input_parser.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Static pattern rules
$(INTEGRATORS) : %.o: %.cc %.h integrator.h guiding_center.h em_fields.h
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

$(MODELS) : %.o: %.cc %.h em_fields.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(FIELDS) : %.o: %.cc %.h em_fields.h
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
