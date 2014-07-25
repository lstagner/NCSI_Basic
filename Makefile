MODELS = guiding_center.o 
INTEGRATORS = runge-kutta4.o noncanonical_symplectic.o 
FIELDS = axisymmetric_tokamak.o

DRIVER_DEPENDS = input_parser.o $(INTEGRATORS) $(MODELS) $(FIELDS)

# CUDA_DRIVER: Compile objects using NVCC instead g++
CUDA_INTEGRATORS = $(INTEGRATORS:.o=_cu.o) 
CUDA_MODELS = $(MODELS:.o=_cu.o)
CUDA_FIELDS = $(FIELDS:.o=_cu.o)
CUDA_DRIVER_DEPENDS = input_parser.o $(CUDA_INTEGRATORS) $(CUDA_MODELS) $(CUDA_FIELDS)

#CXXFLAGS = -g -Wall -Wextra -std=c++0x 
CXXFLAGS = -O3
# Check whether we're on the K20 node which uses sm_35
ifeq ($(HOSTNAME),gpusrv02.pppl.gov)
NVCC_FLAGS=-rdc=true -arch=sm_35 -g -G #-O3 
else 
NVCC_FLAGS=-rdc=true -arch=sm_20 -g -G #-O3 
endif
BOOST_FLAGS=-L$(BOOST_LIBRARY_DIR) -lboost_program_options

all: driver cuda_driver

driver : driver.o $(DRIVER_DEPENDS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOST_FLAGS)

cuda_driver : cuda_driver.o $(CUDA_DRIVER_DEPENDS)
	nvcc $(NVCC_FLAGS) -o $@ $^ $(BOOST_FLAGS)

cuda_driver.o : cuda_driver.cu
	nvcc $(NVCC_FLAGS) -c -o $@ $<

# DONT use -O3 here. Breaks the linking for some reason.
input_parser.o : input_parser.cc input_parser.h
	$(CXX) -c -o $@ $<

# Static pattern rules
$(INTEGRATORS) : %.o: %.cc %.h integrator.h guiding_center.h em_fields.h
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

$(MODELS) : %.o: %.cc %.h em_fields.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(FIELDS) : %.o: %.cc %.h em_fields.h
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

# Rule for making _cu.o objects
%_cu.o : %.cu
	nvcc $(NVCC_FLAGS) -c -o $@ $^


clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
