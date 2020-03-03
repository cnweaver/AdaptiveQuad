all : quadTest

GSL_FLAGS=$(shell if which pkg-config 2>&1 > /dev/null && pkg-config gsl; then pkg-config gsl --cflags --libs; echo -DHAVE_GSL=1 ; fi)

quadTest : test/quadTest.cpp AdaptiveQuad.h test/gg.h test/gsl_interface.h test/romberg.h Makefile
	$(CXX) -std=c++11 -O3 $(GSL_FLAGS) $(CXXFLAGS) test/quadTest.cpp -lgsl -o quadTest

clean : 
	rm -f quadTest