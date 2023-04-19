all: main_test

ifeq ($(benchmark), 1)
    objects = main_p1.cpp
	output = main_test_p1
else ifeq ($(benchmark), 2)
    objects = main_p2.cpp
	output = main_test_p2
endif

selected_config = ches_config_files/config_file_n_exp_$(config).h

set_config:
	cp $(selected_config) ches_config_files/config_file.h

libblst.a:
	./build.sh

main_test: $(objects) libblst.a set_config
	g++ -std=c++17 -o $(output) -g -O2 $(objects) libblst.a



clean:
	rm -f libblst.a
	rm -f main_test_p1
	rm -f main_test_p2