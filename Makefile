all: main_test

# One need to modify the following two directories to the correct paths.
openssl_include_directory = /opt/homebrew/include/
openssl_library_directory = /opt/homebrew/Cellar/openssl@3/3.6.3/lib

ifeq ($(group), 1)
    objects = msm/bench_g1.cpp
	output = main_test_p1
else ifeq ($(group), 2)
    objects = msm/bench_g2.cpp
	output = main_test_p2
endif

selected_config = config/n_exp_$(config).h

bucket_index ?= direct
bench_test_num ?= 5
bench_loop_num ?= 0
bench_deterministic ?= 0
benchmark_flags = -DMSM_BENCH_TEST_NUM=$(bench_test_num) -DMSM_BENCH_LOOP_NUM=$(bench_loop_num) -DMSM_BENCH_DETERMINISTIC=$(bench_deterministic)
ifeq ($(bucket_index), rank)
    bucket_index_flags = -DMSM_BUCKET_INDEX_RANK=1
else ifneq ($(bucket_index), direct)
    $(error bucket_index must be direct or rank)
endif

set_config:
	cp $(selected_config) config/config_file.h

SDKROOT = $(shell xcrun --show-sdk-path 2>/dev/null)

libblst.a:
	SDKROOT=$(SDKROOT) CFLAGS="-O -fno-builtin -fPIC -Wall -Wextra -Werror -isysroot $(SDKROOT)" ./build.sh

main_test: $(objects) libblst.a set_config
	g++ -std=c++17 -o $(output) -g -O2 $(bucket_index_flags) $(benchmark_flags) -I. $(objects) libblst.a -I$(openssl_include_directory) -L$(openssl_library_directory) -lcrypto

clean:
	rm -f libblst.a
	rm -f main_test_p1
	rm -f main_test_p2
