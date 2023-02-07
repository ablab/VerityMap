build_type?="Release"

.PHONY: clean cmake all

all: veritymap

cmake:
	mkdir -p build
	cd build && cmake .. -DCMAKE_BUILD_TYPE="${build_type}"

veritymap: cmake
	$(MAKE) -C build all
	mkdir -p build/bin
	mv $(abspath build/src/projects/veritymap/veritymap) build/bin/veritymap
	rm -rf build/bin/config
	mv $(abspath build/src/projects/veritymap/config) build/bin

test_launch: veritymap
	build/bin/veritymap \
		--target test_dataset/test_target.fasta \
		--queries test_dataset/test_query.fasta -o test_dataset/test_launch
	grep -q "Thank you for using VerityMap!" test_dataset/test_launch/veritymap.log

clean:
	-rm -r build
	-rm -r test_dataset/test_launch

