build_type?="Release"

.PHONY: clean cmake all

all: tandem_mapper

cmake:
	mkdir -p tandemmapper/build
	cd tandemmapper/build && cmake .. -DCMAKE_BUILD_TYPE="${build_type}"

tandem_mapper: cmake
	$(MAKE) -C tandemmapper/build all
	mkdir -p tandemmapper/build/bin && ln -s -f $(abspath tandemmapper/build/src/projects/tandem_mapper/tandem_mapper) tandemmapper/build/bin/tandem_mapper

test_launch: tandem_mapper
	tandemmapper/build/bin/tandem_mapper \
		--target tandemmapper/test_dataset/test_target.fasta \
		--queries tandemmapper/test_dataset/test_query.fasta -o tandemmapper/test_dataset/test_launch
	grep -q "Thank you for using TandemMapper2!" tandemmapper/test_dataset/test_launch/tandem_mapper.log

clean:
	-rm -r tandemmapper/build
	-rm -r tandemmapper/test_dataset/test_launch
