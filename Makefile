build_type?="Release"

.PHONY: clean cmake all

all: tandem_mapper

cmake:
	mkdir -p tandemmapper2/build
	cd tandemmapper2/build && cmake .. -DCMAKE_BUILD_TYPE="${build_type}"

tandem_mapper: cmake
	$(MAKE) -C tandemmapper2/build all
	mkdir -p tandemmapper2/build/bin && ln -s -f $(abspath tandemmapper2/build/src/projects/tandem_mapper/tandem_mapper) tandemmapper2/build/bin/tandem_mapper

test_launch: tandem_mapper
	tandemmapper2/build/bin/tandem_mapper \
		--target tandemmapper2/test_dataset/test_target.fasta \
		--queries tandemmapper2/test_dataset/test_query.fasta -o tandemmapper2/test_dataset/test_launch
	grep -q "Thank you for using TandemMapper2!" tandemmapper2/test_dataset/test_launch/tandem_mapper.log

clean:
	-rm -r tandemmapper2/build
	-rm -r tandemmapper2/test_dataset/test_launch
