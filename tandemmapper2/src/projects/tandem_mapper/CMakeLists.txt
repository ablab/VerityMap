project(tandem_mapper)

add_library(tandem_mapper_library STATIC tandem_mapper.cpp cigar.cpp config/config.cpp)

target_link_libraries(tandem_mapper_library common sequence m ksw2 sketch bloomfilter)

add_executable(tandem_mapper tandem_mapper.cpp ${VERSION} ${CONFIG_DIR_DEF_CPP})

target_link_libraries(tandem_mapper tandem_mapper_library)

add_custom_command(TARGET tandem_mapper POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/config/ $<TARGET_FILE_DIR:tandem_mapper>/config)
