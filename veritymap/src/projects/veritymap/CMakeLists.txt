project(veritymap)

add_library(veritymap_library STATIC cigar.cpp config/config.cpp)

target_link_libraries(veritymap_library common sequence m ksw2 sketch bloomfilter)

add_executable(veritymap veritymap.cpp ${VERSION} ${CONFIG_DIR_DEF_CPP})

target_link_libraries(veritymap veritymap_library)

add_custom_command(TARGET veritymap POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/config/ $<TARGET_FILE_DIR:veritymap>/config)


# add_executable(test_rh test_rh.cpp ${VERSION} ${CONFIG_DIR_DEF_CPP})
# target_link_libraries(test_rh common sequence m veritymap_library)
