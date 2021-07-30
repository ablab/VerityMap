#pragma once

#include <string>
#include <sys/stat.h>
#include <filesystem>

//TODO: throw exception if this is file
inline void ensure_dir_existance(const std::filesystem::path & path) {
    struct stat statbuf{};
    if (not std::filesystem::is_directory(path)) {
        std::filesystem::create_directories(path);
    }
}

//TODO: throw exception if this is file
inline void recreate_dir(const std::filesystem::path & path) {
    struct stat statbuf{};
    if (std::filesystem::is_directory(path)) {
        std::filesystem::remove_all(path);
    }
    std::filesystem::create_directories(path);
}