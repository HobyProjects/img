#pragma once

#if defined(_WIN32) || defined(_WIN64)
#define IMG_PLATFORM_WINDOWS
#elif defined(__linux__)
#define IMG_PLATFORM_LINUX
#elif defined(__APPLE__) || defined(__MACH__)
#define IMG_PLATFORM_APPLE
#endif

#ifdef IMG_PLATFORM_WINDOWS
#define IMG_API __declspec(dllexport)
#else
#define IMG_API __declspec(dllimport)
#endif

#ifndef IMG_PLATFORM_WINDOWS
#define IMG_API
#endif

#if defined(IMG_ASSERTS_ENABLED) || defined(IMG_LOGGING_ENABLED)
#include <stdio.h>
#endif

#ifdef IMG_ASSERTS_ENABLED
#define IMG_ASSERT(expr, ...) if(!(expr)) { printf(__VA_ARGS__); }
#else
#define IMG_ASSERT(expr, ...)
#endif

#ifdef IMG_LOGGING_ENABLED
#define IMG_LOG(...) printf(__VA_ARGS__);
#else
#define IMG_LOG(...)
#endif

#define IMG_TRUE  1
#define IMG_FALSE 0
#define IMG_NULL  0
#define IMG_SUCCESS  0
#define IMG_FAILURE -1
#define IMG_BIT(n) (1 << n)

typedef char int8;
typedef unsigned char uint8;
typedef short int16;
typedef unsigned short uint16;
typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;

#include <filesystem>
#include <string>
#include <fstream>
#include <sstream>
#include <memory>

namespace img {

    struct IMG_API specification {
        uint8* data;
        uint32 width;
        uint32 height;
        uint32 channels;
        std::string extension;
        std::filesystem::path path;

        specification() = default;
        ~specification() = default;
    };

    IMG_API std::shared_ptr<specification> read(const std::filesystem::path& filePath, bool flip = false);
}



