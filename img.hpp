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

    /**
     * @struct specification
     * @brief Holds metadata for an image.
     *
     * This struct contains information about an image, including its data,
     * dimensions, number of channels, file extension, and file path.
     */
    struct IMG_API specification {
        uint8* data;                    /**< Pointer to image data */
        uint32 width;                   /**< Width of the image */
        uint32 height;                  /**< Height of the image */
        uint32 channels;                /**< Number of channels in the image */
        std::string extension;          /**< File extension of the image */
        std::filesystem::path path;     /**< File path of the image */

        specification() = default;
        ~specification() = default;
    };

    /**
     * @brief Reads an image from a file.
     *
     * This function reads an image from the specified file path and returns
     * a shared pointer to a `specification` struct containing the image metadata.
     *
     * @param filePath The path to the image file.
     * @param flip Whether to flip the image vertically.
     * @return A shared pointer to a `specification` struct with the image metadata.
     */
    IMG_API std::shared_ptr<specification> read(const std::filesystem::path& filePath, bool flip = false);
}



