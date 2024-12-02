/**
 * @file img.hpp
 * @brief Image library header.
 *
 * This file provides the API for the image library, including functions for
 * reading and writing images, as well as data structures for representing
 * image metadata.
 *
 * Copyright (c) 2024, HobyProjects
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

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

/**
 * @brief A bit reader for reading bits from a byte array.
 *
 * This class provides methods for reading bits from a byte array. It is
 * used to read bits from image data in the `img::read` function.
 */
class img_bitreader{
    public:
        /**
         * @brief Default constructor.
         *
         * This constructor initializes the bit reader with a nullptr data
         * pointer and a data size of 0.
         */
        img_bitreader() = default;

        /**
         * @brief Constructor that initializes the bit reader with a data
         *        pointer and size.
         *
         * This constructor initializes the bit reader with the specified data
         * pointer and size.
         *
         * @param dptr The data pointer to read bits from.
         * @param dsize The size of the data in bytes.
         */
        img_bitreader(uint8* dptr, uint32 dsize) : data_ptr(dptr), data_size(dsize) {}

        /**
         * @brief Destructor.
         *
         * The destructor is currently empty, as the bit reader does not
         * allocate any memory.
         */
        ~img_bitreader() = default;

        /**
         * @brief Reads a bit field from the data.
         *
         * This function reads a bit field of the size specified by the
         * `SizeT` template argument from the current position in the data.
         * The bit field is read in little-endian order, i.e. the least
         * significant bit of the first byte is read first, followed by the
         * next least significant bit, and so on.
         *
         * @return The bit field read from the data.
         */
        template<typename SizeT>
        SizeT read_bits(){
            SizeT result = 0;
            for(uint32 i = 0; i < sizeof(SizeT) * 8; i++){
                // Check if we have reached the end of the data.
                if(index >= data_size) break;

                // Read the least significant bit of the current byte.
                result |= (SizeT)(data_ptr[index] & 1) << i;

                // Advance to the next byte.
                index++;
            }

            return result;
        }

    private:
        /// The pointer to the data that we are reading from.
        /// This is set by the constructor and is not modified
        /// by the bit reader.
        uint8* data_ptr{nullptr};

        /// The current index into the data. This is incremented
        /// each time we read a byte.
        uint32 index{0};

        /// The size of the data in bytes. This is set by the
        /// constructor and is not modified by the bit reader.
        uint32 data_size{0};
};

/******************************************************************
 *                         PNG IMAGE FILES                        *
 * ****************************************************************/

#define IMG_PNG_SIGNATURE   0x89504E470D0A1A0A  ///< The PNG file signature.
#define IMG_PNG_CHUNK_IHDR  0x52414449          ///< The IHDR chunk type.
#define IMG_PNG_CHUNK_IEND  0x444E4549          ///< The IEND chunk type
#define IMG_PNG_CHUNK_IDAT  0x54414449          ///< The IDAT chunk type
#define IMG_PNG_CHUNK_PLTE  0x504C5445          ///< The PLTE chunk type

/**
 * @brief Structure representing a PNG image chunk.
 *
 * This structure represents a single chunk in a PNG image file.
 * Each chunk contains a length, a type, a pointer to the data
 * and a CRC checksum.
 */
struct IMG_PNGCHUNK {
    uint32 length;              ///< The length of the chunk.
    uint32 type;                ///< The type of the chunk.
    uint8* data;                ///< The data of the chunk.
    uint32 crc;                 ///< The CRC checksum of the chunk.
    IMG_PNGCHUNK* next;         ///< The next chunk in the list.

    /**
     * @brief Default constructor.
     *
     * This is the default constructor for the `IMG_PNGCHUNK` struct.
     * It is not virtual, since this struct is not intended to be
     * inherited from.
     */
    IMG_PNGCHUNK() = default;

    /**
     * @brief Destructor.
     *
     * This is the destructor for the `IMG_PNGCHUNK` struct. It is not
     * virtual, since this struct is not intended to be inherited from.
     */
    ~IMG_PNGCHUNK() = default;
};

/**
 * @brief Structure representing a PNG image header.
 *
 * This structure holds the metadata of a PNG image, including its
 * dimensions, color type, and compression settings.
 */
struct IMG_PNGIHDR {
    
    uint32 width{0};               ///< The width of the image in pixels.
    uint32 height{0};              ///< The height of the image in pixels.
    uint8 channels{0};             ///< The number of channels in the image.
    uint8 color_type{0};           ///< The color type of the image.
    uint8 compression_method{0};   ///< The compression method used.
    uint8 filter_method{0};        ///< The filter method used.
    uint8 interlace_method{0};     ///< The interlace method used.

    /**
     * @brief Default constructor.
     *
     * This constructor initializes an IMG_PNGIHDR instance
     * with default values for its members.
     */
    IMG_PNGIHDR() = default;

    /**
     * @brief Destructor.
     *
     * This destructor cleans up resources used by the IMG_PNGIHDR
     * instance. Currently, it does not perform any actions as
     * no dynamic memory allocation is involved.
     */
    ~IMG_PNGIHDR() = default;
};

/**
 * @enum IMG_PNGCOLORTYPE
 * @brief Enum representing the color types of PNG images.
 *
 * This enum represents the color types of PNG images. The possible values are
 * as follows:
 *
 * - `GRAYSCALE`: A greyscale image with a single channel.
 * 
 * - `TRUECOLOR`: A full-color image with three channels.
 * 
 * - `INDEXED_COLOR`: An indexed image with a single channel.
 * 
 * - `GRAYSCALE_ALPHA`: A greyscale image with alpha channel.
 * 
 * - `TRUECOLOR_ALPHA`: A full-color image with alpha channel.
 */
enum class IMG_PNGCOLORTYPE {
    GRAYSCALE           = 0,  ///< A greyscale image with a single channel.
    TRUECOLOR           = 2,  ///< A full-color image with three channels.
    INDEXED_COLOR       = 3,  ///< An indexed image with a single channel.
    GRAYSCALE_ALPHA     = 4,  ///< A greyscale image with alpha channel.
    TRUECOLOR_ALPHA     = 6   ///< A full-color image with alpha channel.
};


namespace img {

    /**
     * @brief Represents an image image_specification.
     *
     * This struct represents an image image_specification. It contains the following
     * fields:
     *
     * - `data`: A pointer to the image data.
     * 
     * - `width`: The width of the image in pixels.
     * 
     * - `height`: The height of the image in pixels.
     * 
     * - `channels`: The number of channels in the image. This can be 1 (grayscale), 3 (RGB), or 4 (RGBA).
     * 
     * - `extension`: The file extension of the image file (e.g. "png", "jpg", etc.).
     * 
     * - `path`: The file path of the image file.
     *
     * The `image_specification` struct is used to store the results of reading an
     * image from a file. It is returned by the `img::read` function.
     */
    struct IMG_API image_specification {
        /**
         * @brief A pointer to the image data.
         *
         * This field contains a pointer to the image data. The image data is
         * stored in row-major order, meaning that the first row of the image
         * is stored first, followed by the second row, and so on.
         *
         * The memory layout of the image data is as follows:
         *
         * - Each row of the image consists of `width * channels` bytes of data.
         * 
         * - Each pixel of the image consists of `channels` bytes of data.
         * 
         * - The first byte of each pixel represents the red channel, the second
         *   byte represents the green channel, the third byte represents the
         *   blue channel, and the fourth byte represents the alpha channel.
         */
        uint8* data;

        /**
         * @brief The width of the image in pixels.
         *
         * This field contains the width of the image in pixels.
         */
        uint32 width;

        /**
         * @brief The height of the image in pixels.
         *
         * This field contains the height of the image in pixels.
         */
        uint32 height;

        /**
         * @brief The number of channels in the image.
         *
         * This field contains the number of channels in the image. The number of
         * channels can be 1 (grayscale), 3 (RGB), or 4 (RGBA).
         */
        uint32 channels;

        /**
         * @brief The file extension of the image file (e.g. "png", "jpg", etc.).
         *
         * This field contains the file extension of the image file. This can be
         * used to determine the type of the image file.
         */
        std::string extension;

        /**
         * @brief The file path of the image file.
         *
         * This field contains the file path of the image file. This can be used
         * to open the image file and read its contents.
         */
        std::filesystem::path path;

        /**
         * @brief Default constructor.
         */
        image_specification() = default;

        /**
         * @brief Destructor.
         *
         * This is the destructor for the `image_specification` struct. It is not
         * virtual, since this struct is not intended to be inherited from.
         */
        ~image_specification() = default;
    };


    /**
     * @brief Reads an image from a file.
     *
     * This function attempts to read an image from the specified file path
     * and returns a shared pointer to a `image_specification` struct containing
     * the image's metadata.
     *
     * @param filePath The path to the image file.
     * @param flip A boolean indicating whether to flip the image vertically.
     *             Defaults to false.
     * @return A shared pointer to a `image_specification` struct with the image
     *         metadata, or nullptr if the file could not be read.
     */
    IMG_API std::shared_ptr<image_specification> read(const std::filesystem::path& filePath, bool flip = false);
}