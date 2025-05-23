/*
MIT License

Copyright (c) 2025 HobyProjects

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
#pragma once

#include <filesystem>
#include <memory>
#include <string>

#ifdef IMG_DEBUG
#include <stdio.h>
#define IMG_DEBUG_LOG(...) printf(__VA_ARGS__)
#define IMG_ASSERT(exp, ...) \
  if (!(exp)) {              \
    printf(__VA_ARGS__);     \
    std::terminate();        \
  }
#else
#define IMG_DEBUG_LOG(...)
#define IMG_ASSERT(exp, ...)
#endif

namespace img {

/**
 * @struct image_specification
 * @brief This structure defines the specifications of an image.
 */
struct specification {
  int32_t width{0};        ///< Width of the image in pixels.
  int32_t height{0};       ///< Height of the image in pixels.
  int32_t channels{0};     ///< Number of color channels in the image.
  uint8_t* data{nullptr};  ///< Pointer to the image data.
  uint32_t xid{0};         ///< Unique identifier for the image.

  std::filesystem::path filepath;  ///< Filepath where the image is stored.
  std::string filename;            ///< Filename of the image.
};

/**
 * @enum format
 * @brief Enum class representing different image formats.
 */
enum class format {
  png,   ///< Portable Network Graphics format.
  jpeg,  ///< Joint Photographic Experts Group format.
  bmp,   ///< Bitmap format.
  gif,   ///< Graphics Interchange Format.
  webp   ///< Web Picture format.
};

/**
 * @enum mode
 * @brief Enum class representing the image read modes.
 * read_nessessary: Only reads the necessary chunks for the image.
 * read_all: Reads all chunks in the image.
 */
enum class mode { nessessary, fullbreakdown };

/**
 * @brief Import an image from the given filepath with the specified format and read mode.
 *
 * @param filepath The filepath of the image to import.
 * @param format The image format of the file to import.
 * @param mode The image read mode.
 * @param flip Whether to flip the image vertically. Default is true.
 * @return A shared pointer to the image specification if the import is successful, otherwise
 * nullptr.
 * @throws std::runtime_error if the image could not be imported.
 */
[[nodiscard]] std::shared_ptr<specification>
import(const std::filesystem::path& filepath, format format, mode mode, bool flip = true) noexcept;
}  // namespace img