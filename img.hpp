#pragma once

#include <string>
#include <filesystem>
#include <memory>

#ifdef IMG_DEBUG
#include <stdio.h>
#define IMG_DEBUG_LOG(...) printf(__VA_ARGS__)
#else
#define IMG_DEBUG_LOG(...)
#endif

namespace img{

    /**
     * @struct image_specification
     * @brief This structure defines the specifications of an image.
     */
    struct image_specification{
        int32_t width{0};           ///< Width of the image in pixels.
        int32_t height{0};          ///< Height of the image in pixels.
        int32_t channels{0};        ///< Number of color channels in the image.
        uint8_t* data{nullptr};     ///< Pointer to the image data.
        uint32_t xid{0};            ///< Unique identifier for the image.

        std::filesystem::path filepath; ///< Filepath where the image is stored.
        std::string filename;           ///< Filename of the image.
    };

    /**
     * @enum image_format
     * @brief Enum class representing different image formats.
     */
    enum class image_format{
        png,    ///< Portable Network Graphics format.
        jpeg,   ///< Joint Photographic Experts Group format.
        bmp,    ///< Bitmap format.
        gif,    ///< Graphics Interchange Format.
        webp    ///< Web Picture format.
    };

    /**
     * @brief Loads an image from the given filepath.
     * @param filepath Filepath of the image to be loaded.
     * @param flip Whether to flip the image vertically after loading (default is true).
     * @return A std::shared_ptr to the loaded image.
     */
    [[nodiscard]] std::shared_ptr<image_specification> import(const std::filesystem::path& filepath, image_format format, bool flip = true) noexcept;
}