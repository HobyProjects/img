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
#include "img.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

/**
 * @class bit_reader
 * @brief A class for reading bits from a memory block.
 *
 * This class provides a convenient way to read bits from a memory block
 * one by one. It handles the bit shifting and masking for you, and it
 * also handles the case where the bit offset is not a multiple of 8.
 */
class bit_reader {
 public:
  explicit bit_reader(uint8_t* data, size_t size) : m_data(data), m_size(size), m_bit_offset(0) {}
  bit_reader(const bit_reader&) = delete;
  bit_reader& operator=(const bit_reader&) = delete;
  bit_reader(bit_reader&&) = default;
  bit_reader& operator=(bit_reader&&) = default;
  ~bit_reader() = default;

  /**
   * @brief Reads a integral value of type T from the bit stream.
   * @details
   * This function reads sizeof(T) bytes from the bit stream and
   * returns the result as a value of type T. The least significant
   * bit of the first byte is considered the most significant bit of
   * the result, and the most significant bit of the last byte is
   * considered the least significant bit of the result.
   * @tparam T The type of the value to read. Must be an integral type.
   * @return The value read from the bit stream.
   */
  template <typename T>
  T read() {
    static_assert(std::is_integral_v<T>, "T must be an integral type");
    T result = 0;
    for (size_t i = 0; i < sizeof(T); ++i) {
      result = (result << 8) | m_data[m_bit_offset >> 3];
      m_bit_offset += 8;
    }
    return result;
  }

  /**
   * @brief Reads a byte array of size bytes from the bit stream.
   * @details
   * This function reads size bytes from the bit stream and returns
   * a pointer to the byte array. The least significant bit of the
   * first byte is considered the most significant bit of the first
   * byte in the byte array, and the most significant bit of the
   * last byte is considered the least significant bit of the last
   * byte in the byte array.
   * @param size The size of the byte array to read.
   * @return A pointer to the byte array read from the bit stream.
   */
  uint8_t* read(uint32_t size) noexcept {
    if (m_bit_offset + size > m_size * 8) {
      IMG_DEBUG_LOG("An Attempt to read beyond buffer size\n");
      return nullptr;
    }

    if (size == 0) {
      IMG_DEBUG_LOG("An Attempt to read zero bytes\n");
      return nullptr;
    }

    uint8_t* result = new uint8_t[size];
    for (uint32_t i = 0; i < size; ++i) {
      result[i] = read<uint8_t>();
    }
    return result;
  }

 private:
  uint8_t* m_data{nullptr};  // Pointer to the data buffer
  size_t m_size{0};          // Size of the data buffer in bytes
  size_t m_bit_offset{0};    // Current bit offset in the data buffer
};

/**
 * @struct compressed_file_data
 * @brief This structure represents a compressed file.
 * It contains the filepath, filename, data pointer, and size of the file.
 */
struct compressed_file_data {
  std::filesystem::path filepath;  ///< Filepath of the compressed file
  std::string filename;            ///< Filename of the compressed file
  uint8_t* data{nullptr};          ///< Pointer to the compressed file data
  size_t size{0};                  ///< Size of the compressed file data in bytes

  explicit compressed_file_data(const std::filesystem::path& file) {
    filepath = file;
    filename = file.filename().string();
    std::ifstream file_stream(filepath, std::ios::binary | std::ios::ate);
    if (!file_stream.is_open()) {
      IMG_DEBUG_LOG("Failed to open file: %s\n", filepath.string().c_str());
      throw std::runtime_error("Failed to open file: " + filepath.string());
    }
    size = file_stream.tellg();           // Get the size of the file
    file_stream.seekg(0, std::ios::beg);  // Move the cursor to the beginning of the file
    data = new uint8_t[size];             // Allocate memory for the file data
    file_stream.read(reinterpret_cast<char*>(data),
                     size);  // Read the file data into the buffer
    file_stream.close();     // Close the file stream
  };

  compressed_file_data(const compressed_file_data&) = delete;  // Delete copy constructor
  compressed_file_data&
  operator=(const compressed_file_data&) = delete;         // Delete copy assignment operator
  compressed_file_data(compressed_file_data&&) = default;  // Default move constructor
  compressed_file_data&
  operator=(compressed_file_data&&) = default;  // Default move assignment operator

  ~compressed_file_data() {
    delete[] data;  // Free the allocated memory for the file data
  }

  void clear() noexcept {
    delete[] data;   // Free the allocated memory for the file data
    data = nullptr;  // Set the pointer to nullptr
    size = 0;        // Reset the size to 0
  }
};

/**
 * @brief Generates a unique ID for the image.
 * @details This function generates a unique ID for the image by incrementing
 * a static variable. The ID is used to identify the image in the system.
 * @return A unique ID for the image.
 */
static uint32_t unique_xid() noexcept {
  static uint32_t s_xid = 0;  // Static variable to hold the unique ID
  return ++s_xid;             // Increment and return the unique ID
}

/**
 * @struct img_rgb
 * @brief This structure represents an RGB color.
 * It contains the red, green, and blue channel values.
 */
struct img_rgb {
  uint8_t r{0};  // Red channel value
  uint8_t g{0};  // Green channel value
  uint8_t b{0};  // Blue channel value

  img_rgb() = default;  // Default constructor
  img_rgb(uint8_t red, uint8_t green, uint8_t blue) : r(red), g(green), b(blue) {}
  ~img_rgb() = default;  // Default destructor
};

/******************************************************************
 *                         PNG IMAGE FILES                        *
 * ****************************************************************/
static constexpr uint64_t s_IMG_PNG_SIGNATURE = 0x89504E470D0A1A0A;  // PNG file signature
static constexpr uint32_t s_IMG_PNG_IHDR_CHUNK = 0x52414449;         // IHDR chunk type
static constexpr uint32_t s_IMG_PNG_IEND_CHUNK = 0x444E4549;         // IEND chunk type
static constexpr uint32_t s_IMG_PNG_IDAT_CHUNK = 0x54455849;         // IDAT chunk type
static constexpr uint32_t s_IMG_PNG_PLTE_CHUNK = 0x54494C50;         // PLTE chunk type
static constexpr uint32_t s_IMG_PNG_TRNS_CHUNK = 0x74524E53;         // tRNS chunk type
static constexpr uint32_t s_IMG_PNG_GAMA_CHUNK = 0x67414D41;         // gama chunk type
static constexpr uint32_t s_IMG_PNG_CHRM_CHUNK = 0x6348524D;         // chrm chunk type
static constexpr uint32_t s_IMG_PNG_SRGB_CHUNK = 0x73524742;         // srgb chunk type
static constexpr uint32_t s_IMG_PNG_BKGD_CHUNK = 0x624B4744;         // bKGD chunk type
static constexpr uint32_t s_IMG_PNG_PHYS_CHUNK = 0x70485973;         // pHYs chunk type
static constexpr uint32_t s_IMG_PNG_ITXT_CHUNK = 0x69545874;         // iTXt chunk type
static constexpr uint32_t s_IMG_PNG_TEXT_CHUNK = 0x74455874;         // tEXt chunk type
static constexpr uint32_t s_IMG_PNG_ZTXT_CHUNK = 0x7A545874;         // zTXt chunk type
static constexpr uint32_t s_IMG_PNG_SBIT_CHUNK = 0x73424954;         // sBIT chunk type
static constexpr uint32_t s_IMG_PNG_ICCP_CHUNK = 0x69434350;         // iccp chunk type
static constexpr uint32_t s_IMG_PNG_FRAC_CHUNK = 0x66726163;         // fRAc chunk type
static constexpr uint32_t s_IMG_PNG_HIST_CHUNK = 0x68495354;         // hIST chunk type
static constexpr uint32_t s_IMG_PNG_TIME_CHUNK = 0x74494D45;         // tIME chunk type
static constexpr uint32_t s_IMG_PNG_PCAL_CHUNK = 0x7043414C;         // pCAL chunk type
static constexpr uint32_t s_IMG_PNG_SCAL_CHUNK = 0x7343414C;         // sCAL chunk type

static constexpr uint32_t s_CRC32_POLYNOMIAL = 0xEDB88320;  // CRC32 polynomial
static constexpr uint32_t s_CRC_TABLE_SIZE = 256;           // Size of the CRC table

static uint32_t s_CRC_TABLE[s_CRC_TABLE_SIZE];  // CRC table
static std::once_flag s_CRC_TBL_INIT_FLAG;      // Flag to ensure CRC table is initialized only once

/**
 * @brief Generates the CRC table for CRC32 calculation.
 * @details This function generates the CRC table used for CRC32 calculation.
 * The table is generated using the polynomial defined in s_CRC32_POLYNOMIAL.
 */
static void img_png_generate_crc_table() noexcept {
  for (uint32_t i = 0; i < s_CRC_TABLE_SIZE; ++i) {
    uint32_t crc = i;
    for (uint32_t j = 0; j < 8; ++j) {
      if (crc & 1) {
        crc = (crc >> 1) ^ s_CRC32_POLYNOMIAL;
      } else {
        crc >>= 1;
      }
    }
    s_CRC_TABLE[i] = crc;
  }
}

/**
 * @brief Calculates the CRC32 checksum for a given data buffer.
 * @details This function calculates the CRC32 checksum for a given data buffer
 * using the CRC table generated by img_png_generate_crc_table().
 * @param data Pointer to the data buffer.
 * @param length Length of the data buffer in bytes.
 * @return The calculated CRC32 checksum.
 */
static uint32_t img_png_crc32(const uint8_t* data, size_t length) noexcept {
  uint32_t crc = 0xFFFFFFFF;  // Initial CRC value
  for (size_t i = 0; i < length; ++i) {
    uint8_t byte = data[i];
    crc = (crc >> 8) ^ s_CRC_TABLE[(crc ^ byte) & 0xFF];  // Update CRC value using the table
  }
  return ~crc;  // Return the final CRC value
}

/**
 * @struct img_png_chunk
 * @brief This structure represents a PNG chunk.
 * It contains information about the chunk's length, type, data, CRC value,
 * and pointers to the next and previous chunks in the linked list.
 * It also contains a pointer to the parent chunk (if any).
 */
struct img_png_chunk {
  uint32_t length{0};  // Length of the chunk data
  uint32_t type{0};    // Type of the chunk
  uint32_t crc{0};     // CRC value for the chunk

  img_png_chunk() = default;
  virtual ~img_png_chunk() = default;
};

using img_chunk_ptr = std::shared_ptr<img_png_chunk>;  // Pointer to a PNG chunk

/**
 * @struct img_png_type_chunk_map
 * @brief This structure maps a chunk type to its corresponding chunk.
 * It contains the chunk type and a pointer to the chunk.
 */
struct img_png_type_chunk_map {
  uint32_t type{0};               // Type of the chunk
  img_chunk_ptr chunks{nullptr};  // Pointer to the chunk

  img_png_type_chunk_map() = default;
  img_png_type_chunk_map(uint32_t t, img_chunk_ptr&& c) {
    type = t;
    chunks = std::move(c);
  }
  ~img_png_type_chunk_map() = default;
};

/**
 * @struct img_png_ihdr_chunck
 * @brief This structure represents the IHDR chunk of a PNG image.
 * It contains information about the image's width, height, bit depth,
 * color type, compression method, filter method, and interlace method.
 */
struct img_png_ihdr_chunck : img_png_chunk {
  uint32_t width{0};              // Width of the image in pixels
  uint32_t height{0};             // Height of the image in pixels
  uint8_t bit_depth{0};           // Bit depth of the image
  uint8_t color_type{0};          // Color type of the image
  uint8_t compression_method{0};  // Compression method of the image
  uint8_t filter_method{0};       // Filter method of the image
  uint8_t interlace_method{0};    // Interlace method of the image

  img_png_ihdr_chunck() = default;
  virtual ~img_png_ihdr_chunck() = default;
};

struct img_png_iend_chunk : img_png_chunk {
  bool is_iend{false};
  img_png_iend_chunk() = default;
  virtual ~img_png_iend_chunk() = default;
};

/**
 * @struct img_png_idat_chunk
 * @brief This structure represents the IDAT chunk of a PNG image.
 * It contains a pointer to the IDAT chunk data and its length.
 */
struct img_png_idat_chunk : img_png_chunk {
  uint8_t* data{nullptr};                  // Pointer to the IDAT chunk data
  std::weak_ptr<img_png_ihdr_chunck> ihdr;  // Pointer to the IHDR chunk

  img_png_idat_chunk() = default;
  virtual ~img_png_idat_chunk() {
    delete[] data;  // Free the IDAT chunk data
  }
};

/**
 * @struct img_png_plte_chunk
 * @brief This structure represents the PLTE chunk of a PNG image.
 * It contains a vector of RGB colors and the length of the PLTE chunk data.
 */
struct img_png_plte_chunk : img_png_chunk {
  std::vector<img_rgb> palette;  // Vector to hold the palette colors

  img_png_plte_chunk() = default;
  virtual ~img_png_plte_chunk() = default;
};

/**
 * @struct img_png_trns_chunk
 * @brief This structure represents the tRNS chunk of a PNG image.
 * It contains a vector of alpha values and the length of the tRNS chunk data.
 */
struct img_png_trns_chunk : img_png_chunk {
  std::vector<uint8_t> alpha;  // Vector to hold the alpha values

  img_png_trns_chunk() = default;
  virtual ~img_png_trns_chunk() = default;
};

struct img_png_gama_chunk;  // Forward declaration of img_png_gama_chunk
struct img_png_chrm_chunk;  // Forward declaration of img_png_chrm_chunk
struct img_png_iccp_chunk;  // Forward declaration of img_png_iccp_chunk
struct img_png_srgb_chunk;  // Forward declaration of img_png_srgb_chunk

/**
 * @struct img_png_gama_chunk
 * @brief This structure represents the gama chunk of a PNG image.
 * It contains the gamma value for the image.
 */
struct img_png_gama_chunk : img_png_chunk {
  uint32_t gamma{0};      // Gamma value for the image
  bool has_gamma{false};  // Flag to indicate if gama chunk is present

  std::weak_ptr<img_png_chrm_chunk> chrm;  // Pointer to the chrm chunk (if present)
  std::weak_ptr<img_png_iccp_chunk> iccp;  // Pointer to the iccp chunk (if present)
  std::weak_ptr<img_png_srgb_chunk> srgb;  // Pointer to the srgb chunk (if present)

  img_png_gama_chunk() = default;
  virtual ~img_png_gama_chunk() = default;
};

/**
 * @struct img_png_chrm_chunk
 * @brief This structure represents the chrm chunk of a PNG image.
 * It contains the chromaticity values for the image.
 */
struct img_png_chrm_chunk : img_png_chunk {
  float white_x{0.0f}, white_y{0.0f};
  float red_x{0.0f}, red_y{0.0f};
  float green_x{0.0f}, green_y{0.0f};
  float blue_x{0.0f}, blue_y{0.0f};
  bool has_chrm{false};

  std::weak_ptr<img_png_gama_chunk> gama;  // Pointer to the gama chunk (if present)
  std::weak_ptr<img_png_iccp_chunk> iccp;  // Pointer to the iccp chunk (if present)
  std::weak_ptr<img_png_srgb_chunk> srgb;  // Pointer to the srgb chunk (if present)

  img_png_chrm_chunk() = default;
  virtual ~img_png_chrm_chunk() = default;
};

/**
 * @struct img_png_iccp_chunk
 * @brief This structure represents the iccp chunk of a PNG image.
 * It contains the ICC profile data for the image.
 */
struct img_png_iccp_chunk : img_png_chunk {
  std::string profile_name;                 // Name of the ICC profile
  uint8_t compression_method{0};            // Compression method of the ICC profile
  std::vector<uint8_t> compressed_profile;  // Compressed ICC profile data
  bool has_iccp{false};

  std::weak_ptr<img_png_gama_chunk> gama;  // Pointer to the gama chunk (if present)
  std::weak_ptr<img_png_chrm_chunk> chrm;  // Pointer to the chrm chunk (if present)
  std::weak_ptr<img_png_srgb_chunk> srgb;  // Pointer to the srgb chunk (if present)

  img_png_iccp_chunk() = default;
  virtual ~img_png_iccp_chunk() = default;
};

/**
 * @enum img_png_srgb_rendering_intent
 * @brief This enum represents the rendering intent for the srgb chunk of a PNG
 * image.
 */
enum class img_png_srgb_rendering_intent : uint8_t {
  perceptual = 0,
  relative_colorimetic = 1,
  saturation = 2,
  absolute_colorimetric = 3
};

/**
 * @struct img_png_srgb_chunk
 * @brief This structure represents the srgb chunk of a PNG image.
 * It contains the rendering intent for the image.
 */
struct img_png_srgb_chunk : img_png_chunk {
  img_png_srgb_rendering_intent rendering_intent{img_png_srgb_rendering_intent::perceptual};
  bool has_srgb{false};

  std::weak_ptr<img_png_gama_chunk> gama;  // Pointer to the gama chunk (if present)
  std::weak_ptr<img_png_iccp_chunk> iccp;  // Pointer to the iccp chunk (if present)
  std::weak_ptr<img_png_chrm_chunk> chrm;  // Pointer to the chrm chunk (if present)

  img_png_srgb_chunk() = default;
  virtual ~img_png_srgb_chunk() = default;
};

/**
 * @struct img_png_bkgd_chunk
 * @brief This structure represents the bKGD chunk of a PNG image.
 * It contains the background color for the image.
 */
struct img_png_bkgd_chunk : img_png_chunk {
  uint16_t gray{0};           // For grayscale
  uint16_t r{0}, g{0}, b{0};  // For truecolor
  uint8_t palette_index{0};   // For indexed
  uint8_t color_type{255};    // To track what kind of color this represents

  img_png_plte_chunk* plte;   // Pointer to the plte chunk (if present)
  img_png_ihdr_chunck* ihdr;  // Pointer to the IHDR chunk

  img_png_bkgd_chunk() = default;
  virtual ~img_png_bkgd_chunk() = default;
};

/**
 * @brief Verifies the PNG signature.
 * @param signature The PNG signature to verify.
 * @return True if the signature is valid, false otherwise.
 */
static bool img_png_verify_signature(uint64_t signature) noexcept {
  return signature == s_IMG_PNG_SIGNATURE;
}

/**
 * @brief Verifies the CRC32 checksum of a PNG chunk.
 * @details This function calculates the CRC32 checksum for a given PNG chunk
 * consisting of a type and data, and compares it to the provided CRC value.
 * It constructs the input for the CRC calculation by concatenating the chunk
 * type and data, calculates the CRC, and returns the calculated CRC if it matches
 * the provided CRC, otherwise returns 0.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param crc The CRC value to verify against.
 * @param data Pointer to the chunk data.
 * @return The calculated CRC32 value if it matches the provided CRC, 0 otherwise.
 */
static uint32_t img_png_verify_chunk_crc(uint32_t type, uint32_t len, uint32_t crc,
                                         uint8_t* data) noexcept {
  size_t crc_input_size = 4 + len;
  uint8_t* crc_input = new uint8_t[crc_input_size];

  // Write type (big endian) first
  crc_input[0] = (type >> 24) & 0xFF;
  crc_input[1] = (type >> 16) & 0xFF;
  crc_input[2] = (type >> 8) & 0xFF;
  crc_input[3] = type & 0xFF;
  memcpy(crc_input + 4, data, len);

  // Now calculate CRC over type + data
  uint32_t calculated_crc = img_png_crc32(crc_input, crc_input_size);
  delete[] crc_input;

  return (calculated_crc == crc) ? calculated_crc : 0;
}

/**
 * @brief Reads the IHDR chunk of a PNG image.
 * @details This function reads the IHDR chunk of a PNG image and populates an
 * img_png_ihdr_chunck structure with the information. It verifies the CRC32
 * checksum of the chunk, and returns a pointer to the structure if the chunk
 * is read successfully, otherwise returns nullptr.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @param chunck_map The vector to store the read chunk.
 * @return A pointer to the img_png_ihdr_chunck structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_ihdr_chunck>
img_png_read_ihdr_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc,
                        std::vector<std::shared_ptr<img_png_type_chunk_map>>& chunck_map) noexcept {
  std::shared_ptr<img_png_ihdr_chunck> ihdr = std::make_shared<img_png_ihdr_chunck>();

  ihdr->width = (data[0] << 24) | (data[1] << 16) | (data[2] << 8) | data[3];
  ihdr->height = (data[4] << 24) | (data[5] << 16) | (data[6] << 8) | data[7];
  ihdr->bit_depth = data[8];
  ihdr->color_type = data[9];
  ihdr->compression_method = data[10];
  ihdr->filter_method = data[11];
  ihdr->interlace_method = data[12];

  ihdr->length = len;
  ihdr->crc = crc;
  ihdr->type = type;

  std::shared_ptr<img_png_type_chunk_map> ihdr_map =
      std::make_shared<img_png_type_chunk_map>(type, std::move(ihdr));
  chunck_map.emplace_back(std::move(ihdr_map));

  IMG_DEBUG_LOG(
      "IHDR chunk read successfully,\n\t\twidth: %d pxels x height: %d pxels (Total: %d "
      "pxels)\n\t\tbit_depth: %d\n"
      "\t\tcolor_type: %d\n\t\tcompression_method: %d\n\t\tfilter_method: %d\n"
      "\t\tinterlace_method: %d\n",
      ihdr->width, ihdr->height, ihdr->width * ihdr->height, ihdr->bit_depth, ihdr->color_type,
      ihdr->compression_method, ihdr->filter_method, ihdr->interlace_method);

  return ihdr;
}

/**
 * @brief Reads the IDAT chunk of a PNG image.
 * @details This function reads the IDAT chunk of a PNG image and populates an
 * img_png_idat_chunk structure with the information. It verifies the CRC32
 * checksum of the chunk, and returns a pointer to the structure if the chunk
 * is read successfully, otherwise returns nullptr.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @param chunck_map The vector to store the read chunk.
 * @return A pointer to the img_png_idat_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_idat_chunk>
img_png_read_idat_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc,
                        std::vector<std::shared_ptr<img_png_type_chunk_map>>& chunck_map) noexcept {
  std::shared_ptr<img_png_idat_chunk> idat = std::make_shared<img_png_idat_chunk>();
  idat->length = len;
  idat->data = data;

  idat->crc = crc;
  idat->type = type;
  idat->length = len;

  std::shared_ptr<img_png_type_chunk_map> idat_map =
      std::make_shared<img_png_type_chunk_map>(type, std::move(idat));
  chunck_map.emplace_back(std::move(idat_map));

  IMG_DEBUG_LOG("IDAT chunk read successfully\n,\t\tlength: %d\n", len);
  return idat;
}

/**
 * @brief Reads the PLTE chunk of a PNG image.
 * @details This function reads the PLTE chunk of a PNG image and populates an
 * img_png_plte_chunk structure with the information. It verifies the CRC32
 * checksum of the chunk, and returns a pointer to the structure if the chunk
 * is read successfully, otherwise returns nullptr.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @param chunck_map The vector to store the read chunk.
 * @return A pointer to the img_png_plte_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_plte_chunk>
img_png_read_plte_chunck(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc,
                         std::vector<std::shared_ptr<img_png_type_chunk_map>>& chunk_map) noexcept {
  std::shared_ptr<img_png_plte_chunk> plte = std::make_shared<img_png_plte_chunk>();
  plte->palette.reserve(len / 3);
  for (uint32_t i = 0; i < len; i += 3) {
    plte->palette.emplace_back(img_rgb(data[i], data[i + 1], data[i + 2]));
  }

  plte->length = len;
  plte->crc = crc;
  plte->type = type;

  std::shared_ptr<img_png_type_chunk_map> plte_map =
      std::make_shared<img_png_type_chunk_map>(type, std::move(plte));
  chunk_map.emplace_back(std::move(plte_map));

  IMG_DEBUG_LOG("PLTE chunk read successfully\n,\t\tlength: %d\n", len);
  return plte;
}

/**
 * @brief Reads the tRNS chunk of a PNG image.
 * @details This function reads the tRNS chunk of a PNG image and populates an
 * img_png_trns_chunk structure with the information. It verifies the CRC32
 * checksum of the chunk, and returns a pointer to the structure if the chunk
 * is read successfully, otherwise returns nullptr.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @param chunck_map The vector to store the read chunk.
 * @return A pointer to the img_png_trns_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_trns_chunk>
img_png_read_trns_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc,
                        std::vector<std::shared_ptr<img_png_type_chunk_map>>& chunk_map) noexcept {
  std::shared_ptr<img_png_trns_chunk> tRNS = std::make_shared<img_png_trns_chunk>();
  tRNS->alpha.reserve(len);
  for (uint32_t i = 0; i < len; ++i) {
    tRNS->alpha.emplace_back(data[i]);
  }

  tRNS->length = len;
  tRNS->crc = crc;
  tRNS->type = type;

  std::shared_ptr<img_png_type_chunk_map> tRNS_map =
      std::make_shared<img_png_type_chunk_map>(type, std::move(tRNS));
  chunk_map.emplace_back(std::move(tRNS_map));

  IMG_DEBUG_LOG("TRNS chunk read successfully\n,\t\tlength: %d\n", len);
  return tRNS;
}

/**
 * @brief Reads the gAMA chunk of a PNG image.
 * @details This function reads the gAMA chunk of a PNG image and populates an
 * img_png_gama_chunk structure with the information. It verifies the CRC32
 * checksum of the chunk, and returns a pointer to the structure if the chunk
 * is read successfully, otherwise returns nullptr.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @param chunck_map The vector to store the read chunk.
 * @return A pointer to the img_png_gama_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_gama_chunk>
img_png_read_gama_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc,
                        std::vector<std::shared_ptr<img_png_type_chunk_map>>& chunk_map) noexcept {
  std::shared_ptr<img_png_gama_chunk> gama = std::make_shared<img_png_gama_chunk>();

  if (len != 4 || len != 6) {
    gama->has_gamma = false;
    IMG_DEBUG_LOG("Invalid gama chunk length: %x, expected: 4 or 6\n", len);
    return nullptr;
  } else {
    gama->has_gamma = true;
    gama->gamma = (data[0] << 24) | (data[1] << 16) | (data[2] << 8) | data[3];

    gama->length = len;
    gama->crc = crc;
    gama->type = type;
  }

  std::shared_ptr<img_png_type_chunk_map> gama_map =
      std::make_shared<img_png_type_chunk_map>(type, std::move(gama));
  chunk_map.emplace_back(std::move(gama_map));

  IMG_DEBUG_LOG(
      "GAMA chunk read successfully\n,\t\tlength: %d\n\t\thas gamma: %s\n\t\tgamma level: %f\n",
      len, gama->has_gamma ? "true" : "false", (gama->has_gamma) ? gama->gamma / 100000.0f : 0.0f);

  return gama;
}

/**
 * @brief Reads the cHRM chunk of a PNG image.
 * @details This function reads the cHRM chunk of a PNG image and populates an
 * img_png_chrm_chunk structure with the chromaticity values for the image.
 * It verifies the CRC32 checksum of the chunk, checks the length, and returns
 * a pointer to the structure if the chunk is read successfully, otherwise
 * returns nullptr. The function logs the chromaticity values if successful.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @param chunk_map The vector to store the read chunk.
 * @return A pointer to the img_png_chrm_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */

static std::shared_ptr<img_png_chrm_chunk>
img_png_read_chrm_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc,
                        std::vector<std::shared_ptr<img_png_type_chunk_map>>& chunk_map) noexcept {
  std::unique_ptr<img_png_chrm_chunk> chrm = std::make_unique<img_png_chrm_chunk>();

  if (len != 32) {
    chrm->has_chrm = false;
    IMG_DEBUG_LOG("Invalid chrm chunk length: %x, expected: 32\n", len);
    return nullptr;
  } else {
    chrm->has_chrm = true;
    auto read_fixed_point = [&](int offset) {
      return ((data[offset] << 24) | (data[offset + 1] << 16) |
              (data[offset + 2] << 8) | data[offset + 3]) /
              100000.0f;
    };

    chrm->white_x = read_fixed_point(0);
    chrm->white_y = read_fixed_point(4);
    chrm->red_x = read_fixed_point(8);
    chrm->red_y = read_fixed_point(12);
    chrm->green_x = read_fixed_point(16);
    chrm->green_y = read_fixed_point(20);
    chrm->blue_x = read_fixed_point(24);
    chrm->blue_y = read_fixed_point(28);

    chrm->length = len;
    chrm->crc = crc;
    chrm->type = type;
  }

  std::shared_ptr<img_png_type_chunk_map> chrm_map =
      std::make_shared<img_png_type_chunk_map>(type, std::move(chrm));
  chunk_map.emplace_back(std::move(chrm_map));

  IMG_DEBUG_LOG(
      "CHRM chunk read successfully\n,\t\tlength: %d\n\t\thas chrm: %s\n\t\twhite point: (%f, "
      "%f)\n\t\tred point: (%f, %f)\n\t\tgreen point: (%f, %f)\n\t\tblue point: (%f, %f)\n",
      len, chrm->has_chrm ? "true" : "false", chrm->white_x, chrm->white_y,
      chrm->red_x, chrm->red_y, chrm->green_x, chrm->green_y, chrm->blue_x, chrm->blue_y);

  return chrm;
}

/**
 * @brief Reads a PNG file and extracts its chunks.
 * @details This function reads a PNG file and extracts its chunks. It verifies
 * the PNG signature, calculates the CRC for each chunk, and stores the chunks
 * in a vector. It also handles the IHDR, IDAT, and PLTE chunks specifically.
 * @param spec Pointer to the image specification structure.
 * @param data Pointer to the PNG file data.
 * @param size Size of the PNG file data in bytes.
 * @return True if the PNG file was read successfully, false otherwise.
 */
static bool img_png_read(const std::shared_ptr<img::image_specification>& spec, uint8_t* data,
                         size_t size) noexcept {
  bit_reader bit(data, size);                 // Create a bit bit for the data
  uint64_t signature = bit.read<uint64_t>();  // Read the PNG signature
  if (!img_png_verify_signature(signature)) {
    IMG_DEBUG_LOG("Invalid PNG signature: %llx, File:%s is not a PNG file\n", signature,
                  spec->filename.c_str());
    return false;
  }

  IMG_DEBUG_LOG("PNG Signature: %llx matched\n", signature);

  std::call_once(s_CRC_TBL_INIT_FLAG,
                 img_png_generate_crc_table);                       // Generate the CRC table
  std::vector<std::shared_ptr<img_png_type_chunk_map>> png_chunks;  // Vector to hold PNG chunks

  std::weak_ptr<img_png_ihdr_chunck> temp_ihdr_ptr;  // Pointer to IHDR chunk
  std::weak_ptr<img_png_gama_chunk> temp_gamma_ptr;  // Pointer to gama chunk
  std::weak_ptr<img_png_chrm_chunk> temp_chrm_ptr;   // Pointer to chrm chunk
  std::weak_ptr<img_png_iccp_chunk> temp_iccp_ptr;   // Pointer to iccp chunk
  std::weak_ptr<img_png_srgb_chunk> temp_srgb_ptr;   // Pointer to srgb chunk
  std::weak_ptr<img_png_plte_chunk> temp_plte_ptr;   // Pointer to plte chunk

  bool is_file_read_successful = false;

  while (true) {
    img_png_chunk chunk;
    chunk.length = bit.read<uint32_t>();
    chunk.type = bit.read<uint32_t>();

    if (chunk.length == 0) {
      IMG_DEBUG_LOG("Invalid chunk, length: %x, Possible corruption\n", chunk.length);
      break;
    }

    uint8_t* chunk_data = bit.read(chunk.length);
    if (chunk_data == nullptr) {
      IMG_DEBUG_LOG("Failed to read chunk data\n");
      break;
    }

    chunk.crc = bit.read<uint32_t>();  // Read the CRC value for the chunk
    uint32_t crc = img_png_verify_chunk_crc(chunk.type, chunk.length, chunk.crc, chunk_data);
    if (crc) {
      IMG_DEBUG_LOG("Invalid CRC for chunk type: %x, expected: %x\n", chunk.type, crc);
      delete[] chunk_data;  // Free the chunk data
      break;
    }

    IMG_DEBUG_LOG("=====================================================================================\n");
    IMG_DEBUG_LOG("Reading chunk, \n\t\ttype: %x\n\t\tlength: %x\n\t\tcrc: %x\n", chunk.type,
                  chunk.length, chunk.crc);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_IHDR_CHUNK) {
      temp_ihdr_ptr =
          img_png_read_ihdr_chunk(chunk.type, chunk.length, chunk_data, chunk.crc, png_chunks);
      IMG_DEBUG_LOG("Critical error, png chunck point location undefined -> %x",
                    s_IMG_PNG_IHDR_CHUNK);
      continue;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_IDAT_CHUNK) {
      std::weak_ptr<img_png_idat_chunk> idat =
          img_png_read_idat_chunk(chunk.type, chunk.length, chunk_data, chunk.crc, png_chunks);
      if (!idat.expired() && !temp_ihdr_ptr.expired()) {
        std::shared_ptr<img_png_idat_chunk> idat_ptr = idat.lock();
        idat_ptr->ihdr = temp_ihdr_ptr.lock();
      } else {
        IMG_DEBUG_LOG("Critical error, png chunck point location undefined -> %x and %x",
                      s_IMG_PNG_IHDR_CHUNK, s_IMG_PNG_IDAT_CHUNK);
        delete[] chunk_data;
      }

      continue;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_PLTE_CHUNK) {
      temp_plte_ptr =
          img_png_read_plte_chunck(chunk.type, chunk.length, chunk_data, chunk.crc, png_chunks);
      if (temp_plte_ptr.expired()) {
        IMG_DEBUG_LOG("Critical error, png chunck point location undefined -> %x",
                      s_IMG_PNG_PLTE_CHUNK);
        delete[] chunk_data;
      }

      continue;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_GAMA_CHUNK) {
      temp_gamma_ptr =
          img_png_read_gama_chunk(chunk.type, chunk.length, chunk_data, chunk.crc, png_chunks);
      if (!temp_gamma_ptr.expired()) {
        std::shared_ptr<img_png_gama_chunk> gama_ptr = temp_gamma_ptr.lock();
        if (gama_ptr->has_gamma) {
          gama_ptr->chrm = temp_chrm_ptr.lock();
          gama_ptr->iccp = temp_iccp_ptr.lock();
          gama_ptr->srgb = temp_srgb_ptr.lock();
        }
      } else {
        IMG_DEBUG_LOG("Critical error, png chunck point location undefined -> %x",
                      s_IMG_PNG_GAMA_CHUNK);
        delete[] chunk_data;
      }

      continue;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_CHRM_CHUNK) {
      temp_chrm_ptr = img_png_read_chrm_chunk(chunk.type, chunk.length, chunk_data, chunk.crc, png_chunks);
      if(!temp_chrm_ptr.expired()){
        std::shared_ptr<img_png_chrm_chunk> chrm_ptr = temp_chrm_ptr.lock();
        if(chrm_ptr->has_chrm){
          chrm_ptr->gama = temp_gamma_ptr.lock();
          chrm_ptr->iccp = temp_iccp_ptr.lock();
          chrm_ptr->srgb = temp_srgb_ptr.lock();
        }
      }
      else{
        IMG_DEBUG_LOG("Critical error, png chunck point location undefined -> %x",
                      s_IMG_PNG_CHRM_CHUNK);
        delete[] chunk_data;
      }

      continue;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_ICCP_CHUNK) {
      std::unique_ptr<img_png_iccp_chunk> iccp = std::make_unique<img_png_iccp_chunk>();
      temp_iccp_ptr = iccp.get();  // Store the pointer to iccp chunk

      size_t index = 0;
      while (index < chunk.length && chunk_data[index] != '\0') {
        iccp->profile_name += static_cast<char>(chunk_data[index++]);
      }
      ++index;  // skip null terminator

      if (index >= chunk.length) {
        IMG_DEBUG_LOG("Malformed ICCP chunk: missing compression method and data\n");
        break;
      }

      iccp->compression_method = chunk_data[index++];
      if (iccp->compression_method != 0) {
        IMG_DEBUG_LOG(
            "Unsupported ICCP compression method: %u\n (Possible reason, File is not compressed "
            "using standard DEFLATE algorithm)",
            iccp->compression_method);
        break;
      }

      // Remaining is compressed ICC profile
      iccp->compressed_profile.insert(iccp->compressed_profile.end(), chunk_data + index,
                                      chunk_data + chunk.length);

      iccp->length = chunk.length;
      iccp->crc = chunk.crc;
      iccp->type = chunk.type;

      iccp->gama = temp_gamma_ptr;  // Link iccp to gama chunk
      iccp->chrm = temp_chrm_ptr;   // Link iccp to chrm chunk
      iccp->srgb = temp_srgb_ptr;   // Link iccp to srgb chunk

      std::shared_ptr<img_png_type_chunk_map> iccp_map =
          std::make_shared<img_png_type_chunk_map>(chunk.type, std::move(iccp));
      png_chunks.emplace_back(std::move(iccp_map));

      IMG_DEBUG_LOG(
          "ICCP chunk read successfully\n,\t\tlength: %d\n\t\tcompression method: %d\n\t\tprofile "
          "name: "
          "%s\n",
          chunk.length, iccp->compression_method, iccp->profile_name.c_str());
    }

    if (chunk.type == s_IMG_PNG_SRGB_CHUNK) {
      std::unique_ptr<img_png_srgb_chunk> srgb = std::make_unique<img_png_srgb_chunk>();
      temp_srgb_ptr = srgb.get();  // Store the pointer to srgb chunk

      if (chunk.length != 1) {
        srgb->has_srgb = false;
        IMG_DEBUG_LOG("Corrupted SRGB chunk -> length: %x, expected: 1\n", chunk.length);
        break;

      } else {
        srgb->has_srgb = true;
        uint8_t rendering_intent = chunk_data[0];
        switch (rendering_intent) {
          case 0:
            srgb->rendering_intent = img_png_srgb_rendering_intent::perceptual;
            break;
          case 1:
            srgb->rendering_intent = img_png_srgb_rendering_intent::relative_colorimetic;
            break;
          case 2:
            srgb->rendering_intent = img_png_srgb_rendering_intent::saturation;
            break;
          case 3:
            srgb->rendering_intent = img_png_srgb_rendering_intent::absolute_colorimetric;
            break;
          default:
            srgb->rendering_intent = img_png_srgb_rendering_intent::perceptual;
            break;
        };

        srgb->length = chunk.length;
        srgb->crc = chunk.crc;
        srgb->type = chunk.type;

        srgb->gama = temp_gamma_ptr;  // Link srgb to gama chunk
        srgb->chrm = temp_chrm_ptr;   // Link srgb to chrm chunk
        srgb->iccp = temp_iccp_ptr;   // Link srgb to iccp chunk
      }

      std::shared_ptr<img_png_type_chunk_map> srgb_map =
          std::make_shared<img_png_type_chunk_map>(chunk.type, std::move(srgb));
      png_chunks.emplace_back(std::move(srgb_map));

      IMG_DEBUG_LOG("SRGB chunk read successfully\n,\t\tlength: %d\n\t\trendering intent: %d\n",
                    chunk.length, srgb->rendering_intent);
    }

    if (chunk.type == s_IMG_PNG_BKGD_CHUNK) {
      std::unique_ptr<img_png_bkgd_chunk> bkgd = std::make_unique<img_png_bkgd_chunk>();

      bkgd->length = chunk.length;
      bkgd->crc = chunk.crc;
      bkgd->type = chunk.type;

      bkgd->ihdr = temp_ihdr_ptr;  // Link bkgd to ihdr chunk
      bkgd->plte = temp_plte_ptr;  // Link bkgd to plte chunk

      if (!bkgd->ihdr) {
        IMG_DEBUG_LOG(
            "Something went wrong!!!, Parser reaches the BKGD chunk before IHDR chunk (Possible "
            "reason : IHDR chunk corrupted)\n");
        delete[] chunk_data;
        break;
      }

      uint8_t color_type = bkgd->ihdr->color_type;
      bkgd->color_type = color_type;

      switch (color_type) {
        case 0:  // Grayscale
          if (chunk.length != 2) {
            IMG_DEBUG_LOG("Invalid BKGD size for grayscale: %u\n", chunk.length);
            delete[] chunk_data;
            return false;
          }
          bkgd->gray = (chunk_data[0] << 8) | chunk_data[1];
          break;

        case 2:  // Truecolor
          if (chunk.length != 6) {
            IMG_DEBUG_LOG("Invalid BKGD size for truecolor: %u\n", chunk.length);
            delete[] chunk_data;
            return false;
          }
          bkgd->r = (chunk_data[0] << 8) | chunk_data[1];
          bkgd->g = (chunk_data[2] << 8) | chunk_data[3];
          bkgd->b = (chunk_data[4] << 8) | chunk_data[5];
          break;

        case 3:  // Indexed
          if (chunk.length != 1) {
            IMG_DEBUG_LOG("Invalid BKGD size for indexed-color: %u\n", chunk.length);
            delete[] chunk_data;
            return false;
          }
          bkgd->palette_index = chunk_data[0];
          break;

        default:
          IMG_DEBUG_LOG("Unsupported color type for BKGD: %u\n", color_type);
          delete[] chunk_data;
          return false;
      };

      std::shared_ptr<img_png_type_chunk_map> bkgd_map =
          std::make_shared<img_png_type_chunk_map>(chunk.type, std::move(bkgd));
      png_chunks.emplace_back(std::move(bkgd_map));

      IMG_DEBUG_LOG("BKGD chunk read successfully\n,\t\tlength: %d\n\t\tcolor type: %d\n",
                    chunk.length, bkgd->color_type);
    }

    if (chunk.type == s_IMG_PNG_IEND_CHUNK) {
      std::unique_ptr<img_png_iend_chunk> iend = std::make_unique<img_png_iend_chunk>();
      iend->is_iend = true;

      iend->length = chunk.length;
      iend->crc = chunk.crc;
      iend->type = chunk.type;

      std::shared_ptr<img_png_type_chunk_map> iend_map =
          std::make_shared<img_png_type_chunk_map>(chunk.type, std::move(iend));
      png_chunks.emplace_back(std::move(iend_map));
      is_file_read_successful = true;
      delete[] chunk_data;
      break;
    }
  }

  if (!is_file_read_successful) {
    return false;
  }
}

namespace img {
std::shared_ptr<image_specification> import(const std::filesystem::path& filepath,
                                            image_format format, bool flip) noexcept {
  if (format == image_format::png) {
    std::shared_ptr<image_specification> spec = std::make_shared<image_specification>();

    spec->filepath = filepath;
    spec->filename = filepath.filename().string();
    spec->xid = unique_xid();

    compressed_file_data file(filepath);
    if (img_png_read(spec, file.data, file.size)) {
      IMG_DEBUG_LOG("PNG file read successfully: %s\n", filepath.string().c_str());
      return spec;
    } else {
      IMG_DEBUG_LOG("Failed to read PNG file: %s\n", filepath.string().c_str());
      return nullptr;
    }
  }
}
}  // namespace img
