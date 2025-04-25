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
  uint32_t xid{0};                // Unique identifier
  uint32_t type{0};               // Type of the chunk
  img_chunk_ptr chunks{nullptr};  // Pointer to the chunk

  img_png_type_chunk_map() = default;
  img_png_type_chunk_map(uint32_t id, uint32_t t, img_chunk_ptr&& c) {
    xid = id;
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

/**
 * @struct img_png_iend_chunk
 * @brief This structure represents the IEND chunk of a PNG image.
 * It contains a boolean indicating whether the chunk is an IEND chunk or not.
 */
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
  std::vector<uint8_t> data;                // Pointer to the IDAT chunk data
  std::weak_ptr<img_png_ihdr_chunck> ihdr;  // Pointer to the IHDR chunk

  img_png_idat_chunk() = default;
  virtual ~img_png_idat_chunk() = default;
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
struct img_png_scal_chunk;  // Forward declaration of img_png_scal_chunk
struct img_png_pcal_chunk;  // Forward declaration of img_png_pcal_chunk

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

  std::weak_ptr<img_png_plte_chunk> plte;   // Pointer to the plte chunk (if present)
  std::weak_ptr<img_png_ihdr_chunck> ihdr;  // Pointer to the IHDR chunk

  img_png_bkgd_chunk() = default;
  virtual ~img_png_bkgd_chunk() = default;
};

/**
 * @struct img_png_phys_chunk
 * @brief This structure represents the pHYs chunk of a PNG image.
 * It contains the physical pixel dimensions and unit type for the image.
 */
struct img_png_phys_chunk : img_png_chunk {
  uint32_t x_pixels_per_unit{0};  ///< The number of pixels per unit along the x-axis.
  uint32_t y_pixels_per_unit{0};  ///< The number of pixels per unit along the y-axis.
  uint32_t unit_type{0};          ///< The unit type (0 for unknown, 1 for meter).

  std::weak_ptr<img_png_ihdr_chunck> ihdr;  // Pointer to the IHDR chunk
  std::weak_ptr<img_png_scal_chunk> scal;   // Pointer to the scal chunk
  std::weak_ptr<img_png_pcal_chunk> pcal;   // Pointer to the pcal chunk

  img_png_phys_chunk() = default;
  virtual ~img_png_phys_chunk() = default;
};

/**
 * @struct img_png_iTXt_chunk
 * @brief This structure represents the iTXt chunk of a PNG image.
 * It contains the keyword, whether the text is compressed, the language tag,
 * the translated keyword, and the text itself.
 */
struct img_png_itxt_chunk : img_png_chunk {
  std::string keyword;
  uint8_t compression_flag{0};
  uint8_t compression_method{0};
  std::string language_tag;
  std::string translated_keyword;
  std::string text;                      // Uncompressed
  std::vector<uint8_t> compressed_text;  // Store raw compressed for now
  bool is_compressed{false};

  img_png_itxt_chunk() = default;
  virtual ~img_png_itxt_chunk() = default;
};

/******************************************************************
 *                         PNG CHUNK READERS                      *
 * ****************************************************************/

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
 * @return A pointer to the img_png_ihdr_chunck structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_ihdr_chunck>
img_png_read_ihdr_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc) noexcept {
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

  IMG_DEBUG_LOG("IHDR CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("WIDTH = %d PIXELS | HEIGHT = %d PIXELS (TOTAL PIXELS = %d)\n", ihdr->width,
                ihdr->height, ihdr->width * ihdr->height);
  IMG_DEBUG_LOG("BIT DEPTH = %d\n", ihdr->bit_depth);
  IMG_DEBUG_LOG("COLOR TYPE = %d\n", ihdr->color_type);
  IMG_DEBUG_LOG("COMPRESSION METHOD = %d\n", ihdr->compression_method);
  IMG_DEBUG_LOG("FILTER METHOD = %d\n", ihdr->filter_method);
  IMG_DEBUG_LOG("INTERLACE METHOD = %d\n", ihdr->interlace_method);

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
 * @return A pointer to the img_png_idat_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_idat_chunk>
img_png_read_idat_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc) noexcept {
  std::shared_ptr<img_png_idat_chunk> idat = std::make_shared<img_png_idat_chunk>();
  idat->data.reserve(len);
  idat->data.assign(data, data + len);

  idat->length = len;
  idat->crc = crc;
  idat->type = type;

  IMG_DEBUG_LOG("IDAT CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("COMPRESSED DATA LENGTH = %d BYTES\n", idat->length);

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
 * @return A pointer to the img_png_plte_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_plte_chunk>
img_png_read_plte_chunck(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc) noexcept {
  std::shared_ptr<img_png_plte_chunk> plte = std::make_shared<img_png_plte_chunk>();
  plte->palette.reserve(len / 3);
  for (uint32_t i = 0; i < len; i += 3) {
    plte->palette.emplace_back(img_rgb(data[i], data[i + 1], data[i + 2]));
  }

  plte->length = len;
  plte->crc = crc;
  plte->type = type;

  IMG_DEBUG_LOG("PLTE CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("LENGTH = %d BYTES\n", plte->length);
  IMG_DEBUG_LOG("COLORS = %d\n", plte->palette.size());

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
 * @return A pointer to the img_png_trns_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_trns_chunk>
img_png_read_trns_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc) noexcept {
  std::shared_ptr<img_png_trns_chunk> tRNS = std::make_shared<img_png_trns_chunk>();
  tRNS->alpha.reserve(len);
  for (uint32_t i = 0; i < len; ++i) {
    tRNS->alpha.emplace_back(data[i]);
  }

  tRNS->length = len;
  tRNS->crc = crc;
  tRNS->type = type;

  IMG_DEBUG_LOG("TRNS CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("LENGTH = %d BYTES\n", tRNS->length);
  IMG_DEBUG_LOG("ALPHA = %d\n", tRNS->alpha.size());

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
 * @return A pointer to the img_png_gama_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_gama_chunk>
img_png_read_gama_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc) noexcept {
  std::shared_ptr<img_png_gama_chunk> gama = std::make_shared<img_png_gama_chunk>();

  if (len != 4 && len != 6) {
    gama->has_gamma = false;
    IMG_DEBUG_LOG("INVALID GAMA CHUNK LENGHT -> %d, EXPECTED 4 OR 6\n", len);
    return nullptr;
  } else {
    gama->has_gamma = true;
    gama->gamma = (data[0] << 24) | (data[1] << 16) | (data[2] << 8) | data[3];

    gama->length = len;
    gama->crc = crc;
    gama->type = type;
  }

  IMG_DEBUG_LOG("GAMA CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("LENGTH = %d BYTES\n", gama->length);
  IMG_DEBUG_LOG("GAMMA = %d\n", gama->gamma);

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
 * @return A pointer to the img_png_chrm_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */

static std::shared_ptr<img_png_chrm_chunk>
img_png_read_chrm_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc) noexcept {
  std::shared_ptr<img_png_chrm_chunk> chrm = std::make_shared<img_png_chrm_chunk>();

  if (len != 32) {
    chrm->has_chrm = false;
    IMG_DEBUG_LOG("INVALID CHRM CHUNK LENGHT -> %d, EXPECTED 32\n", len);
    return nullptr;
  } else {
    chrm->has_chrm = true;
    auto read_fixed_point = [&](int offset) {
      return ((data[offset] << 24) | (data[offset + 1] << 16) | (data[offset + 2] << 8) |
              data[offset + 3]) /
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

  IMG_DEBUG_LOG("CHRM CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("LENGTH = %d BYTES\n", chrm->length);
  IMG_DEBUG_LOG("WHITE POINTS = (%f, %f)\n", chrm->white_x, chrm->white_y);
  IMG_DEBUG_LOG("RED POINTS = (%f, %f)\n", chrm->red_x, chrm->red_y);
  IMG_DEBUG_LOG("GREEN POINTS = (%f, %f)\n", chrm->green_x, chrm->green_y);
  IMG_DEBUG_LOG("BLUE POINTS = (%f, %f)\n", chrm->blue_x, chrm->blue_y);

  return chrm;
}

/**
 * @brief Reads the iCCP chunk of a PNG image.
 * @details This function reads the iCCP chunk of a PNG image and populates an
 * img_png_iccp_chunk structure with the ICC profile information. It extracts
 * the profile name, verifies the compression method, and stores the compressed
 * ICC profile data. It verifies the CRC32 checksum of the chunk and returns a
 * pointer to the structure if the chunk is read successfully, otherwise returns
 * nullptr. The function logs the profile name, compression method, and any
 * errors encountered during processing.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @return A pointer to the img_png_iccp_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_iccp_chunk>
img_png_read_iccp_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc) noexcept {
  std::shared_ptr<img_png_iccp_chunk> iccp = std::make_shared<img_png_iccp_chunk>();

  size_t index = 0;
  while (index < len && data[index] != '\0') {
    iccp->profile_name += static_cast<char>(data[index++]);
  }
  ++index;  // skip null terminator

  if (index >= len) {
    IMG_DEBUG_LOG("MALFORMED ICCP CHUNK -> NO NULL TERMINATOR\n");
    return nullptr;
  }

  iccp->compression_method = data[index++];
  if (iccp->compression_method != 0) {
    IMG_DEBUG_LOG(
        "UNSUPPORTED ICCP CHUNK -> COMPRESSION METHOD %d NOT SUPPORTED!! (POSSIBLE REASON FILE IS "
        "NOT COMPRESSED USING STANDARD DEFLATE ALGORITHM)\n",
        iccp->compression_method);
    return nullptr;
  }

  // Remaining is compressed ICC profile
  iccp->compressed_profile.insert(iccp->compressed_profile.end(), data + index, data + len);

  iccp->length = len;
  iccp->crc = crc;
  iccp->type = type;

  IMG_DEBUG_LOG("ICCP CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("LENGTH = %d BYTES\n", iccp->length);
  IMG_DEBUG_LOG("PROFILE NAME = %s\n", iccp->profile_name.c_str());
  IMG_DEBUG_LOG("COMPRESSION METHOD = %d\n", iccp->compression_method);

  return iccp;
}

/**
 * @brief Reads the sRGB chunk of a PNG image.
 * @details This function reads the sRGB chunk of a PNG image and populates an
 * img_png_srgb_chunk structure with the rendering intent of the image. It
 * verifies the length of the chunk and the CRC32 checksum of the chunk. It
 * returns a pointer to the structure if the chunk is read successfully,
 * otherwise returns nullptr. The function logs the rendering intent and any
 * errors encountered during processing.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @return A pointer to the img_png_srgb_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_srgb_chunk>
img_png_read_srgb_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc) noexcept {
  std::shared_ptr<img_png_srgb_chunk> srgb = std::make_shared<img_png_srgb_chunk>();

  if (len != 1) {
    srgb->has_srgb = false;
    IMG_DEBUG_LOG("MALFORMED SRGB CHUNK MUST BE 1 BYTE LONG, READ VALUE = %d", len);
    return nullptr;

  } else {
    srgb->has_srgb = true;
    uint8_t rendering_intent = data[0];
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

    srgb->length = len;
    srgb->crc = crc;
    srgb->type = type;
  }

  IMG_DEBUG_LOG("SRGB CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("LENGTH = %d BYTES\n", srgb->length);
  IMG_DEBUG_LOG("RENDERING INTENT = %d\n", srgb->rendering_intent);

  return srgb;
}

/**
 * @brief Reads the BKGD chunk of a PNG image.
 * @details This function reads the BKGD chunk of a PNG image and populates an
 * img_png_bkgd_chunk structure with the background color of the image. It
 * verifies the length of the chunk and the CRC32 checksum of the chunk. It
 * returns a pointer to the structure if the chunk is read successfully,
 * otherwise returns nullptr. The function logs the color type and any errors
 * encountered during processing.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @param ihdr The pointer to the IHDR chunk of the image.
 * @param plte The pointer to the PLTE chunk of the image.
 * @return A pointer to the img_png_bkgd_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_bkgd_chunk>
img_png_read_bkgd_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc,
                        const std::shared_ptr<img_png_ihdr_chunck>& ihdr,
                        const std::shared_ptr<img_png_plte_chunk>& plte) noexcept {
  std::shared_ptr<img_png_bkgd_chunk> bkgd = std::make_shared<img_png_bkgd_chunk>();

  bkgd->length = len;
  bkgd->crc = crc;
  bkgd->type = type;

  bkgd->ihdr = ihdr;
  bkgd->plte = plte;

  if (bkgd->ihdr.expired()) {
    IMG_DEBUG_LOG("SOMETHING WENT WRONG!!, CAN NOT FIND IHDR CHUNK\n");
    return nullptr;
  }

  uint8_t color_type = ihdr->color_type;
  bkgd->color_type = color_type;

  switch (color_type) {
    case 0:  // Grayscale
      if (len != 2) {
        IMG_DEBUG_LOG("INVALID BKGD SIZE FOR GRAYSCALE: %u\n", len);
        return nullptr;
      }
      bkgd->gray = (data[0] << 8) | data[1];
      break;

    case 2:  // Truecolor
      if (len != 6) {
        IMG_DEBUG_LOG("INVALID BKGD SIZE FOR TRUECOLOR: %u\n", len);
        return nullptr;
      }
      bkgd->r = (data[0] << 8) | data[1];
      bkgd->g = (data[2] << 8) | data[3];
      bkgd->b = (data[4] << 8) | data[5];
      break;

    case 3:  // Indexed
      if (len != 1) {
        IMG_DEBUG_LOG("INVALID BKGD SIZE FOR INDEXED-COLOR: %u\n", len);
        return nullptr;
      }
      bkgd->palette_index = data[0];
      break;

    default:
      IMG_DEBUG_LOG("UNSUPPORTED BKGD COLOR TYPE: %u\n", color_type);
      return nullptr;
  };

  IMG_DEBUG_LOG("BKGD CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("LENGTH = %d BYTES\n", bkgd->length);
  IMG_DEBUG_LOG("COLOR TYPE = %d\n", bkgd->color_type);
  IMG_DEBUG_LOG("GRAY = %d\n", bkgd->gray);
  IMG_DEBUG_LOG("R = %d\n", bkgd->r);
  IMG_DEBUG_LOG("G = %d\n", bkgd->g);
  IMG_DEBUG_LOG("B = %d\n", bkgd->b);
  IMG_DEBUG_LOG("PALETTE INDEX = %d\n", bkgd->palette_index);

  return bkgd;
}

/**
 * @brief Reads the pHYs chunk of a PNG image.
 * @details This function reads the pHYs chunk of a PNG image and populates an
 * img_png_phys_chunk structure with the physical pixel dimensions and unit type.
 * It verifies the length of the chunk and the CRC32 checksum, and returns a pointer
 * to the structure if the chunk is read successfully, otherwise returns nullptr.
 * The function logs the physical dimensions and any errors encountered during processing.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @param ihdr A pointer to the IHDR chunk, if present.
 * @param pcal A pointer to the pCAL chunk, if present.
 * @param scal A pointer to the sCAL chunk, if present.
 * @return A pointer to the img_png_phys_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_phys_chunk>
img_png_read_phys_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc,
                        const std::shared_ptr<img_png_ihdr_chunck>& ihdr,
                        const std::shared_ptr<img_png_pcal_chunk>& pcal,
                        const std::shared_ptr<img_png_scal_chunk>& scal) noexcept {
  std::shared_ptr<img_png_phys_chunk> phys = std::make_shared<img_png_phys_chunk>();
  if (len != 9) {
    IMG_DEBUG_LOG("INVALID LENGHT FOR PHYS CHUNK: %u (EXPECTED 9)\n", len);
    return nullptr;
  } else {
    phys->x_pixels_per_unit = (data[0] << 24) | (data[1] << 16) | (data[2] << 8) | data[3];
    phys->y_pixels_per_unit = (data[4] << 24) | (data[5] << 16) | (data[6] << 8) | data[7];
    phys->unit_type = data[8];

    phys->length = len;
    phys->type = type;
    phys->crc = crc;

    phys->ihdr = ihdr;
    phys->pcal = pcal;
    phys->scal = scal;
  }

  IMG_DEBUG_LOG("PHYS CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("LENGTH = %d BYTES\n", phys->length);
  IMG_DEBUG_LOG("X PIXELS PER UNIT = %d\n", phys->x_pixels_per_unit);
  IMG_DEBUG_LOG("Y PIXELS PER UNIT = %d\n", phys->y_pixels_per_unit);
  IMG_DEBUG_LOG("UNIT TYPE = %d\n", phys->unit_type);

  return phys;
}

/**
 * @brief Reads the iTXt chunk of a PNG image.
 * @details This function reads the iTXt chunk of a PNG image and populates an
 * img_png_itxt_chunk structure with the keyword, compression flag, compression
 * method, language tag, translated keyword, and text. It verifies the length of
 * the chunk and the CRC32 checksum, and returns a pointer to the structure if
 * the chunk is read successfully, otherwise returns nullptr. The function logs
 * the keyword, compression flag, compression method, language tag, translated
 * keyword, and text and any errors encountered during processing.
 * @param type The type of the PNG chunk.
 * @param len The length of the chunk data in bytes.
 * @param data Pointer to the chunk data.
 * @param crc The CRC value to verify against.
 * @return A pointer to the img_png_itxt_chunk structure if the chunk is
 * read successfully, otherwise returns nullptr.
 */
static std::shared_ptr<img_png_itxt_chunk>
img_png_read_itxt_chunk(uint32_t type, uint32_t len, uint8_t* data, uint32_t crc) noexcept {
  std::shared_ptr<img_png_itxt_chunk> itxt = std::make_shared<img_png_itxt_chunk>();
  itxt->length = len;
  itxt->type = type;
  itxt->crc = crc;

  size_t offset = 0;
  while (offset < len && data[offset] != '\0') {
    itxt->keyword += static_cast<char>(data[offset++]);
  }
  offset++;

  if (offset >= len)
    return itxt;
  itxt->compression_flag = data[offset++];
  itxt->compression_method = data[offset++];

  while (offset < len && data[offset] != '\0') {
    itxt->language_tag += static_cast<char>(data[offset++]);
  }
  offset++;

  while (offset < len && data[offset] != '\0') {
    itxt->translated_keyword += static_cast<char>(data[offset++]);
  }
  offset++;

  size_t remaining = len - offset;
  if (itxt->compression_flag == 1 && itxt->compression_method == 0) {
    // It's compressed using zlib
    itxt->is_compressed = true;
    itxt->compressed_text.assign(data + offset, data + offset + remaining);
  } else {
    // It's just plain UTF-8 text
    itxt->is_compressed = false;
    itxt->text.assign(reinterpret_cast<const char*>(data + offset), remaining);
  }

  IMG_DEBUG_LOG("ITXT CHUNCK ------------------------------------------------- \n");
  IMG_DEBUG_LOG("LENGTH = %d BYTES\n", itxt->length);
  IMG_DEBUG_LOG("KEYWORD = %s\n", itxt->keyword.c_str());
  IMG_DEBUG_LOG("COMPRESSION FLAG = %d\n", itxt->compression_flag);
  IMG_DEBUG_LOG("COMPRESSION METHOD = %d\n", itxt->compression_method);
  IMG_DEBUG_LOG("LANGUAGE TAG = %s\n", itxt->language_tag.c_str());
  IMG_DEBUG_LOG("TRANSLATED KEYWORD = %s\n", itxt->translated_keyword.c_str());
  IMG_DEBUG_LOG("TEXT = %s\n", itxt->text.c_str());

  return itxt;
}

/**
 * @brief Map of PNG chunks.
 * @details This map stores the PNG chunks, with the chunk type as the key and
 * a pointer to the chunk data as the value. The map is populated by the
 * img_png_read_chunks() function.
 */
static std::unordered_map<uint32_t, std::shared_ptr<img_png_type_chunk_map>> s_PNG_CHUNKS_MAP;
#define DELETE_AND_CONTINUE(arry_ptr) \
  if (arry_ptr != nullptr) {          \
    delete[] arry_ptr;                \
    continue;                         \
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
static bool img_png_read(const std::shared_ptr<img::specification>& spec, img::mode mode,
                         uint8_t* data, size_t size) noexcept {
  bit_reader bit(data, size);                 // Create a bit bit for the data
  uint64_t signature = bit.read<uint64_t>();  // Read the PNG signature
  if (!img_png_verify_signature(signature)) {
    IMG_DEBUG_LOG(
        "PNG SIGNATURE: %llx DOES NOT MATCHED, EXPECTED: %llx\nFILE : %s IS NOT A PNG FILE",
        signature, s_IMG_PNG_SIGNATURE, spec->filename.c_str());
    return false;
  }

  IMG_DEBUG_LOG("PNG SIGNATURE: %llx MATCHED, FILE: %s IS A PNG FILE\n", signature,
                spec->filename.c_str());

  std::call_once(s_CRC_TBL_INIT_FLAG,
                 img_png_generate_crc_table);  // Generate the CRC table

  std::weak_ptr<img_png_ihdr_chunck> temp_ihdr_ptr;  // Pointer to IHDR chunk
  std::weak_ptr<img_png_gama_chunk> temp_gamma_ptr;  // Pointer to gama chunk
  std::weak_ptr<img_png_chrm_chunk> temp_chrm_ptr;   // Pointer to chrm chunk
  std::weak_ptr<img_png_iccp_chunk> temp_iccp_ptr;   // Pointer to iccp chunk
  std::weak_ptr<img_png_srgb_chunk> temp_srgb_ptr;   // Pointer to srgb chunk
  std::weak_ptr<img_png_plte_chunk> temp_plte_ptr;   // Pointer to plte chunk
  std::weak_ptr<img_png_pcal_chunk> temp_pcal_ptr;   // Pointer to pcal chunk
  std::weak_ptr<img_png_scal_chunk> temp_scal_ptr;   // Pointer to scal chunk
  std::weak_ptr<img_png_phys_chunk> temp_phys_ptr;   // Pointer to phys chunk

  bool is_file_read_successful = false;

  while (true) {
    img_png_chunk chunk;
    chunk.length = bit.read<uint32_t>();
    chunk.type = bit.read<uint32_t>();

    if (chunk.length == 0) {
      IMG_DEBUG_LOG(
          "INVALID CHUNK LENGTH (%d) FOR CHUNK: %x. (POSSIBLE REASON: FILE IS CORRUPTED)\n",
          chunk.length, chunk.type);
      break;
    }

    uint8_t* chunk_data = bit.read(chunk.length);
    if (chunk_data == nullptr) {
      IMG_DEBUG_LOG("FAILED TO READ: %x CHUNK DATA\n", chunk.type);
      break;
    }

    chunk.crc = bit.read<uint32_t>();  // Read the CRC value for the chunk
    uint32_t crc = img_png_verify_chunk_crc(chunk.type, chunk.length, chunk.crc, chunk_data);
    if (crc) {
      IMG_DEBUG_LOG("INVALID CRC VALUE FOR CHUNK: %x, EXPECTED: %x\n", chunk.type, chunk.crc);
      delete[] chunk_data;  // Free the chunk data
      break;
    }

    IMG_DEBUG_LOG(
        "=====================================================================================\n");
    IMG_DEBUG_LOG("CHUNK TYPE: %x\n", chunk.type);
    IMG_DEBUG_LOG("CHUNK LENGTH: %d\n", chunk.length);
    IMG_DEBUG_LOG("CHUNK CRC: %x\n", chunk.crc);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_IHDR_CHUNK) {
      temp_ihdr_ptr = img_png_read_ihdr_chunk(chunk.type, chunk.length, chunk_data, chunk.crc);
      if (std::shared_ptr<img_png_ihdr_chunck> ihdr_ptr = temp_ihdr_ptr.lock()) {
        std::shared_ptr<img_png_type_chunk_map> ihdr_map = std::make_shared<img_png_type_chunk_map>(
            spec->xid, s_IMG_PNG_IHDR_CHUNK, std::move(ihdr_ptr));
        s_PNG_CHUNKS_MAP[spec->xid] = ihdr_map;
      } else {
        IMG_DEBUG_LOG("CRITICAL CHUNCK : %x WAS NOT PROPERLY READ\n", s_IMG_PNG_IHDR_CHUNK);
        break;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_IDAT_CHUNK) {
      std::weak_ptr<img_png_idat_chunk> temp_idat_ptr =
          img_png_read_idat_chunk(chunk.type, chunk.length, chunk_data, chunk.crc);
      if (std::shared_ptr<img_png_idat_chunk> idat_ptr = temp_idat_ptr.lock()) {
        std::shared_ptr<img_png_type_chunk_map> idat_map = std::make_shared<img_png_type_chunk_map>(
            spec->xid, s_IMG_PNG_IDAT_CHUNK, std::move(idat_ptr));
        s_PNG_CHUNKS_MAP[spec->xid] = idat_map;
      } else {
        IMG_DEBUG_LOG("CRITICAL CHUNCK WAS NOT PROPERLY READ\n");
        break;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_PLTE_CHUNK &&
        (mode == img::mode::nessessary || mode == img::mode::fullbreakdown)) {
      temp_plte_ptr = img_png_read_plte_chunck(chunk.type, chunk.length, chunk_data, chunk.crc);
      if (std::shared_ptr<img_png_plte_chunk> plte_ptr = temp_plte_ptr.lock()) {
        std::shared_ptr<img_png_type_chunk_map> plte_map = std::make_shared<img_png_type_chunk_map>(
            spec->xid, s_IMG_PNG_PLTE_CHUNK, std::move(plte_ptr));
        s_PNG_CHUNKS_MAP[spec->xid] = plte_map;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_TRNS_CHUNK &&
        (mode == img::mode::nessessary || mode == img::mode::fullbreakdown)) {
      std::weak_ptr<img_png_trns_chunk> temp_trns_ptr =
          img_png_read_trns_chunk(chunk.type, chunk.length, chunk_data, chunk.crc);
      if (std::shared_ptr<img_png_trns_chunk> trns_ptr = temp_trns_ptr.lock()) {
        std::shared_ptr<img_png_type_chunk_map> trns_map = std::make_shared<img_png_type_chunk_map>(
            spec->xid, s_IMG_PNG_TRNS_CHUNK, std::move(trns_ptr));
        s_PNG_CHUNKS_MAP[spec->xid] = trns_map;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_GAMA_CHUNK &&
        (mode == img::mode::nessessary || mode == img::mode::fullbreakdown)) {
      temp_gamma_ptr = img_png_read_gama_chunk(chunk.type, chunk.length, chunk_data, chunk.crc);
      if (std::shared_ptr<img_png_gama_chunk> gama_ptr = temp_gamma_ptr.lock()) {
        if (gama_ptr->has_gamma) {
          gama_ptr->chrm = temp_chrm_ptr.lock();
          gama_ptr->iccp = temp_iccp_ptr.lock();
          gama_ptr->srgb = temp_srgb_ptr.lock();
        }
        std::shared_ptr<img_png_type_chunk_map> gama_map = std::make_shared<img_png_type_chunk_map>(
            spec->xid, s_IMG_PNG_GAMA_CHUNK, std::move(gama_ptr));
        s_PNG_CHUNKS_MAP[spec->xid] = gama_map;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_CHRM_CHUNK &&
        (mode == img::mode::nessessary || mode == img::mode::fullbreakdown)) {
      temp_chrm_ptr = img_png_read_chrm_chunk(chunk.type, chunk.length, chunk_data, chunk.crc);
      if (!temp_chrm_ptr.expired()) {
        std::shared_ptr<img_png_chrm_chunk> chrm_ptr = temp_chrm_ptr.lock();
        if (chrm_ptr->has_chrm) {
          chrm_ptr->gama = temp_gamma_ptr.lock();
          chrm_ptr->iccp = temp_iccp_ptr.lock();
          chrm_ptr->srgb = temp_srgb_ptr.lock();
        }

        std::shared_ptr<img_png_type_chunk_map> chrm_map = std::make_shared<img_png_type_chunk_map>(
            spec->xid, s_IMG_PNG_CHRM_CHUNK, std::move(chrm_ptr));
        s_PNG_CHUNKS_MAP[spec->xid] = chrm_map;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_ICCP_CHUNK &&
        (mode == img::mode::nessessary || mode == img::mode::fullbreakdown)) {
      temp_iccp_ptr = img_png_read_iccp_chunk(chunk.type, chunk.length, chunk_data, chunk.crc);
      if (std::shared_ptr<img_png_iccp_chunk> iccp_ptr = temp_iccp_ptr.lock()) {
        if (iccp_ptr->has_iccp) {
          iccp_ptr->gama = temp_gamma_ptr.lock();
          iccp_ptr->chrm = temp_chrm_ptr.lock();
          iccp_ptr->srgb = temp_srgb_ptr.lock();
        }

        std::shared_ptr<img_png_type_chunk_map> iccp_map = std::make_shared<img_png_type_chunk_map>(
            spec->xid, s_IMG_PNG_ICCP_CHUNK, std::move(iccp_ptr));
        s_PNG_CHUNKS_MAP[spec->xid] = iccp_map;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_SRGB_CHUNK &&
        (mode == img::mode::nessessary || mode == img::mode::fullbreakdown)) {
      temp_srgb_ptr = img_png_read_srgb_chunk(chunk.type, chunk.length, chunk_data, chunk.crc);
      if (std::shared_ptr<img_png_srgb_chunk> srgb_ptr = temp_srgb_ptr.lock()) {
        if (srgb_ptr->has_srgb) {
          srgb_ptr->gama = temp_gamma_ptr.lock();
          srgb_ptr->chrm = temp_chrm_ptr.lock();
          srgb_ptr->iccp = temp_iccp_ptr.lock();
        }

        std::shared_ptr<img_png_type_chunk_map> srgb_map = std::make_shared<img_png_type_chunk_map>(
            spec->xid, s_IMG_PNG_SRGB_CHUNK, std::move(srgb_ptr));
        s_PNG_CHUNKS_MAP[spec->xid] = srgb_map;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_BKGD_CHUNK &&
        (mode == img::mode::nessessary || mode == img::mode::fullbreakdown)) {
      if (!temp_ihdr_ptr.expired() && !temp_plte_ptr.expired()) {
        std::weak_ptr<img_png_bkgd_chunk> temp_bkgd_ptr =
            img_png_read_bkgd_chunk(chunk.type, chunk.length, chunk_data, chunk.crc,
                                    temp_ihdr_ptr.lock(), temp_plte_ptr.lock());

        if ( std::shared_ptr<img_png_bkgd_chunk> bkgd_ptr = temp_bkgd_ptr.lock()) {
          std::shared_ptr<img_png_type_chunk_map> bkgd_map =
              std::make_shared<img_png_type_chunk_map>(spec->xid, s_IMG_PNG_BKGD_CHUNK,
                                                       std::move(bkgd_ptr));
          s_PNG_CHUNKS_MAP[spec->xid] = bkgd_map;
        }

      } else {
        IMG_DEBUG_LOG("CRITICAL ERROR, CHUNCK POINT LOCATION UNDEFINED -> %x and %x",
                      s_IMG_PNG_IHDR_CHUNK, s_IMG_PNG_PLTE_CHUNK);
        break;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_PHYS_CHUNK && mode == img::mode::fullbreakdown) {
      temp_phys_ptr =
          img_png_read_phys_chunk(chunk.type, chunk.length, chunk_data, chunk.crc,
                                  temp_ihdr_ptr.lock(), temp_pcal_ptr.lock(), temp_scal_ptr.lock());
      if ( std::shared_ptr<img_png_phys_chunk> phys_ptr = temp_phys_ptr.lock()) {
        std::shared_ptr<img_png_type_chunk_map> phys_map = std::make_shared<img_png_type_chunk_map>(
            spec->xid, s_IMG_PNG_PHYS_CHUNK, std::move(phys_ptr));
        s_PNG_CHUNKS_MAP[spec->xid] = phys_map;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_ITXT_CHUNK && mode == img::mode::fullbreakdown) {
      std::weak_ptr<img_png_itxt_chunk> temp_itxt_ptr =
          img_png_read_itxt_chunk(chunk.type, chunk.length, chunk_data, chunk.crc);
      if (     std::shared_ptr<img_png_itxt_chunk> itxt_ptr = temp_itxt_ptr.lock()) {
        std::shared_ptr<img_png_type_chunk_map> itxt_map = std::make_shared<img_png_type_chunk_map>(
            spec->xid, s_IMG_PNG_ITXT_CHUNK, std::move(itxt_ptr));
        s_PNG_CHUNKS_MAP[spec->xid] = itxt_map;
      }

      DELETE_AND_CONTINUE(chunk_data);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (chunk.type == s_IMG_PNG_IEND_CHUNK) {
      std::shared_ptr<img_png_iend_chunk> iend = std::make_shared<img_png_iend_chunk>();
      iend->is_iend = true;

      iend->length = chunk.length;
      iend->crc = chunk.crc;
      iend->type = chunk.type;

      std::shared_ptr<img_png_type_chunk_map> iend_map = std::make_shared<img_png_type_chunk_map>(
          spec->xid, s_IMG_PNG_IEND_CHUNK, std::move(iend));
      s_PNG_CHUNKS_MAP[spec->xid] = iend_map;
      is_file_read_successful = true;
      delete[] chunk_data;
      break;
    }
  }
}

namespace img {
std::shared_ptr<specification> import(const std::filesystem::path& filepath, format format,
                                      mode mode, bool flip) noexcept {
  if (format == format::png) {
    std::shared_ptr<specification> spec = std::make_shared<specification>();

    spec->filepath = filepath;
    spec->filename = filepath.filename().string();
    spec->xid = unique_xid();

    compressed_file_data file(filepath);
    if (img_png_read(spec, mode, file.data, file.size)) {
      IMG_DEBUG_LOG("PNG file read successfully: %s\n", filepath.string().c_str());
      return spec;
    } else {
      IMG_DEBUG_LOG("Failed to read PNG file: %s\n", filepath.string().c_str());
      return nullptr;
    }
  }
}
}  // namespace img
