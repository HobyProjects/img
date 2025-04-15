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
#include "img_pch.hpp"

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
  explicit bit_reader(uint8_t* data, size_t size)
      : m_data(data), m_size(size), m_bit_offset(0) {}
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
  size_t size{0};  ///< Size of the compressed file data in bytes

  explicit compressed_file_data(const std::filesystem::path& file) {
    filepath = file;
    filename = file.filename().string();
    std::ifstream file_stream(filepath, std::ios::binary | std::ios::ate);
    if (!file_stream.is_open()) {
      IMG_DEBUG_LOG("Failed to open file: %s\n", filepath.string().c_str());
      throw std::runtime_error("Failed to open file: " + filepath.string());
    }
    size = file_stream.tellg();  // Get the size of the file
    file_stream.seekg(
        0, std::ios::beg);     // Move the cursor to the beginning of the file
    data = new uint8_t[size];  // Allocate memory for the file data
    file_stream.read(reinterpret_cast<char*>(data),
                     size);  // Read the file data into the buffer
    file_stream.close();     // Close the file stream
  };

  compressed_file_data(const compressed_file_data&) =
      delete;  // Delete copy constructor
  compressed_file_data& operator=(const compressed_file_data&) =
      delete;  // Delete copy assignment operator
  compressed_file_data(compressed_file_data&&) =
      default;  // Default move constructor
  compressed_file_data& operator=(compressed_file_data&&) =
      default;  // Default move assignment operator

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
  img_rgb(uint8_t red, uint8_t green, uint8_t blue)
      : r(red), g(green), b(blue) {}
  ~img_rgb() = default;  // Default destructor
};

/**
 * @struct img_rgba
 * @brief This structure represents an RGBA color.
 * It contains the red, green, blue, and alpha channel values.
 */
struct img_rgba {
  uint8_t r{0};  // Red channel value
  uint8_t g{0};  // Green channel value
  uint8_t b{0};  // Blue channel value
  uint8_t a{0};  // Alpha channel value

  img_rgba() = default;  // Default constructor
  img_rgba(uint8_t red, uint8_t green, uint8_t blue, uint8_t alpha)
      : r(red), g(green), b(blue), a(alpha) {}
  ~img_rgba() = default;  // Default destructor
};

/******************************************************************
 *                         PNG IMAGE FILES                        *
 * ****************************************************************/
static constexpr uint64_t s_IMG_PNG_SIGNATURE =
    0x89504E470D0A1A0A;  // PNG file signature
static constexpr uint32_t s_IMG_PNG_IHDR_CHUNK = 0x52414449;  // IHDR chunk type
static constexpr uint32_t s_IMG_PNG_IEND_CHUNK = 0x444E4549;  // IEND chunk type
static constexpr uint32_t s_IMG_PNG_IDAT_CHUNK = 0x54455849;  // IDAT chunk type
static constexpr uint32_t s_IMG_PNG_PLTE_CHUNK = 0x54494C50;  // PLTE chunk type
static constexpr uint32_t s_IMG_PNG_tRNS_CHUNK = 0x74524E53;  // tRNS chunk type
static constexpr uint32_t s_IMG_PNG_gAMA_CHUNK = 0x67414D41;  // gAMA chunk type
static constexpr uint32_t s_IMG_PNG_cHRM_CHUNK = 0x6348524D;  // cHRM chunk type
static constexpr uint32_t s_IMG_PNG_sRGB_CHUNK = 0x73524742;  // sRGB chunk type
static constexpr uint32_t s_IMG_PNG_bKGD_CHUNK = 0x624B4744;  // bKGD chunk type
static constexpr uint32_t s_IMG_PNG_pHYs_CHUNK = 0x70485973;  // pHYs chunk type
static constexpr uint32_t s_IMG_PNG_iTXt_CHUNK = 0x69545874;  // iTXt chunk type
static constexpr uint32_t s_IMG_PNG_tEXt_CHUNK = 0x74455874;  // tEXt chunk type
static constexpr uint32_t s_IMG_PNG_zTXt_CHUNK = 0x7A545874;  // zTXt chunk type
static constexpr uint32_t s_IMG_PNG_sBIT_CHUNK = 0x73424954;  // sBIT chunk type
static constexpr uint32_t s_IMG_PNG_iCCP_CHUNK = 0x69434350;  // iCCP chunk type
static constexpr uint32_t s_IMG_PNG_fRAc_CHUNK = 0x66726163;  // fRAc chunk type
static constexpr uint32_t s_IMG_PNG_hIST_CHUNK = 0x68495354;  // hIST chunk type
static constexpr uint32_t s_IMG_PNG_tIME_CHUNK = 0x74494D45;  // tIME chunk type
static constexpr uint32_t s_IMG_PNG_pCAL_CHUNK = 0x7043414C;  // pCAL chunk type
static constexpr uint32_t s_IMG_PNG_sCAL_CHUNK = 0x7343414C;  // sCAL chunk type

static constexpr uint32_t s_CRC32_POLYNOMIAL = 0xEDB88320;  // CRC32 polynomial
static constexpr uint32_t s_CRC_TABLE_SIZE = 256;  // Size of the CRC table

static uint32_t s_crc_table[s_CRC_TABLE_SIZE];  // CRC table
static std::once_flag
    s_crc_table_init_flag;  // Flag to ensure CRC table is initialized only once

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
    s_crc_table[i] = crc;
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
    crc = (crc >> 8) ^
          s_crc_table[(crc ^ byte) & 0xFF];  // Update CRC value using the table
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

using img_chunk_ptr = std::unique_ptr<img_png_chunk>;  // Pointer to a PNG chunk

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
  img_png_ihdr_chunck* ihdr_ptr{nullptr};  // Pointer to the IHDR chunk

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
 * @struct img_png_tRNS_chunk
 * @brief This structure represents the tRNS chunk of a PNG image.
 * It contains a vector of alpha values and the length of the tRNS chunk data.
 */
struct img_png_tRNS_chunk : img_png_chunk {
  std::vector<uint8_t> alpha;  // Vector to hold the alpha values

  img_png_tRNS_chunk() = default;
  virtual ~img_png_tRNS_chunk() = default;
};

struct img_png_gAMA_chunk;  // Forward declaration of img_png_gAMA_chunk
struct img_png_cHRM_chunk;  // Forward declaration of img_png_cHRM_chunk
struct img_png_iCCP_chunk;  // Forward declaration of img_png_iCCP_chunk
struct img_png_sRGB_chunk;  // Forward declaration of img_png_sRGB_chunk

/**
 * @struct img_png_gAMA_chunk
 * @brief This structure represents the gAMA chunk of a PNG image.
 * It contains the gamma value for the image.
 */
struct img_png_gAMA_chunk : img_png_chunk {
  uint32_t gamma{0};      // Gamma value for the image
  bool has_gamma{false};  // Flag to indicate if gAMA chunk is present

  img_png_cHRM_chunk* cHRM{nullptr};  // Pointer to the cHRM chunk (if present)
  img_png_iCCP_chunk* iCCP{nullptr};  // Pointer to the iCCP chunk (if present)
  img_png_sRGB_chunk* sRGB{nullptr};  // Pointer to the sRGB chunk (if present)

  img_png_gAMA_chunk() = default;
  virtual ~img_png_gAMA_chunk() = default;
};

/**
 * @struct img_png_cHRM_chunk
 * @brief This structure represents the cHRM chunk of a PNG image.
 * It contains the chromaticity values for the image.
 */
struct img_png_cHRM_chunk : img_png_chunk {
  float white_x{0.0f}, white_y{0.0f};
  float red_x{0.0f}, red_y{0.0f};
  float green_x{0.0f}, green_y{0.0f};
  float blue_x{0.0f}, blue_y{0.0f};
  bool has_chrm{false};

  img_png_gAMA_chunk* gAMA{nullptr};  // Pointer to the gAMA chunk (if present)
  img_png_iCCP_chunk* iCCP{nullptr};  // Pointer to the iCCP chunk (if present)
  img_png_sRGB_chunk* sRGB{nullptr};  // Pointer to the sRGB chunk (if present)

  img_png_cHRM_chunk() = default;
  virtual ~img_png_cHRM_chunk() = default;
};

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
static bool img_png_read(const std::shared_ptr<img::image_specification>& spec,
                         uint8_t* data, size_t size) noexcept {
  bit_reader bit(data, size);                 // Create a bit bit for the data
  uint64_t signature = bit.read<uint64_t>();  // Read the PNG signature
  if (signature != s_IMG_PNG_SIGNATURE) {
    IMG_DEBUG_LOG("Invalid PNG signature: %llx, File:%s is not a PNG file\n",
                  signature, spec->filename.c_str());
    return false;
  }

  std::call_once(s_crc_table_init_flag,
                 img_png_generate_crc_table);  // Generate the CRC table
  std::vector<std::shared_ptr<img_png_type_chunk_map>>
      png_chunks;  // Vector to hold PNG chunks

  img_png_ihdr_chunck* temp_ihdr_ptr = nullptr;  // Pointer to IHDR chunk
  img_png_gAMA_chunk* temp_gamma_ptr = nullptr;  // Pointer to gAMA chunk
  img_png_cHRM_chunk* temp_chrm_ptr = nullptr;   // Pointer to cHRM chunk
  img_png_iCCP_chunk* temp_iccp_ptr = nullptr;   // Pointer to iCCP chunk
  img_png_sRGB_chunk* temp_srgb_ptr = nullptr;   // Pointer to sRGB chunk

  while (true) {
    img_png_chunk chunk;
    chunk.length = bit.read<uint32_t>();
    chunk.type = bit.read<uint32_t>();

    if (chunk.length == 0) {
      IMG_DEBUG_LOG("Invalid chunk, length: %x, Possible corruption\n",
                    chunk.length);
      return false;
    }

    uint8_t* chunk_data = bit.read(chunk.length);
    if (chunk_data == nullptr) {
      IMG_DEBUG_LOG("Failed to read chunk data\n");
      return false;
    }

    chunk.crc = bit.read<uint32_t>();  // Read the CRC value for the chunk
    size_t crc_input_size = 4 + chunk.length;
    uint8_t* crc_input = new uint8_t[crc_input_size];

    // Write type (big endian) first
    crc_input[0] = (chunk.type >> 24) & 0xFF;
    crc_input[1] = (chunk.type >> 16) & 0xFF;
    crc_input[2] = (chunk.type >> 8) & 0xFF;
    crc_input[3] = chunk.type & 0xFF;
    memcpy(crc_input + 4, chunk_data, chunk.length);

    // Now calculate CRC over type + data
    uint32_t crc = img_png_crc32(crc_input, crc_input_size);
    delete[] crc_input;

    if (crc != chunk.crc) {
      IMG_DEBUG_LOG("Invalid CRC for chunk type: %x, expected: %x\n",
                    chunk.type, crc);
      delete[] chunk_data;  // Free the chunk data
      return false;
    }

    if (chunk.type == s_IMG_PNG_IHDR_CHUNK) {
      std::unique_ptr<img_png_ihdr_chunck> ihdr =
          std::make_unique<img_png_ihdr_chunck>();
      temp_ihdr_ptr = ihdr.get();

      ihdr->width = (chunk_data[0] << 24) | (chunk_data[1] << 16) |
                    (chunk_data[2] << 8) | chunk_data[3];
      ihdr->height = (chunk_data[4] << 24) | (chunk_data[5] << 16) |
                     (chunk_data[6] << 8) | chunk_data[7];
      ihdr->bit_depth = chunk_data[8];
      ihdr->color_type = chunk_data[9];
      ihdr->compression_method = chunk_data[10];
      ihdr->filter_method = chunk_data[11];
      ihdr->interlace_method = chunk_data[12];

      ihdr->length = chunk.length;
      ihdr->crc = chunk.crc;
      ihdr->type = chunk.type;

      std::shared_ptr<img_png_type_chunk_map> ihdr_map =
          std::make_shared<img_png_type_chunk_map>(chunk.type, std::move(ihdr));
      png_chunks.emplace_back(std::move(ihdr_map));
    }

    if (chunk.type == s_IMG_PNG_IDAT_CHUNK) {
      std::unique_ptr<img_png_idat_chunk> idat =
          std::make_unique<img_png_idat_chunk>();
      idat->ihdr_ptr = temp_ihdr_ptr;

      idat->length = chunk.length;
      idat->data = chunk_data;

      idat->crc = chunk.crc;
      idat->type = chunk.type;
      idat->length = chunk.length;

      std::shared_ptr<img_png_type_chunk_map> idat_map =
          std::make_shared<img_png_type_chunk_map>(chunk.type, std::move(idat));
      png_chunks.emplace_back(std::move(idat_map));
    }

    if (chunk.type == s_IMG_PNG_PLTE_CHUNK) {
      std::unique_ptr<img_png_plte_chunk> plte =
          std::make_unique<img_png_plte_chunk>();
      plte->palette.reserve(chunk.length / 3);
      for (uint32_t i = 0; i < chunk.length; i += 3) {
        plte->palette.emplace_back(
            img_rgb(chunk_data[i], chunk_data[i + 1], chunk_data[i + 2]));
      }

      plte->length = chunk.length;
      plte->crc = chunk.crc;
      plte->type = chunk.type;

      std::shared_ptr<img_png_type_chunk_map> plte_map =
          std::make_shared<img_png_type_chunk_map>(chunk.type, std::move(plte));
      png_chunks.emplace_back(std::move(plte_map));
    }

    if (chunk.type == s_IMG_PNG_tRNS_CHUNK) {
      std::unique_ptr<img_png_tRNS_chunk> tRNS =
          std::make_unique<img_png_tRNS_chunk>();
      tRNS->alpha.reserve(chunk.length);
      for (uint32_t i = 0; i < chunk.length; ++i) {
        tRNS->alpha.emplace_back(chunk_data[i]);
      }

      tRNS->length = chunk.length;
      tRNS->crc = chunk.crc;
      tRNS->type = chunk.type;

      std::shared_ptr<img_png_type_chunk_map> tRNS_map =
          std::make_shared<img_png_type_chunk_map>(chunk.type, std::move(tRNS));
      png_chunks.emplace_back(std::move(tRNS_map));
    }

    if (chunk.type == s_IMG_PNG_gAMA_CHUNK) {
      std::unique_ptr<img_png_gAMA_chunk> gAMA =
          std::make_unique<img_png_gAMA_chunk>();
      temp_gamma_ptr = gAMA.get();  // Store the pointer to gAMA chunk

      if (chunk.length != 4) {
        gAMA->has_gamma = false;
        IMG_DEBUG_LOG("Invalid gAMA chunk length: %x, expected: 4\n",
                      chunk.length);
      } else {
        gAMA->has_gamma = true;
        gAMA->gamma = (chunk_data[0] << 24) | (chunk_data[1] << 16) |
                      (chunk_data[2] << 8) | chunk_data[3];

        gAMA->length = chunk.length;
        gAMA->crc = chunk.crc;
        gAMA->type = chunk.type;

        gAMA->cHRM = temp_chrm_ptr;  // Link gAMA to cHRM chunk
        gAMA->iCCP = temp_iccp_ptr;  // Link gAMA to iCCP chunk
        gAMA->sRGB = temp_srgb_ptr;  // Link gAMA to sRGB chunk
      }
    }

    if (chunk.type == s_IMG_PNG_cHRM_CHUNK) {
      std::unique_ptr<img_png_cHRM_chunk> cHRM =
          std::make_unique<img_png_cHRM_chunk>();
      temp_chrm_ptr = cHRM.get();  // Store the pointer to cHRM chunk

      if (chunk.length != 32) {
        cHRM->has_chrm = false;
        IMG_DEBUG_LOG("Invalid cHRM chunk length: %x, expected: 32\n",
                      chunk.length);
      } else {
        cHRM->has_chrm = true;
        auto read_fixed_point = [&](int offset) {
          return ((chunk_data[offset] << 24) | (chunk_data[offset + 1] << 16) |
                  (chunk_data[offset + 2] << 8) | chunk_data[offset + 3]) /
                 100000.0f;
        };

        cHRM->white_x = read_fixed_point(0);
        cHRM->white_y = read_fixed_point(4);
        cHRM->red_x = read_fixed_point(8);
        cHRM->red_y = read_fixed_point(12);
        cHRM->green_x = read_fixed_point(16);
        cHRM->green_y = read_fixed_point(20);
        cHRM->blue_x = read_fixed_point(24);
        cHRM->blue_y = read_fixed_point(28);

        cHRM->length = chunk.length;
        cHRM->crc = chunk.crc;
        cHRM->type = chunk.type;

        cHRM->gAMA = temp_gamma_ptr;  // Link cHRM to gAMA chunk
        cHRM->iCCP = temp_iccp_ptr;   // Link cHRM to iCCP chunk
        cHRM->sRGB = temp_srgb_ptr;   // Link cHRM to sRGB chunk

        std::shared_ptr<img_png_type_chunk_map> cHRM_map =
            std::make_shared<img_png_type_chunk_map>(chunk.type,
                                                     std::move(cHRM));
        png_chunks.emplace_back(std::move(cHRM_map));
      }
    }

    if (chunk.type == s_IMG_PNG_IEND_CHUNK) {
      std::unique_ptr<img_png_iend_chunk> iend =
          std::make_unique<img_png_iend_chunk>();
      iend->is_iend = true;

      iend->length = chunk.length;
      iend->crc = chunk.crc;
      iend->type = chunk.type;

      std::shared_ptr<img_png_type_chunk_map> iend_map =
          std::make_shared<img_png_type_chunk_map>(chunk.type, std::move(iend));
      png_chunks.emplace_back(std::move(iend_map));
      break;
    }
  }
}

namespace img {
std::shared_ptr<image_specification> import(
    const std::filesystem::path& filepath, image_format format,
    bool flip) noexcept {
  if (format == image_format::png) {
    std::shared_ptr<image_specification> spec =
        std::make_shared<image_specification>();

    spec->filepath = filepath;
    spec->filename = filepath.filename().string();
    spec->xid = unique_xid();

    compressed_file_data file(filepath);
    if (img_png_read(spec, file.data, file.size)) {
      IMG_DEBUG_LOG("PNG file read successfully: %s\n",
                    filepath.string().c_str());
      return spec;
    } else {
      IMG_DEBUG_LOG("Failed to read PNG file: %s\n", filepath.string().c_str());
      return nullptr;
    }
  }
}
}  // namespace img
