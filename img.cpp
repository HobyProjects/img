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
        explicit bit_reader(uint8_t* data, size_t size) : 
            m_data(data), m_size(size), m_bit_offset(0) {}
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
        template<typename T>
        T read() {
            static_assert(std::is_integral_v<T>, "T must be an integral type");
            T result = 0;
            for (size_t i = 0; i < sizeof(T); ++i) {
                result |= (m_data[m_bit_offset >> 3] >> (m_bit_offset & 7)) << (i * 8);
                ++m_bit_offset;
            }
            return result;
        }

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
        uint8_t* m_data{nullptr};   // Pointer to the data buffer
        size_t m_size{0};           // Size of the data buffer in bytes
        size_t m_bit_offset{0};     // Current bit offset in the data buffer
};

/**
 * @struct compressed_file_data
 * @brief This structure represents a compressed file.
 * It contains the filepath, filename, data pointer, and size of the file.
 */
struct compressed_file_data{
    std::filesystem::path filepath; ///< Filepath of the compressed file
    std::string filename;           ///< Filename of the compressed file
    uint8_t* data{nullptr};         ///< Pointer to the compressed file data
    size_t size{0};                 ///< Size of the compressed file data in bytes

    explicit compressed_file_data(const std::filesystem::path& file) {
        filepath = file;
        filename = file.filename().string();
        std::ifstream file_stream(filepath, std::ios::binary | std::ios::ate);
        if (!file_stream.is_open()) {
            IMG_DEBUG_LOG("Failed to open file: %s\n", filepath.string().c_str());
            throw std::runtime_error("Failed to open file: " + filepath.string());
        }
        size = file_stream.tellg();                                     // Get the size of the file
        file_stream.seekg(0, std::ios::beg);                            // Move the cursor to the beginning of the file
        data = new uint8_t[size];                                       // Allocate memory for the file data
        file_stream.read(reinterpret_cast<char*>(data), size);          // Read the file data into the buffer    
        file_stream.close();                                            // Close the file stream
    };

    compressed_file_data(const compressed_file_data&) = delete;                 // Delete copy constructor
    compressed_file_data& operator=(const compressed_file_data&) = delete;      // Delete copy assignment operator
    compressed_file_data(compressed_file_data&&) = default;                     // Default move constructor
    compressed_file_data& operator=(compressed_file_data&&) = default;          // Default move assignment operator

    ~compressed_file_data() {
        delete[] data; // Free the allocated memory for the file data
    }

    void clear() noexcept {
        delete[] data;          // Free the allocated memory for the file data
        data = nullptr;         // Set the pointer to nullptr
        size = 0;               // Reset the size to 0
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

/******************************************************************
 *                         PNG IMAGE FILES                        *
 * ****************************************************************/
static constexpr uint64_t s_IMG_PNG_SIGNATURE  = 0x89504E470D0A1A0A; // PNG file signature
static constexpr uint32_t s_IMG_PNG_IHDR_CHUNK = 0x52414449; // IHDR chunk type
static constexpr uint32_t s_IMG_PNG_IEND_CHUNK = 0x444E4549; // IEND chunk type
static constexpr uint32_t s_IMG_PNG_IDAT_CHUNK = 0x54455849; // IDAT chunk type
static constexpr uint32_t s_IMG_PNG_PLTE_CHUNK = 0x54494C50; // PLTE chunk type
static constexpr uint32_t s_IMG_PNG_tRNS_CHUNK = 0x74524E53; // tRNS chunk type
static constexpr uint32_t s_IMG_PNG_gAMA_CHUNK = 0x67414D41; // gAMA chunk type
static constexpr uint32_t s_IMG_PNG_cHRM_CHUNK = 0x6348524D; // cHRM chunk type
static constexpr uint32_t s_IMG_PNG_sRGB_CHUNK = 0x73524742; // sRGB chunk type
static constexpr uint32_t s_IMG_PNG_bKGD_CHUNK = 0x624B4744; // bKGD chunk type
static constexpr uint32_t s_IMG_PNG_pHYs_CHUNK = 0x70485973; // pHYs chunk type
static constexpr uint32_t s_IMG_PNG_iTXt_CHUNK = 0x69545874; // iTXt chunk type
static constexpr uint32_t s_IMG_PNG_tEXt_CHUNK = 0x74455874; // tEXt chunk type
static constexpr uint32_t s_IMG_PNG_zTXt_CHUNK = 0x7A545874; // zTXt chunk type
static constexpr uint32_t s_IMG_PNG_sBIT_CHUNK = 0x73424954; // sBIT chunk type
static constexpr uint32_t s_IMG_PNG_iCCP_CHUNK = 0x69434350; // iCCP chunk type
static constexpr uint32_t s_IMG_PNG_fRAc_CHUNK = 0x66726163; // fRAc chunk type
static constexpr uint32_t s_IMG_PNG_hIST_CHUNK = 0x68495354; // hIST chunk type
static constexpr uint32_t s_IMG_PNG_tIME_CHUNK = 0x74494D45; // tIME chunk type
static constexpr uint32_t s_IMG_PNG_pCAL_CHUNK = 0x7043414C; // pCAL chunk type
static constexpr uint32_t s_IMG_PNG_sCAL_CHUNK = 0x7343414C; // sCAL chunk type

/**
 * @struct img_png_chunk
 * @brief This structure represents a PNG chunk.
 * It contains information about the chunk's length, type, data, CRC value,
 * and pointers to the next and previous chunks in the linked list.
 * It also contains a pointer to the parent chunk (if any).
 */
struct img_png_chunk{
    uint32_t length{0};         // Length of the chunk data
    uint32_t type{0};           // Type of the chunk
    uint32_t crc{0};            // CRC value for the chunk

    img_png_chunk() = default; 
    virtual ~img_png_chunk() = default;
};

/**
 * @struct img_png_type_chunk_map
 * @brief This structure maps a chunk type to its corresponding chunk.
 * It contains the chunk type and a pointer to the chunk.
 */
struct img_png_type_chunk_map{
    uint32_t type{0};                   // Type of the chunk
    img_png_chunk* chunk{nullptr};      // Pointer to the chunk

    img_png_type_chunk_map() = default;
    img_png_type_chunk_map(uint32_t t, img_png_chunk* c) : type(t), chunk(c) {}
};

/**
 * @struct img_png_ihdr_chunck
 * @brief This structure represents the IHDR chunk of a PNG image.
 * It contains information about the image's width, height, bit depth,
 * color type, compression method, filter method, and interlace method.
 */
struct img_png_ihdr_chunck : img_png_chunk {
    uint32_t width{0};                          // Width of the image in pixels
    uint32_t height{0};                         // Height of the image in pixels
    uint8_t bit_depth{0};                       // Bit depth of the image
    uint8_t color_type{0};                      // Color type of the image
    uint8_t compression_method{0};              // Compression method of the image
    uint8_t filter_method{0};                   // Filter method of the image
    uint8_t interlace_method{0};                // Interlace method of the image
    img_png_chunk* parent_chunck{nullptr};      // Pointer to the parent chunk (IHDR chunk)

    img_png_ihdr_chunck() = default;
    virtual ~img_png_ihdr_chunck() = default;
};

/**
 * @struct img_png_idat_chunk
 * @brief This structure represents the IDAT chunk of a PNG image.
 * It contains a pointer to the IDAT chunk data and its length.
 */
struct img_png_idat_chunk : img_png_chunk {
    uint8_t* data{nullptr};                       // Pointer to the IDAT chunk data
    uint32_t length{0};                           // Length of the IDAT chunk data
    img_png_chunk* parent_chunck{nullptr};        // Pointer to the IDAT chunk

    img_png_idat_chunk() = default;
    virtual ~img_png_idat_chunk() {
        delete[] data;                            // Free the IDAT chunk data
    }
};

static bool img_png_read(const std::shared_ptr<img::image_specification>& spec, uint8_t* data, size_t size) noexcept {
    bit_reader reader(data, size);                  // Create a bit reader for the data
    uint64_t signature = reader.read<uint64_t>();   // Read the PNG signature
    if(signature != s_IMG_PNG_SIGNATURE){
        IMG_DEBUG_LOG("Invalid PNG signature: %llx\nFile:%s is not a PNG file", signature, spec->filename.c_str());
        return false;
    }

    std::vector<img_png_type_chunk_map> png_chunks; // Vector to hold PNG chunks 

    while(true){
        img_png_chunk chunk;                           
        chunk.length = reader.read<uint32_t>();        
        chunk.type = reader.read<uint32_t>();          

        if(chunk.type == s_IMG_PNG_IHDR_CHUNK){
            img_png_ihdr_chunck* ihdr = new img_png_ihdr_chunck();
            ihdr->width = reader.read<uint32_t>();                  
            ihdr->height = reader.read<uint32_t>();                 
            ihdr->bit_depth = reader.read<uint8_t>();               
            ihdr->color_type = reader.read<uint8_t>();              
            ihdr->compression_method = reader.read<uint8_t>();      
            ihdr->filter_method = reader.read<uint8_t>();           
            ihdr->interlace_method = reader.read<uint8_t>();        
            chunk.crc = reader.read<uint32_t>();                    
            png_chunks.emplace_back(img_png_type_chunk_map(chunk.type, ihdr)); 
        }

        if(chunk.type == s_IMG_PNG_IDAT_CHUNK){
            img_png_idat_chunk* idat = new img_png_idat_chunk();
            idat->length = chunk.length;                             
            idat->data = reader.read(chunk.length);                  
            chunk.crc = reader.read<uint32_t>();                     
            png_chunks.emplace_back(img_png_type_chunk_map(chunk.type, idat));
        }

        if(chunk.type == s_IMG_PNG_IEND_CHUNK){
            break;
        }
    }
}

namespace img{
    std::shared_ptr<image_specification> import(const std::filesystem::path & filepath, image_format format, bool flip) noexcept{
        if(format == image_format::png){
            std::shared_ptr<image_specification> spec = std::make_shared<image_specification>();

            spec->filepath = filepath;
            spec->filename = filepath.filename().string();
            spec->xid = unique_xid();

            compressed_file_data file(filepath);
            if(img_png_read(spec, file.data, file.size)){
                IMG_DEBUG_LOG("PNG file read successfully: %s\n", filepath.string().c_str());
                return spec;
            }else{
                IMG_DEBUG_LOG("Failed to read PNG file: %s\n", filepath.string().c_str());
                return nullptr;
            }
        }
    }
}


