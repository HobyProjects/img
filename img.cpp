#include "img.hpp"

namespace img{
    std::shared_ptr<specification> read(const std::filesystem::path & filePath, bool flip){
        if(!std::filesystem::exists(filePath)){
            IMG_ASSERT(false, "File does not exist: ", filePath.string().c_str());
            IMG_LOG("File does not exist: %s", filePath.string().c_str());
            return nullptr;
        }

        std::ifstream file(filePath, std::ios::binary);
        if(!file.is_open()){
            IMG_ASSERT(false, "Failed to open file: ", filePath.string().c_str());
            IMG_LOG("Failed to open file: %s", filePath.string().c_str());
            return nullptr;
        }

        file.seekg(0, std::ios::end);
        std::streampos fileSize = file.tellg();
        file.seekg(0, std::ios::beg);

        uint8* data = new uint8[fileSize];
        file.read((char*)data, fileSize);
        file.close();

        std::shared_ptr<specification> img_spec = std::make_shared<specification>();
    }
}