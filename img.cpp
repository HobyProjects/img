/**
 * @file img.cpp
 * @brief Image library implementation.
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
        img_bitreader redr(data, fileSize);
        uint64 img_signature = redr.read_bits<uint64>();
        
        /******************************************************************
         *                         PNG IMAGE FILES                        *
         ******************************************************************/
        if(img_signature == IMG_PNG_SIGNATURE){
            IMG_PNGCHUNK* root{nullptr};
            IMG_PNGCHUNK* temp_chunk{nullptr};
            IMG_PNGIHDR* img_ihdr{nullptr};
            while(IMG_TRUE){
                IMG_PNGCHUNK* current_chunck = new IMG_PNGCHUNK();
                if(root == nullptr){
                    root = current_chunck;
                }

                current_chunck->length = redr.read_bits<uint32>();
                current_chunck->type = redr.read_bits<uint32>();

                if(current_chunck->type == IMG_PNG_CHUNK_IHDR){
                    img_ihdr = new IMG_PNGIHDR();
                    img_ihdr->width = redr.read_bits<uint32>();
                    img_ihdr->height = redr.read_bits<uint32>();
                    img_ihdr->channels = redr.read_bits<uint8>();
                    img_ihdr->color_type = redr.read_bits<uint8>();
                    img_ihdr->compression_method = redr.read_bits<uint8>();
                    img_ihdr->filter_method = redr.read_bits<uint8>();
                    img_ihdr->interlace_method = redr.read_bits<uint8>();
                }
                
                if(current_chunck->type == IMG_PNG_CHUNK_IEND)
                    break;
            }
        }
    }
}