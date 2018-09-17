//
// Created by Luis Yanes (EI) on 21/09/2018.
//

#ifndef BSG_EXTENDEDKMERBIN_HPP
#define BSG_EXTENDEDKMERBIN_HPP

#include <vector>
#include <cstdint>
#include <mutex>
#include <atomic>
#include <queue>
#include <fstream>
#include <iterator>

/*
 * This class needs to live within a parent that manages it's memory!
 */

class ExtendedKmerBin;
class BinBufferWriterQueue{
    std::vector<std::mutex> bin_mtx;                                    /// Lock per file/queue
    std::vector<std::queue<std::vector<unsigned char>>> bin_buffers;    /// Vector bin buffers *shared*
    std::vector<unsigned int> bin_buffer_size;                          /// Current bin_no size
    std::atomic<unsigned int> max_bin{0};                               /// Largest buffer (0 represents empty or no bin)
    std::atomic<uint64_t> max_bin_size{0};                              /// Size of largest buffer
    std::vector<std::ofstream> bin_files;                               /// Files containing bins *shared*
public:
    std::atomic<uint64_t> num_writers;
    BinBufferWriterQueue(int num_bins) {
        bin_buffers.resize(1);
        for (int i = 0; i < num_bins; i++) {
            std::string name("file"+std::to_string(i)+".bin");
            bin_files.emplace_back(name);
        }
    }
    std::vector<ExtendedKmerBin> bins;                      /// Objects managing access to the bin_buffers
    unsigned int getMaxBin() const { return max_bin; }
    void insert(int bin_no, std::vector<unsigned char> buf, size_t len, uint32_t n_superkmers, uint32_t n_kmers) {
        std::lock_guard<std::mutex> lg(bin_mtx[bin_no]);
        if (bin_buffer_size[bin_no] > max_bin_size) {
            max_bin_size = bin_buffer_size[bin_no];
            max_bin = bin_no;
        }
        bin_buffer_size[bin_no]+=len;
        bin_buffers[bin_no].push(buf);
    }

    void write_bin(int bin_no) {
        std::vector<unsigned char> buffer;
        while (!bin_buffers[bin_no].empty()) {
            const auto &bbref(bin_buffers[bin_no].front());
            std::copy(bbref.begin(), bbref.end(), buffer.end());
            bin_buffer_size[bin_no] -= buffer.size();
            bin_buffers[bin_no].pop();
            std::copy(buffer.begin(), buffer.end(), std::ostream_iterator<unsigned char>(bin_files[bin_no]));
        }

        bin_buffer_size[bin_no] = 0;

        max_bin_size = bin_buffer_size[0];
        max_bin = 0;

        for (int i = 1; i <= bin_buffers.size(); i++) {
            if (bin_buffer_size[i] > max_bin_size) {
                max_bin_size = bin_buffer_size[i];
                max_bin = i;
            }
        }
    }

    bool empty() {
        if (num_writers > 0) return false;
        if (max_bin>0) return false;
        return true;
    }
};

class ExtendedKmerBin {
    const uint8_t K;
    const uint16_t bin_no;
    unsigned int buffer_size;
    unsigned int buffer_pos;

    unsigned int n_super_kmers;
    unsigned int n_kmers;
    std::vector<unsigned char> buffer;
    BinBufferWriterQueue& queue;
public:
    ExtendedKmerBin(uint16_t bin_no, uint8_t k, BinBufferWriterQueue& queue, unsigned int max_size) :
    K(k),
    queue(queue),
    bin_no(bin_no)
    {
        buffer_size = max_size;
        buffer.resize(max_size);
    }

    void store_superkmer(const unsigned char *seq, unsigned int len){

        uint32_t bytes = 1 + (len + 3) / 4;
        if(buffer_pos + bytes > buffer_size)
        {
            flush();
//            pmm_bins->reserve(buffer); // Here synchronizes for memory
            buffer_pos		= 0;
            n_kmers			= 0;
            n_super_kmers		= 0;
        }


        buffer[buffer_pos++] = static_cast<unsigned char>(len - K);
        for(uint32_t i = 0, j = 0 ; i < len / 4 ; ++i,j+=4)
            buffer[buffer_pos++] = (seq[j] << 6) + (seq[j + 1] << 4) + (seq[j + 2] << 2) + seq[j + 3];
        switch (len%4)
        {
            case 1:
                buffer[buffer_pos++] = (seq[len-1] << 6);
                break;
            case 2:
                buffer[buffer_pos++] = (seq[len-2] << 6) + (seq[len-1] << 4);
                break;
            case 3:
                buffer[buffer_pos++] = (seq[len-3] << 6) + (seq[len-2] << 4) + (seq[len-1] << 2);
                break;
            default:break;
        }

        ++n_super_kmers;
        n_kmers += len - K + 1;
    }

    void flush() {
        // Write buffer to output queue (sync)
        queue.insert(bin_no, buffer,buffer_pos, n_super_kmers, n_kmers);
    }
};


#endif //BSG_EXTENDEDKMERBIN_HPP
