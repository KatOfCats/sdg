//
// Created by Luis Yanes (EI) on 21/09/2018.
//

#ifndef BSG_EXTENDEDKMERBIN_HPP
#define BSG_EXTENDEDKMERBIN_HPP

#include <vector>
#include <cstdint>
#include <mutex>
#include <atomic>
#include <list>
#include <fstream>
#include <iterator>
#include "Seq2bit.hpp"

/*
 * This class needs to live within a parent that manages it's memory!
 */

class ExtendedKmerBin;

// TODO: Use this class instead of the vectors in BinBufferWriterQueue
class BinBuffer{
    int bin_id;
    std::mutex bin_mtx;
    std::list<std::vector<unsigned char>> bin_buffer;
    int current_size;
    std::ofstream bin_file;
    uint64_t total_kmers;
    uint64_t total_size;

public:
    BinBuffer(int id) : bin_id(id)
    {
        // Open file init
        // Init buffer
    }
    void insert(std::vector<unsigned char> &buf, size_t len, uint32_t n_superkmers, uint32_t n_kmers) {}
    void write_bin(){}

};

class BinInfo {
public:
    uint64_t size{0};
    uint64_t n_rec{0};
    uint64_t n_plus_x_recs{0};

    const bool operator<(const BinInfo& o) const {
        return std::tie(size,n_rec,n_plus_x_recs) < std::tie(o.size, o.n_rec, o.n_plus_x_recs);
    }

    const bool operator==(const BinInfo& o) const {
        return std::tie(size,n_rec,n_plus_x_recs) == std::tie(o.size, o.n_rec, o.n_plus_x_recs);
    }
};

class BinBufferWriterQueue{

    std::vector<std::mutex> bin_mtx;                                   /// Lock per file/queue
    std::vector<std::list<std::vector<unsigned char>>> bin_buffers;    /// Vector bin buffers *shared*
    std::vector<int> bin_buffer_size;                                  /// Current bin_no size
    std::atomic<unsigned int> max_bin{0};                              /// Largest buffer (0 represents empty or no bin)
    std::atomic<uint64_t> max_bin_size{0};                             /// Size of largest buffer

public:
    std::vector<BinInfo> bin_desc;
    std::vector<std::FILE*> bin_files;                              /// Files containing bins *shared*
    std::atomic<uint64_t> num_writers{0};
    BinBufferWriterQueue(int num_bins);

    unsigned int getMaxBin() const { return max_bin; }
    void insert(int bin_no, std::vector<unsigned char> &buf, size_t len, uint32_t n_superkmers, uint32_t n_kmers);

    void write_bin(int bin_no);

    bool empty();

    void flush() {
        // Make sure all bins are empty before closing the files and destroying all objects
        for (int i = 0; i < bin_buffers.size(); i++) {
            write_bin(i);
            std::fflush(bin_files[i]);
            std::rewind(bin_files[i]);
        }
    }

    ~BinBufferWriterQueue();
};

class ExtendedKmerBin {
    const uint8_t K;
    const uint16_t bin_no{0};
    unsigned int buffer_size;
    unsigned int buffer_pos=0;
    unsigned int n_super_kmers{0};
    unsigned int n_kmers{0};

    std::vector<unsigned char> buffer;
    BinBufferWriterQueue& queue;
public:
    unsigned int total_super_kmers{0};
    unsigned int total_kmers{0};
    ExtendedKmerBin(uint16_t bin_no, uint8_t k, BinBufferWriterQueue& queue, unsigned int max_size);

    void store_superkmer(const Seq2bit &seq, unsigned int offset, unsigned int len);

    void flush();
};


#endif //BSG_EXTENDEDKMERBIN_HPP
