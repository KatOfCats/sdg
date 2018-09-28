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
class BinBufferWriterQueue{

    std::vector<std::mutex> bin_mtx;                                    /// Lock per file/queue
    std::vector<std::list<std::vector<unsigned char>>> bin_buffers;    /// Vector bin buffers *shared*
    std::vector<uint64_t> bin_buffer_size;                          /// Current bin_no size
    std::atomic<unsigned int> max_bin{0};                               /// Largest buffer (0 represents empty or no bin!)
    std::atomic<uint64_t> max_bin_size{0};                              /// Size of largest buffer
    std::vector<std::ofstream> bin_files;                               /// Files containing bins *shared*
public:
    std::atomic<uint64_t> num_writers;
    BinBufferWriterQueue(int num_bins);
    std::vector<ExtendedKmerBin> bins;                      /// Objects managing access to the bin_buffers

    unsigned int getMaxBin() const { return max_bin; }
    void insert(int bin_no, std::vector<unsigned char> &buf, size_t len, uint32_t n_superkmers, uint32_t n_kmers);

    void write_bin(int bin_no);

    bool empty();

    ~BinBufferWriterQueue();
};

class ExtendedKmerBin {
    const uint8_t K;
    const uint16_t bin_no;
    unsigned int buffer_size;
    unsigned int buffer_pos=0;

    unsigned int n_super_kmers;
    unsigned int n_kmers;
    std::vector<unsigned char> buffer;
    BinBufferWriterQueue& queue;
public:
    ExtendedKmerBin(uint16_t bin_no, uint8_t k, BinBufferWriterQueue& queue, unsigned int max_size);

    void store_superkmer(const Seq2bit &seq, unsigned int offset, unsigned int len);

    void flush();
};


#endif //BSG_EXTENDEDKMERBIN_HPP
