//
// Created by Luis Yanes (EI) on 12/09/2018.
//

#ifndef BSG_KMERCOUNTER_HPP
#define BSG_KMERCOUNTER_HPP

#include <vector>
#include <thread>
#include <sglib/readers/FileReader.h>
#include <sglib/datastores/PairedReadsDatastore.hpp>
#include "MPMCQueue.hpp"

/**
 * The members of this class are organised as follows:
 *      Processor
 *      Communication mechanism
 *      Processor
 *      Communication mechanism
 *      .
 *      .
 *      .
 *
 * The idea is to have the code represent the flow of the data through the pipeline,
 * the processors communicate using asynchronous operations over multi-producer multi-consumer
 * queues
 * TODO: Make the queues lock-free and test performance
 *
 * Finally the last step of the pipeline should write the kmers to a database of kmer-counts.
 */

class BinChunk{};
class SortedBin{};
class DiskBin{};

class KDB{};

class KmerCounter {
    using Kmer = std::pair<bool, uint64_t>;
    std::vector<std::thread> file_readers;
    MPMCBoundedQueue<std::vector<Kmer>> file_parts;
    std::vector<std::thread> kmer_producers;
    MPMCBoundedQueue<BinChunk> bin_chunks;
    std::vector<std::thread> disk_writer;

    std::vector<std::ofstream> bin_files;

    std::vector<std::thread> bin_readers;
    MPMCBoundedQueue<DiskBin> bins;
    std::vector<std::thread> bin_sorter;
    MPMCBoundedQueue<SortedBin> sorted_bins;
    std::vector<std::thread> kmer_completer;
    KDB kmer_database;

    uint8_t K = 27;
public:
    explicit KmerCounter(uint8_t k) : K(k), file_parts(2^14), bin_chunks(2^10), bins(2^8), sorted_bins(2^8) {}
    std::function<void()> read_file(PairedReadsDatastore &r);

    void processReads(std::vector<PairedReadsDatastore> &read_files);
    void processBins();
};


#endif //BSG_KMERCOUNTER_HPP
