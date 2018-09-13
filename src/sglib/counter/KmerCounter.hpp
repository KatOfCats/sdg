//
// Created by Luis Yanes (EI) on 12/09/2018.
//

#ifndef BSG_KMERCOUNTER_HPP
#define BSG_KMERCOUNTER_HPP

#include <vector>
#include <thread>
#include <sglib/readers/FileReader.h>
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
    std::vector<std::thread> file_readers;
    MPMCBoundedQueue<std::vector<FastqRecord>> file_parts;
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
};


#endif //BSG_KMERCOUNTER_HPP
