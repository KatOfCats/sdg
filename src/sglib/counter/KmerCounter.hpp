//
// Created by Luis Yanes (EI) on 12/09/2018.
//

#ifndef BSG_KMERCOUNTER_HPP
#define BSG_KMERCOUNTER_HPP

#include <vector>
#include <thread>
#include <list>
#include <numeric>
#include <sglib/readers/FileReader.hpp>
#include <sglib/datastores/PairedReadsDatastore.hpp>
#include <sglib/factories/KMerCountFactory.hpp>
#include "MinimiserTable.hpp"
#include "ExtendedKmerBin.hpp"

class KmerCount;

class KmerCounter {
    const uint8_t K;
    const uint8_t signature_len{5};
    const unsigned int n_bins{512};

    bool loadDist={false};
    bool loadParts={false};

    int num_datastore_readers{8};
    int num_bin_writers{2};

    std::mutex mmer_table_lock;

    std::vector<std::thread> kmer_producers;

    std::vector<std::thread> kmer_counters;
    std::vector<std::thread> kmer_bin_writers;

    std::vector<std::thread> kmer_dist_generators;

    std::atomic<uint64_t> total_superkmers{0};
    std::atomic<uint64_t> total_reads_kmerised{0};
    std::atomic<uint64_t> total_kmers_produced{0};

    MinimiserTable mmer_table;

    std::ofstream output_file;

    void read_file_small_k(const PairedReadsDatastore &r, int thread_id);
    void read_file_large_k(const PairedReadsDatastore &r, BinBufferWriterQueue &wrt_queue, int thread_id);
    void count_small_k(int thread_id, std::vector<uint64_t> &smallK_counts);

    uint64_t expand64(std::vector<unsigned char> &skmers_from_disk, std::vector<KmerCount64> &kmerCounts_to_sort);
    uint64_t expand(std::vector<unsigned char> &skmers_from_disk, std::vector<KmerCount> &kmerCounts_to_sort);
public:
    explicit KmerCounter(uint8_t k, unsigned int mmer_len, unsigned int n_bins, bool loadDist, bool loadParts) :
    K(k),
    n_bins(n_bins),
    signature_len(mmer_len),
    loadDist(loadDist),
    loadParts(loadParts),
    mmer_table(signature_len, n_bins)
    {
    }

    void processReads(std::vector<PairedReadsDatastore> &read_files);
    void loadDistribution(){mmer_table.loadDistribution();};
    void processReadsSmallK(std::vector<PairedReadsDatastore> &read_files);
    void processReadsSmallKSequential(std::vector<PairedReadsDatastore> &read_files);
    void processReadsLargeK(std::vector<PairedReadsDatastore> &read_files);

    void processBins(BinBufferWriterQueue &wrt_queue);

    void mmer_table_initialiser(const PairedReadsDatastore &r, MinimiserTable &minimiserTable);

    void make_minimiser_send_to_bin(std::string sequence, std::vector<ExtendedKmerBin> &bins,
                                    Seq2bit &read_2bit);

    void calculateDistribution(const std::vector<PairedReadsDatastore> &read_files);

    void createPartitions(BinBufferWriterQueue &wrt_queue, std::vector<PairedReadsDatastore> &read_files);

    void sort_and_store_partitions(BinBufferWriterQueue &wrt_queue);

    void loadPartitions(BinBufferWriterQueue &wrt_queue);
};


#endif //BSG_KMERCOUNTER_HPP
