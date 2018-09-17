//
// Created by Luis Yanes (EI) on 12/09/2018.
//

#ifndef BSG_KMERCOUNTER_HPP
#define BSG_KMERCOUNTER_HPP

#include <vector>
#include <thread>
#include <list>
#include <sglib/readers/FileReader.h>
#include <sglib/datastores/PairedReadsDatastore.hpp>
#include <sglib/counter/LBQueue.hpp>
#include "MPMCQueue.hpp"
#include "ExtendedKmerBin.hpp"

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
 *
 * Finally the last step of the pipeline should write the kmers to a database of kmer-counts.
 */
class MinimiserTable{
public:
    unsigned int get_bin_id(unsigned int n) {return n;}
};
class Mmer{
    const uint8_t K;
    uint32_t value;
public:
    Mmer(uint8_t k) : K(k) {}
    uint32_t get(){return value;}
    void insert(const unsigned char *c){

    }
    void set(Mmer c){value = c.value;}

    bool operator<(const Mmer &o) const { return value < o.value; }
    bool operator>(const Mmer &o) const { return value > o.value; }
    bool operator==(const Mmer &o) const {return value == o.value;}
    bool operator<=(const Mmer &o) const {return !(*this < o); }
};


class Seq2bit {
    std::vector<unsigned char> m_data;
    constexpr static std::array<char,256> codes =
            {
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, //32
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, //64
            0,4,1,4,4,4,4,4,2,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4, //96
            0,4,1,4,4,4,4,4,2,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4, //128
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, //160
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, //192
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, //224
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, //256
            };

public:
    Seq2bit(std::string seq) {
        for (const auto &c: seq){
            m_data.emplace_back(codes[c]);
        }
    }
    std::size_t size() {return m_data.size();}
    unsigned char operator[](int i) const { return static_cast<unsigned char>(m_data[i]); }
    void insert_nt(unsigned char nt) {m_data.emplace_back(codes[nt]);}
    unsigned char* data() {return m_data.data();}
};

constexpr std::array<char,256> Seq2bit::codes;





class KmerCounter {
    const uint8_t K;
    const uint8_t signature_len{5};
    const unsigned int n_bins{512};
    std::array<char,256> codes={4};

    int num_datastore_readers{4};
    int num_kmer_splitters{4};
    int num_bin_writers{2};

    int num_smallK_counters{4};
    std::mutex write_lock;
    std::mutex smallK_merge_lock;
    using Kmer = std::pair<bool, uint64_t>;
    using KmerArray = std::array<Kmer, 300>;
    std::vector<std::thread> kmer_producers;
    KmerArrayQueue LBkmer_chunks;
    MPMCBoundedQueue<KmerArray> kmer_chunks;

    std::vector<std::thread> kmer_partitioners;
    std::vector<std::thread> kmer_counters;
    std::vector<std::thread> kmer_bin_writers;

    std::atomic<uint64_t> total_reads_processed{0};
    std::atomic<uint64_t> total_reads_kmerised{0};
    std::atomic<uint64_t> total_kmers_produced{0};

    MinimiserTable mmer_table;

    void read_file_small_k(const PairedReadsDatastore &r, int thread_id);
    void read_file_large_k(const PairedReadsDatastore &r, BinBufferWriterQueue &wrt_queue, int thread_id);

    void count_small_k(int thread_id, std::vector<uint64_t> &smallK_counts);

    void LBread_file(PairedReadsDatastore &r, int thread_id);
    void LBprocess_kmers(int thread_id);
public:
    explicit KmerCounter(uint8_t k) : K(k), signature_len(9), LBkmer_chunks(num_datastore_readers), kmer_chunks(1024) {
        codes['a'] = codes['A'] = 0;
        codes['c'] = codes['C'] = 1;
        codes['g'] = codes['G'] = 2;
        codes['t'] = codes['T'] = 3;
    }

    void processReads(std::vector<PairedReadsDatastore> &read_files, bool sequential);

    void processReadsSmallK(std::vector<PairedReadsDatastore> &read_files);
    void processReadsSmallKSequential(std::vector<PairedReadsDatastore> &read_files);
    void processReadsLargeK(std::vector<PairedReadsDatastore> &read_files);

    void processBins(BinBufferWriterQueue &wrt_queue);

    void make_minimiser_send_to_bin(std::string sequence, std::vector<ExtendedKmerBin> &bins,
                                    std::list<Seq2bit> &reads_2bit);
};


#endif //BSG_KMERCOUNTER_HPP
