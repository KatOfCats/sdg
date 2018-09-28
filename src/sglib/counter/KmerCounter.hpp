//
// Created by Luis Yanes (EI) on 12/09/2018.
//

#ifndef BSG_KMERCOUNTER_HPP
#define BSG_KMERCOUNTER_HPP

#include <vector>
#include <thread>
#include <list>
#include <numeric>
#include <sglib/readers/FileReader.h>
#include <sglib/datastores/PairedReadsDatastore.hpp>
#include <sglib/counter/LBQueue.hpp>
#include "MPMCQueue.hpp"
#include "ExtendedKmerBin.hpp"


class Mmer{
    uint32_t str;
    uint32_t mask;
    uint32_t* norm;
    static uint32_t norm5[1 << 10];
    static uint32_t norm6[1 << 12];
    static uint32_t norm7[1 << 14];
    static uint32_t norm8[1 << 16];
    static uint32_t norm9[1 << 18];
    static uint32_t norm10[1 << 20];
    static uint32_t norm11[1 << 22];

    const uint8_t K;
    uint32_t value;

    friend class CSignatureMapper;
    struct _si
    {
        static uint32_t get_rev(uint32_t mmer, uint32_t len)
        {
            uint32_t rev = 0;
            uint32_t shift = len*2 - 2;
            for(uint32_t i = 0 ; i < len ; ++i)
            {
                rev += (3 - (mmer & 3)) << shift;
                mmer >>= 2;
                shift -= 2;
            }
            return rev;
        }



        static void init_norm(uint32_t* norm, uint32_t len)
        {
            uint32_t special = 1 << len * 2;
            for(uint32_t i = 0 ; i < special ; ++i)
            {
                uint32_t rev = get_rev(i, len);
                uint32_t str_val = is_allowed(i, len) ? i : special;
                uint32_t rev_val = is_allowed(rev, len) ? rev : special;
                norm[i] = std::min(str_val, rev_val);
            }
        }

        _si()
        {
            init_norm(norm5, 5);
            init_norm(norm6, 6);
            init_norm(norm7, 7);
            init_norm(norm8, 8);
            init_norm(norm9, 9);
            init_norm(norm10, 10);
            init_norm(norm11, 11);
        }

    }static _init;

public:
    Mmer(uint8_t k) : K(k) {
        switch (K)
        {
            case 5:
                norm = norm5;
                break;
            case 6:
                norm = norm6;
                break;
            case 7:
                norm = norm7;
                break;
            case 8:
                norm = norm8;
                break;
            case 9:
                norm = norm9;
                break;
            case 10:
                norm = norm10;
                break;
            case 11:
                norm = norm11;
                break;
            default:
                break;
        }
        mask = (1 << K * 2) - 1;
        str = 0;

    }
    inline void insert(unsigned char symb)
    {
        str <<= 2;
        str += symb;
        str &= mask;

        value = norm[str];
    }

    uint32_t get(){return value;}
    void insert(const unsigned char *seq){
        switch (K)
        {
            case 5:
                str = (seq[0] << 8) + (seq[1] << 6) + (seq[2] << 4) + (seq[3] << 2) + (seq[4]);
                break;
            case 6:
                str = (seq[0] << 10) + (seq[1] << 8) + (seq[2] << 6) + (seq[3] << 4) + (seq[4] << 2) + (seq[5]);
                break;
            case 7:
                str = (seq[0] << 12) + (seq[1] << 10) + (seq[2] << 8) + (seq[3] << 6) + (seq[4] << 4 ) + (seq[5] << 2) + (seq[6]);
                break;
            case 8:
                str = (seq[0] << 14) + (seq[1] << 12) + (seq[2] << 10) + (seq[3] << 8) + (seq[4] << 6) + (seq[5] << 4) + (seq[6] << 2) + (seq[7]);
                break;
            case 9:
                str = (seq[0] << 16) + (seq[1] << 14) + (seq[2] << 12) + (seq[3] << 10) + (seq[4] << 8) + (seq[5] << 6) + (seq[6] << 4) + (seq[7] << 2) + (seq[8]);
                break;
            case 10:
                str = (seq[0] << 18) + (seq[1] << 16) + (seq[2] << 14) + (seq[3] << 12) + (seq[4] << 10) + (seq[5] << 8) + (seq[6] << 6) + (seq[7] << 4) + (seq[8] << 2) + (seq[9]);
                break;
            case 11:
                str = (seq[0] << 20) + (seq[1] << 18) + (seq[2] << 16) + (seq[3] << 14) + (seq[4] << 12) + (seq[5] << 10) + (seq[6] << 8) + (seq[7] << 6) + (seq[8] << 4) + (seq[9] << 2) + (seq[10]);
                break;
            default:
                break;
        }

        value = norm[str];
    }
    void set(Mmer c){
        str = c.str;
        value = c.value;
    }

    static bool is_allowed(uint32_t mmer, uint32_t len) {
        if ((mmer & 0x3f) == 0x3f)            // TTT suffix
            return false;
        if ((mmer & 0x3f) == 0x3b)            // TGT suffix
            return false;
        if ((mmer & 0x3c) == 0x3c)            // TG* suffix
            return false;

        for (uint32_t j = 0; j < len - 3; ++j)
            if ((mmer & 0xf) == 0)                // AA inside
                return false;
            else
                mmer >>= 2;

        if (mmer == 0)            // AAA prefix
            return false;
        if (mmer == 0x04)        // ACA prefix
            return false;
        if ((mmer & 0xf) == 0)    // *AA prefix
            return false;

        return true;
    }

    bool operator<(const Mmer &o) const { return value < o.value; }
    bool operator>(const Mmer &o) const { return value > o.value; }
    bool operator==(const Mmer &o) const {return value == o.value;}
    bool operator<=(const Mmer &o) const {return !(*this < o); }
};

class MinimiserTable{
    uint32_t mm_len;
    const uint32_t num_bins;
    std::vector<uint32_t> minimiser_table;
    std::vector<uint32_t> bin_map;
public:
    MinimiserTable(uint32_t mm_len, uint32_t num_bins) : mm_len(mm_len), num_bins(num_bins) {
        minimiser_table.resize(1<<(2*mm_len) + 1);
        std::fill(minimiser_table.begin(), minimiser_table.end(), 0);
    }

    void combine(const MinimiserTable& o) {
        std::transform(o.minimiser_table.cbegin(), o.minimiser_table.cend(), o.minimiser_table.cbegin(), minimiser_table.begin(), std::plus<uint32_t>());
    }
    std::size_t size() const {return minimiser_table.size();}
    unsigned int get_bin_id(unsigned int n) const {return bin_map[n];}
    uint32_t& operator[](const std::size_t i){return minimiser_table[i];};

    void calcBins() {
        bin_map.resize(minimiser_table.size());
        std::fill(bin_map.begin(), bin_map.end(), 1);
        auto n_bins = num_bins+1;
        std::vector<uint32_t > sorted(minimiser_table.size());
        std::iota(sorted.begin(),sorted.end(),0);
        std::sort(sorted.begin(),sorted.end(), [&](uint32_t a, uint32_t b) {return minimiser_table[a]>minimiser_table[b];});
        std::list<std::pair<uint32_t, uint64_t>> _stats;
        for (uint32_t i = 0; i < sorted.size() ; ++i)
        {
            if (Mmer::is_allowed(sorted[i], mm_len))
                _stats.push_back(std::make_pair(sorted[i], minimiser_table[sorted[i]]));
        }

        std::list<std::pair<uint32_t, uint64_t>> group;
        uint32_t bin_no = 0;
        uint32_t rev_bin = 512;
        //counting sum
        double sum = 0.0;
        for (auto &i : _stats)
        {
            i.second += 1000;
            sum += i.second;
        }

        double mean = sum / n_bins;
        double max_bin_size = 1.1 * mean;
        uint32_t n = n_bins-2; //one is needed for disabled signatures
        uint32_t max_bins = n_bins - 2;

        while (_stats.size() > n)
        {
            std::pair<uint32_t, uint64_t>& max = _stats.front();

            if (max.second > mean)
            {
                bin_map[max.first] = (rev_bin - bin_no);
                bin_no++;
                sum -= max.second;
                mean = sum / (max_bins - bin_no);
                max_bin_size = 1.1 * mean;

                _stats.pop_front();
                --n;
            }
            else
            {
                //heuristic
                group.clear();
                double tmp_sum = 0.0;
                uint32_t in_current = 0;
                for (auto it = _stats.begin(); it != _stats.end();)
                {
                    if (tmp_sum + it->second < max_bin_size)
                    {
                        tmp_sum += it->second;
                        group.push_back(*it);
                        it = _stats.erase(it);
                        ++in_current;
                    }
                    else
                        ++it;
                }

                for (auto i = group.begin(); i != group.end(); ++i)
                {
                    bin_map[i->first] = rev_bin-bin_no;
                }
                --n;
                ++bin_no;

                sum -= tmp_sum;
                mean = sum / std::max((max_bins - bin_no),uint32_t(1));
                max_bin_size = 1.1 * mean;
            }
        }
        if (_stats.size() > 0)
        {
            for (auto i = _stats.begin(); i != _stats.end(); ++i)
            {
                bin_map[i->first] = rev_bin-bin_no;
                bin_no++;
                //cout << "rest bin: " << i->second << "\n";
            }
        }
        bin_map[1 << 2 * mm_len] = rev_bin-bin_no;
        std::ofstream bins("bin_map.txt");
        for (int i = 0; i < bin_map.size(); i++) {
            bins << i << " " << bin_map[i] << "\n";
        }
    }
};

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

class KmerCounter {
    const uint8_t K;
    const uint8_t signature_len{5};
    const unsigned int n_bins{512};
    std::array<char,256> codes={4};

    int num_datastore_readers{2};
    int num_kmer_splitters{1};
    int num_bin_writers{6};

    int num_smallK_counters{2};
    std::mutex write_lock;
    std::mutex smallK_merge_lock;
    std::mutex mmer_table_lock;
    using Kmer = std::pair<bool, uint64_t>;
    using KmerArray = std::array<Kmer, 300>;
    std::vector<std::thread> kmer_producers;
    KmerArrayQueue LBkmer_chunks;
    MPMCBoundedQueue<KmerArray> kmer_chunks;

    std::vector<std::thread> kmer_partitioners;
    std::vector<std::thread> kmer_counters;
    std::vector<std::thread> kmer_bin_writers;

    std::vector<std::thread> kmer_dist_generators;

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
    explicit KmerCounter(uint8_t k) : K(k), signature_len(9), mmer_table(signature_len, n_bins), LBkmer_chunks(num_datastore_readers), kmer_chunks(1024) {
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

    void mmer_table_initialiser(const PairedReadsDatastore &r, MinimiserTable &minimiserTable);

    void make_minimiser_send_to_bin(std::string sequence, std::vector<ExtendedKmerBin> &bins,
                                    std::list<Seq2bit> &reads_2bit);

};


#endif //BSG_KMERCOUNTER_HPP
