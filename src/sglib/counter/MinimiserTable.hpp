//
// Created by Luis Yanes (EI) on 2018-12-13.
//

#ifndef BSG_MINIMISERTABLE_HPP
#define BSG_MINIMISERTABLE_HPP

#include <cstdint>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <list>
#include "Mmer.hpp"

class MinimiserTable{
    uint32_t mm_len;                        /// Size of the minimiser, the tables will then be of size 4^mmer_len
    const uint32_t num_bins;                /// Number of parts onto which the kmer space will be subdivided
    std::vector<uint32_t> minimiser_table;  /// Stores the number of times a minimiser has been seen.
public:
    std::vector<uint32_t> bin_map;          /// Indexing a minimiser returns the bin it belongs to
    MinimiserTable(uint32_t mm_len, uint32_t num_bins) : mm_len(mm_len), num_bins(num_bins) {
        minimiser_table.resize((1<<(2*mm_len))+1);
        std::fill(minimiser_table.begin(), minimiser_table.end(), 0);
    }

    void combine(const MinimiserTable& o) {
        std::transform(o.minimiser_table.cbegin(), o.minimiser_table.cend(), o.minimiser_table.cbegin(), minimiser_table.begin(), std::plus<uint32_t>());
    }
    std::size_t size() const {return minimiser_table.size();}
    unsigned int get_bin_id(unsigned int n) const {return bin_map[n];}
    uint32_t& operator[](const std::size_t i){return minimiser_table[i];};

    void loadDistribution(){
        std::ifstream bins("bin_map.txt");
        std::string sep, lkp;
        int pos;
        this->bin_map.resize(this->size());
        std::getline(bins, sep);
        for (int i = 0; i < bin_map.size(); i++) {
            bins >> pos >> sep >> lkp;
            bin_map[pos] = (lkp=="DISABLED_SIGNATURE"?0u:std::stoi(lkp));
//            std::cout << pos << "," << this->bin_map[i] << std::endl;
        }
    }
    void calcBins() {
        // How to redistribute this uniformly?
        // Accumulate the sizes and calculate median
        // Order the bins by size, largest to smallest (larger than median bins can't be split)
        // After the bins have been ordered, take all the bins too large to be grouped and keep them "separate"
        // now, for bins left, take the largest and then all the smallest until median bin size has been achieved

        // 1) print the bins largest first

        bin_map.resize(minimiser_table.size());
        std::fill(bin_map.begin(), bin_map.end(), 0);
        auto n_bins = num_bins;
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
        //counting sum
        double sum = 0.0;
        for (auto &i : _stats)
        {
            i.second += 10000;
            sum += i.second;
        }

        double mean = sum / n_bins;
        double max_bin_size = 1.1 * mean;
        uint32_t n = n_bins-1; //one is needed for disabled signatures
        uint32_t max_bins = n_bins - 1;

        while (_stats.size() > n)
        {
            std::pair<uint32_t, uint64_t>& max = _stats.front();

            if (max.second > mean)
            {
                bin_map[max.first] = (bin_no+1);
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
                    bin_map[i->first] = bin_no+1;
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
                bin_map[i->first] = bin_no+1;
                bin_no++;
            }
        }

        bin_map[1 << 2 * mm_len] = bin_no+1; // Store disabled signatures in a single bin

        {
            std::ofstream bins("bin_map.txt");
            char ACGT[10];
            ACGT[mm_len] = '\0';
            bins << "SIGNATURE\tSEQ\tBIN" << std::endl;
            for (int i = 0; i < bin_map.size(); i++) {
                for (int j = mm_len - 1; j >= 0; --j) ACGT[mm_len- j - 1] = "ACGT"[(i >> 2 * j) & 3];
                bins << i << "\t\t" << ACGT << "\t" << ( (bin_map[i]>0)?std::to_string(bin_map[i]):"DISABLED_SIGNATURE") << "\n";
//                std::cout << i << "," << bin_map[i] << std::endl;
            }
        }
    }
};

#endif //BSG_MINIMISERTABLE_HPP
