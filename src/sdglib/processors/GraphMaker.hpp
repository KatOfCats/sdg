//
// Created by Bernardo Clavijo (EI) on 05/07/2018.
//

// This class IS a re-invention of the wheel. It creates a DBG from a number of different sources and puts it into a
// SequenceDistanceGraph.


#ifndef BSG_GRAPHMAKER_HPP
#define BSG_GRAPHMAKER_HPP

#include <sdglib/graph/SequenceDistanceGraph.hpp>
/**
 * This kmer class can give possible neighbours FW and BW to help build a graph more easily.
 *
 */
class GraphMakerKmer{
public:
    GraphMakerKmer(uint64_t _kmervalue, uint8_t _k):kmervalue(_kmervalue),k(_k){};

    uint64_t kmervalue;
    uint8_t k;
};

typedef union {
    struct {
        unsigned int fw_A:1;
        unsigned int fw_C:1;
        unsigned int fw_T:1;
        unsigned int fw_G:1;
        unsigned int bw_A:1;
        unsigned int bw_C:1;
        unsigned int bw_T:1;
        unsigned int bw_G:1;
    };
    struct {
        unsigned int fw:4;
        unsigned int bw:4;
    };
    unsigned int all = 0;
} connectivity;

class GraphMaker {
public:
    GraphMaker(SequenceDistanceGraph & _sg): sg(_sg){};
    void new_graph_from_kmerset_trivial(const std::unordered_set<uint64_t> & kmerset,uint8_t k);
    void new_graph_from_kmerset_trivial128(const std::unordered_set<__uint128_t, int128_hash> & kmerset,uint8_t k);
    void new_graph_from_kmerlist_trivial128(const std::vector<__uint128_t> & kmerset,uint8_t k);
//    //Minimum cleanup options
    void tip_clipping(int tip_size);
    void remove_small_unconnected(int min_size);
private:
    SequenceDistanceGraph & sg;
};


#endif //BSG_GRAPHMAKER_HPP
