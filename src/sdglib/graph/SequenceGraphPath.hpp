//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#ifndef BSG_SEQUENCEGRAPHPATH_HPP
#define BSG_SEQUENCEGRAPHPATH_HPP

#include <string>
#include <vector>
#include <set>
#include <sdglib/types/GenericTypes.hpp>

class SequenceSubGraph;
class SequenceDistanceGraph;
class SequenceGraphPath {
public:
    explicit SequenceGraphPath(SequenceDistanceGraph & _sg, const std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    explicit SequenceGraphPath(const SequenceDistanceGraph &_sg, const std::vector<sgNodeID_t> _nodes={}) : sg(_sg), nodes(_nodes) {};

    SequenceGraphPath(const SequenceGraphPath& sgp) : nodes(sgp.nodes), sg(sgp.sg) {};

    SequenceGraphPath& operator=(const SequenceGraphPath &other);

    std::string get_fasta_header(bool use_oldnames = false) const;
    std::string get_sequence() const;
    size_t get_sequence_size_fast() const;
    std::vector<Link> get_next_links();
    void reverse();
    bool is_canonical();
    std::set<sgNodeID_t> make_set_of_nodes() const;
    bool operator==(const SequenceGraphPath& rhs) const;
    bool operator<(const SequenceGraphPath& rhs) const;

    bool extend_if_coherent(SequenceGraphPath s){};
    void clear() {
        nodes.clear();
    };

    bool is_unitig();
    std::vector<sgNodeID_t >& getNodes() {return nodes;}
    const std::vector<sgNodeID_t >& getNodes() const {return nodes;}

    std::vector<sgNodeID_t> nodes;

private:
    const SequenceDistanceGraph& sg;
};


#endif //BSG_SEQUENCEGRAPHPATH_HPP
