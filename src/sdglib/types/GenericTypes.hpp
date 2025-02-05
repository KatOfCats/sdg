//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#ifndef BSG_GENERICTYPES_HPP
#define BSG_GENERICTYPES_HPP

#include <iostream>
#include <string>
#include <tuple>
#include <limits>
#include <cstdint>
#include <xxhash.h>
#include "hashing_helper.hpp"

using sgNodeID_t = int64_t; //first node is 1; negatives are RC
using seqID_t = int32_t; //first sequence is 0;

enum class NodeStatus:uint8_t {Active, Deleted};

enum class SupportType:uint8_t {Undefined,Operation,SequenceDistanceGraph,DistanceGraph,PairedRead,LinkedRead,LinkedTag,LongRead,ReadPath,KmerCoverage};

class Support{
public:
    SupportType type=SupportType::Undefined;
    uint16_t index=0;
    uint64_t id=0;
};

/**
 * The Node contains the sequence of a node and its status {Active, Deleted}
 */
class Node{
public:
    Node(const std::string &_seq, NodeStatus _status) : sequence(_seq), status(_status){};
    Node(const std::string &_seq) : sequence(_seq),status(NodeStatus::Active){};
    Node() = default;
    bool operator==(const Node &o) const {
        return std::tie(status,sequence) == std::tie(o.status,o.sequence);
    }

    bool is_canonical();
    void make_rc();

    friend std::ostream &operator<<(std::ostream &os, const Node &node) {
        if (node.sequence.length() > 20) {
            os << node.sequence.substr(0, 20) << " ... " << node.sequence.substr(node.sequence.length() - 20, 20);
        } else {
            os << node.sequence;
        }
        return os;
    }

    std::string sequence = "";
    NodeStatus status = NodeStatus::Active;
    Support support;
};

/**
 * The Link represents a connection between one end of a sequence and the end of another sequence
 *
 * +AAAAAAAAA- connected to +BBBBBBBBBB-
 *
 * would be represented as the -A,B link. Links also contain a distance parameter in case of overlap the distance is
 * negative and in case of "scaffolding" link the distance would be positive.
 */
class Link{
public:
    Link(){};
    Link( sgNodeID_t _src, sgNodeID_t _dst, int32_t _dist, Support _support = {}) : source(_src), dest(_dst), dist(_dist), support(_support) {};

    bool operator==( const  Link) const;
    bool operator<(const Link)const;

    friend std::ostream &operator<<(std::ostream &os, const Link &link) {
        os << link.source << " -> " << link.dest;
        return os;
    }
    sgNodeID_t source = 0;
    sgNodeID_t dest = 0;
    int32_t dist = 0;
    Support support;
};

struct link_hash{
    size_t operator()(const Link& l) const {
        std::tuple<sgNodeID_t , sgNodeID_t , int32_t > tp (l.source,l.dest,l.dist);
        sdglib::hash<std::tuple<sgNodeID_t , sgNodeID_t , int32_t>> h;
        return h (tp);
    }
};


struct int128_hash {
    size_t operator()( const __int128 &x) const
    {
        const void *buffer = (unsigned char *) &x;
        uint64_t tmp_hash = XXH64(buffer, 16, 0);
        return tmp_hash;
    }
};

#endif //BSG_GENERICTYPES_HPP
