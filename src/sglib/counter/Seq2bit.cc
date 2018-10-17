//
// Created by Luis Yanes (EI) on 17/10/2018.
//

#include <sglib/counter/Seq2bit.hpp>

void Seq2bit::insert_nt(unsigned char nt) {m_data.emplace_back(codes[nt]);}

std::size_t Seq2bit::size() {return m_data.size();}

Seq2bit::Seq2bit(std::string seq) {
    for (const auto &c: seq){
        m_data.emplace_back(codes[c]);
    }
}

unsigned char &Seq2bit::operator[](int i) { return m_data[i]; }

unsigned char Seq2bit::operator[](int i) const { return m_data[i]; }

unsigned char *Seq2bit::data() {return m_data.data();}
