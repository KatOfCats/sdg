//
// Created by Luis Yanes (EI) on 17/10/2018.
//

#include <sglib/counter/Seq2bit.hpp>
#include <algorithm>

std::size_t Seq2bit::size() {return m_data.size();}

Seq2bit::Seq2bit(const std::string &seq) {
    std::transform(seq.cbegin(), seq.cend(), std::back_inserter(m_data), [](const char c){return codes[c];});
}

unsigned char &Seq2bit::operator[](int i) { return m_data[i]; }

unsigned char Seq2bit::operator[](int i) const { return m_data[i]; }

unsigned char *Seq2bit::data() {return m_data.data();}
