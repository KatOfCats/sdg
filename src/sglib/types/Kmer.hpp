//
// Created by Luis Yanes (EI) on 06/06/2018.
//

#ifndef BSG_KMER_HPP
#define BSG_KMER_HPP

#include <cstdint>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <sglib/hash/xxhash.h>

class revcomp_DNA2 {
public:
    const unsigned char & operator() (int i) const { return lut_array[i]; }
private:
    unsigned char lut_array[4] {3,2,1,0};
};

class revcomp_2DNA2 {
    // TODO: Make lookup table to optimise and remove comments
    /*
    static const unsigned char lut_array[256] {
            0x0F, 0x0B, 0x07, 0x03, 0x0E
    };
    */
    const unsigned char rev(const unsigned char &x) const {
        unsigned char b=x;
        b = (b & 0x33) <<  2 | (b & 0xCC) >>  2;
        b = (b & 0x0F) <<  4 | (b & 0xF0) >>  4;
        b = ~b;
        return b;
    }
public:
    const unsigned char /*&*/ operator() (unsigned char i) const {  return rev(i); /*return lut_array[i];*/ }
};

// number of bytes for k-mer
constexpr int GET_KMER_BYTES(int k, int bpn) { return k / (bpn * 2) +1; }

// k-mer type, use of 32 Bit structures if b <= 4, else 64 Bit structures
constexpr int GET_KMER_TYPE(int b) { return b <= 4 ? b <= 2 ? 2 : 4 : 8; }

// number of necessary t-bytes-structures to store a single k-mer
constexpr int GET_KMER_NELEMENTS(int b, int t) { return (b - 1) / t ; }




// This is a class to hold the bits of a nucleotide
// This class can represent sequences of any size and any alphabet
// Kmers are represented in bits read from left to right

// TODO: Check for unnecessary copies when creating/updating KMER objects, and remove them!
// FIXME: Check speed of "append" and "prepend" methods for fwd and rc creation vs fwd + RC(fwd)
// FIXME: Add a hashing function for all NTBits size types

template<
        int K,
        int B = GET_KMER_BYTES(K, 2),
        int T = GET_KMER_TYPE(B),
        int NELEMENTS = GET_KMER_NELEMENTS(B, T)>
class KMER {};

template<int K, int B, int NELEMENTS>
class KMER<K, B, 8, NELEMENTS> {
    static constexpr unsigned int BPN = 2;
    using T=uint64_t;
    static constexpr unsigned int MOD = BPN*K % std::numeric_limits<T>::digits;
    static constexpr unsigned int LBP1 = ( (MOD == 0) ? NELEMENTS-1 : NELEMENTS);
    static constexpr unsigned int LBP2 = ( (MOD == 0) ? std::numeric_limits<T>::digits-BPN : MOD - BPN);
    std::array<T,NELEMENTS+1> subparts;
public:
    KMER() {
        std::memset(&subparts, 0, (NELEMENTS+1)*sizeof(T));
    }
    T get() const { return this->subparts[0]; }

    inline void set(T val) {
        this->subparts[0] = val;
        for (int i = 1; i < NELEMENTS; i++) {
            this->subparts[i] = 0;
        }
    }

    inline void set(const KMER &o) {
        for (int i = 0; i < NELEMENTS+1; i++) {
            this->subparts[i] = o.subparts[i];
        }
    }

    static inline const size_t getSize() {
        return static_cast<const size_t>(std::numeric_limits<T>::digits * NELEMENTS);
    }

    inline KMER operator|(const KMER &o) const {
        KMER result;
        for (int i = 0; i < NELEMENTS; i++)
            result.subparts[i] = this->subparts[i] | o.subparts[i];
        return result;
    }

    inline KMER &operator|=(const KMER &o) {
        for (int i = 0; i < NELEMENTS; i++) { this->subparts[i] |= o.subparts[i]; }
        return *this;
    }

    inline KMER &operator|=(const int &o) {
        subparts[0] |= o;
        return *this;
    }

    inline KMER operator&(const int &o) const {
        KMER result;
        result.subparts[0] = subparts[0] & o;
        return result;
    }

    inline KMER operator&(const KMER &o) const {
        KMER result;
        for (int i = 0; i < NELEMENTS; i++)
            result.subparts[i] = this->subparts[i] & o.subparts[i];
        return result;
    }

    inline KMER &operator&=(const KMER &o) {
        for (int i = 0; i < NELEMENTS; i++) { this->subparts[i] &= o.subparts[i]; }
        return *this;
    }

    inline KMER operator~() const {
        KMER result;
        for (int i = 0; i < NELEMENTS; i++)
            result.subparts[i] = ~this->subparts[i];
        return result;
    }

    inline KMER operator<<(const int &dist_to_push) const {
        KMER result;
        result.set(0);

        int large_shift = dist_to_push / std::numeric_limits<T>::digits;
        int small_shift = dist_to_push % std::numeric_limits<T>::digits;

        for (int i = large_shift; i < NELEMENTS; i++) {
            result.subparts[i] = result.subparts[i] | (this->subparts[i - large_shift] << small_shift);

            result.subparts[i + 1] =
                    this->subparts[i - large_shift] >> (std::numeric_limits<T>::digits - small_shift);

        }
        result.subparts[NELEMENTS ] =
                result.subparts[NELEMENTS ] | (this->subparts[NELEMENTS - large_shift] << small_shift);
        return result;
    }

    inline KMER &operator<<=(const int &dist_to_push) {
        if (IsZero() || dist_to_push == 0) return *this;

        *(this) = (*this) << dist_to_push;
        return *this;
    }

    inline KMER operator>>(const int &dist_to_push) const {
        KMER result;
        result.set(0);

        int large_shift (dist_to_push / std::numeric_limits<T>::digits);
        int small_shift (dist_to_push % std::numeric_limits<T>::digits);
        int csmall_shift (std::numeric_limits<T>::digits - small_shift);
        if (small_shift == 0 and large_shift == 0) {
            for (int i = 0; i < NELEMENTS - large_shift; i++) {
                result.subparts[i] = result.subparts[i+large_shift];
            }
            return result;
        }
        result.subparts[NELEMENTS-large_shift] = this->subparts[NELEMENTS-large_shift] >> small_shift;
        for (int i = NELEMENTS-large_shift-1; i >=0; i--){
            result.subparts[i] = (this->subparts[i] >> small_shift) | (this->subparts[i+1] << csmall_shift);
        }
        return result;
    }

    inline KMER &operator>>=(const int &dist_to_push) {
        if (IsZero() || dist_to_push == 0) return *this;

        *(this) = (*this) >> dist_to_push;
        return *this;
    }

    inline bool operator!=(const KMER &o) const {
        for (int i = 0; i < NELEMENTS; i++)
            if (this->subparts[i] != o.subparts[i])
                return true;
        return false;
    }

    inline bool operator==(const KMER &o) const {
        for (int i = 0; i < NELEMENTS; i++)
            if (this->subparts[i] != o.subparts[i])
                return false;
        return true;
    }

    inline bool operator<(const KMER &o) const {
        for (int i = NELEMENTS - 1; i >= 0; --i)
            if (this->subparts[i] != o.subparts[i])
                return this->subparts[i] < o.subparts[i];

        return false;
    }

    inline bool operator<=(const KMER &o) const {
        return operator==(o) || operator<(o);
    }

    inline u_int8_t operator[](size_t idx) const {
        return static_cast<u_int8_t>((this->subparts[idx / (std::numeric_limits<T>::digits / BPN)]
                >> (BPN * (idx % (std::numeric_limits<T>::digits / BPN)))) & 3);
    }

    inline void add(uint8_t c) {
        *this <<= BPN;
        *this |= c;
    }

    inline void prepend(unsigned char c) {
        *this >>= BPN;
        this->subparts[LBP1] |= (static_cast<T>(c) << LBP2);
    }

    inline void append(unsigned char c) {
        add(c);
    }

    inline KMER& RC() {
        *(this) = revcomp(*this);
        return *this;
    }

    inline KMER revcomp(const KMER &x) {
        KMER res;
        res.set(0);
        auto *kmerrev = (unsigned char *) (&(res.subparts[0]));
        auto *kmer = (unsigned char *) (&(x.subparts[0]));
        const auto FORSIZE( std::ceil((float)K*BPN/std::numeric_limits<unsigned char>::digits));
        for (size_t i = 0; i < FORSIZE; ++i) {
            const int POS(FORSIZE - 1 - i);
            kmerrev[POS] = revcomp_2DNA2()(kmer[i]);
        }
        const int DIGITS(std::numeric_limits<T>::digits);
        if (K==32) return res;
        // Cleanup if required.
        constexpr auto NTRest(8-((K*BPN%8) ? K*BPN%8 : 8));
        res >>= NTRest;
        return res;
    }

    inline bool IsZero() const {
        for (int i = NELEMENTS; i >= 0; i--) {
            if (subparts[i] != 0) return false;
        }
        return true;
    }

    inline std::string toString() const {
        std::string result;
        KMER cp;
        cp.set(*this);
        result.resize(K);
        T tmp;
        for (int i = K-1; i >= 0; i--) {
            auto AP(i*BPN/std::numeric_limits<T>::digits);
            auto BP(i*BPN%std::numeric_limits<T>::digits);
            tmp = (cp.subparts[AP] >> BP);
            tmp &= 3;
            result[K - 1 - i] = "ACGT"[tmp];
        }
        return result;
    }

    inline uint64_t hash() const {
        return subparts[0];
    }
};

template<int K, int B>
class KMER<K, B, 8, 0> {
    static const unsigned int BPN = 2;
    using T=uint64_t;
    static constexpr unsigned int LBP = BPN*K % std::numeric_limits<T>::digits - BPN;
    T subparts = 0;
public:
    T get() const { return this->subparts; }

    inline void set(T val) {
        this->subparts = val;
    }

    inline void set(const KMER &o) {
        this->subparts = o.subparts;
    }

    static inline const size_t getSize() { return static_cast<const size_t>(std::numeric_limits<T>::digits); }

    inline KMER operator|(const KMER &o) const {
        KMER result;
        result.subparts = this->subparts | o.subparts;
        return result;
    }

    inline KMER &operator|=(const KMER &o) {
        this->subparts |= o.subparts;
        return *this;
    }

    inline KMER operator&(const KMER &o) const {
        KMER result;
        result.subparts = this->subparts & o.subparts;
        return result;
    }

    inline KMER &operator&=(const KMER &o) {
        this->subparts &= o.subparts;
        return *this;
    }

    inline KMER operator~() const {
        KMER result;
        result.subparts = ~this->subparts;
        return result;
    }

    inline KMER operator<<(const int &dist_to_push) const {
        KMER result;
        result.subparts = result.subparts | (this->subparts << dist_to_push);
        return result;
    }

    inline KMER &operator<<=(const int &dist_to_push) {
        *(this) = (*this) << dist_to_push;
        return *this;
    }

    inline KMER operator>>(const int &dist_to_push) const {
        KMER result;
        result.set(0);

        result.subparts = this->subparts >> dist_to_push;
        return result;
    }

    inline KMER &operator>>=(const int &dist_to_push) {
        *(this) = (*this) >> dist_to_push;
        return *this;
    }

    inline bool operator!=(const KMER &o) const {
        if (this->subparts != o.subparts)
            return true;
        return false;
    }

    inline bool operator==(const KMER &o) const {
        if (this->subparts != o.subparts)
            return false;
        return true;
    }

    inline bool operator<(const KMER &o) const {
        return this->subparts < o.subparts;
    }

    inline bool operator<=(const KMER &o) const {
        return this->subparts <= o.subparts;
    }

    inline u_int8_t operator[](size_t idx) const {
        return static_cast<u_int8_t>((this->subparts >> (BPN * (idx % (std::numeric_limits<T>::digits / BPN)))) &
                                     3);
    }

    inline void add(unsigned char c) {
        subparts <<= BPN;
        subparts |= c;
    }

    inline void prepend(unsigned char c) {
        subparts >>= BPN;
        subparts |= static_cast<T>(c) << LBP;
    }

    inline void append(unsigned char c) {
        add(c);
    }

    inline KMER& RC() {
        *this = revcomp(*this);
        return *this;
    }

    inline KMER revcomp(const KMER &o) {
        KMER res;
        res.set(0);
        T x = o.subparts;

        x = (x & 0x3333333333333333) << 2  | (x & 0xCCCCCCCCCCCCCCCC) >> 2;
        x = (x & 0x0F0F0F0F0F0F0F0F) << 4  | (x & 0xF0F0F0F0F0F0F0F0) >> 4;
        x = (x & 0x00FF00FF00FF00FF) << 8  | (x & 0xFF00FF00FF00FF00) >> 8;
        x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16;
        x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32;
        res.subparts = ~x;
        res.subparts = res.subparts >> BPN * ((std::numeric_limits<T>::digits / BPN) - K);
        return res;
        //return (res >> (BPN * ((std::numeric_limits<T>::digits / BPN) - K)));
    }

    inline bool IsZero() const {
        if (this->subparts != 0)
            return false;
        return true;
    }

    inline std::string toString() const {
        std::string result;
        T cp(subparts);
        result.resize(K);
        for (int i = K - 1; i >= 0; i--) {
            result[K - 1 - i] = "ACGT"[cp >> (i * BPN) & 3];
        }
        return result;
    }

    inline uint64_t hash() const {
        return subparts;
    }

};

template<int K, int B>
class KMER<K, B, 4, 0> {
    static const unsigned int BPN = 2;
    using T=uint32_t;
    static constexpr unsigned int LBP = BPN*K % std::numeric_limits<T>::digits - BPN;
    T subparts = 0;
public:
    T get() const { return this->subparts; }

    inline void set(T val) {
        this->subparts = val;
    }

    inline void set(const KMER &o) {
        this->subparts = o.subparts;
    }

    static inline const size_t getSize() { return static_cast<const size_t>(std::numeric_limits<T>::digits); }

    inline KMER operator|(const KMER &o) const {
        KMER result;
        result.subparts = this->subparts | o.subparts;
        return result;
    }

    inline KMER &operator|=(const KMER &o) {
        this->subparts |= o.subparts;
        return *this;
    }

    inline KMER operator&(const KMER &o) const {
        KMER result;
        result.subparts = this->subparts & o.subparts;
        return result;
    }

    inline KMER &operator&=(const KMER &o) {
        this->subparts &= o.subparts;
        return *this;
    }

    inline KMER operator~() const {
        KMER result;
        result.subparts = ~this->subparts;
        return result;
    }

    inline KMER operator<<(const int &dist_to_push) const {
        KMER result;
        result.subparts = result.subparts | (this->subparts << dist_to_push);
        return result;
    }

    inline KMER &operator<<=(const int &dist_to_push) {
        *(this) = (*this) << dist_to_push;
        return *this;
    }

    inline KMER operator>>(const int &dist_to_push) const {
        KMER result;
        result.set(0);

        result.subparts = this->subparts >> dist_to_push;
        return result;
    }

    inline KMER &operator>>=(const int &dist_to_push) {
        *(this) = (*this) >> dist_to_push;
        return *this;
    }

    inline bool operator!=(const KMER &o) const {
        if (this->subparts != o.subparts)
            return true;
        return false;
    }

    inline bool operator==(const KMER &o) const {
        if (this->subparts != o.subparts)
            return false;
        return true;
    }

    inline bool operator<(const KMER &o) const {
        return this->subparts < o.subparts;
    }

    inline bool operator<=(const KMER &o) const {
        return this->subparts <= o.subparts;
    }

    inline u_int8_t operator[](size_t idx) const {
        return static_cast<u_int8_t>((this->subparts >> (BPN * (idx % (std::numeric_limits<T>::digits / BPN)))) &
                                     3);
    }

    inline void add(unsigned char c) {
        subparts <<= BPN;
        subparts |= c;
    }

    inline void prepend(unsigned char c) {
        subparts >>= BPN;
        subparts |= static_cast<T>(c) << LBP;

    }

    inline void append(unsigned char c) {
        add(c);
    }


    inline KMER& RC() {
        *this = revcomp(*this);
        return *this;
    }

    inline KMER revcomp(const KMER &o) {
        KMER res;
        res.set(0);
        T x = o.subparts;

        x = (x & 0x33333333) << 2 | (x & 0xCCCCCCCC) >> 2;
        x = (x & 0x0F0F0F0F) << 4 | (x & 0xF0F0F0F0) >> 4;
        x = (x & 0x00FF00FF) << 8 | (x & 0xFF00FF00) >> 8;
        x = (x & 0x0000FFFF) << 16 | (x & 0xFFFF0000) >> 16;
        res.subparts = ~x;
        res.subparts = res.subparts >> BPN * ((std::numeric_limits<T>::digits / BPN) - K);
        return res;
    }

    inline bool IsZero() const {
        if (subparts != 0)
            return false;
        return true;
    }

    inline std::string toString() const {
        std::string result;
        T cp(subparts);
        result.resize(K);
        for (int i = K - 1; i >= 0; i--) {
            result[K - 1 - i] = "ACGT"[cp >> (i * BPN) & 3];
        }
        return result;
    }

    inline uint64_t hash() const {
        return subparts;
    }

};

#endif //BSG_KMER_HPP
