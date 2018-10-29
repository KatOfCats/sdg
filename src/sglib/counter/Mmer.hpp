//
// Created by Luis Yanes (EI) on 2018-12-13.
//

#ifndef BSG_MMER_HPP
#define BSG_MMER_HPP

#include <cstdint>
#include <algorithm>

class Mmer{

    uint32_t str;
    uint32_t mask;

    // INFO: This pointer appears once per mmer and points to the init'ed one on the normN vars
    uint32_t* norm;         /// Canonical mmer value for each mmer or special value for non-allowed mmers

    // INFO: These appear a single time in the program
    static uint32_t norm5[1 << 10];
    static uint32_t norm6[1 << 12];
    static uint32_t norm7[1 << 14];
    static uint32_t norm8[1 << 16];
    static uint32_t norm9[1 << 18];
    static uint32_t norm10[1 << 20];
    static uint32_t norm11[1 << 22];

    const uint8_t K;    /// Size in nts of the mmer
    uint32_t value;     /// Corresponding canonical mmer value

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

    }
    static _init;

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

    /// Inserts a single nt and populates the mmer value using the table of valid mmers (norm)
    /// \param symb NT to insert
    inline void insert(unsigned char symb)
    {
        str <<= 2;
        str += symb;
        str &= mask;

        value = norm[str];
    }

    uint32_t get(){return value;}

    /// Inserts K nts to the minimiser and populates it's value using the table of valid mmers (norm)
    /// \param seq 2bit sequence of nts with size >= K
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

    /// Copies the values important to evaluate a mmer against the bins
    /// \param c
    void set(Mmer c){
        str = c.str;
        value = c.value;
    }

    /// Defines which minimisers will be used, is used when calculating which mmer goes to which bin
    /// \param mmer mmer to be evaluated
    /// \param len length in nts of this minimiser
    /// \return returns whether the minimiser is valid or not
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
    bool operator<=(const Mmer &o) const {return this->value <= o.value; }
};

#endif //BSG_MMER_HPP
