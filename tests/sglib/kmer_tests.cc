//
// Created by Luis Yanes (EI) on 26/10/2018.
//

#include "catch.hpp"

#include <string>
#include <sglib/types/Kmer.hpp>
#include <math.h>
using namespace std;

TEST_CASE("Array based kmer", "[KMER][KMERVEC]") {
    const int K=32;
    unsigned char b2f[255]{4};b2f['a'] = b2f['A'] = 0;b2f['c'] = b2f['C'] = 1;b2f['g'] = b2f['G'] = 2;b2f['t'] = b2f['T'] = 3;
    unsigned char b2r[255]{4};b2r['a'] = b2r['A'] = 3;b2r['c'] = b2r['C'] = 2;b2r['g'] = b2r['G'] = 1;b2r['t'] = b2r['T'] = 0;

    SECTION("Validate create from string") {
        const string stringKmer   ("GAAAGGGGTTTTTTTTTTTTTTTTCCCCGGTG");
        NTBits<K> kmer;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
        }
        const string kmerString2(kmer.toString());
        REQUIRE(stringKmer == kmerString2);
    }

    SECTION("Validate reverse complement") {
        const string stringKmer   ("GAAAGGGGTTTTTTTTTTTTTTTTCCCCGGTG");
        NTBits<K> kmer;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
        }
        const string kmerString2RC("CACCGGGGAAAAAAAAAAAAAAAACCCCTTTC");
        kmer = kmer.RC();
        string rc = kmer.toString();
        REQUIRE(kmerString2RC == rc);
    }

    SECTION("Validate != comparison") {
        const string stringKmer   ("GAAAGGGGTTTTTTTTTTTTTTTTCCCCGGTG");
        NTBits<K> kmer, kmer2;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
        }
        REQUIRE(kmer != kmer2);
    }

    SECTION("Validate < comparison") {
        const string stringKmer   ("GAAAGGGGTTTTTTTTTTTTTTTTCCCCGGTG");
        const string stringKmer2  ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        NTBits<K> kmer, kmer2;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
            kmer2.add(b2f[stringKmer2[i]]);
        }
        bool isLess (kmer2 < kmer);
        REQUIRE(true == isLess);
    }

}

TEST_CASE("uint64 kmer", "[KMER][UINT64 KMER]") {
    unsigned char b2f[255]{4};b2f['a'] = b2f['A'] = 0;b2f['c'] = b2f['C'] = 1;b2f['g'] = b2f['G'] = 2;b2f['t'] = b2f['T'] = 3;
    unsigned char b2r[255]{4};b2r['a'] = b2r['A'] = 3;b2r['c'] = b2r['C'] = 2;b2r['g'] = b2r['G'] = 1;b2r['t'] = b2r['T'] = 0;

    SECTION("Validate create from string") {
        const string stringKmer   ("GAAAGGGGTTTTTTTTTTT");
        const int K=19;
        NTBits<K> kmer;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
        }
        const string kmerString2(kmer.toString());
        REQUIRE(stringKmer == kmerString2);
    }

    SECTION("Validate reverse complement") {
        const string stringKmer   ("TTTTTTTTTTTCCCCGGTG");
        const string kmerString2RC("CACCGGGGAAAAAAAAAAA");
        const int K=19;
        NTBits<K> kmer;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
        }
        kmer = kmer.RC();
        string rc = kmer.toString();
        REQUIRE(kmerString2RC == rc);
    }

    SECTION("Validate != comparison") {
        const string stringKmer   ("GAAAGGGGTTTTTTTTTTT");
        const int K=14;
        NTBits<K> kmer, kmer2;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
        }
        REQUIRE(kmer != kmer2);
    }

    SECTION("Validate < comparison") {
        const string stringKmer   ("CACCGGGGAAAAAAAAAAA");
        const string stringKmer2  ("AAAAAAAAAAAAAAAAAAA");
        const int K=19;
        NTBits<K> kmer, kmer2;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
            kmer2.add(b2f[stringKmer2[i]]);
        }
        REQUIRE(kmer2 < kmer);
    }

}

TEST_CASE("uint32 kmer", "[KMER][UINT32 KMER]") {
    unsigned char b2f[255]{4};b2f['a'] = b2f['A'] = 0;b2f['c'] = b2f['C'] = 1;b2f['g'] = b2f['G'] = 2;b2f['t'] = b2f['T'] = 3;
    unsigned char b2r[255]{4};b2r['a'] = b2r['A'] = 3;b2r['c'] = b2r['C'] = 2;b2r['g'] = b2r['G'] = 1;b2r['t'] = b2r['T'] = 0;

    SECTION("Validate create from string") {
        const string stringKmer   ("GAAAGGGGTTTTTT");
        const int K=14;
        NTBits<K> kmer;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
        }
        const string kmerString2(kmer.toString());
        REQUIRE(stringKmer == kmerString2);
    }

    SECTION("Validate reverse complement") {
        const string stringKmer   ("TTTTTTCCCCGGTG");
        const int K=14;
        NTBits<K> kmer;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
        }
        const string kmerString2RC("CACCGGGGAAAAAA");
        kmer = kmer.RC();
        string rc = kmer.toString();
        REQUIRE(kmerString2RC == rc);
    }

    SECTION("Validate != comparison") {
        const string stringKmer   ("TTTTTTCCCCGGTG");
        const int K=14;
        NTBits<K> kmer, kmer2;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
        }
        REQUIRE(kmer != kmer2);
    }

    SECTION("Validate < comparison") {
        const string stringKmer   ("TTTTTTCCCCGGTG");
        const string stringKmer2   ("AAAAAAAAAAAAAA");
        const int K=14;
        NTBits<K> kmer, kmer2;
        for (int i = 0; i < K; i++) {
            kmer.add(b2f[stringKmer[i]]);
            kmer2.add(b2f[stringKmer2[i]]);
        }
        REQUIRE(kmer2 < kmer);
    }

}

TEST_CASE("uint32 kmer,rc", "[KMER][UINT32 KMER]") {
    unsigned char b2f[255]{4};b2f['a'] = b2f['A'] = 0;b2f['c'] = b2f['C'] = 1;b2f['g'] = b2f['G'] = 2;b2f['t'] = b2f['T'] = 3;
    unsigned char b2r[255]{4};b2r['a'] = b2r['A'] = 3;b2r['c'] = b2r['C'] = 2;b2r['g'] = b2r['G'] = 1;b2r['t'] = b2r['T'] = 0;
    const int K=14;
    NTBits<K> km,rc;
    const string stringKmer   ("TTTTTTCCCCGGTG");
    const string kmerString2RC("CACCGGGGAAAAAA");
    for (int i = 0; i < K; i++) {
        km.append(b2f[stringKmer[i]]);
        rc.prepend(b2r[stringKmer[i]]);
    }

    SECTION("Validate create from string") {
        REQUIRE(stringKmer == km.toString());
    }
    SECTION("Validate create RC from string") {
        REQUIRE(kmerString2RC == rc.toString());
    }
}

TEST_CASE("uint64 kmer,rc", "[KMER][UINT64 KMER]"){
    unsigned char b2f[255]{4};b2f['a'] = b2f['A'] = 0;b2f['c'] = b2f['C'] = 1;b2f['g'] = b2f['G'] = 2;b2f['t'] = b2f['T'] = 3;
    unsigned char b2r[255]{4};b2r['a'] = b2r['A'] = 3;b2r['c'] = b2r['C'] = 2;b2r['g'] = b2r['G'] = 1;b2r['t'] = b2r['T'] = 0;
    const string stringKmer   ("TTTTTTTTTTTCCCCGGTG");
    const string kmerString2RC("CACCGGGGAAAAAAAAAAA");
    const int K=19;
    NTBits<K> km,rc;
    for (int i = 0; i < K; i++) {
        km.append(b2f[stringKmer[i]]);
        rc.prepend(b2r[stringKmer[i]]);
    }

    SECTION("Validate create from string") {
        REQUIRE(stringKmer == km.toString());
    }
    SECTION("Validate create RC from string") {
        REQUIRE(kmerString2RC == rc.toString());
    }

}

TEST_CASE("36mer,rc", "[KMER][ARRAY KMER]") {
    unsigned char b2f[255]{4};b2f['a'] = b2f['A'] = 0;b2f['c'] = b2f['C'] = 1;b2f['g'] = b2f['G'] = 2;b2f['t'] = b2f['T'] = 3;
    unsigned char b2r[255]{4};b2r['a'] = b2r['A'] = 3;b2r['c'] = b2r['C'] = 2;b2r['g'] = b2r['G'] = 1;b2r['t'] = b2r['T'] = 0;
    const int K=36;
    NTBits<K> km,rc;
    const string stringKmer   ("GGAAAGGGGTTTTTTTTTTTTTTTTCCCCGGTGAAA");
    const string kmerString2RC("TTTCACCGGGGAAAAAAAAAAAAAAAACCCCTTTCC");
    for (int i = 0; i < K; i++) {
        km.append(b2f[stringKmer[i]]);
        rc.prepend(b2r[stringKmer[i]]);
    }

    SECTION("Validate create from string") {
        std::string result(km.toString());
        REQUIRE(stringKmer == result);
    }

    SECTION("Validate create RC from string") {
        std::string result (rc.toString());
        REQUIRE(kmerString2RC == result);
    }

    SECTION("Create RC from kmer") {
        km.RC();
        std::string result(km.toString());
        REQUIRE(kmerString2RC == result);
    }
}

TEST_CASE("35mer,rc", "[KMER][ARRAY KMER]") {
    unsigned char b2f[255]{4};b2f['a'] = b2f['A'] = 0;b2f['c'] = b2f['C'] = 1;b2f['g'] = b2f['G'] = 2;b2f['t'] = b2f['T'] = 3;
    unsigned char b2r[255]{4};b2r['a'] = b2r['A'] = 3;b2r['c'] = b2r['C'] = 2;b2r['g'] = b2r['G'] = 1;b2r['t'] = b2r['T'] = 0;
    const int K=35;
    NTBits<K> km,rc;
    const string stringKmer   ("GGAAAGGGGTTTTTTTTTTTTTTTTCCCCGGTGAA");
    const string kmerString2RC("TTCACCGGGGAAAAAAAAAAAAAAAACCCCTTTCC");
    for (int i = 0; i < K; i++) {
        km.append(b2f[stringKmer[i]]);
        rc.prepend(b2r[stringKmer[i]]);
    }

    SECTION("Validate create from string") {
        std::string result(km.toString());
        REQUIRE(stringKmer == result);
    }

    SECTION("Validate create RC from string") {
        std::string result (rc.toString());
        REQUIRE(kmerString2RC == result);
    }

    SECTION("Create RC from kmer") {
        km.RC();
        std::string result(km.toString());
        REQUIRE(kmerString2RC == result);
    }
}

TEST_CASE("63mer,rc", "[KMER][ARRAY KMER]") {
    unsigned char b2f[255]{4};b2f['a'] = b2f['A'] = 0;b2f['c'] = b2f['C'] = 1;b2f['g'] = b2f['G'] = 2;b2f['t'] = b2f['T'] = 3;
    unsigned char b2r[255]{4};b2r['a'] = b2r['A'] = 3;b2r['c'] = b2r['C'] = 2;b2r['g'] = b2r['G'] = 1;b2r['t'] = b2r['T'] = 0;
    const int K=63;
    NTBits<K> km,rc;
    const string kmerStringFW("CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC");
    const string kmerStringRC("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG");
    for (int i = 0; i < K; i++) {
        km.append(b2f[kmerStringFW[i]]);
        rc.prepend(b2r[kmerStringFW[i]]);
    }

    SECTION("Validate create from string") {
        std::string result(km.toString());
        REQUIRE(kmerStringFW == result);
    }

    SECTION("Validate create RC from string") {
        std::string result (rc.toString());
        REQUIRE(kmerStringRC == result);
    }

    SECTION("Create RC from kmer") {
        km.RC();
        std::string result(km.toString());
        REQUIRE(kmerStringRC == result);
    }
}

TEST_CASE("64mer,rc", "[KMER][ARRAY KMER]") {
    unsigned char b2f[255]{4};b2f['a'] = b2f['A'] = 0;b2f['c'] = b2f['C'] = 1;b2f['g'] = b2f['G'] = 2;b2f['t'] = b2f['T'] = 3;
    unsigned char b2r[255]{4};b2r['a'] = b2r['A'] = 3;b2r['c'] = b2r['C'] = 2;b2r['g'] = b2r['G'] = 1;b2r['t'] = b2r['T'] = 0;
    const int K=64;
    NTBits<K> km,rc;
    const string kmerStringFW("CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG");
    const string kmerStringRC("CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG");
    for (int i = 0; i < K; i++) {
        km.append(b2f[kmerStringFW[i]]);
        rc.prepend(b2r[kmerStringFW[i]]);
    }

    SECTION("Validate create from string") {
        std::string result(km.toString());
        REQUIRE(kmerStringFW == result);
    }

    SECTION("Validate create RC from string") {
        std::string result (rc.toString());
        REQUIRE(kmerStringRC == result);
    }

    SECTION("Create RC from kmer") {
        km.RC();
        std::string result(km.toString());
        REQUIRE(kmerStringRC == result);
    }
}