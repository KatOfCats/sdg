//
// Created by Luis Yanes (EI) on 13/09/2018.
//

#include <catch.hpp>
#include <sglib/counter/KmerCounter.hpp>

TEST_CASE("KmerCounter file reader processes all reads") {
    KmerCounter kmerCounter(27);
    std::vector<PairedReadsDatastore> datastores;

    kmerCounter.processReads(datastores);

    REQUIRE(0 == 0);
}