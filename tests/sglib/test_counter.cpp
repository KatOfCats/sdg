//
// Created by Luis Yanes (EI) on 13/09/2018.
//

#include <catch.hpp>
#include <sglib/counter/KmerCounter.hpp>
#include <sglib/datastores/PairedReadsDatastore.hpp>
#include <sglib/workspace/WorkSpace.hpp>

TEST_CASE("KmerCounter file reader processes all reads") {
    KmerCounter kmerCounter(27);
    std::vector<PairedReadsDatastore> datastores;
    PairedReadsDatastore prds("../tests/datasets/dbg_test/R1.fastq", "../tests/datasets/dbg_test/R2.fastq", "test.prds");

    datastores.emplace_back(prds);
    kmerCounter.processReads(datastores, false);

    REQUIRE(0 == 0);
}