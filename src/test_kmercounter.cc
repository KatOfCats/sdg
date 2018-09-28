#include <iostream>
#include <fstream>
#include <sglib/counter/KmerCounter.hpp>
#include "sglib/logger/OutputLog.h"
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    KmerCounter kmerCounter(21);
    std::vector<PairedReadsDatastore> datastores;
    //PE_R1_1M
//    PairedReadsDatastore prds("../tests/datasets/dbg_test/R1.fastq", "../tests/datasets/dbg_test/R2.fastq", "test.prds");
    PairedReadsDatastore prds("/Users/yanesl/rundir/yellow_rust/PE_R1_1M.fastq", "/Users/yanesl/rundir/yellow_rust/PE_R2_1M.fastq", "test.prds", 0, 310);

    datastores.emplace_back(prds);
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    kmerCounter.processReads(datastores, false);
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double > time_span = std::chrono::duration_cast<std::chrono::duration<double >>(t2 - t1);

    std::cout <<"Total time elapsed= " << time_span.count() << std::endl;
    return 0;
}

