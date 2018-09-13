//
// Created by Luis Yanes (EI) on 12/09/2018.
//

#include <sglib/processors/KmerCompressionIndex.hpp>
#include "KmerCounter.hpp"

void KmerCounter::processReads(std::vector<PairedReadsDatastore> &read_files) {

    for (const auto &r: read_files) {
        for (int i = 0; i < 2; i++) {
            std::thread th(&KmerCounter::read_file, r);
            file_readers.push_back(std::move(th));
        }
    }
}

void KmerCounter::processBins() {

}

std::function<void()> KmerCounter::read_file(PairedReadsDatastore &ds) {
    BufferedPairedSequenceGetter bpsg(ds,100000,1000);
    const size_t local_kmers_size = 2000000;
    std::vector<Kmer> readkmers(local_kmers_size);
    StringKMerFactory skf(K);

    for (uint64_t rid = 1; rid <= ds.size(); ++rid) {
        readkmers.clear();
        skf.create_kmers(bpsg.get_read_sequence(rid), readkmers);
        file_parts.enqueue(readkmers);
    }

    return 0;
}
