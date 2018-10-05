//
// Created by Luis Yanes (EI) on 23/08/2018.
//

#include <catch.hpp>
#include <sglib/readers/FileReader.h>
#include <sglib/workspace/WorkSpace.hpp>
#include <random>

TEST_CASE("Workspace create, read, write") {
    REQUIRE(0==0);
}

TEST_CASE("Long reads datastore create, read, write") {
/*    {
        std::string lr_filepath("../tests/datasets/workspace/strwb_lr.fastq");
        std::string lrds_output_path("long_reads.loseq");
        LongReadsDatastore ds(lr_filepath, lrds_output_path);
    }
*/
    LongReadsDatastore ds("long_reads.loseq");
    BufferedSequenceGetter bufferedSequenceGetter(ds);

    //random number engine
    std::mt19937 gen(10); // Always using same seed to get same results
    std::uniform_int_distribution<> dis(1, ds.size());

    std::array<uint64_t, 300> reads_to_check;
    for (auto &r:reads_to_check){
        r = dis(gen);
    }
    const char* seq_ptr;
    std::unordered_map<uint32_t, size_t> read_sz1;
    StringKMerFactory skf(15);
    std::vector<std::pair<bool, uint64_t>> read_kmers;
    for (const auto &r: reads_to_check){
        read_kmers.clear();
        seq_ptr = bufferedSequenceGetter.get_read_sequence(r);
        skf.create_kmers(seq_ptr, read_kmers);
        read_sz1[r]=read_kmers.size();
    }

    bool equal=true;
    uint32_t it=0;
    int i=0;
    while(equal and it < 100) {
        for (i = 0; i < reads_to_check.size() and equal; i++) {
            seq_ptr = bufferedSequenceGetter.get_read_sequence(reads_to_check[i]);
            read_kmers.clear();
            skf.create_kmers(seq_ptr, read_kmers);
            if (read_kmers.size() != read_sz1[reads_to_check[i]]) equal = false;
        }
        it++;
    }

    if (!equal) {
        std::cout << "Failed during iteration " << it << ", on read " << reads_to_check[i];
    }

    REQUIRE(equal);
}

TEST_CASE("Fastq file reader") {
    FastqReader<FastqRecord> fastqReader({0}, "../tests/datasets/test.fastq");
    FastqRecord read;
    bool c;
    c = fastqReader.next_record(read);
    uint num_reads(0);
    while (c) {
        num_reads++;
        c = fastqReader.next_record(read);
    }

    REQUIRE(num_reads==10);
}

TEST_CASE("Fasta file reader") {
    FastaReader<FastaRecord> fastaReader({0}, "../tests/datasets/test.fasta");
    FastaRecord read;
    bool c;
    c = fastaReader.next_record(read);
    uint num_reads(0);
    while (c) {
        num_reads++;
        c = fastaReader.next_record(read);
    }

    REQUIRE(num_reads==3);
}