#include <iostream>
#include <fstream>
#include <sglib/counter/KmerCounter.hpp>
#include "sglib/logger/OutputLog.hpp"
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    unsigned int n_bins{256};
    bool loadDist=false;
    bool loadParts=false;
    std::string sample;
    uint8_t K=21;
    uint8_t mmer=5;
    uint num_threads=8;
    {
        cxxopts::Options options("bsg-kcount", "Count kmers from a sample");
        options.add_options()
                ("help", "Print help")
                ("s,sample", "Sample dataset to count", cxxopts::value(sample))
                ("d,distribution", "Load distributions from file", cxxopts::value(loadDist))
                ("p,partitions", "Load partitions from directory", cxxopts::value(loadParts))
                ("b,bins", "Number of bins used for distribution", cxxopts::value(n_bins))
                ("k", "Kmer size to count", cxxopts::value(K))
                ("m,minimiser", "Size of the minimisers", cxxopts::value(mmer))
                ("t,threads", "Number of threads", cxxopts::value(num_threads));
        try {
            auto result(options.parse(argc, argv));
            if (result.count("help")) {
                std::cout << options.help({""}) << std::endl;
                exit(0);
            }
            if (result.count("sample") == 0){
                throw cxxopts::missing_argument_exception("sample");
            }
        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl;
            std::cout << options.help({""}) << std::endl;
            exit(1);
        }
    }


    KmerCounter kmerCounter(K, mmer, n_bins, loadDist, loadParts);
    std::vector<PairedReadsDatastore> datastores;

    std::ifstream t(sample.data());
    if (t.is_open()) {
        std::cout << "Loading from disk" << std::endl;
        PairedReadsDatastore prds(sample.data());
        datastores.emplace_back(prds);
    } else {
        std::cout << "Creating from PE data" << std::endl;
        PairedReadsDatastore prds(argv[1], argv[2], sample.data(), 0, 310);
        datastores.emplace_back(prds);
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    kmerCounter.processReads(datastores);

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double >>(t2 - t1);

    std::cout << "Total time elapsed= " << time_span.count() << std::endl;
    return 0;
}