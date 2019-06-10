#include <iostream>
#include <fstream>
#include <sglib/workspace/WorkSpace.hpp>
#include "sglib/logger/OutputLog.hpp"
#include "cxxopts.hpp"

int main(int argc, char * argv[]) {
    std::cout << "Welcome to bsg-lrmapper"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    bool sat_kmer_index=false;
    unsigned int k=15;
    unsigned int max_filter=200;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string workspace_file,output_prefix;
    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;
    try
    {
        cxxopts::Options options("bsg-lrmapper", "Long_reads-to-graph mapper for bsg worskpaces");

        options.add_options()
                ("h,help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value(workspace_file))
                ("o,output", "output file prefix", cxxopts::value(output_prefix))
                ("k","long read indexing/mapping kmer size", cxxopts::value(k)->default_value("15"))
                ("m,max_kmer_repeat", "maximum number of times a kmer appears (LongReadMapper)", cxxopts::value(max_filter)->default_value("200"))
                ("s,use_sat-index", "Use saturated small-k index", cxxopts::value(sat_kmer_index)->default_value("false"));



        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("w")!=1 or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input workspace and output prefix");
        }



    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    std::cout<<std::endl;
    WorkSpace ws;
    sglib::OutputLog()<<"Loading Workspace..."<<std::endl;
    ws.load_from_disk(workspace_file);
    ws.add_log_entry("bsg-mapper run started");
    sglib::OutputLog()<<"Loading Workspace DONE"<<std::endl;
    sglib::OutputLog()<<"Mapping reads..."<<std::endl;
    auto pri=0;
    for (uint32_t lrds_idx=0; lrds_idx < ws.long_read_datastores.size(); lrds_idx++) {
        sglib::OutputLog()<<"Mapping reads from long reads library..."<<std::endl;
        ws.long_read_mappers[lrds_idx] = LongReadMapper(ws.sg, ws.long_read_datastores[lrds_idx], k, sat_kmer_index);
        m.map_reads(max_filter);
        ws.add_log_entry("reads from "+m.datastore.filename+" re-mapped to current graph");
        sglib::OutputLog()<<"Mapping reads from long reads library DONE."<<std::endl;
    }
    ws.add_log_entry("bsg-mapper run finished");
    ws.dump_to_disk(output_prefix+".bsgws");
    sglib::OutputLog()<<"Mapping reads DONE."<<std::endl;
    return 0;
}

