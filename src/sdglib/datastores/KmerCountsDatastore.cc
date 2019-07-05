//
// Created by Bernardo Clavijo (EI) on 2019-06-25.
//

#include "KmerCountsDatastore.hpp"
#include <sdglib/workspace/WorkSpace.hpp>

KmerCountsDatastore::KmerCountsDatastore(const WorkSpace &_ws, std::ifstream &infile): ws(_ws) {
    read(infile);
}

void KmerCountsDatastore::index_sdg() {
    //add all k-mers from SDG
    counts.clear();
    count_names.clear();
    uint64_t t=0;
    for(auto &n:ws.sdg.nodes) if (n.sequence.size()>=k) t+=n.sequence.size()+1-k;
    kindex.reserve(t);
    StringKMerFactory skf(k);
    for(auto &n:ws.sdg.nodes) if (n.sequence.size()>=k) skf.create_kmers(n.sequence,kindex);
    //sort
    std::sort(kindex.begin(),kindex.end());

    //create a first count
    counts.emplace_back();
    count_names.emplace_back("sdg");
    auto &c=counts.back();
    c.reserve(kindex.size());

    //collapse, but save coverage to the first count
    auto wi=kindex.begin();
    auto ri=kindex.begin();
    for (;ri<kindex.end();++wi){
        *wi=*ri;
        c.emplace_back(1);
        while(++ri<kindex.end() and *ri==*wi) ++(c.back());
    }
    kindex.resize(c.size());
}

void KmerCountsDatastore::add_count(const std::string &count_name, const std::vector<std::string> &filenames) {
    count_names.emplace_back(count_name);
    counts.emplace_back(kindex.size());
    uint64_t present(0), absent(0), rp(0);
    sdglib::OutputLog(sdglib::INFO)<<"Populating lookup map"<<std::endl;
    std::unordered_map<uint64_t,uint64_t> kmer_map;
    kmer_map.reserve(kindex.size());
    for (uint64_t i=0;i<kindex.size();++i) kmer_map[kindex[i]]=i;
    sdglib::OutputLog(sdglib::INFO)<<"Map populated with "<<kmer_map.size()<<" entries"<< std::endl;
    for (auto filename:filenames) {
        sdglib::OutputLog(sdglib::INFO) << "Counting from file: " << filename << std::endl;
        FastqReader<FastqRecord> fastqReader({0}, filename);
#pragma omp parallel shared(fastqReader)
        {
            uint64_t thread_present(0), thread_absent(0), thread_rp(0);
            const size_t local_kmers_size = 2000000;
            std::vector<uint64_t> found_kmers;
            found_kmers.reserve(local_kmers_size);
            FastqRecord read;
            std::vector<uint64_t> readkmers;
            StringKMerFactory skf(k);

            bool c;
#pragma omp critical(fastqreader)
            {
                c = fastqReader.next_record(read);
            }
            while (c) {
                readkmers.clear();
                skf.create_kmers(read.seq,readkmers);


                for (auto &rk:readkmers) {
                    auto findk = kmer_map.find(rk);
                    if (kmer_map.end() != findk) {
                        //++thread_counts[findk->second];
                        found_kmers.emplace_back(findk->second);
                        if (found_kmers.size() == local_kmers_size) {
#pragma omp critical(results_merge)
                            {
                                auto &arc = counts.back();
                                for (auto &x:found_kmers) if (arc[x] < UINT16_MAX) ++arc[x];
                                if (rp / 100000 != (rp + thread_rp) / 100000)
                                    sdglib::OutputLog(sdglib::INFO) << rp << " reads processed " << present << " / "
                                                                    << present + absent << " kmers found" << std::endl;
                                rp += thread_rp;
                                present += thread_present;
                                absent += thread_absent;
                            }
                            found_kmers.clear();
                            thread_absent = 0;
                            thread_present = 0;
                            thread_rp = 0;

                        }
                        ++thread_present;
                    } else ++thread_absent;
                }
                ++thread_rp;
#pragma omp critical(fastqreader)
                c = fastqReader.next_record(read);
            }
#pragma omp critical(results_merge)
            {
                auto &arc = counts.back();
                for (auto &x:found_kmers) if (arc[x] < UINT16_MAX) ++arc[x];
                rp += thread_rp;
                present += thread_present;
                absent += thread_absent;
            }
        }
        sdglib::OutputLog(sdglib::INFO) << rp << " reads processed " << present << " / " << present + absent
                                        << " kmers found" << std::endl;
    }
}

/** This template is used to do the counts from the datastores, it is templatised here rather than on the header **/
template<class T>
void add_count_to_kds( KmerCountsDatastore & kds, const std::string & count_name, const T & datastore){
    kds.count_names.emplace_back(count_name);
    kds.counts.emplace_back(kds.kindex.size());
    uint64_t present(0), absent(0), rp(0);
    sdglib::OutputLog(sdglib::INFO)<<"Populating lookup map"<<std::endl;
    std::unordered_map<uint64_t,uint64_t> kmer_map;
    kmer_map.reserve(kds.kindex.size());
    for (uint64_t i=0;i<kds.kindex.size();++i) kmer_map[kds.kindex[i]]=i;
    sdglib::OutputLog(sdglib::INFO)<<"Map populated with "<<kmer_map.size()<<" entries, counting from datastore: " << datastore.filename << std::endl;
#pragma omp parallel
    {
        ReadSequenceBuffer bpsg(datastore);
        uint64_t thread_present(0), thread_absent(0), thread_rp(0);
        const size_t local_kmers_size = 2000000;
        std::vector<uint64_t> found_kmers; // kmer index of found kmers is saved here, increments are done in the critical
        found_kmers.reserve(local_kmers_size);
        std::vector<uint64_t> readkmers;
        CStringKMerFactory cskf(kds.get_k());
#pragma omp for schedule(static,10000)
        for (uint64_t rid = 1; rid <= datastore.size(); ++rid) {
            readkmers.clear();
            cskf.create_kmers(readkmers, bpsg.get_read_sequence(rid));

            for (auto &rk:readkmers) {
                auto findk = kmer_map.find(rk);
                if (kmer_map.end() != findk) {
                    //++thread_counts[findk->second];
                    found_kmers.emplace_back(findk->second);
                    if (found_kmers.size() == local_kmers_size) {
#pragma omp critical(results_merge)
                        {
                            auto &arc=kds.counts.back();
                            for (auto &x:found_kmers) if (arc[x] < UINT16_MAX) ++arc[x];
                            if (rp / 100000 != (rp + thread_rp) / 100000)
                                sdglib::OutputLog(sdglib::INFO) << rp << " reads processed " << present << " / "
                                                                << present + absent << " kmers found" << std::endl;
                            rp += thread_rp;
                            present += thread_present;
                            absent += thread_absent;
                        }
                        found_kmers.clear();
                        thread_absent = 0;
                        thread_present = 0;
                        thread_rp = 0;

                    }
                    ++thread_present;
                } else ++thread_absent;
            }
            ++thread_rp;
        }
#pragma omp critical(results_merge)
        {
            auto &arc=kds.counts.back();
            for (auto &x:found_kmers) if (arc[x] < UINT16_MAX) ++arc[x];
            rp += thread_rp;
            present += thread_present;
            absent += thread_absent;
        }
    }
    sdglib::OutputLog(sdglib::INFO) << rp << " reads processed " << present << " / " << present + absent
                                    << " kmers found" << std::endl;
}

void KmerCountsDatastore::add_count(const std::string & count_name, const PairedReadsDatastore & datastore){
    add_count_to_kds(*this,count_name,datastore);
}
void KmerCountsDatastore::add_count(const std::string & count_name, const LinkedReadsDatastore & datastore){
    add_count_to_kds(*this,count_name,datastore);
}
void KmerCountsDatastore::add_count(const std::string & count_name, const LongReadsDatastore & datastore){
    add_count_to_kds(*this,count_name,datastore);
}

std::vector<uint16_t> KmerCountsDatastore::project_count(const uint16_t count_idx, const std::string &s) {
    std::vector<uint64_t> skmers;
    StringKMerFactory skf(k);
    skf.create_kmers(s,skmers);
    std::vector<uint16_t> kcov;
    for (auto &kmer: skmers){
        auto nk = std::lower_bound(kindex.begin(), kindex.end(), kmer);

        if (nk!=kindex.end() and *nk == kmer) {
            kcov.push_back(counts[count_idx][nk-kindex.begin()]);

        } else {
            kcov.push_back(0);
        }
    }
    return kcov;
}

std::vector<uint16_t> KmerCountsDatastore::project_count(const std::string &count_name, const std::string &s) {
    auto cnitr=std::find(count_names.begin(),count_names.end(),count_name);
    if (cnitr!=count_names.end()){
        return project_count(cnitr-count_names.begin(),s);
    }
    return {};
}

void KmerCountsDatastore::read(std::ifstream &input_file) {
    input_file.read((char *) &k, sizeof(k));
    sdglib::read_string(input_file,name);
    sdglib::read_stringvector(input_file,count_names);
    sdglib::read_flat_vectorvector(input_file,counts);

}

void KmerCountsDatastore::write(std::ofstream &output_file) {
    output_file.write((char *) &k, sizeof(k));
    sdglib::write_string(output_file,name);
    sdglib::write_stringvector(output_file,count_names);
    sdglib::write_flat_vectorvector(output_file,counts);
}
