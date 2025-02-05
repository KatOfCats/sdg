//
// Created by Bernardo Clavijo (EI) on 2019-06-25.
//

#include "KmerCounts.hpp"
#include <sdglib/workspace/WorkSpace.hpp>

KmerCounts::KmerCounts(const WorkSpace &_ws, std::ifstream &infile): ws(_ws) {
    read(infile);
}

void KmerCounts::index_sdg() {
    //add all k-mers from SDG
    counts.clear();
    count_names.clear();
    uint64_t t=0;
    for(auto &n:ws.sdg.nodes) if (n.sequence.size()>=k) t+=n.sequence.size()+1-k;
    kindex.reserve(t);
    if (count_mode==Canonical) {
        StringKMerFactory skf(k);
        for (auto &n:ws.sdg.nodes) if (n.sequence.size() >= k) skf.create_kmers(n.sequence, kindex);
    } else if (count_mode==NonCanonical) {
        StringKMerFactoryNC skf(k);
        for (auto &n:ws.sdg.nodes) if (n.sequence.size() >= k) skf.create_kmers(n.sequence, kindex);
    }
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

void KmerCounts::add_count(const std::string &count_name, const std::vector<std::string> &filenames, bool fastq) {
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
        FastaReader<FastaRecord> fastaReader({0}, filename);
#pragma omp parallel shared(fastqReader)
        {
            uint64_t thread_present(0), thread_absent(0), thread_rp(0);
            const size_t local_kmers_size = 2000000;
            std::vector<uint64_t> found_kmers;
            found_kmers.reserve(local_kmers_size);
            FastqRecord read;
            FastaRecord reada;

            std::vector<uint64_t> readkmers;
            StringKMerFactory skf(k);
            StringKMerFactoryNC skfnc(k);

            bool c;
#pragma omp critical(fastqreader)
            {
                if (fastq) c = fastqReader.next_record(read);
                else c= fastaReader.next_record(reada);
            }
            while (c) {
                readkmers.clear();
                if (count_mode==Canonical) {
                    skf.create_kmers((fastq?read.seq:reada.seq),readkmers);
                } else if (count_mode==NonCanonical) {
                    skfnc.create_kmers((fastq?read.seq:reada.seq),readkmers);
                }



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
                {
                    if (fastq) c = fastqReader.next_record(read);
                    else c= fastaReader.next_record(reada);
                }
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
void add_count_to_kds( KmerCounts & kds, const std::string & count_name, const T & datastore) {
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

void KmerCounts::add_count(const std::string & count_name, const PairedReadsDatastore & datastore){
    add_count_to_kds(*this,count_name,datastore);
}
void KmerCounts::add_count(const std::string & count_name, const LinkedReadsDatastore & datastore){
    add_count_to_kds(*this,count_name,datastore);
}
void KmerCounts::add_count(const std::string & count_name, const LongReadsDatastore & datastore){
    add_count_to_kds(*this,count_name,datastore);
}

std::vector<uint16_t> KmerCounts::project_count(const uint16_t count_idx, const std::string &s) {
    std::vector<uint64_t> skmers;

    //StringKMerFactory skf(k);
    //skf.create_kmers(s,skmers);
    if (count_mode==Canonical) {
        StringKMerFactory skf(k);
        skf.create_kmers(s,skmers);
    } else if (count_mode==NonCanonical) {
        StringKMerFactoryNC skf(k);
        skf.create_kmers(s,skmers);
    }
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

std::vector<uint16_t> KmerCounts::project_count(const std::string &count_name, const std::string &s) {
    auto cnitr=std::find(count_names.begin(),count_names.end(),count_name);
    if (cnitr!=count_names.end()){
        return project_count(cnitr-count_names.begin(),s);
    }
    return {};
}

void KmerCounts::read(std::ifstream &ws_file) {
    std::string filepath;
    sdglib::read_string(ws_file, filepath);
    std::ifstream count_file(filepath);
    read_counts(count_file);
}

void KmerCounts::read_counts(std::ifstream &count_file) {
    count_file.read((char *) &k, sizeof(k));
    count_file.read((char *) &count_mode, sizeof(count_mode));
    sdglib::read_string(count_file,name);
    sdglib::read_stringvector(count_file,count_names);
    sdglib::read_flat_vector(count_file,kindex);
    sdglib::read_flat_vectorvector(count_file,counts);

}

void KmerCounts::write(std::ofstream &output_file) const {
    sdglib::write_string(output_file, name+".count");
    std::ofstream count_file(name+".count");
    write_counts(count_file);
}
void KmerCounts::write(std::fstream &output_file) const {
    sdglib::write_string(output_file, name+".count");
    std::ofstream count_file(name+".count");
    write_counts(count_file);
}

void KmerCounts::write_counts(std::ofstream &count_file) const {
    count_file.write((char *) &k, sizeof(k));
    count_file.write((char *) &count_mode, sizeof(count_mode));
    sdglib::write_string(count_file,name);
    sdglib::write_stringvector(count_file,count_names);
    sdglib::write_flat_vector(count_file,kindex);
    sdglib::write_flat_vectorvector(count_file,counts);
}

std::vector<std::string> KmerCounts::get_count_names() {
    return count_names;
}

const std::vector<uint16_t> &KmerCounts::get_count_by_name(const std::string &name) const {
    for (int i = 0; i < count_names.size(); i++) {
        if (count_names[i] == name) {
            return counts[i];
        }
    }

    throw std::runtime_error("Couldn't find a count named: "+name);

}
