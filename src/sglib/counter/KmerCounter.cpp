//
// Created by Luis Yanes (EI) on 12/09/2018.
//

#include <map>
#include <sglib/processors/KmerCompressionIndex.hpp>
#include "KmerCounter.hpp"

void KmerCounter::processReads(std::vector<PairedReadsDatastore> &read_files, bool sequential) {
 if (K>12) {
     processReadsLargeK(read_files);
 } else {
     if (sequential) {
         processReadsSmallKSequential(read_files);
     } else {
         processReadsSmallK(read_files);
     }
 }

    std::cout << "\n\n";
    std::cout << "Total no. of reads: " << total_reads_kmerised << "\n";
    std::cout << "Total no. of reads processed: " << total_reads_processed << "\n";
    std::cout << "Total no. of k-mers: " << total_kmers_produced << "\n";
    std::cout << "Done!" << std::endl;
}

void KmerCounter::processBins(BinBufferWriterQueue &wrt_queue) {
    wrt_queue.num_writers++;
    while (!wrt_queue.empty()) {
        // Get largest bin
        auto max_bin = wrt_queue.getMaxBin();
        // Write and pop
        if (max_bin) wrt_queue.write_bin(max_bin);
    }
    wrt_queue.num_writers--;
}


void KmerCounter::read_file_small_k(const PairedReadsDatastore &ds, int thread_id) {

    auto bpsg=BufferedPairedSequenceGetter(ds,128*1024,ds.readsize*2+2);
    std::vector<Kmer> readkmers;
    StringKMerFactory skf(K);
    uint64_t reads_kmerised(0);
    uint64_t kmers_produced(0);
    KmerArray insertKmers;
    auto first_read = 1+(thread_id*ds.size()/num_datastore_readers);
    auto last_read = 1+( (thread_id+1)*ds.size()/num_datastore_readers);
    uint64_t rid;
    kmer_chunks.num_writers++;
    for (rid = first_read; rid < last_read and rid <= ds.size(); ++rid) {
        readkmers.clear();
        auto read_sequence = bpsg.get_read_sequence(rid);
        skf.create_kmers(read_sequence, readkmers);
        kmers_produced+=readkmers.size();
        insertKmers[0].second = readkmers.size();
        for (int i = 0; i < readkmers.size(); i++) {
            insertKmers[i+1] = readkmers[i];
        }
        reads_kmerised++;
        while(!kmer_chunks.enqueue(insertKmers)) {
            std::this_thread::yield();
        }
    }
    total_reads_kmerised+=reads_kmerised;
    total_kmers_produced+=kmers_produced;
    kmer_chunks.num_writers--;
}

void KmerCounter::count_small_k(int thread_id, std::vector<uint64_t> &smallK_counts) {
    KmerArray insertKmers;
    uint chunks_processed(0);
    std::vector<uint64_t> thread_counts(smallK_counts.size());

    for (;;){
        if (kmer_chunks.dequeue(insertKmers)) {
            chunks_processed++;
            for (int i = 0; i < insertKmers[0].second; i++) {
                thread_counts[insertKmers[i+1].second] += 1;
            }
        } else {
            if (!kmer_chunks.empty()) {
                std::this_thread::yield();
            } else {
                break;
            }
        }
    }

    {
        std::lock_guard<std::mutex> lg(smallK_merge_lock);
        for (unsigned int i = 0; i < smallK_counts.size(); i++) {
            if (smallK_counts[i] + thread_counts[i] > 2500) {
                smallK_counts[i] = 2500;
            } else {
                smallK_counts[i] += thread_counts[i];
            }
        }
    }
    total_reads_processed += chunks_processed;
}

void KmerCounter::processReadsSmallK(std::vector<PairedReadsDatastore> &read_files) {
    std::vector<uint64_t> smallK_counts(1 << (2 * K));

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    for (const auto &r: read_files) {
        for (int i = 0; i < num_datastore_readers; i++) {
            std::thread th(&KmerCounter::read_file_small_k, this, std::ref(r), i);
            kmer_producers.push_back(std::move(th));
        }

        usleep(250);

        for (int i = 0; i < num_smallK_counters; i++) {
            std::thread th(&KmerCounter::count_small_k, this, i+1, std::ref(smallK_counts));
            kmer_counters.push_back(std::move(th));
        }

        for (auto &f: kmer_producers){
            if (f.joinable()) f.join();
        }

        for (auto &f: kmer_counters){
            if (f.joinable()) f.join();
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double > read_processing = std::chrono::duration_cast<std::chrono::duration<double >>(t2 - t1);
    t1 = t2;

    std::vector<uint64_t> histogram(1024);
    std::ofstream hist_file("/tmp/bsg.hist");
    for (unsigned int i = 0; i < smallK_counts.size(); i++) {
        histogram[smallK_counts[i]]++;
    }

    for (unsigned int i = 1; i < histogram.size(); i++) {
        hist_file << i << " " <<  histogram[i] << "\n";
    }

    t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double > histogram_making = std::chrono::duration_cast<std::chrono::duration<double >>(t2 - t1);
    t1 = t2;
    std::cout <<"Read processing time = " << read_processing.count() << std::endl;
    std::cout <<"Histogram making time elapsed = " << histogram_making.count() << std::endl;
}

void KmerCounter::processReadsSmallKSequential(std::vector<PairedReadsDatastore> &read_files) {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::vector<uint64_t> smallK_counts(1 << (2 * K));
    for (const auto &ds: read_files) {
        auto bpsg = BufferedPairedSequenceGetter(ds, 128 * 1024, ds.readsize * 2 + 2);
        std::vector<Kmer> readkmers;
        StringKMerFactory skf(K);
        uint64_t reads_kmerised(0);
        uint64_t kmers_produced(0);
        KmerArray insertKmers;
        uint64_t rid;
        for (rid = 1; rid <= ds.size(); ++rid) {
            readkmers.clear();
            auto read_sequence = bpsg.get_read_sequence(rid);
            skf.create_kmers(read_sequence, readkmers);

            kmers_produced+=readkmers.size();

            for (const auto &rk:readkmers) {
                smallK_counts[rk.second] += 1;
            }

            reads_kmerised++;
        }
        total_reads_processed+=reads_kmerised;
        total_kmers_produced+=kmers_produced;
        total_reads_kmerised+=reads_kmerised;
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double > read_processing = std::chrono::duration_cast<std::chrono::duration<double >>(t2 - t1);
    t1 = t2;


    std::vector<uint64_t> histogram(1024);
    std::ofstream hist_file("/tmp/bsg.hist");
    for (unsigned int i = 0; i < smallK_counts.size(); i++) {
        histogram[smallK_counts[i]]++;
    }

    for (unsigned int i = 1; i < histogram.size(); i++) {
        hist_file << i << " " <<  histogram[i] << "\n";
    }
    t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double > histogram_making = std::chrono::duration_cast<std::chrono::duration<double >>(t2 - t1);
    t1 = t2;
    std::cout <<"Read processing time = " << read_processing.count() << std::endl;
    std::cout <<"Histogram making time elapsed = " << histogram_making.count() << std::endl;
    std::cout << std::endl;
}




void KmerCounter::read_file_large_k(const PairedReadsDatastore &ds, BinBufferWriterQueue &wrt_queue, int thread_id) {

    auto bpsg=BufferedPairedSequenceGetter(ds,128*1024,ds.readsize*2+2);
    std::vector<Kmer> readkmers;
    StringKMerFactory skf(K);
    uint64_t reads_kmerised(0);
    uint64_t kmers_produced(0);
    KmerArray insertKmers;
    std::list<Seq2bit> this_thread_reads_2bit;
    auto first_read = 1+(thread_id*ds.size()/num_datastore_readers);
    auto last_read = 1+( (thread_id+1)*ds.size()/num_datastore_readers);
    uint64_t rid;
    kmer_chunks.num_writers++;

    for (rid = first_read; rid < last_read and rid <= ds.size(); ++rid) {
        readkmers.clear();
        auto read_sequence = bpsg.get_read_sequence(rid);

        // Produce the minimiser based super-kmers
        make_minimiser_send_to_bin(read_sequence, wrt_queue.bins, this_thread_reads_2bit);

        // Insert in the min bin queue
        while(!kmer_chunks.enqueue(insertKmers)) {
            std::this_thread::yield();
        }
    }
    total_reads_kmerised+=reads_kmerised;
    total_kmers_produced+=kmers_produced;
    kmer_chunks.num_writers--;
}

void KmerCounter::processReadsLargeK(std::vector<PairedReadsDatastore> &read_files) {

    BinBufferWriterQueue wrt_queue(512);
    for (int i = 0; i < n_bins; i++) {
        wrt_queue.bins.emplace_back(i, K, wrt_queue, 32*1024);
    }

    for (const auto &r: read_files) {
        for (int i = 0; i < num_datastore_readers; i++) {
            std::thread th(&KmerCounter::read_file_large_k, this, std::ref(r), std::ref(wrt_queue), i);
            kmer_producers.push_back(std::move(th));
        }

        usleep(250);

        for (int i = 0; i < num_bin_writers; i++) {
            std::thread th(&KmerCounter::processBins, this, std::ref(wrt_queue));
            kmer_bin_writers.push_back(std::move(th));
        }

        for (auto &f: kmer_producers){
            if (f.joinable()) f.join();
        }

        for (auto &f: kmer_bin_writers){
            if (f.joinable()) f.join();
        }

    }
}

void KmerCounter::make_minimiser_send_to_bin(std::string in_seq, std::vector<ExtendedKmerBin> &bins,
                                             std::list<Seq2bit> &reads_2bit) {
    int i = 0;
    int len = 0;
    unsigned int bin_no = 0;
    Mmer current_mmer(signature_len), end_mmer(signature_len);
    reads_2bit.emplace_back(Seq2bit(in_seq));
    auto& seq = reads_2bit.back();
    while (in_seq[i]!='\0')
    {
        bool contains_N = false;
        //building first signature after 'N' or at the read begining
        for (uint32_t j = 0; j < signature_len; ++j, ++i)
            if (seq[i] == 4)//'N'
            {
                contains_N = true;
                break;
            }
        //signature must be shorter than k-mer so if signature contains 'N', k-mer will contains it also
        if (contains_N)
        {
            ++i;
            continue;
        }
        len = signature_len;
        auto signature_start_pos = i - signature_len;
        current_mmer.insert(reads_2bit.back().data()+signature_start_pos);
        end_mmer.set(current_mmer);
        for (; in_seq[i]!='\0'; ++i)
        {
            if (seq[i] > 3)//'N'
            {
                if (len >= K)
                {
                    bin_no = mmer_table.get_bin_id(current_mmer.get());
                    bins[bin_no].store_superkmer(reads_2bit.back().data() + i - len, len);
                }
                len = 0;
                ++i;
                break;
            }
            end_mmer.insert(reads_2bit.back().data()+i);
            if (end_mmer < current_mmer)//signature at the end of current k-mer is lower than current
            {
                if (len >= K)
                {
                    bin_no = mmer_table.get_bin_id(current_mmer.get());
                    bins[bin_no].store_superkmer(reads_2bit.back().data() + i - len, len);
                    len = K - 1;
                }
                current_mmer.set(end_mmer);
                signature_start_pos = i - signature_len + 1;
            }
            else if (end_mmer == current_mmer)
            {
                current_mmer.set(end_mmer);
                signature_start_pos = i - signature_len + 1;
            }
            else if (signature_start_pos + K - 1 < i)//need to find new signature
            {
                bin_no = mmer_table.get_bin_id(current_mmer.get());
                bins[bin_no].store_superkmer(reads_2bit.back().data() + i - len, len);
                len = K - 1;
                //looking for new signature
                ++signature_start_pos;
                //building first signature in current k-mer
                end_mmer.insert(reads_2bit.back().data()+signature_start_pos);
                current_mmer.set(end_mmer);
                for (uint32_t j = signature_start_pos + signature_len; j <= i; ++j)
                {
                    end_mmer.insert(reads_2bit.back().data()+j);
                    if (end_mmer <= current_mmer)
                    {
                        current_mmer.set(end_mmer);
                        signature_start_pos = j - signature_len + 1;
                    }
                }
            }
            ++len;
            if (len == K + 255) //one byte is used to store counter of additional symbols in extended k-mer
            {
                bin_no = mmer_table.get_bin_id(current_mmer.get());
                bins[bin_no].store_superkmer(reads_2bit.back().data() + i + 1 - len, len);
                i -= K - 2;
                len = 0;
                break;
            }

        }
    }
    if (len >= K)//last one in read
    {
        bin_no = mmer_table.get_bin_id(current_mmer.get());
        bins[bin_no].store_superkmer(seq.data() + i - len, len);
    }

}





















void KmerCounter::LBread_file(PairedReadsDatastore &ds, int thread_id) {

    BufferedPairedSequenceGetter bpsg(ds,100000,1000);
    std::vector<std::pair<bool, uint64_t>> readkmers(300);
    StringKMerFactory skf(K);
    uint64_t reads_kmerised(0);
    KmerArray insertKmers;
    auto first_read = 1+(thread_id*ds.size()/4);
    auto last_read = 1+( (thread_id+1)*ds.size()/4);
    std::cout <<"Reader " << thread_id <<" started" <<std::endl;
    uint64_t rid;

    for (rid = first_read; rid < last_read and rid <= ds.size(); ++rid) {
        readkmers.clear();
        skf.create_kmers(bpsg.get_read_sequence(rid), readkmers);
        for (int i = 0; i < readkmers.size(); i++) {
            insertKmers[i] = readkmers[i];
        }
        reads_kmerised++;
        LBkmer_chunks.push(insertKmers);
    }
    total_reads_kmerised+=reads_kmerised;
    LBkmer_chunks.mark_completed();
    std::cout << "Reader thread " << thread_id << " done!" << std::endl;
}

// Get kmers from the queue and try to insert them in the appropriate minimiser bin
void KmerCounter::LBprocess_kmers(int thread_id) {
    KmerArray data, proc_data;
    uint chunks_processed(0);

    std::cout << "Process thread started"<<std::endl;
    usleep(thread_id*200);

    while (!LBkmer_chunks.completed()) {
        LBkmer_chunks.pop(data);
        chunks_processed++;
        if (chunks_processed % 100000 == 0){
            std::lock_guard<std::mutex> lg(write_lock);
            std::cout << "Processor " << thread_id << ": ";
            for (const auto &d: data) {
                std::cout << d.second << ", ";
            }
            std::cout << std::endl;
        }
    }
    {
        std::cout << "Thread " << thread_id << " processed " << chunks_processed << " reads" << std::endl;
    }
    total_reads_processed += chunks_processed;
}
