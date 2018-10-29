//
// Created by Luis Yanes (EI) on 12/09/2018.
//

#include <map>
#include <random>
#include <iomanip>
#include <sglib/utilities/omp_safe.hpp>
#include <sglib/processors/KmerCompressionIndex.hpp>
#include "KmerCounter.hpp"

void KmerCounter::processReads(std::vector<PairedReadsDatastore> &read_files) {
 if (K>12) {
     if (K<64) {
         processReadsLargeK(read_files);
     }
     else {
         throw std::logic_error("Function not yet implemented");
     }
 } else {
     throw std::logic_error("Function not yet implemented");
 }

    std::cout << "\n\n";
    std::cout << "Total no. of reads: " << total_reads_kmerised << "\n";
    std::cout << "Total no. of k-mers: " << total_kmers_produced << "\n";
    std::cout << "Total no. of superkmers: " << total_superkmers<< "\n";
    std::cout << "Done!" << std::endl;
}


void KmerCounter::processBins(BinBufferWriterQueue &wrt_queue) {
    while (!wrt_queue.empty()) {
        // Get largest bin
        auto max_bin = wrt_queue.getMaxBin();
        // Write and pop
        if (max_bin) wrt_queue.write_bin(max_bin);
    }
}


/**
 * Calculates the mmer distribution a sample of the input files
 *
 * in increasing number of kmers and not uniformly
 * @param read_files
 */
void KmerCounter::calculateDistribution(const std::vector<PairedReadsDatastore> &read_files) {
    {
        if (loadDist) {
            mmer_table.loadDistribution();
            return;
        }

        for (const auto &r: read_files) {
            std::thread th(&KmerCounter::mmer_table_initialiser, this, std::ref(r), std::ref(mmer_table));
            kmer_dist_generators.push_back(std::move(th));
        }

        for (auto &f: kmer_dist_generators) {
            if (f.joinable()) f.join();
        }
    }

    mmer_table.calcBins();
    // Dump this to disk, so that it can be used/loaded for dbg purposes
}


void KmerCounter::loadPartitions(BinBufferWriterQueue &wrt_queue) {
    std::ifstream partDesc("partDesc.txt");
    std::cout << "Loading partitions" << std::endl;
    std::string sep;
    // Make sure all bins are empty before closing the files and destroying all objects
    for (int i = 0; i < wrt_queue.bin_desc.size(); i++) {
        uint64_t bin, n_rec, n_plus_x_recs, size, dummy;
        partDesc >> bin >> n_rec >> n_plus_x_recs >> size >> dummy >> dummy;
        wrt_queue.bin_desc[bin].n_rec = n_rec;
        wrt_queue.bin_desc[bin].n_plus_x_recs = n_plus_x_recs;
        wrt_queue.bin_desc[bin].size = size;
        std::cout << bin+1 << "\t\t" << n_rec << "\t" << n_plus_x_recs << "\t" << size << std::endl;
    }
}

/**
 * Sends the kmers to partitions based on the distribution calculated on the previous step
 * This function distributes the superkmers in the reads to the bin files.
 * @param wrt_queue
 * @param read_files
 */
void KmerCounter::createPartitions(BinBufferWriterQueue &wrt_queue, std::vector<PairedReadsDatastore> &read_files) {
    if (loadParts) {
        loadPartitions(wrt_queue);
        return;
    }
    for (const auto &r: read_files) {
        for (int i = 0; i < num_datastore_readers; i++) {
            std::thread th(&KmerCounter::read_file_large_k, this, std::ref(r), std::ref(wrt_queue), i);
            kmer_producers.push_back(std::move(th));
        }

        for (int i = 0; i < num_bin_writers; i++) {
            std::thread th(&KmerCounter::processBins, this, std::ref(wrt_queue));
            kmer_bin_writers.push_back(std::move(th));
        }
    }

    for (auto &f: kmer_producers){
        if (f.joinable()) f.join();
    }

    for (auto &f: kmer_bin_writers){
        if (f.joinable()) f.join();
    }

    wrt_queue.flush();

    // Do some reporting
    uint64_t total_kmers{0};
    uint64_t total_superkmers{0};
    for (int i = 0; i < wrt_queue.bin_desc.size(); i++) {
        std::cout << std::setw(5) << std::setfill('0') << i << "\t" << wrt_queue.bin_desc[i].n_rec << "\t" << wrt_queue.bin_desc[i].n_plus_x_recs << std::endl;
        total_kmers +=wrt_queue.bin_desc[i].n_rec;
        total_superkmers += wrt_queue.bin_desc[i].n_plus_x_recs;
    }
    std::cout << "Total kmers in bin_desc: " << total_kmers << std::endl;

}

std::string asString(int K, uint64_t kc){
    std::string result;
    result.resize(K);
    for (int j = 0; j < K; j++) result[j] = "ACGT"[(kc >> 2 * j) & 3];
    return result;
};

std::string asString64(int K, __uint128_t kc){
    std::string result;
    result.resize(K);
    for (int j = 0; j < K; j++) result[j] = "ACGT"[(kc >> 2 * j) & 3];
    return result;
};

uint64_t rc(int K, uint64_t kc) {
    uint64_t result{0};
    for (int j = K - 1; j >= 0; --j) {
        unsigned char symb = (unsigned char) 3 & (unsigned char)(kc>>2*j);
        result |= (uint64_t)(3-symb)<<(K-1-j)*2;
    }
    return result;
}

/**
 * Load the bins one by one, expand the total number of kmers from the superkmers stored
 * Sort the expanded kmers from the bins and then use a counting unique to store the total number of appearances.
 *
 * This function should also store the kmers in a database (which hasn't been designed yet)
 * @param wrt_queue
 */
void KmerCounter::sort_and_store_partitions(BinBufferWriterQueue &wrt_queue) {
    // Load bins and sort them
    output_file.open("count_dump.txt");
    uint64_t total_counter_kmers{0};
    uint64_t total_kmers_in_bins{0};
    uint64_t total_kmers_expanded{0};
#pragma omp parallel
#pragma omp for schedule(static,1)
    for (int i = 0; i <= n_bins; i++) {
        if (wrt_queue.bin_desc[i].size == 0) continue;
        std::vector<unsigned char> buf;
        // Read the whole bin
        buf.resize(wrt_queue.bin_desc[i].size);
        auto read_sz = std::fread(&buf[0], sizeof(unsigned char), buf.size(), wrt_queue.bin_files[i]);

        if (read_sz < buf.size()) {
            throw std::runtime_error("Error reading bin "+std::to_string(i));
        }

        // Check number of kmers expanded from bins corresponds with kmers in bins
        total_kmers_in_bins+=wrt_queue.bin_desc[i].n_rec;

        uint64_t expanded_kmers{0};
        uint64_t n_bin_kmers{0};
        // Expand the bin, into new memory
        if (K < 31) {
            std::vector<KmerCount> bin_kmers;
            bin_kmers.resize(wrt_queue.bin_desc[i].n_rec);
            expanded_kmers = expand(buf, bin_kmers);

            total_kmers_expanded += expanded_kmers;

            if (wrt_queue.bin_desc[i].n_rec != expanded_kmers) {
                throw std::runtime_error("Error expanding bin " + std::to_string(i) +
                                         ", the expected number of kmers did not match the number of kmers expanded "
                                         + std::to_string(wrt_queue.bin_desc[i].n_rec) +
                                         std::to_string(expanded_kmers));
            }

            // Sort the bin
            std::sort(bin_kmers.begin(), bin_kmers.end());

            auto ritr = bin_kmers.cbegin();
            auto witr = bin_kmers.begin();
            for (; ritr != bin_kmers.cend();) {
                auto bitr = ritr;
                while (ritr != bin_kmers.cend() and bitr->kmer == ritr->kmer) {
                    ++ritr;
                }
                if (ritr - bitr > 0 /*min_k_count*/) {
                    witr->kmer = bitr->kmer;
                    if (ritr - bitr <= 255 /*max_k_count*/) {
                        witr->count = static_cast<uint8_t>(ritr - bitr);
                    } else {
                        witr->count = std::numeric_limits<uint8_t>::max()/*max_k_count*/;
                    }
                    ++witr;
                }
            }
            bin_kmers.resize(witr - bin_kmers.begin());
            n_bin_kmers = bin_kmers.size();

            // TODO: Store the bin (Design a DB for it, partition based or hashtable based?)
#pragma omp critical (fout)
            {
                for (const auto &kc: bin_kmers) {
                    output_file << ">" << (int) kc.count << "\n" << asString(K, kc.kmer) << "\n";
                }
            }
        }
        if (K>=31) {
            std::vector<KmerCount64> bin_kmers;
            bin_kmers.resize(wrt_queue.bin_desc[i].n_rec);
            expanded_kmers = expand64(buf, bin_kmers);

            total_kmers_expanded += expanded_kmers;

            if (wrt_queue.bin_desc[i].n_rec != expanded_kmers) {
                throw std::runtime_error("Error expanding bin " + std::to_string(i) +
                                         ", the expected number of kmers did not match the number of kmers expanded "
                                         + std::to_string(wrt_queue.bin_desc[i].n_rec) +
                                         std::to_string(expanded_kmers));
            }

            // Sort the bin
            std::sort(bin_kmers.begin(), bin_kmers.end());

            auto ritr = bin_kmers.cbegin();
            auto witr = bin_kmers.begin();
            for (; ritr != bin_kmers.cend();) {
                auto bitr = ritr;
                while (ritr != bin_kmers.cend() and bitr->kmer == ritr->kmer) {
                    ++ritr;
                }
                if (ritr - bitr > 0 /*min_k_count*/) {
                    witr->kmer = bitr->kmer;
                    if (ritr - bitr <= 255 /*max_k_count*/) {
                        witr->count = static_cast<uint8_t>(ritr - bitr);
                    } else {
                        witr->count = std::numeric_limits<uint8_t>::max()/*max_k_count*/;
                    }
                    ++witr;
                }
            }
            bin_kmers.resize(witr - bin_kmers.begin());
            n_bin_kmers = bin_kmers.size();

            // TODO: Store the bin (Design a DB for it, partition based or hashtable based?)
#pragma omp critical (fout)
            {
                for (const auto &kc: bin_kmers) {
                    output_file << ">" << (int) kc.count << "\n" << asString64(K, kc.kmer) << "\n";
                }
            }
        }



#pragma omp critical (cout)
        {
            std::cout << i << "\t" << expanded_kmers << "\t" << n_bin_kmers << std::endl;
        }
        // Total kmers counted after filtering
        total_counter_kmers+=n_bin_kmers;

    }
    std::cout << "Total counted kmers = " << total_counter_kmers << std::endl;
    std::cout << "Total kmers in bins = " << total_kmers_in_bins << std::endl;
    std::cout << "Total kmers expanded = " << total_kmers_expanded << std::endl;
}

void KmerCounter::read_file_large_k(const PairedReadsDatastore &ds, BinBufferWriterQueue &wrt_queue, int thread_id) {
    wrt_queue.num_writers++;
    std::vector<ExtendedKmerBin> bins;
    for (int bin=0; bin <= n_bins; bin++) {
        bins.emplace_back(bin, K, wrt_queue, 64*1024);
    }
    auto bpsg=BufferedPairedSequenceGetter(ds,128*1024,ds.readsize*2+2);
    StringKMerFactory skf(K);
    std::vector<Seq2bit> this_thread_reads_2bit;
    auto first_read = 1+(thread_id*ds.size()/num_datastore_readers);
    auto last_read = 1+( (thread_id+1)*ds.size()/num_datastore_readers);
    uint64_t rid;
    for (rid = first_read; rid < last_read and rid <= ds.size(); ++rid) {
        auto read_sequence = bpsg.get_read_sequence(rid);
        Seq2bit read2bit(read_sequence);
        // Produce the minimiser based super-kmers
        make_minimiser_send_to_bin(read_sequence, bins, read2bit);
        if (rid-first_read%100000) this_thread_reads_2bit.clear();
    }
    // FLUSH the buffers too!!
    for (int i = 0; i < bins.size(); i++) {
        if (!bins[i].empty())
            bins[i].flush();
        total_kmers_produced += bins[i].total_kmers;
        total_superkmers += bins[i].total_super_kmers;
    }
    total_reads_kmerised+= last_read-first_read;
    wrt_queue.num_writers--;
}

void KmerCounter::processReadsLargeK(std::vector<PairedReadsDatastore> &read_files) {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    calculateDistribution(read_files);

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double > time_span = std::chrono::duration_cast<std::chrono::duration<double >>(t2 - t1);
    std::cout <<"Estimating distribution time elapsed= " << time_span.count() << std::endl;

    t1=t2;

    BinBufferWriterQueue wrt_queue(n_bins + 1, loadParts);

    createPartitions(wrt_queue, read_files); // Add debugging option to keep the files that were created

    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double >>(t2 - t1);
    std::cout <<"Bins to disk time elapsed= " << time_span.count() << std::endl;

    t1=t2;
    sort_and_store_partitions(wrt_queue);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double >>(t2 - t1);
    std::cout <<"Load/Sort bins time elapsed= " << time_span.count() << std::endl;

    std::cout << "DONE FOR NOW!" << std::endl;
}

uint64_t KmerCounter::expand64(std::vector<unsigned char> &skmers_from_disk, std::vector<KmerCount64> &kmerCounts_to_sort) {
    uint64_t pos = 0;
    // TODO: Adjust for variable K size KMER<SIZE> kmer_can;
    typeof(kmerCounts_to_sort[0]) kmer_can;

    uint32_t kmer_bytes = (K + 3) / 4;
    uint32_t kmer_len_shift = sizeof(kmer_can.kmer)*8 - K*2;
    typeof(kmer_can.kmer) kmer_mask = ( ((typeof(kmer_can.kmer))1<<(2*K)) - 1);
    uint64_t input_pos = 0;

    unsigned char additional_symbols;

    unsigned char symb;
    while (pos < skmers_from_disk.size())
    {
        typeof(kmer_can) kmer(0);
        typeof(kmer_can) rev_kmer(0);
        additional_symbols = skmers_from_disk[pos++];

        for (uint32_t i = 0; i < kmer_bytes; ++i) // Read fw and rv
        {
            symb =  skmers_from_disk[pos+i];
            kmer.kmer |= (typeof(kmer_can.kmer)) symb << (((sizeof(kmer_can.kmer)-1)-i)*8);      // Clear some space for 4 nts
            rev_kmer.kmer |= (typeof(kmer_can.kmer)) revcomp_2DNA2()(symb) << (i*8); // Add 4 nts on the top
        }
        if (kmer_len_shift)
            kmer.kmer >>= kmer_len_shift;

        rev_kmer.kmer &= kmer_mask; // Clean up the head
        kmer.kmer &= kmer_mask;     // Clean up the head

        pos += kmer_bytes;
        unsigned char byte_shift = 6 - (K % 4) * 2;
        if (byte_shift != 6)
            --pos;

        kmer_can = kmer < rev_kmer ? kmer : rev_kmer;
        kmerCounts_to_sort[input_pos] = kmer_can;
        input_pos++;

        for (int i = 0; i < additional_symbols; ++i)
        {
            symb = (skmers_from_disk[pos] >> byte_shift) & 3;
            if (byte_shift == 0)
            {
                ++pos;
                byte_shift = 6;
            }
            else
                byte_shift -= 2;

            kmer.kmer = kmer.kmer << 2;
            rev_kmer.kmer = rev_kmer.kmer >> 2;

            kmer.kmer |= symb;
            kmer.kmer &= kmer_mask;

            rev_kmer.kmer |= (typeof(kmer_can.kmer))(3-symb)<<((2*K)-2);
            rev_kmer.kmer &= kmer_mask;

            kmer_can = kmer < rev_kmer ? kmer : rev_kmer;
            kmerCounts_to_sort[input_pos] = kmer_can;
            input_pos++;
        }
        if (byte_shift != 6)
            ++pos;
    }
    return input_pos;
}

uint64_t KmerCounter::expand(std::vector<unsigned char> &skmers_from_disk, std::vector<KmerCount> &kmerCounts_to_sort) {
    uint64_t pos = 0;
    // TODO: Adjust for variable K size KMER<SIZE> kmer_can;
    typeof(kmerCounts_to_sort[0]) kmer_can;

    uint32_t kmer_bytes = (K + 3) / 4;
    uint32_t kmer_len_shift = sizeof(kmer_can.kmer)*8 - K*2;
    typeof(kmer_can.kmer) kmer_mask = ( ((typeof(kmer_can.kmer))1<<(2*K)) - 1);
    uint64_t input_pos = 0;

    unsigned char additional_symbols;

    unsigned char symb;
    while (pos < skmers_from_disk.size())
    {
        typeof(kmer_can) kmer(0);
        typeof(kmer_can) rev_kmer(0);
        additional_symbols = skmers_from_disk[pos++];

        for (uint32_t i = 0; i < kmer_bytes; ++i) // Read fw and rv
        {
            symb =  skmers_from_disk[pos+i];
            kmer.kmer |= (typeof(kmer_can.kmer)) symb << (((sizeof(kmer_can.kmer)-1)-i)*8);      // Clear some space for 4 nts
            rev_kmer.kmer |= (typeof(kmer_can.kmer)) revcomp_2DNA2()(symb) << (i*8); // Add 4 nts on the top
        }
        if (kmer_len_shift)
            kmer.kmer >>= kmer_len_shift;

        rev_kmer.kmer &= kmer_mask; // Clean up the head
        kmer.kmer &= kmer_mask;     // Clean up the head

        pos += kmer_bytes;
        unsigned char byte_shift = 6 - (K % 4) * 2;
        if (byte_shift != 6)
            --pos;

        kmer_can = kmer < rev_kmer ? kmer : rev_kmer;
        kmerCounts_to_sort[input_pos] = kmer_can;
        input_pos++;

        for (int i = 0; i < additional_symbols; ++i)
        {
            symb = (skmers_from_disk[pos] >> byte_shift) & 3;
            if (byte_shift == 0)
            {
                ++pos;
                byte_shift = 6;
            }
            else
                byte_shift -= 2;

            kmer.kmer = kmer.kmer << 2;
            rev_kmer.kmer = rev_kmer.kmer >> 2;

            kmer.kmer |= symb;
            kmer.kmer &= kmer_mask;

            rev_kmer.kmer |= (typeof(kmer_can.kmer))(3-symb)<<((2*K)-2);
            rev_kmer.kmer &= kmer_mask;

            kmer_can = kmer < rev_kmer ? kmer : rev_kmer;
            kmerCounts_to_sort[input_pos] = kmer_can;
            input_pos++;
        }
        if (byte_shift != 6)
            ++pos;
    }
    return input_pos;
}

void KmerCounter::mmer_table_initialiser(const PairedReadsDatastore &ds, MinimiserTable &minimiserTable) {
    auto bpsg=BufferedPairedSequenceGetter(ds,4*1024,ds.readsize*2+2);
    std::default_random_engine generator;
    std::uniform_int_distribution<uint32_t> distribution(1,ds.size());
    MinimiserTable local_minimiser_table(signature_len, n_bins);
    uint64_t total_mmers_seen{0};
    for (unsigned int r = 1; r < std::min(ds.size()*0.01,100000.); r++) {
        auto seq = Seq2bit(std::string(bpsg.get_read_sequence(r)));

        Mmer current_mmer(signature_len), end_mmer(signature_len);

        uint32_t i = 0;
        uint32_t len = 0;//length of extended kmer

        while (i + K - 1 < seq.size())
        {
            bool contains_N = false;
            //building first signature after 'N' or at the read start
            for (uint32_t j = 0; j < signature_len; ++j, ++i)
                if (seq[i] > 3)//'N'
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
            current_mmer.insert(seq.data() + signature_start_pos);
            end_mmer.set(current_mmer);
            for (; i < seq.size(); ++i)
            {
                if (seq[i] > 3)//'N'
                {
                    if (len >= K) {
                        total_mmers_seen+= 1+len-K;
                        local_minimiser_table[current_mmer.get()] += 1 + len - K;
                    }
                    len = 0;
                    ++i;
                    break;
                }
                end_mmer.insert(seq[i]);
                if (end_mmer < current_mmer)//signature at the end of current k-mer is lower than current
                {
                    if (len >= K)
                    {
                        total_mmers_seen+= 1+len-K;
                        local_minimiser_table[current_mmer.get()] += 1 + len - K;
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
                    total_mmers_seen+= 1+len-K;
                    local_minimiser_table[current_mmer.get()] += 1 + len - K;
                    len = K - 1;
                    //looking for new signature
                    ++signature_start_pos;
                    //building first signature in current k-mer
                    end_mmer.insert(seq.data() + signature_start_pos);
                    current_mmer.set(end_mmer);
                    for (uint32_t j = signature_start_pos + signature_len; j <= i; ++j)
                    {
                        end_mmer.insert(seq[j]);
                        if (end_mmer <= current_mmer)
                        {
                            current_mmer.set(end_mmer);
                            signature_start_pos = j - signature_len + 1;
                        }
                    }
                }
                ++len;
            }
        }
        if (len >= K) {//last one in read
            total_mmers_seen+= 1+len-K;
            local_minimiser_table[current_mmer.get()] += 1 + len - K;
        }
    }

    {
        std::lock_guard<std::mutex> lg(mmer_table_lock);
        minimiserTable.combine(local_minimiser_table);
    }

}

void KmerCounter::make_minimiser_send_to_bin(std::string in_seq, std::vector<ExtendedKmerBin> &bins,
                                             Seq2bit &read_2bit) {
    int i = 0;
    int len = 0;
    unsigned int bin_no = 0;
    Mmer current_mmer(signature_len), end_mmer(signature_len);
    auto& seq = read_2bit;
    while (in_seq[i]!='\0' and i + K - 1 < in_seq.size())
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
        current_mmer.insert(seq.data()+signature_start_pos);
        end_mmer.set(current_mmer);
        for (; i < in_seq.size(); ++i)
        {
            if (seq[i] > 3)//'N'
            {
                if (len >= K)
                {
                    bin_no = mmer_table.get_bin_id(current_mmer.get());

                    bins[bin_no].store_superkmer(seq,i - len, len);
                }
                len = 0;
                ++i;
                break;
            }
            end_mmer.insert(seq.data()[i]);
            if (end_mmer < current_mmer)//signature at the end of current k-mer is lower than current
            {
                if (len >= K)
                {
                    bin_no = mmer_table.get_bin_id(current_mmer.get());

                    bins[bin_no].store_superkmer(seq, i - len, len);

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
                bins[bin_no].store_superkmer(seq, i - len, len);

                len = K - 1;
                //looking for new signature
                ++signature_start_pos;
                //building first signature in current k-mer
                end_mmer.insert(seq.data()+signature_start_pos);
                current_mmer.set(end_mmer);
                for (uint32_t j = signature_start_pos + signature_len; j <= i; ++j)
                {
                    end_mmer.insert(seq.data()[j]);
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

                bins[bin_no].store_superkmer(seq, i + 1 - len, len);

                i -= K - 2;
                len = 0;
                break;
            }

        }
    }
    if (len >= K)//last one in read
    {
        bin_no = mmer_table.get_bin_id(current_mmer.get());

        bins[bin_no].store_superkmer(seq, i - len, len);
    }
}

constexpr std::array<char,256> Seq2bit::codes;
uint32_t Mmer::norm5[];
uint32_t Mmer::norm6[];
uint32_t Mmer::norm7[];
uint32_t Mmer::norm8[];
uint32_t Mmer::norm9[];
uint32_t Mmer::norm10[];
uint32_t Mmer::norm11[];

Mmer::_si Mmer::_init;
