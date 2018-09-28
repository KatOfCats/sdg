//
// Created by Luis Yanes (EI) on 21/09/2018.
//

#include "ExtendedKmerBin.hpp"

void BinBufferWriterQueue::insert(int bin_no, std::vector<unsigned char> &buf, size_t len, uint32_t n_superkmers,
                                  uint32_t n_kmers) {
    std::lock_guard<std::mutex> lg(bin_mtx[bin_no]);
    if (bin_buffer_size[bin_no] > max_bin_size) {
        max_bin_size = bin_buffer_size[bin_no];
        max_bin = bin_no;
    }
    bin_buffer_size[bin_no]+=len;
    std::vector<unsigned char> copy(buf);
    bin_buffers[bin_no].push_back(copy);
}

void BinBufferWriterQueue::write_bin(int bin_no) {
    std::lock_guard<std::mutex> lg(bin_mtx[bin_no]);
    std::vector<unsigned char> buffer;
    while (bin_no != 0 and !bin_buffers[bin_no].empty()) {
        const auto &bbref(bin_buffers[bin_no].front());
        std::copy(bbref.begin(), bbref.end(), std::back_inserter(buffer));
        bin_buffer_size[bin_no] -= buffer.size();
        bin_buffers[bin_no].pop_front();
        std::copy(buffer.begin(), buffer.end(), std::ostream_iterator<unsigned char>(bin_files[bin_no]));
    }

    bin_buffer_size[bin_no] = 0;

    max_bin_size = bin_buffer_size[0];
    max_bin = 0;

    for (int i = 1; i < bin_buffers.size(); i++) {
        if (bin_buffer_size[i] > max_bin_size) {
            max_bin_size = bin_buffer_size[i];
            max_bin = i;
        }
    }
}

bool BinBufferWriterQueue::empty() {
    if (num_writers > 0) return false;
    if (max_bin>0) return false;
    return true;
}

BinBufferWriterQueue::BinBufferWriterQueue(int num_bins) {
    std::vector<std::mutex> list(num_bins);
    bin_mtx.swap(list);
    bin_buffers.resize(num_bins);
    for (int i = 0; i < num_bins; i++) {
        std::string name("file"+std::to_string(i)+".bin");
        bin_files.emplace_back(name);
        bin_buffer_size.emplace_back(0);
    }
}

BinBufferWriterQueue::~BinBufferWriterQueue() {
    // Make sure all bins are empty before closing the files and destroying all objects
    for (int i = 0; i < bin_buffers.size(); i++) {
        write_bin(i);
    }
}

ExtendedKmerBin::ExtendedKmerBin(uint16_t bin_no, uint8_t k, BinBufferWriterQueue &queue, unsigned int max_size) :
        K(k),
        queue(queue),
        bin_no(bin_no)
{
    buffer_size = max_size;
    buffer.resize(max_size);
}

void ExtendedKmerBin::store_superkmer(const Seq2bit &seq, unsigned int offset, unsigned int len) {

    uint32_t bytes = 1 + (len + 3) / 4;
    if(buffer_pos + bytes > buffer_size)
    {
        flush();
//            pmm_bins->reserve(buffer); // Here synchronizes for memory
        buffer_pos		= 0;
        n_kmers			= 0;
        n_super_kmers		= 0;
    }


    buffer[buffer_pos++] = static_cast<unsigned char>(len - K); // TODO: Check allocation
    for(uint32_t i = 0, j = 0 ; i < len / 4 ; ++i,j+=4)
        buffer[buffer_pos++] = (seq[j] << 6) + (seq[j + 1] << 4) + (seq[j + 2] << 2) + seq[j + 3];
    switch (len%4)
    {
        case 1:
            buffer[buffer_pos++] = (seq[len-1] << 6);
            break;
        case 2:
            buffer[buffer_pos++] = (seq[len-2] << 6) + (seq[len-1] << 4);
            break;
        case 3:
            buffer[buffer_pos++] = (seq[len-3] << 6) + (seq[len-2] << 4) + (seq[len-1] << 2);
            break;
        default:break;
    }

    ++n_super_kmers;
    n_kmers += len - K + 1;
}

void ExtendedKmerBin::flush() {
    // Write buffer to output queue (sync)
    queue.insert(bin_no, buffer,buffer_pos, n_super_kmers, n_kmers);
}
