//
// Created by Luis Yanes (EI) on 13/09/2018.
//

#include <catch.hpp>
#include <sglib/counter/MPMCQueue.hpp>
#include <xmmintrin.h>
#include <thread>
#include <iostream>
#include <unistd.h>

size_t const thread_count = 4;
size_t const batch_size = 1;
size_t const iter_count = 20000000;

bool volatile g_start = 0;

typedef MPMCBoundedQueue<int> queue_t;

std::function<unsigned()> enqueue_thread(void* ctx, unsigned int id)
{
    queue_t& queue = *(queue_t*)ctx;
    int data;

    srand((unsigned)time(0) + id);
    size_t pause = (size_t)rand() % 1000;

    while (g_start == 0) {
        std::this_thread::yield();
    }

    for (size_t i = 0; i != pause; i += 1) {
        _mm_pause();
    }

    for (int iter = 0; iter != iter_count; ++iter) {
        for (size_t i = 0; i != batch_size; i += 1) {
            while (!queue.enqueue((int) i))
                std::this_thread::yield();
        }
    }
    return 0;
}
std::function<unsigned()> dequeue_thread(void* ctx, unsigned int id) {
    queue_t& queue = *(queue_t*)ctx;
    int data;

    srand((unsigned)time(0) + id);
    size_t pause = (size_t)rand() % 1000;

    while (g_start == 0) {
        std::this_thread::yield();
    }

    for (size_t i = 0; i != pause; i += 1) {
        _mm_pause();
    }

    for (int iter = 0; iter != iter_count; ++iter) {
        for (size_t i = 0; i != batch_size; i += 1) {
            while (!queue.dequeue(data))
                std::this_thread::yield();
        }
    }

        return 0;
}


TEST_CASE("Test MPMCBoundedQueue1024") {
    queue_t queue (1024);

    std::vector<std::thread> threads;
    for (int i = 0; i != thread_count/2; ++i)
    {
        std::thread th(enqueue_thread, &queue, i+1);
        threads.push_back(std::move(th));
    }
    for (int i = 0; i != thread_count/2; ++i) {
        std::thread th(dequeue_thread, &queue, (i*2)+1);
        threads.push_back(std::move(th));
    }
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    usleep(10000);
    g_start=1;
    for (auto & th : threads)
    {
        if (th.joinable())
            th.join();
    }


    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double > time_span = std::chrono::duration_cast<std::chrono::duration<double >>(t2 - t1);

    std::cout <<"Total time elapsed= " << time_span.count() << "\ncycles/op=" << time_span.count()/(batch_size * iter_count * 2 * thread_count) << std::endl;
}
