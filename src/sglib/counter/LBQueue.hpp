//
// Created by Luis Yanes (EI) on 17/09/2018.
//

#ifndef BSG_LBQUEUE_HPP
#define BSG_LBQUEUE_HPP

#include <list>
#include <tuple>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <array>
#include <atomic>

class KmerArrayQueue {
    typedef std::array<std::pair<bool, uint64_t>, 300> elem_t;
    typedef std::queue <elem_t, std::list<elem_t>> queue_t;
    queue_t q;
    std::atomic<int> n_writers{0};
    mutable std::mutex mtx;                        // The mutex to synchronise on
    std::condition_variable cv_queue_empty;

public:
    explicit KmerArrayQueue(int _n_writers) {
        std::lock_guard <std::mutex> lck(mtx);

        n_writers = _n_writers;
    }

    ~KmerArrayQueue() {}

    bool empty() {
        std::lock_guard<std::mutex> lck(mtx);
        return q.empty();
    }
    bool completed() {
        std::lock_guard<std::mutex> lck(mtx);
        return q.empty() && !n_writers;
    }
    void mark_completed() {
        std::lock_guard<std::mutex> lck(mtx);
        n_writers--;
        if(!n_writers)
            cv_queue_empty.notify_all();
    }
    void push(elem_t &elem) {
        std::unique_lock<std::mutex> lck(mtx);

        bool was_empty = q.empty();
        q.push(std::move(elem));

        if(was_empty)
            cv_queue_empty.notify_all();
    }
    bool pop(elem_t &elem) {
        std::unique_lock<std::mutex> lck(mtx);
        cv_queue_empty.wait(lck, [this]{return !q.empty() || !n_writers;});

        if(q.empty())
            return false;
        elem = std::move(q.front());
        q.pop();

        return true;
    }
};
#endif //BSG_LBQUEUE_HPP
