#include <queue>
#include <mutex>
#include <condition_variable>

template<typename T>
class thread_safe_queue {
private:
    std::queue<T> queue;
    std::mutex m;
    std::condition_variable cv;

public:
    void push(T value) {
        std::lock_guard<std::mutex> lock(m);
        queue.push(std::move(value));
        cv.notify_one();
    }

    bool try_pop(T& result) {
        std::lock_guard<std::mutex> lock(m);
        if (queue.empty()) return false;
        result = std::move(queue.front());
        queue.pop();
        return true;
    }

    bool empty() const {
        std::lock_guard<std::mutex> lock(m);
        return queue.empty();
    }
};
