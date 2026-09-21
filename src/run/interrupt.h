#pragma once

#include <atomic>
#include <csignal>

struct Interrupted {};

inline std::atomic<bool> &interrupt_flag()
{
    static std::atomic<bool> flag{false};
    return flag;
}

inline bool interrupted() { return interrupt_flag().load(std::memory_order_relaxed); }

class InterruptScope
{
public:
    InterruptScope()
    {
        interrupt_flag().store(false);
        prev_ = std::signal(SIGINT, [](int) { interrupt_flag().store(true, std::memory_order_relaxed); });
    }
    ~InterruptScope() { std::signal(SIGINT, prev_); }
    InterruptScope(const InterruptScope &) = delete;
    InterruptScope &operator=(const InterruptScope &) = delete;

private:
    void (*prev_)(int);
};
