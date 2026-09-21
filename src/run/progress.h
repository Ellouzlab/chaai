#pragma once

#include <cstdio>
#include <chrono>
#include <string>

struct ProgressBar
{
    int total;
    int done = 0;
    std::string label;
    std::chrono::steady_clock::time_point start;
    static constexpr int BAR_WIDTH = 30;

    ProgressBar(int total_, const std::string &label_ = "")
        : total(total_), label(label_), start(std::chrono::steady_clock::now()) {}

    void update(int completed)
    {
        done = completed;
        auto now = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(now - start).count();

        const float frac = total > 0 ? static_cast<float>(done) / total : 1.0f;
        const int filled = static_cast<int>(frac * BAR_WIDTH);

        char eta_buf[32] = "??:??";
        if (done > 0 && done < total) {
            const double remaining = elapsed * (total - done) / done;
            const int m = static_cast<int>(remaining / 60);
            const int s = static_cast<int>(remaining) % 60;
            snprintf(eta_buf, sizeof(eta_buf), "%02d:%02d", m, s);
        } else if (done >= total) {
            const int m = static_cast<int>(elapsed / 60);
            const int s = static_cast<int>(elapsed) % 60;
            snprintf(eta_buf, sizeof(eta_buf), "%02d:%02d", m, s);
        }

        char bar[64];
        for (int i = 0; i < BAR_WIDTH; i++) {
            bar[i] = (i < filled) ? '#' : '-';
        }
        bar[BAR_WIDTH] = '\0';

        fprintf(stderr, "\r[%s] [%s] %d/%d (%d%%) ETA %s  ",
                label.c_str(), bar, done, total,
                static_cast<int>(frac * 100), eta_buf);
        fflush(stderr);
    }

    void finish()
    {
        update(total);
        fprintf(stderr, "\n");
    }
};
