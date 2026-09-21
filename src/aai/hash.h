#pragma once

#include <cstdint>

static inline uint64_t hash64(uint64_t k)
{
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

static inline uint64_t fmh_threshold(uint32_t scaled) {
    return UINT64_MAX / scaled;
}
