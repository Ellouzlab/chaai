#pragma once

#include <cstdint>

static constexpr uint16_t MIN_SEG_LEN_DEFAULT = 30;
static constexpr int RUN_HIST_MAX = 4096;
static constexpr double PSI_UNSCORED = -100.0;

struct AASegment
{
    uint64_t nuc_start;
    uint16_t aa_len;
    uint8_t frame;
};
static_assert(sizeof(AASegment) == 16, "AASegment layout");
