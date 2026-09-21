#pragma once

#include <cstdint>
#include <cstddef>
#include <algorithm>
#if defined(__AVX2__)
#include <immintrin.h>
#endif

struct NucLookup
{
    uint8_t encode[256];
    uint8_t complement[256];

    constexpr NucLookup() : encode{}, complement{}
    {
        for (int i = 0; i < 256; i++) {
            encode[i] = 255;
            complement[i] = 255;
        }
        encode['A'] = encode['a'] = 0;
        encode['C'] = encode['c'] = 1;
        encode['G'] = encode['g'] = 2;
        encode['T'] = encode['t'] = encode['U'] = encode['u'] = 3;

        complement['A'] = complement['a'] = 3;
        complement['C'] = complement['c'] = 2;
        complement['G'] = complement['g'] = 1;
        complement['T'] = complement['t'] = complement['U'] = complement['u'] = 0;
    }
};

static constexpr NucLookup NUC = NucLookup();

static constexpr uint8_t ENC_N = 4;
static constexpr uint8_t ENC_COMP[5] = {3, 2, 1, 0, 4};
static inline uint8_t enc_of(char c) { const uint8_t e = NUC.encode[static_cast<uint8_t>(c)]; return e <= 3 ? e : ENC_N; }

static inline void encode_bytes(const char *src, char *dst, size_t n)
{
    size_t i = 0;
#if defined(__AVX2__)
    const __m256i up = _mm256_set1_epi8(static_cast<char>(0xDF));
    const __m256i cA = _mm256_set1_epi8('A'), cC = _mm256_set1_epi8('C'), cG = _mm256_set1_epi8('G');
    const __m256i cT = _mm256_set1_epi8('T'), cU = _mm256_set1_epi8('U');
    const __m256i m3 = _mm256_set1_epi8(0x03), m1 = _mm256_set1_epi8(0x01), four = _mm256_set1_epi8(4);
    for (; i + 32 <= n; i += 32) {
        const __m256i v = _mm256_and_si256(_mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + i)), up);
        const __m256i ok = _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi8(v, cA), _mm256_cmpeq_epi8(v, cC)),
                           _mm256_or_si256(_mm256_cmpeq_epi8(v, cG), _mm256_or_si256(_mm256_cmpeq_epi8(v, cT), _mm256_cmpeq_epi8(v, cU))));
        const __m256i t1 = _mm256_and_si256(_mm256_srli_epi16(v, 1), m3);
        const __m256i t2 = _mm256_and_si256(_mm256_srli_epi16(v, 2), m1);
        const __m256i code = _mm256_xor_si256(t1, t2);
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(dst + i), _mm256_blendv_epi8(four, code, ok));
    }
#endif
    for (; i < n; i++) dst[i] = static_cast<char>(enc_of(src[i]));
}

static constexpr uint8_t CODON_TABLE[64] = {
    11,  2, 11,  2,   16, 16, 16, 16,    1, 15,  1, 15,    9,  9, 12,  9,
     5,  8,  5,  8,   14, 14, 14, 14,    1,  1,  1,  1,   10, 10, 10, 10,
     6,  3,  6,  3,    0,  0,  0,  0,    7,  7,  7,  7,   19, 19, 19, 19,
    20, 18, 20, 18,   15, 15, 15, 15,   20,  4, 17,  4,   10, 13, 10, 13
};
