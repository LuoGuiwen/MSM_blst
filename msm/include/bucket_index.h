#ifndef MSM_BUCKET_INDEX_H
#define MSM_BUCKET_INDEX_H

#include <cstddef>
#include <cstdint>

inline size_t bucket_rank_table_word_count(size_t max_bucket_value)
{
    return (max_bucket_value + 64) / 64;
}

inline uint32_t bucket_rank_index(uint32_t bucket_value,
                                  const bucket_rank_word table[])
{
    const bucket_rank_word& entry = table[bucket_value >> 6];
    const uint32_t bit = bucket_value & 63;
    const uint64_t lower_mask = (uint64_t(1) << bit) - 1;
    return entry.rank_before
         + uint32_t(__builtin_popcountll(entry.bits & lower_mask));
}

inline void construct_bucket_rank_table(bucket_rank_word table[],
                                        size_t word_count,
                                        const int bucket_set[],
                                        size_t bucket_set_size)
{
    for (size_t i = 0; i < bucket_set_size; ++i) {
        const uint32_t value = uint32_t(bucket_set[i]);
        table[value >> 6].bits |= uint64_t(1) << (value & 63);
    }

    uint32_t rank = 0;
    for (size_t i = 0; i < word_count; ++i) {
        table[i].rank_before = rank;
        rank += uint32_t(__builtin_popcountll(table[i].bits));
    }
}

inline bool check_bucket_rank_table(const bucket_rank_word table[],
                                    const int bucket_set[],
                                    size_t bucket_set_size)
{
    for (size_t i = 0; i < bucket_set_size; ++i) {
        if (bucket_rank_index(uint32_t(bucket_set[i]), table) != i)
            return false;
    }
    return true;
}

#endif /* MSM_BUCKET_INDEX_H */
