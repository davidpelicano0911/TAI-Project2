/*
 * A Fenwick tree (Binary Indexed Tree) gives O(log N) for every operation:
 *   increment(s)  → one point-update  : O(log N)
 *   getLow(s)     → one prefix query  : O(log N)
 *   getHigh(s)    → one prefix query  : O(log N)
 *   getTotal()    → O(1) (cached)
 *   get(s)        → two prefix queries: O(log N)
 *
 * Layout:
 *   tree[] is 1-indexed, size = symbolLimit + 1.
 *   Symbol s occupies Fenwick position s+1.
 *   query(k) = sum of freqs[0 .. k-1]  (prefix sum of the first k symbols)
 *   getLow(s)  = query(s)
 *   getHigh(s) = query(s+1)
 */

#pragma once

#include <cstdint>
#include <vector>
#include "FrequencyTable.hpp"


class FenwickFrequencyTable final : public FrequencyTable {

private:
    std::uint32_t symbolLimit;
    std::uint32_t highBit_;   
    std::vector<std::uint32_t> tree;       
    std::uint32_t total;      

public:
    explicit FenwickFrequencyTable(std::uint32_t numSyms);

public:
    std::uint32_t getSymbolLimit() const override;
    std::uint32_t get(std::uint32_t symbol) const override;
    void set(std::uint32_t symbol, std::uint32_t freq) override;
    void increment(std::uint32_t symbol) override;
    std::uint32_t getTotal() const override;
    std::uint32_t getLow (std::uint32_t symbol) const override;
    std::uint32_t getHigh(std::uint32_t symbol) const override;
    std::uint32_t findSymbol(std::uint32_t value) const override;


private:
    void update(std::uint32_t i, std::int32_t delta);

    std::uint32_t query(std::uint32_t i) const;

    void checkSymbol(std::uint32_t symbol) const;
};
