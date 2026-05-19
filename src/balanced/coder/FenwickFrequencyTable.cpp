#include <stdexcept>
#include "FenwickFrequencyTable.hpp"

using std::uint32_t;
using std::int32_t;

FenwickFrequencyTable::FenwickFrequencyTable(uint32_t numSyms)
    : symbolLimit(numSyms),
      highBit_(1u),
      tree(numSyms + 1, 0u),   // 1-indexed; index 0 is unused sentinel
      total(0u)
{
    if (numSyms < 1)
        throw std::domain_error("Number of symbols must be positive");
    // Compute highest power of 2 that is <= symbolLimit
    while (highBit_ * 2u <= symbolLimit)
        highBit_ *= 2u;
}

uint32_t FenwickFrequencyTable::getSymbolLimit() const {
    return symbolLimit;
}

uint32_t FenwickFrequencyTable::get(uint32_t symbol) const {
    checkSymbol(symbol);
    // Point query as difference of two consecutive prefix sums.
    return query(symbol + 1) - query(symbol);
}

void FenwickFrequencyTable::set(uint32_t symbol, uint32_t freq) {
    checkSymbol(symbol);
    uint32_t curr = get(symbol);
    total = total - curr + freq;
    // delta may be negative when reducing a frequency.
    update(symbol + 1, static_cast<int32_t>(freq) - static_cast<int32_t>(curr));
}

void FenwickFrequencyTable::increment(uint32_t symbol) {
    checkSymbol(symbol);
    ++total;
    update(symbol + 1, 1);
}

uint32_t FenwickFrequencyTable::getTotal() const {
    return total;
}

uint32_t FenwickFrequencyTable::getLow(uint32_t symbol) const {
    checkSymbol(symbol);
    return query(symbol);        // sum of freqs[0 .. symbol-1]
}

uint32_t FenwickFrequencyTable::getHigh(uint32_t symbol) const {
    checkSymbol(symbol);
    return query(symbol + 1);    // sum of freqs[0 .. symbol]
}


// Add `delta` to 1-indexed Fenwick position `i` and propagate upward.
// Uses signed int indices so that `x & -x` (lowest-set-bit) is well-defined.
void FenwickFrequencyTable::update(uint32_t i, int32_t delta) {
    for (int x = static_cast<int>(i); x <= static_cast<int>(symbolLimit); x += x & -x)
        tree[static_cast<uint32_t>(x)] += static_cast<uint32_t>(delta);
}

// Prefix sum of Fenwick positions 1..i (returns 0 when i == 0).
uint32_t FenwickFrequencyTable::query(uint32_t i) const {
    uint32_t sum = 0;
    for (int x = static_cast<int>(i); x > 0; x -= x & -x)
        sum += tree[static_cast<uint32_t>(x)];
    return sum;
}

void FenwickFrequencyTable::checkSymbol(uint32_t symbol) const {
    if (symbol >= symbolLimit)
        throw std::domain_error("Symbol out of range");
}


// O(log N) Fenwick descent: finds symbol s with getLow(s) <= value < getHigh(s).
//
// The BIT stores prefix sums; we descend bit by bit from the highest bit,
// greedily extending our position whenever the subtree sum fits within `rem`.
// After the loop, `pos` equals the 0-indexed target symbol.
uint32_t FenwickFrequencyTable::findSymbol(uint32_t value) const {
    uint32_t pos = 0;
    uint32_t rem = value;
    for (uint32_t mask = highBit_; mask > 0u; mask >>= 1) {
        uint32_t next = pos + mask;
        if (next <= symbolLimit && tree[next] <= rem) {
            pos  = next;
            rem -= tree[pos];
        }
    }
    return pos;
}
