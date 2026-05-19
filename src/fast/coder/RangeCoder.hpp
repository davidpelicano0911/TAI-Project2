/*
 * Byte-level range coder (LZMA-style carry propagation).
 *
 * Same interface: write(FrequencyTable&, symbol) / read(FrequencyTable&).
 * Renormalises every 8 bits.
 *
 * Invariant: range_ always in [RC_TOP, 0xFFFFFFFF] after each renorm.
 * RC_TOP = 2^24 = 16 777 216.  Works correctly as long as FrequencyTable
 * totals stay below RC_TOP
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <istream>
#include <ostream>
#include "FrequencyTable.hpp"

static constexpr uint32_t RC_TOP = 1u << 24;  // renorm threshold

// ── Range Encoder ─────────────────────────────────────────────────────────────
//
// Uses LZMA-style carry propagation:
//   low_    — 64-bit; lower 32 bits are the actual range-coder state,
//             bit 32 is the carry flag set by  low_ += sym_low * r.
//   cache_  — last byte held back (always in 0x00..0xFE; carry can safely
//             increment it to 0xFF without further propagation).
//   pending_— count of 0xFF bytes buffered after cache_ (become 0x00 on carry).
//
class RangeEncoder {
    std::ostream& out_;
    uint64_t low_;
    uint32_t range_;
    uint8_t cache_;
    uint32_t pending_;

    void shift_low() {
        bool carry = (low_ >> 32) != 0;
        uint8_t top = static_cast<uint8_t>(static_cast<uint32_t>(low_) >> 24);

        if (top < 0xFFu || carry) {
            out_.put(static_cast<char>(cache_ + (carry ? 1u : 0u)));
            const uint8_t fill = carry ? 0x00u : 0xFFu;
            for (uint32_t i = 0; i < pending_; ++i)
                out_.put(static_cast<char>(fill));
            pending_ = 0;
            cache_ = top;
        } else {
            ++pending_;
        }
        // Shift left by 8, discarding carry bit (keeps only lower 32 bits).
        low_ = static_cast<uint64_t>(static_cast<uint32_t>(low_) << 8);
    }

public:
    explicit RangeEncoder(std::ostream& o)
        : out_(o), low_(0), range_(0xFFFFFFFFu), cache_(0), pending_(0) {}

    void write(const FrequencyTable& freqs, uint32_t symbol) {
        uint32_t total = freqs.getTotal();
        uint32_t sym_low = freqs.getLow(symbol);
        uint32_t sym_high = freqs.getHigh(symbol);

        uint32_t r = range_ / total;
        low_ += static_cast<uint64_t>(sym_low) * r;
        range_ = (sym_high < total) ? (sym_high - sym_low) * r : range_ - sym_low * r;

        while (range_ < RC_TOP) {
            range_ <<= 8;
            shift_low();
        }
    }

    // Flush all buffered state. Must be called after the last symbol.
    void finish() {
        for (int i = 0; i < 5; ++i)
            shift_low();
    }
};


// ── Range Decoder ─────────────────────────────────────────────────────────────
//
// Pre-reads 5 bytes to match the encoder's initial cache_ byte + 4-byte window.
//
class RangeDecoder {
    std::istream& in_;
    uint32_t range_;
    uint32_t code_;

    uint8_t read_byte() {
        int b = in_.get();
        return b >= 0 ? static_cast<uint8_t>(b) : 0u;
    }

public:
    explicit RangeDecoder(std::istream& i) : in_(i), range_(0xFFFFFFFFu), code_(0) {
        for (int k = 0; k < 5; ++k)
            code_ = (code_ << 8) | read_byte();
    }

    uint32_t read(const FrequencyTable& freqs) {
        uint32_t total = freqs.getTotal();
        uint32_t r = range_ / total;
        uint32_t value = std::min(code_ / r, total - 1u);

        uint32_t symbol = freqs.findSymbol(value);
        uint32_t sym_low = freqs.getLow(symbol);
        uint32_t sym_high = freqs.getHigh(symbol);

        code_ -= sym_low * r;
        range_ = (sym_high < total) ? (sym_high - sym_low) * r : range_ - sym_low * r;

        while (range_ < RC_TOP) {
            code_ = (code_ << 8) | read_byte();
            range_ <<= 8;
        }
        return symbol;
    }
};
