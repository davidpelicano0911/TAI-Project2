/*
 * Reference arithmetic coding
 *
 * Copyright (c) Project Nayuki
 * MIT License. See readme file.
 * https://www.nayuki.io/page/reference-arithmetic-coding
 */

#pragma once

#include <cstdint>
#include <vector>

/*
 * A table of symbol frequencies indexed from 0 to getSymbolLimit()-1.
 * Supports any alphabet size set at construction time.
 */
class FrequencyTable {

	public: virtual ~FrequencyTable() = 0;

	public: virtual std::uint32_t getSymbolLimit() const = 0;
	public: virtual std::uint32_t get(std::uint32_t symbol) const = 0;
	public: virtual void set(std::uint32_t symbol, std::uint32_t freq) = 0;
	public: virtual void increment(std::uint32_t symbol) = 0;
	public: virtual std::uint32_t getTotal() const = 0;
	public: virtual std::uint32_t getLow(std::uint32_t symbol) const = 0;
	public: virtual std::uint32_t getHigh(std::uint32_t symbol) const = 0;

	// Find symbol s such that getLow(s) <= value < getHigh(s).
	public: virtual std::uint32_t findSymbol(std::uint32_t value) const;

};

/*
 * Flat (uniform) frequency table — each symbol has frequency 1.
 * Useful as a starting point for adaptive coding.
 */
class FlatFrequencyTable final : public FrequencyTable {

	private: std::uint32_t numSymbols;

	public: explicit FlatFrequencyTable(std::uint32_t numSyms);

	public: std::uint32_t getSymbolLimit() const override;
	public: std::uint32_t get(std::uint32_t symbol) const override;
	public: std::uint32_t getTotal() const override;
	public: std::uint32_t getLow(std::uint32_t symbol) const override;
	public: std::uint32_t getHigh(std::uint32_t symbol) const override;
	public: void set(std::uint32_t symbol, std::uint32_t freq) override;
	public: void increment(std::uint32_t symbol) override;
	private: void checkSymbol(std::uint32_t symbol) const;

};

/*
 * Mutable frequency table backed by a vector. Cumulative sums are computed lazily.
 * Pass any std::vector<uint32_t> of frequencies to the constructor — the alphabet
 * size equals the vector length, which can be anything from 1 to 2^32-2.
 */
class SimpleFrequencyTable final : public FrequencyTable {

	private: std::vector<std::uint32_t> frequencies;
	private: mutable std::vector<std::uint32_t> cumulative;
	private: std::uint32_t total;

	public: explicit SimpleFrequencyTable(const std::vector<std::uint32_t> &freqs);
	public: explicit SimpleFrequencyTable(const FrequencyTable &freqs);

	public: std::uint32_t getSymbolLimit() const override;
	public: std::uint32_t get(std::uint32_t symbol) const override;
	public: void set(std::uint32_t symbol, std::uint32_t freq) override;
	public: void increment(std::uint32_t symbol) override;
	public: std::uint32_t getTotal() const override;
	public: std::uint32_t getLow(std::uint32_t symbol) const override;
	public: std::uint32_t getHigh(std::uint32_t symbol) const override;

	private: void initCumulative(bool checkTotal = true) const;
	private: static std::uint32_t checkedAdd(std::uint32_t x, std::uint32_t y);

};
