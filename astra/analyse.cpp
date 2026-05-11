//
// ASTRA — step 1: understand the image + auto-select threshold T
//
// For each image we want to know:
//   - what is the background value (the most common pixel value = mode)
//   - how many pixels are "near" the background (within a small delta)
//   - how many pixels are "anomalies" (far from the background = stars, artifacts)
//   - what T to use automatically, chosen via MAD (median absolute deviation)
//
// MAD is robust: stars are outliers and don't affect it, unlike std deviation.
// We pick T = 4 * MAD, then round up to the nearest power of two for simplicity.
//
// Usage: ./analyse <image_file>
//

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

static constexpr int WIDTH  = 1500;
static constexpr int HEIGHT = 1500;
static constexpr int NPIX   = WIDTH * HEIGHT;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <image_file>\n", argv[0]);
        return 1;
    }

    // ------------------------------------------------------------------
    // 1. Read the raw image (big-endian uint16)
    // ------------------------------------------------------------------
    FILE* f = fopen(argv[1], "rb");
    if (!f) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    std::vector<uint16_t> px(NPIX);
    for (int i = 0; i < NPIX; i++) {
        uint8_t b[2];
        if (fread(b, 1, 2, f) != 2) { fprintf(stderr, "Short read\n"); return 1; }
        px[i] = (uint16_t)((b[0] << 8) | b[1]);
    }
    fclose(f);

    // ------------------------------------------------------------------
    // 2. Build a histogram of all 65536 possible pixel values
    // ------------------------------------------------------------------
    uint32_t hist[65536] = {};
    for (int i = 0; i < NPIX; i++)
        hist[px[i]]++;

    // ------------------------------------------------------------------
    // 3. Find the background = the most frequent pixel value (mode)
    // ------------------------------------------------------------------
    uint16_t background = 0;
    uint32_t background_count = 0;
    for (int v = 0; v < 65536; v++) {
        if (hist[v] > background_count) {
            background_count = hist[v];
            background = (uint16_t)v;
        }
    }

    // ------------------------------------------------------------------
    // 4. Compute MAD (Median Absolute Deviation) from the background
    //
    //    MAD = median( |pixel - background| )  over all pixels
    //
    //    We use the histogram to compute this efficiently without sorting
    //    all 2.25M pixels — we walk the histogram and find the median
    //    of the absolute deviations.
    // ------------------------------------------------------------------

    // Build a histogram of absolute deviations from background
    // Max possible deviation = 65535
    std::vector<uint32_t> dev_hist(65536, 0);
    for (int v = 0; v < 65536; v++) {
        if (hist[v] == 0) continue;
        int dev = std::abs(v - (int)background);
        dev_hist[dev] += hist[v];
    }

    // Find median of absolute deviations = value where cumulative count >= NPIX/2
    uint32_t mad = 0;
    uint64_t cumsum = 0;
    for (int d = 0; d < 65536; d++) {
        cumsum += dev_hist[d];
        if (cumsum >= (uint64_t)NPIX / 2) {
            mad = (uint32_t)d;
            break;
        }
    }

    // Auto-select T = 4 * MAD, rounded up to next power of two
    // (power of two makes it easy to store in 1 byte as an exponent)
    uint32_t T_raw = 4 * mad;
    uint32_t T_auto = 1;
    while (T_auto < T_raw) T_auto <<= 1;
    // Clamp to at least 1
    if (T_auto < 1) T_auto = 1;

    // ------------------------------------------------------------------
    // 5. Count anomalies at the auto-selected threshold
    // ------------------------------------------------------------------
    uint64_t within_auto = 0;
    for (int v = 0; v < 65536; v++) {
        int dev = std::abs(v - (int)background);
        if ((uint32_t)dev <= T_auto)
            within_auto += hist[v];
    }
    uint64_t anomalies_auto = NPIX - within_auto;

    // ------------------------------------------------------------------
    // 6. Print results
    // ------------------------------------------------------------------
    printf("File      : %s\n", argv[1]);
    printf("Background: %u  (appears %u times = %.2f%% of pixels)\n",
           background, background_count, background_count * 100.0 / NPIX);
    printf("MAD       : %u\n", mad);
    printf("Auto T    : %u  (= 4 * MAD, rounded to power of 2)\n", T_auto);
    printf("At auto T : %.2f%% background,  %llu anomalies\n",
           within_auto * 100.0 / NPIX, (unsigned long long)anomalies_auto);
    printf("\n");
    printf("Threshold  |  Pixels within range  |  %% of image  |  Anomalies left\n");
    printf("-----------|----------------------|--------------|----------------\n");

    int thresholds[] = {0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    for (int T : thresholds) {
        uint64_t within = 0;
        for (int v = 0; v < 65536; v++) {
            if (std::abs(v - (int)background) <= T)
                within += hist[v];
        }
        double pct = within * 100.0 / NPIX;
        // mark the auto-selected threshold
        const char* marker = ((uint32_t)T == T_auto) ? " <-- auto" : "";
        printf("  T = %4d  |  %20llu  |  %10.2f%%  |  %llu%s\n",
               T, (unsigned long long)within, pct,
               (unsigned long long)(NPIX - within), marker);
    }

    return 0;
}
