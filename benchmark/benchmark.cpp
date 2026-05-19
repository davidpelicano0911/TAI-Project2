#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <chrono>
#include <vector>
#include <string>
#include <filesystem>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

namespace fs = std::filesystem;

static constexpr int WIDTH  = 1500;
static constexpr int HEIGHT = 1500;
static constexpr int N_PIXELS = WIDTH * HEIGHT;
static constexpr int RAW_BYTES = N_PIXELS * 2;

// ---------------------------------------------------------------------------
// Timing helpers
// ---------------------------------------------------------------------------

using Clock = std::chrono::steady_clock;
using Ms    = std::chrono::duration<double, std::milli>;

static double now_ms() {
    return std::chrono::duration_cast<Ms>(Clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------------------
// Shell helpers
// ---------------------------------------------------------------------------

struct RunResult {
    int    exit_code;
    double wall_ms;
};

static RunResult run(const std::string& cmd) {
    double t0 = now_ms();
    int ret = system(cmd.c_str());
    double t1 = now_ms();
    return { WEXITSTATUS(ret), t1 - t0 };
}

// ---------------------------------------------------------------------------
// Verification: byte-exact comparison of decompressed vs original
// ---------------------------------------------------------------------------

static bool verify(const std::string& original, const std::string& decompressed) {
    FILE* fa = fopen(original.c_str(),     "rb");
    FILE* fb = fopen(decompressed.c_str(), "rb");
    if (!fa || !fb) { if (fa) fclose(fa); if (fb) fclose(fb); return false; }

    const size_t BUF = 65536;
    std::vector<uint8_t> ba(BUF), bb(BUF);
    bool ok = true;
    while (ok) {
        size_t ra = fread(ba.data(), 1, BUF, fa);
        size_t rb = fread(bb.data(), 1, BUF, fb);
        if (ra != rb || memcmp(ba.data(), bb.data(), ra) != 0) { ok = false; break; }
        if (ra == 0) break;
    }
    fclose(fa); fclose(fb);
    return ok;
}

// ---------------------------------------------------------------------------
// File size helper
// ---------------------------------------------------------------------------

static long file_size(const std::string& path) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return -1;
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fclose(f);
    return sz;
}

// ---------------------------------------------------------------------------
// Compressor descriptors
// ---------------------------------------------------------------------------

struct Compressor {
    std::string name;
    // {src} and {dst} are placeholders replaced at runtime
    std::string compress_cmd;
    std::string decompress_cmd;
    std::string compressed_ext; // extension added to compressed file
};

static std::vector<Compressor> make_compressors(const std::string& /*tmp*/) {
    return {
        // --- gzip ---
        { "gzip-1",  "gzip -1 -k -f {src} && mv {src}.gz {dst}", "gzip -dk -f {src} && mv {srcnoext} {dst}", ".gz" },
        { "gzip-6",  "gzip -6 -k -f {src} && mv {src}.gz {dst}", "gzip -dk -f {src} && mv {srcnoext} {dst}", ".gz" },
        { "gzip-9",  "gzip -9 -k -f {src} && mv {src}.gz {dst}", "gzip -dk -f {src} && mv {srcnoext} {dst}", ".gz" },

        // --- bzip2 ---
        { "bzip2-1", "bzip2 -1 -k -f {src} && mv {src}.bz2 {dst}", "bzip2 -dk -f {src} && mv {srcnoext} {dst}", ".bz2" },
        { "bzip2-9", "bzip2 -9 -k -f {src} && mv {src}.bz2 {dst}", "bzip2 -dk -f {src} && mv {srcnoext} {dst}", ".bz2" },

        // --- xz / lzma ---
        { "xz-1",   "xz -1 -k -f {src} && mv {src}.xz {dst}", "xz -dk -f {src} && mv {srcnoext} {dst}", ".xz" },
        { "xz-6",   "xz -6 -k -f {src} && mv {src}.xz {dst}", "xz -dk -f {src} && mv {srcnoext} {dst}", ".xz" },
        { "xz-9",   "xz -9 -k -f {src} && mv {src}.xz {dst}", "xz -dk -f {src} && mv {srcnoext} {dst}", ".xz" },

        // --- zstd ---
        { "zstd-1",  "zstd -1  -q -f {src} -o {dst}", "zstd -dq -f {src} -o {dst}", ".zst" },
        { "zstd-3",  "zstd -3  -q -f {src} -o {dst}", "zstd -dq -f {src} -o {dst}", ".zst" },
        { "zstd-9",  "zstd -9  -q -f {src} -o {dst}", "zstd -dq -f {src} -o {dst}", ".zst" },
        { "zstd-19", "zstd -19 -q -f {src} -o {dst}", "zstd -dq -f {src} -o {dst}", ".zst" },

        // --- our compressors ---
        { "rais",    "../rais/compress {src} {dst} 2>/dev/null",   "../rais/decompress {src} {dst}", ".rais" },
        { "astra",   "../astra/compress {src} {dst} 2>/dev/null",  "../astra/decompress {src} {dst}", ".astr" },
    };
}

// ---------------------------------------------------------------------------
// Simple string replace
// ---------------------------------------------------------------------------

static std::string sreplace(std::string s, const std::string& from, const std::string& to) {
    size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
        s.replace(pos, from.size(), to);
        pos += to.size();
    }
    return s;
}

// Build compress command: src=original file, dst=compressed output path
static std::string build_compress_cmd(const Compressor& c,
                                      const std::string& src,
                                      const std::string& dst) {
    std::string cmd = c.compress_cmd;
    cmd = sreplace(cmd, "{src}", src);
    cmd = sreplace(cmd, "{dst}", dst);
    return cmd;
}

// Build decompress command: src=compressed file, dst=decompressed output path
static std::string build_decompress_cmd(const Compressor& c,
                                        const std::string& src,
                                        const std::string& dst) {
    // For tools that decompress to a fixed name we need {srcnoext}
    fs::path p(src);
    std::string srcnoext = (p.parent_path() / p.stem()).string();

    std::string cmd = c.decompress_cmd;
    cmd = sreplace(cmd, "{srcnoext}", srcnoext);
    cmd = sreplace(cmd, "{src}", src);
    cmd = sreplace(cmd, "{dst}", dst);
    return cmd;
}

// ---------------------------------------------------------------------------
// Per-file result
// ---------------------------------------------------------------------------

struct FileResult {
    std::string file;
    long        original_bytes;
    long        compressed_bytes;
    double      ratio;          // compressed / original
    double      bits_per_pixel; // compressed_bytes*8 / N_PIXELS
    double      compress_ms;
    double      decompress_ms;
    bool        lossless;
};

// ---------------------------------------------------------------------------
// Run one compressor against one file
// ---------------------------------------------------------------------------

static FileResult benchmark_one(const Compressor& c,
                                const std::string& data_path,
                                const std::string& tmp_dir) {
    fs::path src(data_path);
    std::string stem = src.filename().string();

    std::string compressed   = tmp_dir + "/" + stem + c.compressed_ext;
    std::string decompressed = tmp_dir + "/" + stem + ".dec";

    // Compress
    std::string ccmd = build_compress_cmd(c, data_path, compressed);
    RunResult cr = run(ccmd);

    // Decompress
    std::string dcmd = build_decompress_cmd(c, compressed, decompressed);
    RunResult dr = run(dcmd);

    long orig_sz = RAW_BYTES;
    long comp_sz = file_size(compressed);

    bool lossless = (cr.exit_code == 0 && dr.exit_code == 0)
                    && verify(data_path, decompressed);

    // Cleanup temp files
    remove(compressed.c_str());
    remove(decompressed.c_str());

    double ratio = (comp_sz > 0) ? (double)comp_sz / orig_sz : -1.0;
    // bits per BYTE of original (matches teacher's leaderboard metric)
    double bpp   = (comp_sz > 0) ? ((double)comp_sz * 8.0) / orig_sz : -1.0;

    return { stem, orig_sz, comp_sz, ratio, bpp, cr.wall_ms, dr.wall_ms, lossless };
}

// ---------------------------------------------------------------------------
// Print helpers
// ---------------------------------------------------------------------------

static void print_header() {
    printf("%-10s %-8s %12s %12s %8s %10s %12s %12s %10s %8s\n",
           "Compressor", "File",
           "OrigBytes", "CompBytes",
           "Ratio", "b/B",
           "CompressMs", "DecompressMs", "TotalMs", "Lossless");
    printf("%s\n", std::string(120, '-').c_str());
}

static void print_row(const std::string& cname, const FileResult& r) {
    printf("%-10s %-8s %12ld %12ld %8.4f %10.4f %12.1f %12.1f %10.1f %8s\n",
           cname.c_str(), r.file.c_str(),
           r.original_bytes, r.compressed_bytes,
           r.ratio, r.bits_per_pixel,
           r.compress_ms, r.decompress_ms,
           r.compress_ms + r.decompress_ms,
           r.lossless ? "YES" : "NO");
}

static void print_summary(const std::string& cname,
                          const std::vector<FileResult>& results) {
    long  sum_comp  = 0, sum_orig = 0;
    double sum_cms  = 0, sum_dms = 0;
    int n = 0;
    for (const auto& r : results) {
        if (r.compressed_bytes > 0) {
            sum_orig += r.original_bytes;
            sum_comp += r.compressed_bytes;
            sum_cms  += r.compress_ms;
            sum_dms  += r.decompress_ms;
            ++n;
        }
    }
    if (n == 0) return;
    double ratio = (double)sum_comp / sum_orig;
    double bpb   = (double)sum_comp * 8.0 / sum_orig;
    printf("%-10s %-8s %12ld %12ld %8.4f %10.4f %12.1f %12.1f %10.1f\n",
           cname.c_str(), "(total)",
           sum_orig, sum_comp,
           ratio, bpb,
           sum_cms, sum_dms, sum_cms + sum_dms);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    // Usage: benchmark <data_dir> [tmp_dir] [-o results.csv]
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <data_dir> [tmp_dir] [-o out.csv]\n", argv[0]);
        return 1;
    }

    std::string data_dir  = argv[1];
    std::string tmp_dir   = "/tmp/tai_bench";
    std::string csv_path;

    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "-o" && i + 1 < argc) {
            csv_path = argv[++i];
        } else {
            tmp_dir = argv[i];
        }
    }

    fs::create_directories(tmp_dir);

    // Collect data files A–H (or whatever is present)
    std::vector<std::string> data_files;
    for (const auto& entry : fs::directory_iterator(data_dir)) {
        if (entry.is_regular_file()) {
            long sz = fs::file_size(entry.path());
            if (sz == RAW_BYTES) { // only files matching expected size
                data_files.push_back(entry.path().string());
            }
        }
    }
    std::sort(data_files.begin(), data_files.end());

    if (data_files.empty()) {
        fprintf(stderr, "No data files found in %s (expected size %d bytes)\n",
                data_dir.c_str(), RAW_BYTES);
        return 1;
    }

    printf("\nTAI Project 2 — Compression Benchmark\n");
    printf("Data dir : %s\n", data_dir.c_str());
    printf("Tmp dir  : %s\n", tmp_dir.c_str());
    printf("Files    : %zu   (%d bytes each = %dx%d × 16-bit)\n\n",
           data_files.size(), RAW_BYTES, WIDTH, HEIGHT);

    auto compressors = make_compressors(tmp_dir);

    // Open CSV output if requested
    FILE* csv = nullptr;
    if (!csv_path.empty()) {
        csv = fopen(csv_path.c_str(), "w");
        if (!csv) { fprintf(stderr, "Cannot open %s for writing\n", csv_path.c_str()); return 1; }
        fcntl(fileno(csv), F_SETFD, FD_CLOEXEC);
        fprintf(csv, "compressor,file,orig_bytes,comp_bytes,ratio,bpp,compress_ms,decompress_ms,lossless\n");
    }

    print_header();

    for (const auto& c : compressors) {
        std::vector<FileResult> results;
        for (const auto& f : data_files) {
            FileResult r = benchmark_one(c, f, tmp_dir);
            print_row(c.name, r);
            if (csv) {
                fprintf(csv, "%s,%s,%ld,%ld,%.6f,%.6f,%.3f,%.3f,%s\n",
                        c.name.c_str(), r.file.c_str(),
                        r.original_bytes, r.compressed_bytes,
                        r.ratio, r.bits_per_pixel,
                        r.compress_ms, r.decompress_ms,
                        r.lossless ? "true" : "false");
            }
            results.push_back(r);
        }
        printf("%s\n", std::string(110, '-').c_str());
        print_summary(c.name, results);
        printf("\n");
    }

    if (csv) { fclose(csv); printf("Results saved to %s\n", csv_path.c_str()); }

    return 0;
}
