#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <thread>
#include <filesystem>
#include <unistd.h>

namespace fs = std::filesystem;

// Pull in range compression logic as a library (suppresses its main).
#define RANGE_COMPRESS_LIB
#include "range_compress.cpp"
#undef RANGE_COMPRESS_LIB

static constexpr uint8_t OUTER_MAGIC[4] = {'A','P','X','2'};

static bool read_file(const fs::path& path, std::vector<uint8_t>& data) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    rewind(f);
    if (n < 0) { fclose(f); return false; }
    data.resize((size_t)n);
    bool ok = data.empty() || fread(data.data(), 1, data.size(), f) == data.size();
    fclose(f);
    return ok;
}

static bool write_wrapped(const char* out_path, const std::vector<uint8_t>& payload) {
    FILE* f = fopen(out_path, "wb");
    if (!f) return false;
    fwrite(OUTER_MAGIC, 1, 4, f);
    fwrite(payload.data(), 1, payload.size(), f);
    fclose(f);
    return true;
}

static void remove_if_exists(const fs::path& p) {
    std::error_code ec;
    fs::remove(p, ec);
}

int main(int argc, char* argv[]) {
    if (argc != 3 && argc != 5) {
        fprintf(stderr, "Usage: %s [n_rows n_cols] <input_raw> <output>\n", argv[0]);
        return 1;
    }

    int H = (argc == 5) ? std::atoi(argv[1]) : 1500;
    int W = (argc == 5) ? std::atoi(argv[2]) : 1500;
    const char* src_path = (argc == 5) ? argv[3] : argv[1];
    const char* dst_path = (argc == 5) ? argv[4] : argv[2];

    fs::path tmp_base = fs::path("/tmp") / ("apex_" + std::to_string((long long)getpid()));

    static const int block_sizes[] = {250, 300, 400, 500, 750, 1000, 1500};
    static const int n_bs = (int)(sizeof(block_sizes) / sizeof(block_sizes[0]));

    auto tmp_path = [&](int bs) -> fs::path {
        fs::path p = tmp_base;
        p += "_" + std::to_string(bs) + ".nvr";
        return p;
    };

    // Remove any stale temp files.
    for (int i = 0; i < n_bs; i++) remove_if_exists(tmp_path(block_sizes[i]));

    // Run all block sizes in parallel.
    std::vector<std::thread> threads;
    threads.reserve(n_bs);
    for (int i = 0; i < n_bs; i++) {
        int bs = block_sizes[i];
        std::string dst = tmp_path(bs).string();
        threads.emplace_back([=]() {
            range_compress_run(H, W, src_path, dst.c_str(), bs);
        });
    }
    for (auto& t : threads) t.join();

    // Pick the smallest result.
    std::vector<uint8_t> best;
    for (int i = 0; i < n_bs; i++) {
        fs::path p = tmp_path(block_sizes[i]);
        std::vector<uint8_t> payload;
        if (read_file(p, payload)) {
            if (best.empty() || payload.size() < best.size())
                best.swap(payload);
        }
        remove_if_exists(p);
    }

    if (best.empty()) {
        fprintf(stderr, "apex: all candidate compressors failed\n");
        return 1;
    }
    if (!write_wrapped(dst_path, best)) {
        fprintf(stderr, "apex: cannot write %s\n", dst_path);
        return 1;
    }
    return 0;
}
