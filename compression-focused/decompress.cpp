#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <filesystem>
#include <unistd.h>

namespace fs = std::filesystem;

// Pull in range decompression logic as a library (suppresses its main).
#define RANGE_DECOMPRESS_LIB
#include "range_decompress.cpp"
#undef RANGE_DECOMPRESS_LIB

static constexpr uint8_t OUTER_MAGIC[4] = {'A','P','X','2'};

static bool write_file(const fs::path& path, const std::vector<uint8_t>& data) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return false;
    bool ok = data.empty() || fwrite(data.data(), 1, data.size(), f) == data.size();
    fclose(f);
    return ok;
}

static void remove_if_exists(const fs::path& p) {
    std::error_code ec;
    fs::remove(p, ec);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input> <output_raw>\n", argv[0]);
        return 1;
    }

    FILE* f = fopen(argv[1], "rb");
    if (!f) { fprintf(stderr, "apex: cannot open %s\n", argv[1]); return 1; }

    uint8_t magic[4];
    if (fread(magic, 1, 4, f) != 4 || memcmp(magic, OUTER_MAGIC, 4) != 0) {
        fprintf(stderr, "apex: bad magic\n");
        fclose(f);
        return 1;
    }

    std::vector<uint8_t> payload;
    uint8_t buf[65536];
    while (true) {
        size_t n = fread(buf, 1, sizeof(buf), f);
        if (n) payload.insert(payload.end(), buf, buf + n);
        if (n < sizeof(buf)) break;
    }
    fclose(f);

    fs::path tmp_payload = fs::path("/tmp") / ("apex_dec_" + std::to_string((long long)getpid()) + ".nvr");

    if (!write_file(tmp_payload, payload)) {
        fprintf(stderr, "apex: cannot write temporary payload\n");
        return 1;
    }

    int ret = range_decompress_run(tmp_payload.c_str(), argv[2]);
    remove_if_exists(tmp_payload);
    return ret;
}
