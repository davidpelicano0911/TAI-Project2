#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <filesystem>
#include <unistd.h>
#include <sys/wait.h>

namespace fs = std::filesystem;

static constexpr uint8_t MAGIC[4] = {'A','P','X','2'};

static std::string shell_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out += c;
    }
    out += "'";
    return out;
}

static int run_cmd(const std::string& cmd) {
    int ret = system(cmd.c_str());
    if (ret == -1) return 1;
    if (WIFEXITED(ret)) return WEXITSTATUS(ret);
    return 1;
}

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

static bool write_wrapped(const char* out_path, uint8_t codec_id,
                          const std::vector<uint8_t>& payload) {
    FILE* f = fopen(out_path, "wb");
    if (!f) return false;
    fwrite(MAGIC, 1, 4, f);
    fputc(codec_id, f);
    fwrite(payload.data(), 1, payload.size(), f);
    fclose(f);
    return true;
}

static void remove_if_exists(const fs::path& p) {
    std::error_code ec;
    fs::remove(p, ec);
}

int main(int argc, char* argv[]) {
    // Suporta duas formas:
    //   compress n_rows n_cols input output   (interface do professor)
    //   compress input output                 (compatibilidade com benchmark)
    if (argc != 3 && argc != 5) {
        fprintf(stderr, "Usage: %s [n_rows n_cols] <input_raw> <output.apex>\n", argv[0]);
        return 1;
    }

    const char* src_path = (argc == 5) ? argv[3] : argv[1];
    const char* dst_path = (argc == 5) ? argv[4] : argv[2];
    std::string dims     = (argc == 5) ? (std::string(argv[1]) + " " + std::string(argv[2]) + " ") : "";

    fs::path exe = fs::absolute(argv[0]).parent_path();
    fs::path tmp_base = fs::path("/tmp") / ("apex_" + std::to_string((long long)getpid()));
    fs::path range_path = tmp_base; range_path += ".nvr";

    struct Candidate {
        uint8_t id;
        fs::path path;
        std::string command;
    };

    std::string src          = shell_quote(src_path);
    std::string range_compress = shell_quote((exe / "range_compress").string());
    std::vector<Candidate> candidates = {
        {9, range_path, range_compress + " " + dims + src + " " + shell_quote(range_path.string()) + " 500  >/dev/null 2>/dev/null"},
        {9, range_path, range_compress + " " + dims + src + " " + shell_quote(range_path.string()) + " 750  >/dev/null 2>/dev/null"},
        {9, range_path, range_compress + " " + dims + src + " " + shell_quote(range_path.string()) + " 1000 >/dev/null 2>/dev/null"},
        {9, range_path, range_compress + " " + dims + src + " " + shell_quote(range_path.string()) + " 1500 >/dev/null 2>/dev/null"},
    };

    bool have_best = false;
    uint8_t best_id = 0;
    std::vector<uint8_t> best_payload;

    for (const auto& c : candidates) {
        remove_if_exists(c.path);
        if (run_cmd(c.command) != 0) continue;
        std::vector<uint8_t> payload;
        if (!read_file(c.path, payload)) continue;
        if (!have_best || payload.size() < best_payload.size()) {
            have_best = true;
            best_id = c.id;
            best_payload.swap(payload);
        }
    }

    remove_if_exists(range_path);

    if (!have_best) {
        fprintf(stderr, "apex: all candidate compressors failed\n");
        return 1;
    }
    if (!write_wrapped(dst_path, best_id, best_payload)) {
        fprintf(stderr, "apex: cannot write %s\n", dst_path);
        return 1;
    }
    return 0;
}
