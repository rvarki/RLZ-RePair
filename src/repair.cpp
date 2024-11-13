#include "repair.h"

#include <CLI11.hpp>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

#ifdef _WIN32
    #include <windows.h>
#elif __linux__
    #include <fstream>
#elif __APPLE__
    #include <sys/types.h>
    #include <sys/sysctl.h>

#endif

extern "C" {
#include "despair_c.h"
#include "repair_c.h"
}

unsigned long long getTotalSystemMemory() {
#ifdef _WIN32
    // Windows implementation using GlobalMemoryStatusEx
    MEMORYSTATUSEX statex;
    statex.dwLength = sizeof(statex);
    if (GlobalMemoryStatusEx(&statex)) {
        return statex.ullTotalPhys;  // Total physical memory in bytes
    } else {
        std::cerr << "Failed to get memory status on Windows." << std::endl;
        return 0;
    }
#elif __linux__
    // Linux implementation reading from /proc/meminfo
    std::ifstream meminfo("/proc/meminfo");
    std::string line;
    unsigned long long totalMem = 0;

    while (std::getline(meminfo, line)) {
        if (line.find("MemTotal:") == 0) {
            sscanf(line.c_str(), "MemTotal: %llu kB", &totalMem);
            totalMem *= 1024;  // Convert kB to bytes
            break;
        }
    }
    
    if (totalMem == 0) {
        std::cerr << "Failed to read /proc/meminfo on Linux." << std::endl;
    }
    
    return totalMem;

#elif __APPLE__
    // macOS implementation using sysctl
    int mib[2] = {CTL_HW, HW_MEMSIZE};
    unsigned long long totalMem = 0;
    size_t len = sizeof(totalMem);
    
    if (sysctl(mib, 2, &totalMem, &len, NULL, 0) == 0) {
        return totalMem;  // Total physical memory in bytes
    } else {
        std::cerr << "Failed to get memory size on macOS." << std::endl;
        return 0;
    }
#else
    std::cerr << "Unsupported platform." << std::endl;
    return 0;
#endif
}

int main(int argc, char *argv[]) {
    CLI::App app(
        "rlz - Run RePair with the RLZ parse.\n\nImplemented by Rahul Varki");
    std::string ref_file;
    std::string rlz_parse;
    std::string seq_file;
    float factor = 0.75;
    bool verbose = false;
    std::string version = "Version: 1.0.0";
    bool decompress = false;
    int mb = 0;

    auto total_mem = getTotalSystemMemory();
    mb = std::max(1, int(0.95 * total_mem / 1024 / 1024));

    app.add_option("-r,--ref", ref_file,
                   "The reference file used to create the RLZ parse")
        ->configurable()
        ->required();
    app.add_option("-p,--parse", rlz_parse,
                   "The RLZ parse of the sequence file (.rlz)")
        ->configurable()
        ->required();
    app.add_flag("-d,--decompress", decompress,
                 "Decompress the RLZ parse instead of compressing it");
    app.add_option("-f,--factor", factor,
                   "The factor used for the repair algorithm (default: 0.75)")
        ->configurable();
    app.add_option("-m,--mem", mb,
                   "The memory budget in MB (default: 95\% of system memory)")
        ->configurable();
    app.add_option("-s, --seq", seq_file, "The sequence file to use for repair")
        ->configurable();
    app.add_flag("--verbose", verbose, "Verbose output")->configurable();
    app.set_version_flag("-v,--version", version);
    app.footer(
        "Example usage:\n"
        "  Compress: ./repair --ref reference.fasta --parse sequence.rlz --seq sequence.fasta\n");
    app.description("Run RePair on RLZ parse");
    CLI11_PARSE(app, argc, argv);
    if (verbose) {
        spdlog::set_level(spdlog::level::debug);
    }

    spdlog::info("Starting to RePair the RLZ parse");
    spdlog::info("The reference file provided: {}", ref_file);
    spdlog::info("The RLZ parse file provided: {}", rlz_parse);
    spdlog::info("The sequence file provided: {}", seq_file);
    spdlog::info("The factor provided: {}", factor);

    if (decompress) {
        spdlog::debug("Decompressing the RLZ parse");
        auto sw_decompress = spdlog::stopwatch();
        int x = despair(seq_file.c_str());
        if (x != 0) {
            spdlog::error("Error in decompression function");
            return 1;
        }
        return 0;
    }

    spdlog::stopwatch sw_repair;

    spdlog::debug("Starting repair function");


    int x = repair_main(seq_file.c_str(), rlz_parse.c_str(),
    ref_file.c_str(), factor, mb);

    if (x != 0) {
        spdlog::error("Error in repair function");
        return 1;
    }

    auto sw_elapsed = sw_repair.elapsed();
    spdlog::debug("RePair completed in {:.6} seconds", sw_elapsed.count());

    return 0;
}