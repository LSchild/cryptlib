//
// Created by leonard on 3/10/25.
//

#include "pirouette.h"
#include "argparse.hpp"

int main(int argc, char** argv) {

    argparse::ArgumentParser program("pirouette");

    program.add_argument("-r", "--number-of-runs")
    .required()
    .help("Specify the number of runs")
    .scan<'i', uint64_t>();

    program.add_argument("-dbN", "--database-entries")
        .required()
        .help("Specify the number of database entries in *bits*")
        .scan<'i', uint64_t>();

    program.add_argument("-size", "--record-size")
        .required()
        .help("Specify the record size in bits").scan<'i', uint64_t>();

    program.add_argument("-mu1", "--dimension-1")
        .default_value(11)
        .required()
        .help("Specify the first folding dimension").scan<'i', uint64_t>();

    program.add_argument("-mu2", "--dimension-2")
            .default_value(11)
            .required()
            .help("Specify the second folding dimension").scan<'i', uint64_t>();

    program.add_argument("-mu3", "--dimension-3")
        .default_value(11)
        .required()
        .help("Specify the second folding dimension").scan<'i', uint64_t>();

    program.add_argument("-crt", "--use-fast")
        .flag()
        .help("Use CRT version for Phase 1. Note: this is for demonstration purposes, the current implementation contains a bug.");

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    auto n_runs = program.get<uint64_t>("-r");
    auto entries_in_bits = program.get<uint64_t>("-dbN");
    auto record_size_in_bits = program.get<uint64_t>("-size");
    auto mu1 = program.get<uint64_t>("-mu1");
    auto mu2 = program.get<uint64_t>("-mu2");
    auto mu3 = program.get<uint64_t>("-mu3");
    auto mode = program["--use-fast"] == true ? CRT : PRIME;

    auto default_config = PirouettePIR::GenerateDefaultConfig(mode);
    default_config->SetMetaParams(mu1,mu2,mu3,entries_in_bits);

    auto db = GenerateTestDatabase(1ull << entries_in_bits, record_size_in_bits);
    default_config->PrepareDB(entries_in_bits, record_size_in_bits, db);

    std::srand(time(nullptr));
    uint64_t n_records = 1u << entries_in_bits;
    std::vector<std::pair<uint64_t, uint64_t>> queries;
    std::vector<uint64_t> queried_indices;

    uint64_t total_query_sizes = 0;
    for (uint32_t i = 0; i < n_runs; i++) {
        auto idx = std::rand() % n_records;
        queried_indices.push_back(idx);
        auto q = default_config->CreateQuery(idx, std::rand());
        queries.push_back(q);

        auto bits_seed = GetMSB(q.first) - 1;
        auto bits_msg = GetMSB(q.second) - 1;
        total_query_sizes += bits_msg + bits_seed;
    }
    std::vector<RLWECiphertext> results;

    std::cerr << "[!] Starting measurements for :" << std::endl;
    std::cerr << "[!] N_ROUNDS = " << n_runs << std::endl;
    std::cerr << "[!] N_ENTRIES = " << n_records << std::endl;
    std::cerr << "[!] RECORD_SIZE_BITS = " << record_size_in_bits << std::endl;
    std::cerr << "[!] QUERY_SIZE (byte) = " << total_query_sizes / (8 * n_runs) << std::endl;

    auto measurement_start = TIC_HD;
    for (auto& v : queries) {
        results.push_back(default_config->Query(v));
    }
    auto measurement_stop = TIC_HD;
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(measurement_stop - measurement_start).count();

    std::cerr << "[!!!] Evaluation took " << elapsed_time << "ms. or " << elapsed_time / n_runs << "ms. / query" << std::endl;
}