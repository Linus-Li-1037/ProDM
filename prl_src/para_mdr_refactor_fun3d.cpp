#include <mpi.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include "FUN3DPMGARD.hpp"
#include "MDR/Refactor/Refactor.hpp"

template <class T>
int run(int argc, char** argv, MPI_Comm comm) {
    int rank = 0;
    int np = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);

    if (argc != 10) {
        if (rank == 0) {
            std::fprintf(stderr,
                "Usage: %s data_root output_dir variables nt partition_prefix "
                "target_level num_bitplanes write_mode f|d\n", argv[0]);
        }
        return 1;
    }
    if (np < 72 || np % 72 != 0) {
        if (rank == 0) std::fprintf(stderr, "ERROR: ranks must be a multiple of 72\n");
        return 1;
    }

    std::string data_root = argv[1];
    std::string output_root = argv[2];
    if (data_root.back() != '/') data_root += '/';
    if (output_root.back() != '/') output_root += '/';
    const int num_timesteps = std::atoi(argv[4]);
    const int target_level = std::atoi(argv[6]);
    int num_bitplanes = std::atoi(argv[7]);
    const int write_mode = std::atoi(argv[8]);
    if (num_timesteps < 1 || target_level < 0 || num_bitplanes < 1 ||
        (write_mode != 0 && write_mode != 1)) return 1;
    if (num_bitplanes % 2 != 0) ++num_bitplanes;
    const int maximum_bitplanes = 8 * sizeof(T);
    if (num_bitplanes > maximum_bitplanes) num_bitplanes = maximum_bitplanes;

    const int subdomain = rank % 72;
    const int local_partition = rank / 72;
    const int partitions_per_subdomain = np / 72;
    const std::string data_dir =
        PMGARDFUN3D::subdomain_dir(data_root, subdomain);
    const auto variables = PMGARDFUN3D::variables(data_dir, argv[3]);
    if (variables.empty()) return 1;

    std::string partition_prefix = argv[5];
    if (!partition_prefix.empty() && partition_prefix.front() != '/') {
        partition_prefix = data_dir + partition_prefix;
    }
    if (write_mode == 1 &&
        (!PMGARDFUN3D::ensure_directory(output_root) ||
         !PMGARDFUN3D::copy_metadata(
             data_dir, output_root, subdomain, local_partition))) return 1;

    // One file per rank holding every (timestep, field) it owns, instead of a metadata,
    // info and per-level file for each frame.  write_mode 0 keeps it a dry run: blocks
    // are assembled and measured, but no file is created.
    auto archive = std::make_shared<PMGARDFUN3D::SingleFileWriter>(
        write_mode == 1
            ? PMGARDFUN3D::single_file_archive_name(output_root, np, rank)
            : std::string(),
        static_cast<uint64_t>(num_timesteps),
        static_cast<uint64_t>(variables.size()), sizeof(T));
    if (!archive->open()) return 1;

    unsigned long long local_refactored_size = archive->prologue_size();
    unsigned long long local_num_elements = 0;
    double local_refactor_time = 0;
    for (int timestep = 0; timestep < num_timesteps; ++timestep) {
        for (size_t field = 0; field < variables.size(); ++field) {
            const std::string& variable = variables[field];
            std::vector<T> local;
            const std::string input = data_dir + variable + ".dat." +
                                      std::to_string(timestep);
            if (!PMGARDFUN3D::read_local_field(
                    input, partition_prefix, local_partition,
                    partitions_per_subdomain, local)) {
                std::fprintf(stderr, "ERROR: rank %d cannot read %s\n",
                             rank, input.c_str());
                return 1;
            }
            const double range = PMGARDFUN3D::global_range(local, comm);
            auto io = std::make_shared<PMGARDFUN3D::IOState>();
            archive->begin_block();
            PMGARDFUN3D::ArchiveWriter writer(archive, io);

            auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
            auto interleaver = MDR::DirectInterleaver<T>();
            using Stream = typename std::conditional<
                sizeof(T) == 8, uint64_t, uint32_t>::type;
            auto encoder = MDR::PerBitBPEncoder_old<T, Stream>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto collector = MDR::SquaredErrorCollector<T>();
            auto refactor = MDR::ComposedRefactor<
                T, decltype(decomposer), decltype(interleaver),
                decltype(encoder), decltype(compressor), decltype(collector),
                PMGARDFUN3D::ArchiveWriter>(
                    decomposer, interleaver, encoder, compressor, collector,
                    writer);
            refactor.negabinary = false;

            const std::vector<uint32_t> dimensions{
                static_cast<uint32_t>(local.size())};
            MPI_Barrier(comm);
            const double start = MPI_Wtime();
            refactor.refactor(local.data(), dimensions, target_level,
                              num_bitplanes);
            local_refactor_time += MPI_Wtime() - start - io->seconds;

            PMGARDFUN3D::FrameInfo info;
            info.num_elements = local.size();
            info.value_range = range;
            info.target_level = target_level;
            info.num_bitplanes = num_bitplanes;

            // FrameInfo rides at the front of the block's metadata, so a frame is one
            // self-contained object and the reader needs no sidecar file.
            const auto metadata = PMGARDFUN3D::pack_block_metadata(
                info, archive->metadata());
            size_t block_size = 0;
            if (!archive->commit_block(
                    PMGARDFUN3D::frame_block_index(timestep, field, variables.size()),
                    metadata.data(), metadata.size(), block_size)) {
                return 1;
            }
            local_refactored_size += block_size;
            local_num_elements += local.size();
        }
    }
    if (!archive->close()) return 1;

    unsigned long long total_refactored_size = 0;
    unsigned long long total_num_elements = 0;
    double max_refactor_time = 0;
    double local_write_time = archive->io_time();
    double max_write_time = 0;
    MPI_Reduce(&local_refactored_size, &total_refactored_size, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(&local_num_elements, &total_num_elements, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(&local_refactor_time, &max_refactor_time, 1,
               MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&local_write_time, &max_write_time, 1,
               MPI_DOUBLE, MPI_MAX, 0, comm);
    if (rank == 0) {
        std::printf("PMGARD ranks=%d timesteps=%d fields=%zu "
                    "preprocessing=0.000000 refactor=%.6f write=%.6f "
                    "files=%d blocks_per_file=%d "
                    "total_refactored_size=%llu full_ratio=%.4f\n",
                    np, num_timesteps, variables.size(),
                    max_refactor_time, max_write_time,
                    write_mode == 1 ? np : 0,
                    num_timesteps * static_cast<int>(variables.size()),
                    total_refactored_size,
                    static_cast<double>(total_num_elements) * sizeof(T) /
                        total_refactored_size);
    }
    return 0;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    const char dtype = argc > 1 ? argv[argc - 1][0] : 'f';
    int result = 1;
    if (dtype == 'f') result = run<float>(argc, argv, MPI_COMM_WORLD);
    if (dtype == 'd') result = run<double>(argc, argv, MPI_COMM_WORLD);
    MPI_Finalize();
    return result;
}
