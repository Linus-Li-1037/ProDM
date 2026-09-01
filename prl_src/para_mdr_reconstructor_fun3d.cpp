#include <mpi.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include "FUN3DPMGARD.hpp"
#include "MDR/Reconstructor/Reconstructor.hpp"

template <class T>
int run(int argc, char** argv, MPI_Comm comm) {
    int rank = 0;
    int np = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);
    if (argc < 7) {
        if (rank == 0) {
            std::fprintf(stderr, "Usage: %s input_dir variables nt "
                         "num_tolerances tolerance... f|d\n", argv[0]);
        }
        return 1;
    }
    const int num_timesteps = std::atoi(argv[3]);
    const int num_tolerances = std::atoi(argv[4]);
    if (np < 72 || np % 72 != 0 || num_timesteps < 1 ||
        num_tolerances < 1 || argc != 6 + num_tolerances) return 1;

    std::vector<double> tolerances(num_tolerances);
    for (int i = 0; i < num_tolerances; ++i) {
        tolerances[i] = std::atof(argv[5 + i]);
    }
    std::string input_root = argv[1];
    if (input_root.back() != '/') input_root += '/';
    const int subdomain = rank % 72;
    const auto variables = PMGARDFUN3D::variables(
        PMGARDFUN3D::subdomain_dir(input_root, subdomain), argv[2]);
    if (variables.empty()) return 1;

    // One file per rank: the block for (timestep, field) is found through the directory
    // at the front of it, so the read path needs no MPI at all.
    PMGARDFUN3D::SingleFileArchive archive(
        PMGARDFUN3D::single_file_archive_name(input_root, np, rank));
    if (!archive.open()) return 1;
    if (archive.num_fields() != variables.size() ||
        archive.num_timesteps() < static_cast<uint64_t>(num_timesteps)) {
        std::fprintf(stderr,
                     "ERROR: rank %d asked for %d timesteps of %zu fields, %s holds "
                     "%llu timesteps of %llu fields\n",
                     rank, num_timesteps, variables.size(), archive.path().c_str(),
                     static_cast<unsigned long long>(archive.num_timesteps()),
                     static_cast<unsigned long long>(archive.num_fields()));
        return 1;
    }
    if (archive.header().element_size != sizeof(T)) {
        std::fprintf(stderr,
                     "ERROR: %s was refactored from %u-byte values, this run uses %zu "
                     "(f / d mismatch)\n",
                     archive.path().c_str(), archive.header().element_size, sizeof(T));
        return 1;
    }

    std::vector<unsigned long long> local_retrieved(num_tolerances, 0);
    std::vector<unsigned long long> total_retrieved(num_tolerances, 0);
    std::vector<double> local_reconstruct_time(num_tolerances, 0);
    std::vector<double> max_reconstruct_time(num_tolerances, 0);
    unsigned long long local_num_elements = 0;

    for (int timestep = 0; timestep < num_timesteps; ++timestep) {
        for (size_t field = 0; field < variables.size(); ++field) {
            auto io = std::make_shared<PMGARDFUN3D::IOState>();
            PMGARDFUN3D::ArchiveRetriever retriever(&archive, io);
            if (!retriever.select_block(PMGARDFUN3D::frame_block_index(
                    timestep, field, variables.size()))) {
                return 1;
            }
            const PMGARDFUN3D::FrameInfo info = retriever.info();

            auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
            auto interleaver = MDR::DirectInterleaver<T>();
            using Stream = typename std::conditional<
                sizeof(T) == 8, uint64_t, uint32_t>::type;
            auto encoder = MDR::PerBitBPEncoder_old<T, Stream>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto estimator = MDR::MaxErrorEstimatorHB<T>();
            auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<
                MDR::MaxErrorEstimatorHB<T>>(estimator);
            auto reconstructor = MDR::ComposedReconstructor<
                T, decltype(decomposer), decltype(interleaver),
                decltype(encoder), decltype(compressor), decltype(interpreter),
                decltype(estimator), PMGARDFUN3D::ArchiveRetriever>(
                    decomposer, interleaver, encoder, compressor, interpreter,
                    retriever);
            reconstructor.load_metadata();

            for (int i = 0; i < num_tolerances; ++i) {
                const double io_before = io->seconds;
                MPI_Barrier(comm);
                const double start = MPI_Wtime();
                reconstructor.progressive_reconstruct(
                    tolerances[i] * info.value_range, -1);
                local_reconstruct_time[i] += MPI_Wtime() - start -
                    (io->seconds - io_before);
                local_retrieved[i] += reconstructor.get_retrieved_size() +
                                      sizeof(PMGARDFUN3D::FrameInfo);
            }

            // A short read inside MDR's retriever protocol cannot fail the
            // reconstruction, so it is caught here rather than showing up as a bad
            // reconstruction.
            if (!archive.good()) {
                std::fprintf(stderr, "ERROR: rank %d failed to read block %llu of %s\n",
                             rank,
                             static_cast<unsigned long long>(retriever.block()),
                             archive.path().c_str());
                return 1;
            }
            local_num_elements += info.num_elements;
        }
    }

    unsigned long long total_num_elements = 0;
    MPI_Reduce(&local_num_elements, &total_num_elements, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(local_retrieved.data(), total_retrieved.data(), num_tolerances,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(local_reconstruct_time.data(), max_reconstruct_time.data(),
               num_tolerances, MPI_DOUBLE, MPI_MAX, 0, comm);
    if (rank == 0) {
        std::printf("PMGARD ranks=%d preprocessing=0.000000\n", np);
        for (int i = 0; i < num_tolerances; ++i) {
            std::printf("  rel_tol=%.3g retrieved=%llu bitrate=%.4f "
                        "ratio=%.4f reconstruct=%.6f\n",
                        tolerances[i], total_retrieved[i],
                        total_retrieved[i] * 8.0 / total_num_elements,
                        static_cast<double>(total_num_elements) * sizeof(T) /
                            total_retrieved[i], max_reconstruct_time[i]);
        }
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
