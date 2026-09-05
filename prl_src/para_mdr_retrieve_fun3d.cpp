// Retrieve stage of a retrieve-then-transfer workflow for PMGARD.
//
// On the cluster that holds the refactored archives, this runs the half of a progressive
// reconstruction that decides WHAT to fetch -- metadata, error estimator, size
// interpreter -- fetches exactly that, and writes it out as a bundle: one file per rank in
// the same archive layout, every level cut down to the bitplanes the tolerance needs,
// the decision recorded beside the metadata.  Nothing is decoded here.  PMGARD needs no
// mesh, so only metadata.json is copied beside the bundle, which can then be moved to
// another cluster and reconstructed there with para_mdr_reconstruct_retrieved_fun3d.

#include <mpi.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "FUN3DPMGARD.hpp"
#include "MDR/Reconstructor/Reconstructor.hpp"

template <class T>
int run(int argc, char** argv, MPI_Comm comm) {
    int rank = 0;
    int np = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);
    if (argc != 7) {
        if (rank == 0) {
            std::fprintf(stderr, "Usage: %s input_dir output_dir variables nt rel_tol f|d\n",
                         argv[0]);
        }
        return 1;
    }
    if (np < 72 || np % 72 != 0) {
        if (rank == 0) std::fprintf(stderr, "ERROR: ranks must be a multiple of 72\n");
        return 1;
    }

    std::string input_root = argv[1];
    std::string output_root = argv[2];
    if (input_root.back() != '/') input_root += '/';
    if (output_root.back() != '/') output_root += '/';
    const int num_timesteps = std::atoi(argv[4]);
    const double rel_tolerance = std::atof(argv[5]);
    if (num_timesteps < 1 || !(rel_tolerance > 0)) return 1;

    const int subdomain = rank % 72;
    const int local_partition = rank / 72;
    const std::string input_subdomain = PMGARDFUN3D::subdomain_dir(input_root, subdomain);
    const auto variables = PMGARDFUN3D::variables(input_subdomain, argv[3]);
    if (variables.empty()) return 1;
    if (!PMGARDFUN3D::copy_metadata(input_subdomain, output_root, subdomain,
                                    local_partition)) {
        std::fprintf(stderr, "ERROR: rank %d cannot copy metadata.json into %s\n",
                     rank, output_root.c_str());
        return 1;
    }

    PMGARDFUN3D::SingleFileArchive archive(
        PMGARDFUN3D::single_file_archive_name(input_root, np, rank));
    if (!archive.open()) return 1;
    if (archive.num_fields() != variables.size() ||
        archive.num_timesteps() < static_cast<uint64_t>(num_timesteps) ||
        archive.header().element_size != sizeof(T)) {
        std::fprintf(stderr, "ERROR: rank %d: %s does not match this run "
                             "(fields, timesteps or f/d)\n",
                     rank, archive.path().c_str());
        return 1;
    }

    PMGARDFUN3D::SingleFileWriter bundle(
        PMGARDFUN3D::retrieved_archive_name(output_root, np, rank),
        static_cast<uint64_t>(num_timesteps),
        static_cast<uint64_t>(variables.size()), sizeof(T));
    if (!bundle.open()) return 1;

    unsigned long long local_bundle_bytes = bundle.prologue_size();
    unsigned long long local_payload_bytes = 0;   // what the tolerance actually needs
    double local_retrieve_time = 0;               // decide + fetch, excludes the bundle write

    for (int timestep = 0; timestep < num_timesteps; ++timestep) {
        for (size_t field = 0; field < variables.size(); ++field) {
            const uint64_t index =
                PMGARDFUN3D::frame_block_index(timestep, field, variables.size());

            MPI_Barrier(comm);
            const double start = MPI_Wtime();

            auto io = std::make_shared<PMGARDFUN3D::IOState>();
            PMGARDFUN3D::ArchiveRetriever retriever(&archive, io);
            if (!retriever.select_block(index)) return 1;
            const PMGARDFUN3D::FrameInfo info = retriever.info();

            // ---- decide: the same estimator + interpreter the reconstructor uses ----
            auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
            auto interleaver = MDR::DirectInterleaver<T>();
            using Stream = typename std::conditional<
                sizeof(T) == 8, uint64_t, uint32_t>::type;
            auto encoder = MDR::PerBitBPEncoder_old<T, Stream>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto estimator = MDR::MaxErrorEstimatorHB<T>();
            auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<
                MDR::MaxErrorEstimatorHB<T>>(estimator);
            auto planner = MDR::ComposedReconstructor<
                T, decltype(decomposer), decltype(interleaver),
                decltype(encoder), decltype(compressor), decltype(interpreter),
                decltype(estimator), PMGARDFUN3D::ArchiveRetriever>(
                    decomposer, interleaver, encoder, compressor, interpreter,
                    retriever);
            planner.load_metadata();
            const double abs_tolerance = rel_tolerance * info.value_range;
            std::vector<uint8_t> planes;
            const auto level_bytes = planner.plan_retrieval(abs_tolerance, planes);

            // ---- fetch: each level up to the planes taken ----
            bundle.begin_block();
            const auto& layout = retriever.layout();
            std::vector<uint8_t> level;
            for (size_t L = 0; L < level_bytes.size(); ++L) {
                level.assign(level_bytes[L], 0);
                if (level_bytes[L] && (L >= layout.component_offsets.size() ||
                    !archive.read_at(layout.component_offsets[L], level.data(),
                                     level.size()))) {
                    return 1;
                }
                bundle.add_component(level.data(), level.size());
                local_payload_bytes += level.size();
            }
            local_payload_bytes += retriever.block_metadata().size();
            local_retrieve_time += MPI_Wtime() - start;

            // ---- record the decision and write the bundle block ----
            PMGARDFUN3D::RetrievalRecord record;
            record.rel_tolerance = rel_tolerance;
            record.abs_tolerance = abs_tolerance;
            record.taken.assign(planes.begin(), planes.end());
            const auto bundle_metadata =
                PMGARDFUN3D::append_record(retriever.block_metadata(), record);
            size_t block_size = 0;
            if (!bundle.commit_block(index, bundle_metadata.data(),
                                     bundle_metadata.size(), block_size)) {
                return 1;
            }
            local_bundle_bytes += block_size;
        }
    }
    if (!bundle.close()) return 1;

    unsigned long long total_bundle_bytes = 0;
    unsigned long long total_payload_bytes = 0;
    double max_retrieve_time = 0;
    double local_read_time = archive.io_time();
    double max_read_time = 0;
    double local_write_time = bundle.io_time();
    double max_write_time = 0;
    MPI_Reduce(&local_bundle_bytes, &total_bundle_bytes, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(&local_payload_bytes, &total_payload_bytes, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(&local_retrieve_time, &max_retrieve_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&local_read_time, &max_read_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&local_write_time, &max_write_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    if (rank == 0) {
        std::printf("PMGARD-Retrieve ranks=%d timesteps=%d fields=%zu rel_tol=%.3g "
                    "retrieve=%.6f read_io=%.6f write_io=%.6f "
                    "retrieved_payload=%llu bundle_bytes=%llu files=%d\n",
                    np, num_timesteps, variables.size(), rel_tolerance,
                    max_retrieve_time, max_read_time, max_write_time,
                    total_payload_bytes, total_bundle_bytes, np);
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
