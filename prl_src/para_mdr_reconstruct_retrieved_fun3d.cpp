// Reconstruct stage of a retrieve-then-transfer workflow for PMGARD.
//
// Reads the bundle para_mdr_retrieve_fun3d wrote -- one file per rank, each level already
// cut down to the bitplanes one tolerance needs, the decision recorded beside the
// metadata -- and decodes it.  Nothing is decided here: the bitplane counts come out of
// the record through FixedPlanInterpreter, so what is decoded is exactly what was
// shipped and the size interpreter never runs on this side.  The tolerance is the one
// the bundle was retrieved for; it is not a command-line argument.
//
// verify=1 reads the original fields back and reports, per timestep and field, the
// largest absolute deviation against the bound the tolerance asked for.  PMGARD works on
// the flat array of a rank's nodes, so the original is gathered the way the refactor
// gathered it: whole subdomain file for one partition per subdomain, else filtered by
// the .part file under partition_prefix.

#include <mpi.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "FUN3DPMGARD.hpp"
#include "MDR/Reconstructor/Reconstructor.hpp"

// Largest absolute deviation over this rank's nodes; -1 if the lengths disagree.
template <class T>
double max_abs_error(const std::vector<T>& original, const T* reconstructed, size_t count) {
    if (original.size() != count) return -1;
    double worst = 0;
    for (size_t i = 0; i < count; ++i) {
        worst = std::max(worst, std::fabs(static_cast<double>(original[i]) -
                                          static_cast<double>(reconstructed[i])));
    }
    return worst;
}

template <class T>
int run(int argc, char** argv, MPI_Comm comm) {
    int rank = 0;
    int np = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);

    const char* usage =
        "Usage: %s bundle_dir variables nt verify [original_root partition_prefix] f|d\n"
        "  The tolerance is the one the bundle was retrieved for.\n"
        "  verify=1 checks each reconstruction against the original fields under "
        "original_root (partition_prefix selects this rank's nodes when there is more "
        "than one partition per subdomain); verify=0 omits both.\n";
    if (argc != 6 && argc != 8) {
        if (rank == 0) std::fprintf(stderr, usage, argv[0]);
        return 1;
    }
    const int num_timesteps = std::atoi(argv[3]);
    const int verify = std::atoi(argv[4]);
    if (np < 72 || np % 72 != 0 || num_timesteps < 1 || (verify != 0 && verify != 1) ||
        argc != 6 + 2 * verify) {
        if (rank == 0) std::fprintf(stderr, usage, argv[0]);
        return 1;
    }

    std::string bundle_root = argv[1];
    if (bundle_root.back() != '/') bundle_root += '/';
    const int subdomain = rank % 72;
    const int local_partition = rank / 72;
    const int partitions_per_subdomain = np / 72;
    const auto variables = PMGARDFUN3D::variables(
        PMGARDFUN3D::subdomain_dir(bundle_root, subdomain), argv[2]);
    if (variables.empty()) return 1;

    std::string original_dir;
    std::string partition_prefix;
    if (verify == 1) {
        std::string original_root = argv[5];
        if (original_root.back() != '/') original_root += '/';
        original_dir = PMGARDFUN3D::subdomain_dir(original_root, subdomain);
        partition_prefix = argv[6];
        if (!partition_prefix.empty() && partition_prefix.front() != '/') {
            partition_prefix = original_dir + partition_prefix;
        }
    }

    PMGARDFUN3D::SingleFileArchive archive(
        PMGARDFUN3D::retrieved_archive_name(bundle_root, np, rank));
    if (!archive.open()) return 1;
    if (archive.num_fields() != variables.size() ||
        archive.num_timesteps() < static_cast<uint64_t>(num_timesteps) ||
        archive.header().element_size != sizeof(T)) {
        std::fprintf(stderr, "ERROR: rank %d: %s does not match this run "
                             "(fields, timesteps or f/d)\n",
                     rank, archive.path().c_str());
        return 1;
    }

    const size_t num_frames = static_cast<size_t>(num_timesteps) * variables.size();
    unsigned long long local_retrieved = 0;
    double local_reconstruct_time = 0;
    unsigned long long local_num_elements = 0;
    double rel_tolerance = 0;
    std::vector<double> local_errors(verify == 1 ? num_frames : 0, 0);
    std::vector<double> max_errors(local_errors.size(), 0);
    std::vector<double> frame_ranges(verify == 1 ? num_frames : 0, 0);
    std::vector<double> frame_bounds(verify == 1 ? num_frames : 0, 0);

    for (int timestep = 0; timestep < num_timesteps; ++timestep) {
        for (size_t field = 0; field < variables.size(); ++field) {
            const size_t frame = static_cast<size_t>(timestep) * variables.size() + field;
            auto io = std::make_shared<PMGARDFUN3D::IOState>();
            PMGARDFUN3D::ArchiveRetriever retriever(&archive, io);
            if (!retriever.select_block(PMGARDFUN3D::frame_block_index(
                    timestep, field, variables.size()))) {
                return 1;
            }
            const PMGARDFUN3D::FrameInfo info = retriever.info();
            PMGARDFUN3D::RetrievalRecord record;
            if (!PMGARDFUN3D::split_record(retriever.block_metadata(), record)) {
                std::fprintf(stderr, "ERROR: block %llu of %s has no retrieval record; "
                                     "was it written by para_mdr_retrieve_fun3d?\n",
                             static_cast<unsigned long long>(retriever.block()),
                             archive.path().c_str());
                return 1;
            }
            rel_tolerance = record.rel_tolerance;

            std::vector<T> original;
            if (verify == 1) {
                const std::string input = original_dir + variables[field] + ".dat." +
                                          std::to_string(timestep);
                if (!PMGARDFUN3D::read_local_field(input, partition_prefix, local_partition,
                                                   partitions_per_subdomain, original)) {
                    std::fprintf(stderr, "ERROR: rank %d cannot read the original %s\n",
                                 rank, input.c_str());
                    return 1;
                }
                frame_ranges[frame] = info.value_range;
                frame_bounds[frame] = record.abs_tolerance;
            }

            // The planes to decode are the planes the bundle holds.
            std::vector<uint8_t> planes(record.taken.begin(), record.taken.end());
            auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
            auto interleaver = MDR::DirectInterleaver<T>();
            using Stream = typename std::conditional<
                sizeof(T) == 8, uint64_t, uint32_t>::type;
            auto encoder = MDR::PerBitBPEncoder_old<T, Stream>();
            auto compressor = MDR::AdaptiveLevelCompressor(64);
            auto estimator = MDR::MaxErrorEstimatorHB<T>();
            PMGARDFUN3D::FixedPlanInterpreter interpreter(planes);
            auto reconstructor = MDR::ComposedReconstructor<
                T, decltype(decomposer), decltype(interleaver),
                decltype(encoder), decltype(compressor), decltype(interpreter),
                decltype(estimator), PMGARDFUN3D::ArchiveRetriever>(
                    decomposer, interleaver, encoder, compressor, interpreter,
                    retriever);
            reconstructor.load_metadata();

            const double io_before = io->seconds;
            MPI_Barrier(comm);
            const double start = MPI_Wtime();
            const T* decoded = reconstructor.progressive_reconstruct(record.abs_tolerance, -1);
            local_reconstruct_time += MPI_Wtime() - start - (io->seconds - io_before);
            local_retrieved += reconstructor.get_retrieved_size() +
                               sizeof(PMGARDFUN3D::FrameInfo);

            if (!archive.good()) {
                std::fprintf(stderr, "ERROR: rank %d failed to read block %llu of %s\n",
                             rank, static_cast<unsigned long long>(retriever.block()),
                             archive.path().c_str());
                return 1;
            }
            if (verify == 1) {
                local_errors[frame] = max_abs_error(original, decoded, info.num_elements);
            }
            local_num_elements += info.num_elements;
        }
    }

    unsigned long long total_num_elements = 0;
    unsigned long long total_retrieved = 0;
    double max_reconstruct_time = 0;
    MPI_Reduce(&local_num_elements, &total_num_elements, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(&local_retrieved, &total_retrieved, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(&local_reconstruct_time, &max_reconstruct_time, 1,
               MPI_DOUBLE, MPI_MAX, 0, comm);
    if (verify == 1) {
        MPI_Reduce(local_errors.data(), max_errors.data(),
                   static_cast<int>(local_errors.size()), MPI_DOUBLE, MPI_MAX, 0, comm);
    }
    if (rank == 0) {
        std::printf("PMGARD-Retrieved ranks=%d preprocessing=0.000000\n", np);
        std::printf("  rel_tol=%.3g retrieved=%llu bitrate=%.4f ratio=%.4f reconstruct=%.6f\n",
                    rel_tolerance, total_retrieved,
                    total_retrieved * 8.0 / total_num_elements,
                    static_cast<double>(total_num_elements) * sizeof(T) / total_retrieved,
                    max_reconstruct_time);
        if (verify == 1) {
            int violations = 0;
            for (size_t frame = 0; frame < num_frames; ++frame) {
                const int timestep = static_cast<int>(frame / variables.size());
                const std::string& variable = variables[frame % variables.size()];
                const bool violated = !(max_errors[frame] <= frame_bounds[frame]);
                violations += violated ? 1 : 0;
                std::printf("  verify timestep=%d field=%s rel_tol=%.3g range=%.6g bound=%.6g "
                            "max_error=%.6g %s\n",
                            timestep, variable.c_str(), rel_tolerance, frame_ranges[frame],
                            frame_bounds[frame], max_errors[frame],
                            violated ? "VIOLATED" : "ok");
            }
            std::printf("  verify frames=%zu violations=%d\n", num_frames, violations);
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
