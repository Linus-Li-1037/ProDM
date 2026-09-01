#ifndef PRODM_FUN3D_PMGARD_HPP
#define PRODM_FUN3D_PMGARD_HPP

#include <mpi.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <sys/stat.h>

#include "utils.hpp"
#include "MDR/Writer/WriterInterface.hpp"
#include "MDR/Retriever/RetrieverInterface.hpp"
#include "SingleFileArchive.hpp"

namespace PMGARDFUN3D {

struct FrameInfo {
    uint64_t num_elements = 0;
    double value_range = 0;
    int32_t target_level = 0;
    int32_t num_bitplanes = 0;
};

struct IOState {
    double seconds = 0;
    size_t bytes = 0;
    size_t metadata_size = 0;
    size_t retrieved_size = 0;
    std::vector<uint32_t> offsets;
};

inline bool ensure_directory(const std::string& path) {
    std::string current;
    for (char c : path) {
        current.push_back(c);
        if (c == '/' && current != "/" &&
            ::mkdir(current.c_str(), 0775) != 0 && errno != EEXIST) {
            return false;
        }
    }
    return current.empty() || current.back() == '/' ||
           ::mkdir(current.c_str(), 0775) == 0 || errno == EEXIST;
}

inline std::string subdomain_dir(const std::string& root, int subdomain) {
    std::string directory = root;
    if (directory.back() != '/') directory += '/';
    return directory + "subdomain" + std::to_string(subdomain + 1) + '/';
}

inline std::vector<std::string> variables(const std::string& directory,
                                          const std::string& specification) {
    if (specification != "all") {
        std::vector<std::string> result;
        size_t begin = 0;
        while (begin <= specification.size()) {
            const size_t end = specification.find(',', begin);
            result.push_back(specification.substr(begin, end - begin));
            if (end == std::string::npos) break;
            begin = end + 1;
        }
        return result;
    }

    std::ifstream input(directory + "metadata.json");
    std::string text((std::istreambuf_iterator<char>(input)),
                     std::istreambuf_iterator<char>());
    const size_t key = text.find("\"variables\"");
    const size_t left = text.find('[', key);
    const size_t right = text.find(']', left);
    std::vector<std::string> result;
    if (key == std::string::npos || left == std::string::npos ||
        right == std::string::npos) return result;
    for (size_t position = left; position < right;) {
        const size_t first = text.find('"', position + 1);
        if (first == std::string::npos || first >= right) break;
        const size_t second = text.find('"', first + 1);
        result.push_back(text.substr(first + 1, second - first - 1));
        position = second;
    }
    return result;
}

inline bool copy_metadata(const std::string& source_directory,
                          const std::string& output_root,
                          int subdomain, int local_partition) {
    if (local_partition != 0) return true;
    const std::string destination = subdomain_dir(output_root, subdomain);
    if (!ensure_directory(destination)) return false;
    std::ifstream input(source_directory + "metadata.json", std::ios::binary);
    std::ofstream output(destination + "metadata.json",
                         std::ios::binary | std::ios::trunc);
    output << input.rdbuf();
    return input.good() || input.eof();
}

inline std::string frame_base(const std::string& root,
                              const std::string& variable, int timestep,
                              int np, int rank) {
    return root + variable + ".dat." + std::to_string(timestep) +
           ".pmgard.p" + std::to_string(np) + ".rank" +
           std::to_string(rank);
}

inline std::vector<std::string> level_files(const std::string& base,
                                            int target_level) {
    std::vector<std::string> files;
    for (int level = 0; level <= target_level; ++level) {
        files.push_back(base + ".level" + std::to_string(level));
    }
    return files;
}

template <class T>
bool read_local_field(const std::string& data_file,
                      const std::string& partition_prefix,
                      int local_partition, int partitions_per_subdomain,
                      std::vector<T>& local) {
    size_t count = 0;
    auto full = MGARD::readfile<T>(data_file.c_str(), count);
    if (count == 0) return false;
    if (partitions_per_subdomain == 1) {
        local.assign(full.begin(), full.end());
        return true;
    }
    size_t partition_count = 0;
    auto partition = MGARD::readfile<int32_t>(
        (partition_prefix + ".part").c_str(), partition_count);
    if (partition_count != count) return false;
    for (size_t i = 0; i < count; ++i) {
        if (partition[i] == local_partition) local.push_back(full[i]);
    }
    return !local.empty();
}

template <class T>
double global_range(const std::vector<T>& local, MPI_Comm comm) {
    double local_min = std::numeric_limits<double>::max();
    double local_max = -std::numeric_limits<double>::max();
    for (T value : local) {
        local_min = std::min(local_min, static_cast<double>(value));
        local_max = std::max(local_max, static_cast<double>(value));
    }
    double global_min = 0;
    double global_max = 0;
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    return global_max - global_min;
}

class StreamWriter : public MDR::concepts::WriterInterface {
public:
    StreamWriter(std::string metadata_file, std::vector<std::string> level_files,
                bool enabled, std::shared_ptr<IOState> state)
        : metadata_file_(std::move(metadata_file)),
          level_files_(std::move(level_files)), enabled_(enabled),
          state_(std::move(state)) {}

    std::vector<uint32_t> write_level_components(
        const std::vector<std::vector<uint8_t*>>& components,
        const std::vector<std::vector<uint32_t>>& sizes) const override {
        const double start = MPI_Wtime();
        std::vector<uint32_t> level_num;
        for (size_t level = 0; level < components.size(); ++level) {
            size_t level_size = 0;
            for (uint32_t size : sizes[level]) level_size += size;
            state_->bytes += level_size;
            level_num.push_back(1);
            if (!enabled_) continue;
            FILE* file = std::fopen(level_files_[level].c_str(), "wb");
            if (!file) std::abort();
            for (size_t plane = 0; plane < components[level].size(); ++plane) {
                std::fwrite(components[level][plane], 1, sizes[level][plane], file);
            }
            std::fclose(file);
        }
        state_->seconds += MPI_Wtime() - start;
        return level_num;
    }

    void write_metadata(const uint8_t* metadata, uint32_t size) const override {
        const double start = MPI_Wtime();
        state_->bytes += size;
        state_->metadata_size = size;
        if (enabled_) {
            FILE* file = std::fopen(metadata_file_.c_str(), "wb");
            if (!file) std::abort();
            std::fwrite(metadata, 1, size, file);
            std::fclose(file);
        }
        state_->seconds += MPI_Wtime() - start;
    }

    void print() const override {}

private:
    std::string metadata_file_;
    std::vector<std::string> level_files_;
    bool enabled_;
    std::shared_ptr<IOState> state_;
};

class StreamRetriever : public MDR::concepts::RetrieverInterface {
public:
    StreamRetriever(std::string metadata_file, std::vector<std::string> level_files,
                   std::shared_ptr<IOState> state)
        : metadata_file_(std::move(metadata_file)),
          level_files_(std::move(level_files)), state_(std::move(state)) {
        state_->offsets.assign(level_files_.size(), 0);
    }

    std::vector<std::vector<const uint8_t*>> retrieve_level_components(
        const std::vector<std::vector<uint32_t>>& level_sizes,
        const std::vector<uint32_t>& retrieve_sizes,
        const std::vector<uint8_t>& previous_planes,
        const std::vector<uint8_t>& current_planes) override {
        release();
        const double start = MPI_Wtime();
        for (size_t level = 0; level < retrieve_sizes.size(); ++level) {
            uint8_t* buffer = static_cast<uint8_t*>(std::malloc(retrieve_sizes[level]));
            FILE* file = std::fopen(level_files_[level].c_str(), "rb");
            if (!file) std::abort();
            std::fseek(file, state_->offsets[level], SEEK_SET);
            std::fread(buffer, 1, retrieve_sizes[level], file);
            std::fclose(file);
            buffers_.push_back(buffer);
            state_->offsets[level] += retrieve_sizes[level];
        }
        state_->retrieved_size = state_->metadata_size;
        for (uint32_t offset : state_->offsets) state_->retrieved_size += offset;
        state_->seconds += MPI_Wtime() - start;

        std::vector<std::vector<const uint8_t*>> result;
        for (size_t level = 0; level < current_planes.size(); ++level) {
            const uint8_t* position = buffers_[level];
            std::vector<const uint8_t*> planes;
            for (int plane = previous_planes[level]; plane < current_planes[level]; ++plane) {
                planes.push_back(position);
                position += level_sizes[level][plane];
            }
            result.push_back(std::move(planes));
        }
        return result;
    }

    uint8_t* load_metadata() const override {
        const double start = MPI_Wtime();
        FILE* file = std::fopen(metadata_file_.c_str(), "rb");
        if (!file) std::abort();
        std::fseek(file, 0, SEEK_END);
        state_->metadata_size = std::ftell(file);
        std::rewind(file);
        auto* metadata = static_cast<uint8_t*>(std::malloc(state_->metadata_size));
        std::fread(metadata, 1, state_->metadata_size, file);
        std::fclose(file);
        state_->seconds += MPI_Wtime() - start;
        return metadata;
    }

    void release() override {
        for (uint8_t* buffer : buffers_) std::free(buffer);
        buffers_.clear();
    }
    void print() const override {}
    size_t get_retrieved_size() const { return state_->retrieved_size; }
    size_t get_metadata_size() const { return state_->metadata_size; }
    std::vector<uint32_t> get_offsets() const { return state_->offsets; }

private:
    std::string metadata_file_;
    std::vector<std::string> level_files_;
    std::shared_ptr<IOState> state_;
    std::vector<uint8_t*> buffers_;
};

// ---------------------------------------------------------------------------
// Single-file archive: the same MDR writer/retriever protocol, backed by one file per
// rank instead of a metadata, info and per-level file for every frame.  See
// SingleFileArchive.hpp for the layout.
// ---------------------------------------------------------------------------

inline uint64_t frame_block_index(int timestep, size_t field, size_t num_fields) {
    return single_file_block_index(static_cast<uint64_t>(timestep),
                                  static_cast<uint64_t>(field),
                                  static_cast<uint64_t>(num_fields));
}

// Block metadata is FrameInfo followed by MDR's own metadata: the driver needs the
// value range before it can turn a relative tolerance into an absolute one, and MDR
// needs its own bytes back untouched.
inline std::vector<uint8_t> pack_block_metadata(const FrameInfo& info,
                                                const std::vector<uint8_t>& mdr) {
    std::vector<uint8_t> packed(sizeof(FrameInfo) + mdr.size());
    std::memcpy(packed.data(), &info, sizeof(FrameInfo));
    if (!mdr.empty()) {
        std::memcpy(packed.data() + sizeof(FrameInfo), mdr.data(), mdr.size());
    }
    return packed;
}

class ArchiveWriter : public MDR::concepts::WriterInterface {
public:
    ArchiveWriter(std::shared_ptr<SingleFileWriter> archive,
                  std::shared_ptr<IOState> state)
        : archive_(std::move(archive)), state_(std::move(state)) {}

    // Buffers each level as one component.  The concatenation is counted as output
    // work, like the fwrite it replaces, so the driver's refactor time stays the codec
    // alone; the single write of the whole block is the archive's own io_time.
    std::vector<uint32_t> write_level_components(
        const std::vector<std::vector<uint8_t*>>& components,
        const std::vector<std::vector<uint32_t>>& sizes) const override {
        const double start = MPI_Wtime();
        std::vector<uint32_t> level_num;
        for (size_t level = 0; level < components.size(); ++level) {
            size_t level_size = 0;
            for (uint32_t size : sizes[level]) level_size += size;
            state_->bytes += level_size;
            archive_->add_component_parts(components[level], sizes[level]);
            level_num.push_back(1);
        }
        state_->seconds += MPI_Wtime() - start;
        return level_num;
    }

    void write_metadata(const uint8_t* metadata, uint32_t size) const override {
        const double start = MPI_Wtime();
        state_->bytes += size;
        state_->metadata_size = size;
        archive_->set_metadata(metadata, size);
        state_->seconds += MPI_Wtime() - start;
    }

    void print() const override {}

private:
    std::shared_ptr<SingleFileWriter> archive_;
    std::shared_ptr<IOState> state_;
};

class ArchiveRetriever : public MDR::concepts::RetrieverInterface {
public:
    ArchiveRetriever(const SingleFileArchive* archive,
                     std::shared_ptr<IOState> state)
        : archive_(archive), state_(std::move(state)) {}

    // Points the retriever at one (timestep, field) and loads that block's metadata.
    bool select_block(uint64_t index) {
        if (!archive_ || !archive_->block_layout(index, layout_)) return false;
        block_metadata_.resize(layout_.metadata_size);
        if (!block_metadata_.empty() &&
            !archive_->read_at(layout_.metadata_offset, block_metadata_.data(),
                               block_metadata_.size())) {
            return false;
        }
        if (block_metadata_.size() < sizeof(FrameInfo)) {
            std::cerr << "ERROR: block " << index << " of " << archive_->path()
                      << " has no FrameInfo" << std::endl;
            return false;
        }
        std::memcpy(&info_, block_metadata_.data(), sizeof(FrameInfo));
        state_->offsets.assign(layout_.component_sizes.size(), 0);
        state_->metadata_size = block_metadata_.size() - sizeof(FrameInfo);
        state_->retrieved_size = 0;
        block_ = index;
        return true;
    }

    const FrameInfo& info() const { return info_; }
    uint64_t block() const { return block_; }

    std::vector<std::vector<const uint8_t*>> retrieve_level_components(
        const std::vector<std::vector<uint32_t>>& level_sizes,
        const std::vector<uint32_t>& retrieve_sizes,
        const std::vector<uint8_t>& previous_planes,
        const std::vector<uint8_t>& current_planes) override {
        release();
        const double start = MPI_Wtime();
        buffers_.resize(retrieve_sizes.size());
        // Grown rather than assumed: the level list MDR passes here need not match the
        // block's component count (see the bounds check below).
        if (state_->offsets.size() < retrieve_sizes.size()) {
            state_->offsets.resize(retrieve_sizes.size(), 0);
        }
        for (size_t level = 0; level < retrieve_sizes.size(); ++level) {
            // Kept allocated even if the read fails: the pointers handed back must stay
            // valid, and archive_->good() is what reports the failure.
            buffers_[level].assign(retrieve_sizes[level], 0);
            // ComposedReconstructor also calls this with a level list shorter than the
            // block's when it reconstructs below the finest level, so a level index past
            // the components is expected rather than an error -- but reading at
            // component_offsets[level] would be out of bounds.
            if (retrieve_sizes[level] && level < layout_.component_offsets.size()) {
                archive_->read_at(
                    layout_.component_offsets[level] + state_->offsets[level],
                    buffers_[level].data(), retrieve_sizes[level]);
                state_->offsets[level] += retrieve_sizes[level];
            }
        }
        state_->retrieved_size = state_->metadata_size;
        for (uint32_t offset : state_->offsets) state_->retrieved_size += offset;
        state_->seconds += MPI_Wtime() - start;

        std::vector<std::vector<const uint8_t*>> result;
        for (size_t level = 0; level < current_planes.size(); ++level) {
            const uint8_t* position = buffers_[level].data();
            std::vector<const uint8_t*> planes;
            for (int plane = previous_planes[level]; plane < current_planes[level];
                 ++plane) {
                planes.push_back(position);
                position += level_sizes[level][plane];
            }
            result.push_back(std::move(planes));
        }
        return result;
    }

    // MDR's own metadata: the block's bytes past the FrameInfo the driver read.
    uint8_t* load_metadata() const override {
        const size_t size = block_metadata_.size() - sizeof(FrameInfo);
        state_->metadata_size = size;
        auto* metadata = static_cast<uint8_t*>(std::malloc(size ? size : 1));
        if (size) {
            std::memcpy(metadata, block_metadata_.data() + sizeof(FrameInfo), size);
        }
        return metadata;
    }

    void release() override { buffers_.clear(); }
    void print() const override {}
    size_t get_retrieved_size() const { return state_->retrieved_size; }
    size_t get_metadata_size() const { return state_->metadata_size; }
    std::vector<uint32_t> get_offsets() const { return state_->offsets; }

private:
    const SingleFileArchive* archive_ = nullptr;
    std::shared_ptr<IOState> state_;
    SingleFileBlockLayout layout_;
    std::vector<uint8_t> block_metadata_;
    std::vector<std::vector<uint8_t>> buffers_;
    FrameInfo info_;
    uint64_t block_ = 0;
};

}  // namespace PMGARDFUN3D

#endif
