#ifndef _MDR_DECOMPOSER_INTERFACE_HPP
#define _MDR_DECOMPOSER_INTERFACE_HPP
#include <cstdint>

namespace MDR {
    namespace concepts {

        // inplace data decomposer: de-correlates and overwrites original data
        template<class T>
        class DecomposerInterface {
        public:

            virtual ~DecomposerInterface() = default;

            virtual void decompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides) const = 0;

            virtual void recompose(T * data, const std::vector<uint32_t>& dimensions, uint32_t target_level, std::vector<uint32_t> strides) const = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
