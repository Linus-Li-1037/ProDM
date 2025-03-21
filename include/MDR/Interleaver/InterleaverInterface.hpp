#ifndef _MDR_INTERLEAVER_INTERFACE_HPP
#define _MDR_INTERLEAVER_INTERFACE_HPP
#include <cstdint>

namespace MDR {
    namespace concepts {

        // level data interleaver: interleave level coefficients
        template<class T>
        class InterleaverInterface {
        public:

            virtual ~InterleaverInterface() = default;

            virtual void interleave(T const * data, const std::vector<uint32_t>& dims, const std::vector<uint32_t>& dims_fine, const std::vector<uint32_t>& dims_coasre, T * buffer, std::vector<uint32_t> strides=std::vector<uint32_t>()) const = 0;

            virtual void reposition(T const * buffer, const std::vector<uint32_t>& dims, const std::vector<uint32_t>& dims_fine, const std::vector<uint32_t>& dims_coasre, T * data, std::vector<uint32_t> strides=std::vector<uint32_t>()) const = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
