//
// Created by leonard on 7/27/26.
//

#ifndef TOOTHPASTE_DEFINITIONS_H
#define TOOTHPASTE_DEFINITIONS_H

#ifdef __clang__

#define RESTRICTED __restrict

#elif defined(__GNUC__)

#define RESTRICTED __restrict__

#elif defined(_MSC_VER)

#define RESTRICTED __restrict

#elif defined(__INTEL_COMPILER)

#define RESTRICTED __restrict

#else

#define RESTRICTED

#endif

#endif //TOOTHPASTE_DEFINITIONS_H
