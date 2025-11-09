// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_COMPILER_H_
#define MAIA_COMPILER_H_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file deefines compiler-specific behaviour.
////////////////////////////////////////////////////////////////////////////////

/// \brief Detects the compiler
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__) && !defined(__NVCOMPILER)                   \
    && !defined(__PGI) && !defined(_CRAYC) && !defined(_SX) && !defined(__INTEL_LLVM_COMPILER)
#define MAIA_GCC_COMPILER
#elif defined(__clang__)
#define MAIA_CLANG_COMPILER
#elif defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#define MAIA_INTEL_COMPILER
#elif defined(__IBMCPP__)
#define MAIA_IBM_COMPILER
#elif defined(__NVCOMPILER)
#define MAIA_NVHPC_COMPILER
#define MAIA_PGI_COMPILER
#elif defined(__PGI)
#define MAIA_PGI_COMPILER
#elif defined(_CRAYC)
#define MAIA_CRAY_COMPILER
#elif defined(_MSC_VER)
#define MAIA_MS_COMPILER
#elif defined(_SX)
#define MAIA_SX_COMPILER
#else
#error Unsupported compiler
#endif

/// \brief Defines the non-standard macro FUN_, which is replaced by the
/// compiler with the current function name.
#if defined(MAIA_GCC_COMPILER) || defined(MAIA_CLANG_COMPILER) || defined(MAIA_INTEL_COMPILER)                         \
    || defined(MAIA_NVHPC_COMPILER) || defined(MAIA_PGI_COMPILER)
#define FUN_ __PRETTY_FUNCTION__
#elif defined(MAIA_IBM_COMPILER)
#define FUN_ __FUNCTION__
#elif defined(MAIA_CRAY_COMPILER) || defined(MAIA_SX_COMPILER)
#define FUN_ __func__
#elif defined(MAIA_MS_COMPILER)
#define FUN_ __FUNCDNAME__
#else
#error Function name macro is NECESSARY and undefined for current compiler.
#endif

/// \brief Defines the RESTRICT macro (see USE_RESTRICT)
#if defined(USE_RESTRICT)                                                                                              \
    && (defined(MAIA_SX_COMPILER) || defined(MAIA_GCC_COMPILER) || defined(MAIA_CLANG_COMPILER)                        \
        || defined(MAIA_INTEL_COMPILER) || defined(MAIA_MS_COMPILER) || defined(MAIA_NVHPC_COMPILER)                   \
        || defined(MAIA_PGI_COMPILER))
#define RESTRICT __restrict
#elif defined(USE_RESTRICT) && defined(MAIA_IBM_COMPILER)
#error RESTRICT not implemented for the ibm compiler
#elif defined(USE_RESTRICT) && defined(MAIA_CRAY_COMPILER)
#error RESTRICT not implemented for the cray compiler
#elif defined(USE_RESTRICT)
#error RESTRICT is defined but compiler is unknown
#else
#define RESTRICT
#endif


// use default c++11 attribute specifier sequence
#if defined(COMPILER_ATTRIBUTES)
#define ATTRIBUTES1(a) [[a]]
#define ATTRIBUTES2(a, b) [[a, b]]
#define ATTRIBUTES3(a, b, c) [[a, b, c]]
#define ATTRIBUTES4(a, b, c, d) [[a, b, c, d]]
/// \brief no return tells the compiler that a given function will never return
/// (e.g. because it causes the program to terminate)
#define ATTRIBUTE_NORETURN noreturn
/// \brief indicates that function has memory order dependency
#define ATTRIBUTE_CARRIES_DEPENDENCY carries_dependency
#endif

/// \brief Defines the compiler attributes macros
#if defined(COMPILER_ATTRIBUTES)                                                                                       \
    && (defined(MAIA_GCC_COMPILER) || defined(MAIA_CLANG_COMPILER) || defined(MAIA_INTEL_COMPILER)                     \
        || defined(MAIA_NVHPC_COMPILER) || defined(MAIA_PGI_COMPILER))
/// \brief flatten encourages all function calls inside a function definition to be inlined
#define ATTRIBUTE_FLATTEN gnu::flatten
/// \brief  always inline encourages the compiler to inline the function if its
/// marked with the inline keyword
#define ATTRIBUTE_ALWAYS_INLINE gnu::always_inline
/// \brief hot indicates the compiler that the function is a hotspot, the compiler
/// will try to optimize this function more heavily
#define ATTRIBUTE_HOT gnu::hot
/// \brief pure indicates purity in the functional sense
#define ATTRIBUTE_PURE gnu::pure
/// \brief deactivate autovectorization
//#define ATTRIBUTE_NO_AUTOVEC gnu::optimize("no-tree-vectorize")
#define ATTRIBUTE_NO_AUTOVEC
// todo labels:totest this should now be also supported in GCC and Clang
#elif defined(MAIA_MS_COMPILER)
/// \brief flatten encourages all function calls inside a function definition to be inlined
#define ATTRIBUTE_FLATTEN
/// \brief  always inline encourages the compiler to inline the function if its
/// marked with the inline keyword
#define ATTRIBUTE_ALWAYS_INLINE
/// \brief hot indicates the compiler that the function is a hotspot, the compiler
/// will try to optimize this function more heavily
#define ATTRIBUTE_HOT hot
/// \brief pure indicates purity in the functional sense
#define ATTRIBUTE_PURE
/// \brief deactivate autovectorization
#define ATTRIBUTE_NO_AUTOVEC
#elif defined(COMPILER_ATTRIBUTES) && defined(MAIA_IBM_COMPILER)
#pragma message("WARNING: ATTRIBUTES are not implemented for the ibm compiler")
#define ATTRIBUTES1(a)
#define ATTRIBUTES2(a, b)
#define ATTRIBUTES3(a, b, c)
#define ATTRIBUTES4(a, b, c, d)
#define ATTRIBUTE_FLATTEN
#define ATTRIBUTE_NORETURN
#define ATTRIBUTE_ALWAYS_INLINE
#define ATTRIBUTE_HOT
#define ATTRIBUTE_PURE
/// \brief deactivate autovectorization
#define ATTRIBUTE_NO_AUTOVEC
#elif defined(COMPILER_ATTRIBUTES) && defined(MAIA_CRAY_COMPILER)
#pragma message("WARNING: ATTRIBUTES are not implemented for the cray compiler")
#define ATTRIBUTES1(a)
#define ATTRIBUTES2(a, b)
#define ATTRIBUTES3(a, b, c)
#define ATTRIBUTES4(a, b, c, d)
#define ATTRIBUTE_FLATTEN
#define ATTRIBUTE_NORETURN
#define ATTRIBUTE_ALWAYS_INLINE
#define ATTRIBUTE_HOT
#define ATTRIBUTE_PURE
/// \brief deactivate autovectorization
#define ATTRIBUTE_NO_AUTOVEC
#elif defined(MAIA_IBM_COMPILER)
// This is a fix for the IBM compiler until it supports the 'message' pragma
#define ATTRIBUTES1(a)
#define ATTRIBUTES2(a, b)
#define ATTRIBUTES3(a, b, c)
#define ATTRIBUTES4(a, b, c, d)
#define ATTRIBUTE_FLATTEN
#define ATTRIBUTE_NORETURN
#define ATTRIBUTE_ALWAYS_INLINE
#define ATTRIBUTE_HOT
#define ATTRIBUTE_PURE
/// \brief deactivate autovectorization
#define ATTRIBUTE_NO_AUTOVEC
#else
#pragma message("WARNING: ATTRIBUTES are disabled")
#define ATTRIBUTES1(a)
#define ATTRIBUTES2(a, b)
#define ATTRIBUTES3(a, b, c)
#define ATTRIBUTES4(a, b, c, d)
#define ATTRIBUTE_FLATTEN
#define ATTRIBUTE_NORETURN
#define ATTRIBUTE_ALWAYS_INLINE
#define ATTRIBUTE_HOT
#define ATTRIBUTE_PURE
/// \brief deactivate autovectorization
#define ATTRIBUTE_NO_AUTOVEC
#endif

// todo labels:toenhance,totest allignas is now fully supported by both GCC and CLANG. Change and test!!!
#if defined(USE_ALIGNMENT) && (defined(MAIA_GCC_COMPILER) || defined(MAIA_CLANG_COMPILER))
#define ALIGNED_F(a) static_cast<const MFloat*>(__builtin_assume_aligned((a), sizeof(MFloat)))
#define ALIGNED_I(a) static_cast<const MInt*>(__builtin_assume_aligned((a), sizeof(MInt)))
#define ALIGNED_L(a) static_cast<const MLong*>(__builtin_assume_aligned((a), sizeof(MLong)))
#define ALIGNED_MF(a) static_cast<MFloat*>(__builtin_assume_aligned((a), sizeof(MFloat)))
#define ALIGNED_MI(a) static_cast<MInt*>(__builtin_assume_aligned((a), sizeof(MInt)))
#define ALIGNED_ML(a) static_cast<MLong*>(__builtin_assume_aligned((a), sizeof(MLong)))
#elif defined(USE_ALIGNMENT) && defined(MAIA_MS_COMPILER)
#define ALIGNED_F(a) alignas(sizeof(MFloat)) a
#define ALIGNED_I(a) alignas(sizeof(MInt)) a
#define ALIGNED_L(a) alignas(sizeof(MLong)) a
#define ALIGNED_MF(a) alignas(sizeof(MFloat)) a
#define ALIGNED_MI(a) alignas(sizeof(MInt)) a
#define ALIGNED_ML(a) alignas(sizeof(MLong)) a
#elif defined(USE_ALIGNMENT) && defined(MAIA_INTEL_COMPILER)
#error ALIGNMENT not implemented for intel compiler
#elif defined(USE_ALIGNMENT) && defined(MAIA_IBM_COMPILER)
#error ALIGNMENT not implemented for ibm compiler
#elif defined(USE_ALIGNMENT) && defined(MAIA_PGI_COMPILER)
#error ALIGNMENT not implemented for pgi compiler
#elif defined(USE_ALIGNMENT) && defined(MAIA_NVHPC_COMPILER)
#error ALIGNMENT not implemented for nvhpc compiler
#elif defined(USE_ALIGNMENT) && defined(MAIA_CRAY_COMPILER)
#error ALIGNMENT not implemented for cray compiler
#elif defined(USE_ALIGNMENT) && defined(MAIA_SX_COMPILER)
#error ALIGNMENT not implemented for NEC SX compiler
#elif defined(USE_ALIGNMENT)
#error ALIGNMENT is defined but compiler is unknown
#else
#define ALIGNED_F(a) (a)
#define ALIGNED_I(a) (a)
#define ALIGNED_F(a) (a)
#define ALIGNED_MF(a) (a)
#define ALIGNED_MI(a) (a)
#define ALIGNED_MF(a) (a)


#endif

#if defined(MAIA_MS_COMPILER)
#define M_PI 3.14159265359
#define NOMINMAX
#define MAIA_DISABLE_DG // unsure but not needed
#define MAIA_DISABLE_LB // problems with postprocessing solver
#endif

// Enable certain settings to ensure Windows compatibility
#if defined(_WIN64)
#define MAIA_WINDOWS
#define WITH_HDF5
#endif

// Replace 'IF_CONSTEXPR' by 'if constexpr' after switching to C++17
#if __cplusplus >= 201703
#define IF_CONSTEXPR if constexpr
#else
#define IF_CONSTEXPR if
#endif

// TODO labels:gpu This is a work around label used at certain location to
// handle nvhpc's bug of static constexpr
#if defined(MAIA_NVHPC_COMPILER) && defined(MAIA_PSTL)
#define WAR_NVHPC_PSTL
#endif

#endif // MAIA_COMPILER_H_
