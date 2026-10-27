#pragma once
// Portability shim: the project was written for MSVC. On other compilers,
// provide the few MSVC-only functions it uses.
#ifdef _MSC_VER
#include <intrin.h>
#include <io.h>
#include <fcntl.h>
#else
#include <cerrno>
#include <cstdio>
inline int fopen_s(FILE** file, const char* path, const char* mode) {
	*file = std::fopen(path, mode);
	return *file ? 0 : errno;
}
#define fprintf_s fprintf
#define sprintf_s snprintf
inline unsigned int __lzcnt(unsigned int x) { return x ? __builtin_clz(x) : 32; }
inline unsigned int _lzcnt_u64(unsigned long long x) { return x ? __builtin_clzll(x) : 64; }
#define _setmode(fd, mode) ((void)(fd), (void)(mode), 0)
#define _fileno fileno
#define _O_U16TEXT 0
#define _O_TEXT 0
#endif
