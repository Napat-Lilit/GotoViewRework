
#pragma once

#include <cstdarg>
#include <string>
#include <mutex>

void Warning(const char* format, ...);
void Error(const char* format, ...);

#define Assert(expr) \
    ((expr) ? (void)0 : \
        Error("Assertion \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))

#define WarningAssert(expr) \
    ((expr) ? (void)0 : \
        Warning("Warning \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))

#define CHECK_LE(a, b) \
    WarningAssert((a) <= (b))

