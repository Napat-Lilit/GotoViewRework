
#include "error.hpp"

// Error Reporting Functions
static std::string StringVaprintf(const std::string &fmt, va_list args) {
    // Figure out how much space we need to allocate; add an extra character for '\0'.
    va_list argsCopy;
    va_copy(argsCopy, args);
    size_t size = vsnprintf(nullptr, 0, fmt.c_str(), args) + 1;
    std::string str;
    str.resize(size);
    vsnprintf(&str[0], size, fmt.c_str(), argsCopy);
    str.pop_back();  // remove trailing NUL
    return str;
}
static void processError(const char *format, va_list args, const char *errorType) {
    // Build up an entire formatted error string and print it all at once;
    // this way, if multiple threads are printing messages at once, they
    // don't get jumbled up...
    std::string errorString;
    errorString += errorType;
    errorString += ": ";
    errorString += StringVaprintf(format, args);

    // Print the error message (but not more than one time).
    static std::string lastError;
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);
    if (errorString != lastError) {
        fprintf(stderr, "%s\n", errorString.c_str());
        lastError = errorString;
    }
}

// For small problems
void Warning(const char* format, ...) {
    va_list args;
    va_start(args, format);
    processError(format, args, "Warning");
    va_end(args);
}
// For large problems
void Error(const char* format, ...) {
    va_list args;
    va_start(args, format);
    processError(format, args, "Error");
    va_end(args);

    std::abort();
}
