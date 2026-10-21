// A deliberately tiny test framework: one registry, a handful of macros, no
// external dependency. Each test file registers its own cases at static
// initialisation time and test_main.cpp runs them all.
#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace testing {

struct TestCase {
    std::string name;
    std::function<void()> body;
};

// Thrown by a failed assertion; caught by the runner.
struct Failure {
    std::string message;
};

std::vector<TestCase>& registry();

struct Registrar {
    Registrar(const char* name, std::function<void()> body);
};

inline int& failure_count() {
    static int count = 0;
    return count;
}

inline int& check_count() {
    static int count = 0;
    return count;
}

std::string format_double(double value);

template <typename A, typename B>
void check_equal(const A& actual, const B& expected, const char* actual_text,
                 const char* expected_text, const char* file, int line) {
    ++check_count();
    if (!(actual == expected)) {
        std::ostringstream message;
        message << file << ':' << line << ": " << actual_text << " == " << expected_text
                << " failed\n  actual:   " << actual << "\n  expected: " << expected;
        throw Failure{message.str()};
    }
}

inline void check_near(double actual, double expected, double tolerance, const char* actual_text,
                       const char* expected_text, const char* file, int line) {
    ++check_count();
    if (!(std::abs(actual - expected) <= tolerance)) {
        std::ostringstream message;
        message << file << ':' << line << ": " << actual_text << " is not within " << tolerance
                << " of " << expected_text << "\n  actual:   " << format_double(actual)
                << "\n  expected: " << format_double(expected)
                << "\n  delta:    " << format_double(std::abs(actual - expected));
        throw Failure{message.str()};
    }
}

inline void check_true(bool value, const char* text, const char* file, int line) {
    ++check_count();
    if (!value) {
        std::ostringstream message;
        message << file << ':' << line << ": expected true: " << text;
        throw Failure{message.str()};
    }
}

inline void check_false(bool value, const char* text, const char* file, int line) {
    ++check_count();
    if (value) {
        std::ostringstream message;
        message << file << ':' << line << ": expected false: " << text;
        throw Failure{message.str()};
    }
}

// Runs `body` and requires that it throws, optionally with a message that
// contains `needle`.
template <typename Callable>
void check_throws(Callable&& body, const std::string& needle, const char* text, const char* file,
                  int line) {
    ++check_count();
    try {
        body();
    } catch (const Failure&) {
        throw;
    } catch (const std::exception& error) {
        if (!needle.empty() && std::string(error.what()).find(needle) == std::string::npos) {
            std::ostringstream message;
            message << file << ':' << line << ": " << text << " threw the wrong message\n"
                    << "  what():   " << error.what() << "\n  expected to contain: " << needle;
            throw Failure{message.str()};
        }
        return;
    }
    std::ostringstream message;
    message << file << ':' << line << ": " << text << " did not throw";
    throw Failure{message.str()};
}

}  // namespace testing

#define TEST(suite, name)                                                           \
    static void suite##_##name##_body();                                            \
    static ::testing::Registrar suite##_##name##_registrar(#suite "." #name,        \
                                                           suite##_##name##_body);  \
    static void suite##_##name##_body()

#define CHECK(expression) \
    ::testing::check_true((expression), #expression, __FILE__, __LINE__)

#define CHECK_FALSE(expression) \
    ::testing::check_false((expression), #expression, __FILE__, __LINE__)

#define CHECK_EQ(actual, expected) \
    ::testing::check_equal((actual), (expected), #actual, #expected, __FILE__, __LINE__)

#define CHECK_NEAR(actual, expected, tolerance) \
    ::testing::check_near((actual), (expected), (tolerance), #actual, #expected, __FILE__, __LINE__)

#define CHECK_THROWS_AS(statement, needle) \
    ::testing::check_throws([&]() { statement; }, (needle), #statement, __FILE__, __LINE__)
