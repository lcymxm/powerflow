#include <exception>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "test_framework.hpp"

namespace testing {

std::vector<TestCase>& registry() {
    static std::vector<TestCase> tests;
    return tests;
}

Registrar::Registrar(const char* name, std::function<void()> body) {
    registry().push_back(TestCase{name, std::move(body)});
}

std::string format_double(double value) {
    std::ostringstream out;
    out.precision(12);
    out << value;
    return out.str();
}

}  // namespace testing

int main() {
    const auto& tests = testing::registry();
    std::cout << "running " << tests.size() << " test cases\n\n";

    int failed = 0;
    for (const auto& test : tests) {
        try {
            test.body();
            std::cout << "  [ ok ]   " << test.name << '\n';
        } catch (const testing::Failure& failure) {
            ++failed;
            std::cout << "  [ FAIL ] " << test.name << '\n' << failure.message << '\n';
        } catch (const std::exception& error) {
            ++failed;
            std::cout << "  [ FAIL ] " << test.name << "\n    unexpected exception: "
                      << error.what() << '\n';
        } catch (...) {
            ++failed;
            std::cout << "  [ FAIL ] " << test.name << "\n    unexpected non-standard exception\n";
        }
    }

    const int passed = static_cast<int>(tests.size()) - failed;
    std::cout << '\n'
              << passed << " passed, " << failed << " failed, " << testing::check_count()
              << " assertions\n";
    return failed == 0 ? 0 : 1;
}
