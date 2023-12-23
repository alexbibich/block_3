#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

#include <iostream>
#include <fstream>
#include <filesystem>

#include "problems/test_problems.h"

inline std::string get_test_string() {
    auto test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    auto test_string = std::string(test_info->test_case_name()) + "." + std::string(test_info->name());
    return test_string;
}

inline std::string prepare_test_folder()
{
    std::string path = std::string("output/");
    std::filesystem::create_directories(path);
    return path;
}

int main(int argc, char** argv)
{
    prepare_test_folder();
    ::testing::InitGoogleTest(&argc, argv);
#ifdef _WIN32
    std::wcout.imbue(std::locale("rus_rus.866"));
#endif
    int res = RUN_ALL_TESTS();
    return res;
}

