#include <gtest/gtest.h>

#include <valik/shared.hpp>

TEST(adjust_test, try)
{
    std::vector<uint32_t> bin_count_in{8, 63, 64, 65, 96, 97, 159};
    for (auto b : bin_count_in)
    {
        uint32_t adjusted_count = valik::adjust_bin_count(b);

        if (b == std::numeric_limits<uint32_t>::max())
            EXPECT_EQ(adjusted_count, 64UL);    //default
        else if (b < 97UL)
            EXPECT_EQ(adjusted_count, 64UL);
        else
            EXPECT_EQ(adjusted_count, 128UL);
    }
}
