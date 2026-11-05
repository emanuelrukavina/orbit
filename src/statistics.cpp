/**
 * Author: Emanuel Rukavina <emanuelrukavina@icloud.com>
 * 
 * This file is part of the Start Fresh mod for Kerbal Space Program (KSP).
 */

#include "statistics.hpp"
#include <random>
#include <cmath>
#include <limits>
namespace liborbit::detail
{
    double bootstrap_p_value(const Problem &observed, double lo, double hi,
                             double improvement, std::size_t samples, std::uint64_t seed)
    {
        if (!samples)
            return std::numeric_limits<double>::quiet_NaN();
        std::mt19937_64 rng(seed);
        std::normal_distribution<double> standard(0.0, 1.0);
        std::size_t exceed = 0;
        for (std::size_t b = 0; b < samples; ++b)
        {
            Problem null_sample = observed;
            null_sample.null_chi2 = 0;
            for (int j = 0; j < 6; ++j)
            {
                null_sample.residual[j] = standard(rng) / null_sample.inv_sigma[j];
                null_sample.null_chi2 += std::pow(null_sample.residual[j] * null_sample.inv_sigma[j], 2);
            }
            const Fit fit = fit_segment(null_sample, lo, hi);
            if (null_sample.null_chi2 - fit.chi2 >= improvement - 1e-10)
                ++exceed;
        }
        return static_cast<double>(exceed + 1) / static_cast<double>(samples + 1);
    }
}
