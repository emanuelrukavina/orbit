/**
 * Author: Emanuel Rukavina <emanuelrukavina@icloud.com>
 * 
 * This file is part of the Start Fresh mod for Kerbal Space Program (KSP).
 */

#pragma once
#include "math.hpp"
#include <array>
namespace liborbit::detail
{
    struct Fit
    {
        double fraction = 0, mass = 0, chi2 = 0, information = 0;
    };
    struct Problem
    {
        Vec3 d{}, n{};
        double distance = 0;
        std::array<double, 6> residual{};
        std::array<double, 6> inv_sigma{};
        double null_chi2 = 0;
    };
    Problem make_problem(const Input &input);
    Fit fit_segment(const Problem &problem, double lo, double hi);
}
