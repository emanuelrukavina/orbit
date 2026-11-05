/**
 * Author: Emanuel Rukavina <emanuelrukavina@icloud.com>
 * 
 * This file is part of the Start Fresh mod for Kerbal Space Program (KSP).
 */

#pragma once
#include "fit.hpp"
#include <cstdint>
#include <cstddef>
namespace liborbit::detail
{
    double bootstrap_p_value(const Problem &observed, double lo, double hi,
                             double improvement, std::size_t samples, std::uint64_t seed);
}
