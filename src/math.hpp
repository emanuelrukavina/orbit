/**
 * Author: Emanuel Rukavina <emanuelrukavina@icloud.com>
 * 
 * This file is part of the Start Fresh mod for Kerbal Space Program (KSP).
 */

#pragma once
#include <liborbit.hpp>

namespace liborbit::detail
{
    constexpr double G = 6.67430e-11;
    Vec3 add(Vec3 a, Vec3 b);
    Vec3 sub(Vec3 a, Vec3 b);
    Vec3 scale(Vec3 a, double s);
    double dot(Vec3 a, Vec3 b);
    double norm(Vec3 a);
    void validate(const Input &in);
}
