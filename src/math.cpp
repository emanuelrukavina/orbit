/**
 * Author: Emanuel Rukavina <emanuelrukavina@icloud.com>
 * 
 * This file is part of the Start Fresh mod for Kerbal Space Program (KSP).
 */

#include "math.hpp"
#include <cmath>
#include <stdexcept>

namespace liborbit::detail
{
    Vec3 add(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
    Vec3 sub(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
    Vec3 scale(Vec3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
    double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    double norm(Vec3 a) { return std::sqrt(dot(a, a)); }
    static bool finite(Vec3 v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }
    void validate(const Input &in)
    {
        for (const auto *b : {&in.a, &in.b})
        {
            if (!(std::isfinite(b->mass_kg) && b->mass_kg > 0 &&
                  std::isfinite(b->acceleration_sigma_m_s2) && b->acceleration_sigma_m_s2 > 0 &&
                  finite(b->position_m) && finite(b->measured_acceleration_m_s2) &&
                  finite(b->known_external_acceleration_m_s2)))
                throw std::invalid_argument("Body data must be finite; mass and acceleration sigma must be positive");
        }
        if (!(norm(sub(in.a.position_m, in.b.position_m)) > 0))
            throw std::invalid_argument("Body positions must differ");
        if (!(std::isfinite(in.min_fraction) && std::isfinite(in.max_fraction) &&
              in.min_fraction > 0 && in.min_fraction < in.max_fraction && in.max_fraction < 1))
            throw std::invalid_argument("Require 0 < min_fraction < max_fraction < 1");
    }
}
