/**
 * Author: Emanuel Rukavina <emanuelrukavina@icloud.com>
 * 
 * This file is part of the Start Fresh mod for Kerbal Space Program (KSP).
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

namespace liborbit
{

    struct Vec3
    {
        double x = 0.0, y = 0.0, z = 0.0;
    };

    struct BodyObservation
    {
        double mass_kg = 0.0;
        Vec3 position_m{};
        Vec3 measured_acceleration_m_s2{};
        Vec3 known_external_acceleration_m_s2{};

        // Standard deviation of EACH acceleration component; independent Gaussian errors.
        double acceleration_sigma_m_s2 = 0.0;
    };

    struct Input
    {
        BodyObservation a;
        BodyObservation b;

        // Search candidate position A + fraction*(B-A), within these bounds.
        double min_fraction = 0.01;
        double max_fraction = 0.99;

        // Bootstrap draws under no-hidden-body hypothesis; 0 disables p-value.
        std::size_t bootstrap_samples = 999;
        std::uint64_t random_seed = 42;
    };

    enum class Classification
    {
        candidate_supported,
        no_significant_evidence,
        inconclusive
    };

    struct Result
    {
        Classification classification = Classification::inconclusive;
        std::string explanation;
        double estimated_mass_kg = 0.0;
        Vec3 estimated_position_m{};
        double fraction_from_a = 0.0;
        double distance_from_a_m = 0.0;
        double distance_from_b_m = 0.0;
        double null_chi_squared = 0.0;
        double candidate_chi_squared = 0.0;
        double improvement_chi_squared = 0.0;

        // Monte Carlo p-value for the complete one-dimensional position search.
        // NaN if bootstrap_samples == 0.
        double global_p_value = 0.0;

        // Approximate 5-sigma mass sensitivity at fitted location, position fixed.
        double fixed_position_5sigma_mass_kg = 0.0;
        
        // Conditional 1-sigma mass error, holding fitted position fixed.
        double conditional_mass_sigma_kg = 0.0;
    };

    // Throws std::invalid_argument for invalid inputs.
    // This analyzes ONE epoch and ONE positive point mass on the A-B segment.
    // It does not fit a time-dependent orbit or prove that a real object exists.
    Result analyze(const Input &input);

} // namespace liborbit
