/**
 * Author: Emanuel Rukavina <emanuelrukavina@icloud.com>
 * 
 * This file is part of the Start Fresh mod for Kerbal Space Program (KSP).
 */

#include <liborbit.hpp>
#include "fit.hpp"
#include "statistics.hpp"
#include <cmath>

namespace liborbit
{
    Result analyze(const Input &input)
    {
        detail::validate(input);
        const auto p = detail::make_problem(input);
        const auto fit = detail::fit_segment(p, input.min_fraction, input.max_fraction);
        Result r;
        r.estimated_mass_kg = fit.mass;
        r.fraction_from_a = fit.fraction;
        r.estimated_position_m = detail::add(input.a.position_m, detail::scale(p.d, fit.fraction));
        r.distance_from_a_m = fit.fraction * p.distance;
        r.distance_from_b_m = (1 - fit.fraction) * p.distance;
        r.null_chi_squared = p.null_chi2;
        r.candidate_chi_squared = fit.chi2;
        r.improvement_chi_squared = std::max(0.0, p.null_chi2 - fit.chi2);
        r.conditional_mass_sigma_kg = 1.0 / std::sqrt(fit.information);
        r.fixed_position_5sigma_mass_kg = 5 * r.conditional_mass_sigma_kg;
        r.global_p_value = detail::bootstrap_p_value(p, input.min_fraction, input.max_fraction,
                                                     r.improvement_chi_squared,
                                                     input.bootstrap_samples, input.random_seed);
        if (!input.bootstrap_samples)
        {
            r.classification = Classification::inconclusive;
            r.explanation = "A candidate was fitted, but significance was not calibrated (bootstrap disabled).";
        }
        else if (r.global_p_value <= 0.05 && fit.mass > 0)
        {
            r.classification = Classification::candidate_supported;
            r.explanation = "Residual accelerations favor a positive point mass on the segment under the stated noise and background assumptions; this is not proof of a body.";
        }
        else
        {
            r.classification = Classification::no_significant_evidence;
            r.explanation = "No statistically significant evidence for a single positive point mass on the searched segment; this does not rule out other bodies.";
        }
        return r;
    }
}
