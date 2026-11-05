/**
 * Author: Emanuel Rukavina <emanuelrukavina@icloud.com>
 * 
 * This file is part of the Start Fresh mod for Kerbal Space Program (KSP).
 */

#include "fit.hpp"
#include <cmath>
#include <limits>

namespace liborbit::detail
{
    Problem make_problem(const Input &in)
    {
        Problem p;
        p.d = sub(in.b.position_m, in.a.position_m);
        p.distance = norm(p.d);
        p.n = scale(p.d, 1.0 / p.distance);
        const double accA = G * in.b.mass_kg / (p.distance * p.distance);
        const double accB = G * in.a.mass_kg / (p.distance * p.distance);
        const Vec3 rA = sub(sub(in.a.measured_acceleration_m_s2,
                                scale(p.n, accA)),
                            in.a.known_external_acceleration_m_s2);
        const Vec3 rB = sub(add(in.b.measured_acceleration_m_s2,
                                scale(p.n, accB)),
                            in.b.known_external_acceleration_m_s2);
        p.residual = {rA.x, rA.y, rA.z, rB.x, rB.y, rB.z};
        const double wa = 1.0 / in.a.acceleration_sigma_m_s2;
        const double wb = 1.0 / in.b.acceleration_sigma_m_s2;
        p.inv_sigma = {wa, wa, wa, wb, wb, wb};
        for (int i = 0; i < 6; ++i)
            p.null_chi2 += std::pow(p.residual[i] * p.inv_sigma[i], 2);
        return p;
    }

    static Fit fit_at(const Problem &p, double f)
    {
        const double a = G / (p.distance * p.distance * f * f);
        const double b = -G / (p.distance * p.distance * (1 - f) * (1 - f));
        const std::array<double, 6> shape = {a * p.n.x, a * p.n.y, a * p.n.z, b * p.n.x, b * p.n.y, b * p.n.z};
        double information = 0, cross = 0;
        for (int i = 0; i < 6; ++i)
        {
            const double weighted = shape[i] * p.inv_sigma[i];
            information += weighted * weighted;
            cross += weighted * (p.residual[i] * p.inv_sigma[i]);
        }
        const double mass = std::max(0.0, cross / information);
        const double improvement = mass > 0 ? cross * cross / information : 0;
        return {f, mass, std::max(0.0, p.null_chi2 - improvement), information};
    }
    Fit fit_segment(const Problem &p, double lo, double hi)
    {
        // Dense scan plus local golden-section refinement. Avoids assuming unimodality globally.
        constexpr int grid = 400;
        Fit best = fit_at(p, lo);
        int best_idx = 0;
        for (int i = 1; i <= grid; ++i)
        {
            Fit f = fit_at(p, lo + (hi - lo) * i / grid);
            if (f.chi2 < best.chi2)
            {
                best = f;
                best_idx = i;
            }
        }
        double left = lo + (hi - lo) * std::max(0, best_idx - 1) / grid;
        double right = lo + (hi - lo) * std::min(grid, best_idx + 1) / grid;
        constexpr double phi = 0.6180339887498948482;
        double x1 = right - phi * (right - left), x2 = left + phi * (right - left);
        Fit f1 = fit_at(p, x1), f2 = fit_at(p, x2);
        for (int i = 0; i < 65; ++i)
        {
            if (f1.chi2 < f2.chi2)
            {
                right = x2;
                x2 = x1;
                f2 = f1;
                x1 = right - phi * (right - left);
                f1 = fit_at(p, x1);
            }
            else
            {
                left = x1;
                x1 = x2;
                f1 = f2;
                x2 = left + phi * (right - left);
                f2 = fit_at(p, x2);
            }
        }
        if (f1.chi2 < best.chi2)
            best = f1;
        if (f2.chi2 < best.chi2)
            best = f2;
        return best;
    }
}
