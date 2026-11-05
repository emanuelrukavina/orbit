#include <liborbit.hpp>
#include <cassert>
#include <cmath>
#include <stdexcept>

int main()
{
    using namespace liborbit;
    constexpr double G = 6.67430e-11, D = 1e11, mA = 5e24, mB = 6e24, mC = 1e24, f = 0.6;
    Input in;
    in.a.mass_kg = mA;
    in.b.mass_kg = mB;
    in.a.position_m = {0, 0, 0};
    in.b.position_m = {D, 0, 0};
    in.a.acceleration_sigma_m_s2 = 1e-10;
    in.b.acceleration_sigma_m_s2 = 1e-10;
    in.a.measured_acceleration_m_s2 = {G * mB / (D * D) + G * mC / (f * f * D * D), 0, 0};
    in.b.measured_acceleration_m_s2 = {-G * mA / (D * D) - G * mC / ((1 - f) * (1 - f) * D * D), 0, 0};
    in.bootstrap_samples = 99;
    auto r = analyze(in);
    assert(std::abs(r.fraction_from_a - f) < 1e-5);
    assert(std::abs(r.estimated_mass_kg / mC - 1) < 1e-5);
    assert(r.classification == Classification::candidate_supported);
    // A no-hidden-body case should not trigger a detection.
    in.a.measured_acceleration_m_s2 = {G * mB / (D * D), 0, 0};
    in.b.measured_acceleration_m_s2 = {-G * mA / (D * D), 0, 0};
    r = analyze(in);
    assert(r.classification == Classification::no_significant_evidence);
    assert(r.estimated_mass_kg < 1e10);
    // Invalid uncertainty must be rejected.
    in.a.acceleration_sigma_m_s2 = 0;
    bool threw = false;
    try
    {
        (void)analyze(in);
    }
    catch (const std::invalid_argument &)
    {
        threw = true;
    }
    assert(threw);
}
