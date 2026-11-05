#include <liborbit.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

// Use the liborbit namespace so we can write Input, Result,
// and analyze() without putting liborbit:: in front of them.
using namespace liborbit;

// Newton's gravitational constant in SI units.
constexpr double G = 6.67430e-11;

// Distance between the two known bodies in meters.
constexpr double D = 1.0e11;

// Mass of body A in kilograms.
constexpr double mA = 5.0e24;

// Mass of body B in kilograms.
constexpr double mB = 6.0e24;

// Print the results of any example in the same format.
void print_result(const Result& result)
{
    // Print the classification as its numeric enum value.
    std::cout << "Classification: "
              << static_cast<int>(result.classification) << '\n';

    // Print the estimated mass of the unseen body.
    std::cout << "Mass (kg): "
              << result.estimated_mass_kg << '\n';

    // Print the estimated position between A and B.
    std::cout << "Position fraction: "
              << result.fraction_from_a << '\n';

    // Print the bootstrap p-value.
    std::cout << "Bootstrap p-value: "
              << result.global_p_value << '\n';

    // Print how much the hidden-body model improved the fit.
    std::cout << "Fit improvement: "
              << result.improvement_chi_squared << '\n';

    // Print liborbit's explanation of the result.
    std::cout << result.explanation << "\n\n";
}

// Create the common input used by all examples.
Input make_input()
{
    // Create an empty input structure.
    Input in;

    // Set the masses of the known bodies.
    in.a.mass_kg = mA;
    in.b.mass_kg = mB;

    // Place body A at the origin.
    in.a.position_m = {0, 0, 0};

    // Place body B along the positive X axis.
    in.b.position_m = {D, 0, 0};

    // Set the measurement uncertainty for each body.
    in.a.acceleration_sigma_m_s2 = 1e-10;
    in.b.acceleration_sigma_m_s2 = 1e-10;

    // Run 999 bootstrap simulations for statistical testing.
    in.bootstrap_samples = 999;

    // Return the configured input.
    return in;
}

// EXAMPLE 1: HIDDEN BODY
// Detect a synthetic hidden body between A and B.
void example_hidden_body()
{
    // Print the example title.
    std::cout << "=== Example 1: Hidden body ===\n";

    // Define the hidden body's mass.
    constexpr double hidden_mass = 1.0e24;

    // Place the hidden body 60% of the distance from A to B.
    constexpr double f = 0.6;

    // Create the standard input.
    Input in = make_input();

    // Synthetic measurements: mutual gravity plus one unseen body at f*D.
    in.a.measured_acceleration_m_s2 = {
        G * mB / (D * D) +
        G * hidden_mass / (f * f * D * D),
        0, 0
    };

    // Body B is pulled toward A and toward the hidden body.
    in.b.measured_acceleration_m_s2 = {
        -G * mA / (D * D) -
        G * hidden_mass / ((1 - f) * (1 - f) * D * D),
        0, 0
    };

    // Analyze the measurements.
    const Result result = analyze(in);

    // Print the estimated mass, position, and significance.
    print_result(result);
}

// EXAMPLE 2: NO HIDDEN BODY
// Both bodies behave exactly as predicted by mutual gravity.
void example_no_hidden_body()
{
    // Print the example title.
    std::cout << "=== Example 2: No hidden body ===\n";

    // Create the standard input.
    Input in = make_input();

    // Body A accelerates only because of body B.
    in.a.measured_acceleration_m_s2 = {
        G * mB / (D * D),
        0, 0
    };

    // Body B accelerates only because of body A.
    in.b.measured_acceleration_m_s2 = {
        -G * mA / (D * D),
        0, 0
    };

    // Analyze the measurements.
    const Result result = analyze(in);

    // Print the result, which should show no improvement.
    print_result(result);
}

// EXAMPLE 3: SMALL HIDDEN BODY
// Test how liborbit behaves when the signal is weaker.
void example_small_hidden_body()
{
    // Print the example title.
    std::cout << "=== Example 3: Small hidden body ===\n";

    // Use a hidden mass much smaller than in example 1.
    constexpr double hidden_mass = 1.0e20;

    // Place the hidden mass at the midpoint.
    constexpr double f = 0.5;

    // Create the standard input.
    Input in = make_input();

    // Calculate body A's acceleration with the small hidden mass.
    in.a.measured_acceleration_m_s2 = {
        G * mB / (D * D) +
        G * hidden_mass / (f * f * D * D),
        0, 0
    };

    // Calculate body B's acceleration with the small hidden mass.
    in.b.measured_acceleration_m_s2 = {
        -G * mA / (D * D) -
        G * hidden_mass / ((1 - f) * (1 - f) * D * D),
        0, 0
    };

    // Analyze the measurements.
    const Result result = analyze(in);

    // Print the result.
    print_result(result);
}

// EXAMPLE 4: HIDDEN BODY CLOSER TO A
// See how the inferred position changes.
void example_body_near_a()
{
    // Print the example title.
    std::cout << "=== Example 4: Hidden body near A ===\n";

    // Define the hidden body's mass.
    constexpr double hidden_mass = 1.0e24;

    // Place it 20% of the distance from A toward B.
    constexpr double f = 0.2;

    // Create the standard input.
    Input in = make_input();

    // Body A experiences stronger hidden-body gravity
    // because the unseen mass is closer to A.
    in.a.measured_acceleration_m_s2 = {
        G * mB / (D * D) +
        G * hidden_mass / (f * f * D * D),
        0, 0
    };

    // Body B experiences a weaker hidden-body contribution.
    in.b.measured_acceleration_m_s2 = {
        -G * mA / (D * D) -
        G * hidden_mass / ((1 - f) * (1 - f) * D * D),
        0, 0
    };

    // Analyze the measurements.
    const Result result = analyze(in);

    // Print the estimated position.
    print_result(result);
}

// EXAMPLE 5: MEASUREMENT NOISE
// Add small errors to the ideal accelerations.
void example_noisy_measurements()
{
    // Print the example title.
    std::cout << "=== Example 5: Noisy measurements ===\n";

    // Define the hidden body's mass.
    constexpr double hidden_mass = 1.0e24;

    // Place the hidden body 60% of the distance from A.
    constexpr double f = 0.6;

    // Create the standard input.
    Input in = make_input();

    // Add a small fixed measurement error to A's acceleration.
    in.a.measured_acceleration_m_s2 = {
        G * mB / (D * D) +
        G * hidden_mass / (f * f * D * D) + 2e-10,
        -1e-10,
        1e-10
    };

    // Add different measurement errors to B's acceleration.
    in.b.measured_acceleration_m_s2 = {
        -G * mA / (D * D) -
        G * hidden_mass / ((1 - f) * (1 - f) * D * D) - 1e-10,
        2e-10,
        -1e-10
    };

    // Analyze the imperfect measurements.
    const Result result = analyze(in);

    // Print the estimated mass and position.
    print_result(result);
}

// EXAMPLE 6: KNOWN EXTERNAL ACCELERATION
// Account for another known gravitational influence.
void example_external_gravity()
{
    // Print the example title.
    std::cout << "=== Example 6: External gravity ===\n";

    // Define the hidden body's mass.
    constexpr double hidden_mass = 1.0e24;

    // Define its position between A and B.
    constexpr double f = 0.6;

    // Create the standard input.
    Input in = make_input();

    // Define a known external acceleration acting on A.
    in.a.known_external_acceleration_m_s2 = {
        1e-6, 2e-7, 0
    };

    // Define a different known external acceleration acting on B.
    in.b.known_external_acceleration_m_s2 = {
        8e-7, -1e-7, 0
    };

    // The observed acceleration of A includes mutual gravity,
    // hidden-body gravity, and known external acceleration.
    in.a.measured_acceleration_m_s2 = {
        G * mB / (D * D) +
        G * hidden_mass / (f * f * D * D) + 1e-6,
        2e-7,
        0
    };

    // The observed acceleration of B includes the same
    // categories of forces, but with its own directions.
    in.b.measured_acceleration_m_s2 = {
        -G * mA / (D * D) -
        G * hidden_mass / ((1 - f) * (1 - f) * D * D) + 8e-7,
        -1e-7,
        0
    };

    // Analyze the measurements after accounting for known forces.
    const Result result = analyze(in);

    // Print the result.
    print_result(result);
}

// EXAMPLE 7: DISABLE BOOTSTRAP
// Fit the mass and position without statistical calibration.
void example_without_bootstrap()
{
    // Print the example title.
    std::cout << "=== Example 7: Without bootstrap ===\n";

    // Define the hidden body's mass.
    constexpr double hidden_mass = 1.0e24;

    // Define its position.
    constexpr double f = 0.6;

    // Create the standard input.
    Input in = make_input();

    // Set the measurements for body A.
    in.a.measured_acceleration_m_s2 = {
        G * mB / (D * D) +
        G * hidden_mass / (f * f * D * D),
        0, 0
    };

    // Set the measurements for body B.
    in.b.measured_acceleration_m_s2 = {
        -G * mA / (D * D) -
        G * hidden_mass / ((1 - f) * (1 - f) * D * D),
        0, 0
    };

    // Disable bootstrap simulations.
    in.bootstrap_samples = 0;

    // Fit the candidate without computing a calibrated p-value.
    const Result result = analyze(in);

    // Print the result; classification should be inconclusive.
    print_result(result);
}

// EXAMPLE 8: RESTRICTED POSITION SEARCH
// Only search the middle part of the A-B segment.
void example_restricted_search()
{
    // Print the example title.
    std::cout << "=== Example 8: Restricted search ===\n";

    // Define the hidden body's mass.
    constexpr double hidden_mass = 1.0e24;

    // Place the hidden body at 60% of the segment.
    constexpr double f = 0.6;

    // Create the standard input.
    Input in = make_input();

    // Search only from 40% to 80% of the A-B distance.
    in.min_fraction = 0.4;
    in.max_fraction = 0.8;

    // Calculate the observed acceleration of A.
    in.a.measured_acceleration_m_s2 = {
        G * mB / (D * D) +
        G * hidden_mass / (f * f * D * D),
        0, 0
    };

    // Calculate the observed acceleration of B.
    in.b.measured_acceleration_m_s2 = {
        -G * mA / (D * D) -
        G * hidden_mass / ((1 - f) * (1 - f) * D * D),
        0, 0
    };

    // Analyze only the configured search interval.
    const Result result = analyze(in);

    // Print the result.
    print_result(result);
}

// Run all examples from a single executable.
int main()
{
    // Use scientific notation for large and small values.
    std::cout << std::scientific << std::setprecision(6);

    // Run the original synthetic hidden-body example.
    example_hidden_body();

    // Run the no-hidden-body example.
    example_no_hidden_body();

    // Run the small-hidden-body example.
    example_small_hidden_body();

    // Run the example with a hidden body closer to A.
    example_body_near_a();

    // Run the example with imperfect acceleration measurements.
    example_noisy_measurements();

    // Run the example with known external gravitational influences.
    example_external_gravity();

    // Run the example without bootstrap significance testing.
    example_without_bootstrap();

    // Run the example with a restricted position search.
    example_restricted_search();

    // Exit successfully.
    return 0;
}
