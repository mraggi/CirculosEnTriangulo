#pragma once

#include "IntegerInterval.hpp"
#include "Projectors.hpp"
#include "Random.hpp"
#include <array>
#include <iostream>

// This class MINIMIZES! Remember that!

template <class Vector, class Cost>
class DifferentialEvolver
{
public:
    using Population = std::vector<Vector>;
    using idx = std::ptrdiff_t;

    DifferentialEvolver(Population InitialPopulation,
                        Cost f,
                        int num_untouchable = 0)
        : f(f), X(InitialPopulation), n_untouchable(num_untouchable)
    {
        if (population_size() < 4)
        {
            std::cerr << "Error: Population must have size at least 4!"
                      << std::endl;
            std::abort();
        }

        costs.reserve(X.size());

        std::transform(X.begin(),
                       X.end(),
                       std::back_inserter(costs),
                       [f](const Vector& x) { return f(x); });

        idx best_idx = std::min_element(costs.begin(), costs.end()) -
          costs.begin();

        best = X[best_idx];
        best_cost = costs[best_idx];
    };

    template <class Projector = Identity<Vector>>
    void step(double change_prob = 0.5,
              double force = 0.1,
              Projector project_to_manifold = Identity<Vector>())
    {
        long xsz = X.size() - n_untouchable;
        for (long i = 0; i < xsz; ++i) // leave someone alone
        {
            auto& x = X[i];
            double fx = costs[i];

            Vector y = morphed(x, change_prob, force);
            project_to_manifold(y);
            double fy = f(y);

            if (fy < fx)
            {
                costs[i] = fy;
                X[i] = y;
                if (fy < best_cost)
                {
                    best_cost = fy;
                    best = y;
                }
            }
        }
    }

    const Population& population() const { return X; }

    //     std::pair<Vector,double> best

    idx population_size() const { return X.size(); }

    void InsertIndividual(const Vector& x)
    {
        auto fx = f(x);
        if (fx < best_cost)
        {
            best_cost = fx;
            best = x;
        }
        X.emplace_back(x);
        costs.emplace_back(fx);
    }

private:
    std::array<idx, 3> pick_three() const
    {
        idx n = population_size();

        idx a = random.random_int<idx>(0, n);
        idx b = random.random_int<idx>(0, n);
        while (b == a)
        {
            b = random.random_int<idx>(0, n);
        }

        idx c = random.random_int<idx>(0, n);
        while (c == a || c == b)
        {
            c = random.random_int<idx>(0, n);
        }
        return {a, b, c};
    }

    Vector morphed(const Vector& x, double change_prob, double force)
    {
        auto [a, b, c] = pick_three();
        auto& A = X[a];
        auto& B = X[b];
        auto& C = X[c];

        idx R = random.random_int(0, population_size());
        Vector y = x;
        for (auto i : indices(x))
        {
            if (i == R || random.probability_of_true(change_prob))
            {
                y[i] = A[i] + force*(B[i] - C[i]);
            }
        }
        return y;
    }

private:
    Cost f;
    std::vector<Vector> X;
    std::vector<double> costs;
    int n_untouchable{0};

public:
    Vector best;
    double best_cost{10000000.0};

private:
    mutable Random random;
};
