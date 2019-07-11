#pragma once

#include <CGAL/Simple_cartesian.h>
#include <algorithm>
#include <discreture.hpp>
#include <iostream>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Vector_2 = Kernel::Vector_2;
using Line_2 = Kernel::Line_2;
using FT = Kernel::FT;

Point_2 Origin(0, 0);

/*
Funciones importantes
--------------------------------
Para caso afín:
rescale_TK(F,k) - Uniformly rescales circles so that they have T(k), use when
you want same radii greedy_rescale_TK(F,k) - Rescales circles so that they have
T(k), use if you don't care about same radii min_inflate_sq(F) - The square of
the minimum inflation so that the circles have T
--------------------------------
Para caso centralmente simétrico:
cs_greedy_rescale_TK(F,k) - Rescales circles so that they have T_O(k)
cs_min_inflate_sq(F) - The square of the minimum inflation so that the circles
have T_O
--------------------------------
<< imprime Circle y Sphere en formato de GeoGebra
*/

struct Circle
{
    Point_2 c;
    FT r;
};

std::ostream& operator<<(std::ostream& os, const Circle& C)
{
    os << "Circle((" << C.c.x() << "," << C.c.y() << ")," << C.r << ")";
    return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& C)
{
    os << '\n';
    for (auto& c : C)
        os << c << '\n';
    return os;
}

Line_2 transversal_candidate_bc(const Circle& A, const Circle& B, const Circle& C)
{
    return {A.c + A.r*(B.c - A.c)/(A.r + B.r),
            A.c + A.r*(C.c - A.c)/(A.r + C.r)}; // Los puntos de los
                                                    // segmentos
}

std::array<Line_2, 3> transversal_candidates(const Circle& A,
                                             const Circle& B,
                                             const Circle& C)
{
    return {transversal_candidate_bc(A, B, C),
            transversal_candidate_bc(B, A, C),
            transversal_candidate_bc(C, A, B)};
}

FT min_inflate_sq(const std::vector<Circle>& F, Line_2 L)
{
    FT worst = 0;
    for (const Circle& circle : F)
    {
        auto d2 = squared_distance(circle.c, L);
        auto ratio = d2/(circle.r*circle.r);
        if (ratio > worst)
            worst = ratio;
    }

    return worst;
}

FT min_inflate_sq(const std::vector<Circle>& F)
{
    FT best = 99999999999.0;
    for (auto&& vc : discreture::combinations(F, 3))
    {
        auto L = transversal_candidates(vc[0], vc[1], vc[2]);
        for (auto&& l : L)
        {
            FT candidate = min_inflate_sq(F, l);
            if (candidate < best)
                best = candidate;
        }
    }
    return best;
}

FT min_inflate_TK_sq(const std::vector<Circle>& F, int k)
{
    FT worst = 0.0;
    for (auto&& kada : discreture::combinations(F, k))
    {
        std::vector<Circle> T(kada.begin(), kada.end());
        FT lll = min_inflate_sq(T);
        if (lll > worst)
            worst = lll;
    }
    return worst;
}

void rescale_TK(std::vector<Circle>& F, int k)
{
    FT l = sqrt(min_inflate_TK_sq(F, k));
    for (auto& circle : F)
        circle.r *= l;
}

void greedy_rescale_TK(std::vector<Circle>& F, int k)
{
    int N = F.size();
    std::vector<FT> worst(F.size(), 0.0);
    for (auto&& kada : discreture::combinations(N, k))
    {
        std::vector<Circle> T = discreture::compose(F, kada);
        FT lll = min_inflate_sq(T);
        for (auto i : kada)
        {
            if (worst[i] < lll)
                worst[i] = lll;
        }
    }

    for (int i = 0; i < N; ++i)
    {
        F[i].r *= sqrt(worst[i]);
    }
}

// ***** Centrally Symmetric Part *****

Line_2 cs_transversal_candidate(Point_2 a, Point_2 b, FT r)
{
    return {Origin, a + r*(b - a)};
}

Line_2 cs_transversal_candidate(Circle A, Circle B)
{
    return cs_transversal_candidate(A.c, B.c, A.r/(A.r + B.r));
}

FT cs_min_inflate_sq(const std::vector<Circle>& F)
{
    FT best = 99999999999.0;
    for (auto&& vc : discreture::combinations(F, 2))
    {
        Line_2 L = cs_transversal_candidate(vc[0], vc[1]);
        FT candidate = min_inflate_sq(F, L); // Ta bien sin cs
        if (candidate < best)
            best = candidate;
        Circle reflejado = {{-vc[1].c.x(), -vc[1].c.y()}, vc[1].r};
        L = cs_transversal_candidate(vc[0], reflejado);
        candidate = min_inflate_sq(F, L);
        if (candidate < best)
            best = candidate;
    }
    return best;
}

FT cs_min_inflate_TK_sq(const std::vector<Circle>& F, int k)
{
    FT worst = 0.0;

    for (auto&& kada : discreture::combinations(F, k))
    {
        std::vector<Circle> T(kada.begin(), kada.end());
        FT lll = cs_min_inflate_sq(T);
        if (lll > worst)
            worst = lll;
    }

    return worst;
}

void cs_rescale_TK(std::vector<Circle>& F, int k)
{
    FT l = sqrt(cs_min_inflate_TK_sq(F, k));
    for (auto& circle : F)
        circle.r *= l;
}

void cs_greedy_rescale_TK(std::vector<Circle>& F, int k)
{
    int N = F.size();
    std::vector<FT> worst(F.size(), 0.0);
    for (auto&& kada : discreture::combinations(N, k))
    {
        std::vector<Circle> T = discreture::compose(F, kada);
        FT lll = cs_min_inflate_sq(T);
        for (auto i : kada)
        {
            if (worst[i] < lll)
                worst[i] = lll;
        }
    }

    for (int i = 0; i < N; ++i)
    {
        F[i].r *= sqrt(worst[i]);
    }
}