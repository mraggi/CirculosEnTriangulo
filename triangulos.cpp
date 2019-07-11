#include "DifferentialEvolution.hpp"

#include "include/Point.hpp"
#include "include/TimeHelpers.hpp"
#include <cmath>
// #include <discreture.hpp>
#include <iomanip>

using Individual = std::vector<Point>;

std::array<Point,3> T;
Point REFLECT;
Random RANDOM;


bool is_in_triangle(const Point& P)
{
    if (!P.IsToTheLeftOfLine(T[0], T[1])) return false;
    if (!P.IsToTheLeftOfLine(T[1], T[2])) return false;
    if (!P.IsToTheLeftOfLine(T[2], T[0])) return false;

    return true;
}

Point random_point_in_triangle()
{
    Point P = {RANDOM.random_real(0., 1.), RANDOM.random_real(0., 1.)};
    Point W = P.x*T[1] + P.y*T[2];
	if (P.x + P.y > 1.)
		return REFLECT-W;
	else
		return W;
}

Individual random_points_in_triangle(int n)
{
    Individual I(n);
    std::generate(I.begin(), I.end(), []() { return random_point_in_triangle(); });
    return I;
}

void GetInTriangle(Individual& X)
{
    for (auto& x : X)
    {
        if (!is_in_triangle(x))
        {
            x = random_point_in_triangle();
        }
    }
}

double get_radius(const Individual& X)
{
    int n = X.size();
    double to_others = 1.;
    for (int j = 1; j < n; ++j)
    {
        for (int i = 0; i < j; ++i)
        {
            double d2 = distancesq(X[i], X[j]);
            if (d2 < to_others)
                to_others = d2;
        }
    }

	double to_sides = 1.;
    for (auto& x : X)
    {
        for (int i = 0; i < 3; ++i)
        {
            int j = (i + 1)%3;
            double d = x.DistanceSqToLine(T[i], T[j]);
            if (d < to_sides)
                to_sides = d;
        }
    }
    
    double best = std::min(to_others, 4*to_sides);

    return std::sqrt(best)/2;
}

int main()
{
    Chronometer C;
    //     std::cout << "Different radii\nk, num_circles, pop_size,
    //     num_epochs\n"; std::cin >> k >> family_size >> pop_size >> num_epochs;
	T[0] = {0.,0.};
	T[1] = {1., 0.};
	T[2] = {0.5, std::sqrt(3.0)/2};
	REFLECT = (T[1]+T[2]);
	
    int k = 13;
    int num_epochs = 10000;
	int meta_epochs = 30;
	
	std::cout << "k, num_epochs, meta_epochs\n";
    std::cout << k << ' ' << num_epochs << ' ' << meta_epochs << std::endl;
	
	Individual best_overall_individual = random_points_in_triangle(k);
	double best_overall_cost = 0.;
	
	for (int meta_epoch : NN(meta_epochs))
	{
		std::cout << "STARTING META EPOCH " << meta_epoch << std::endl;
		int pop_size = RANDOM.random_int(20, 100);
		
		std::vector<Individual> Population(pop_size);
		std::generate(Population.begin(), Population.end(), [k]() {
			return random_points_in_triangle(k);
		});
		
		Population[RANDOM.random_int(0,pop_size)] = best_overall_individual;
		
		DifferentialEvolver D(Population, [](const Individual& X) {
			return -get_radius(X); // quiero MAXIMIZAR la distancia mínima.
		});

		for (int epoch : NN(num_epochs))
		{
			double prob = RANDOM.random_real(0.1, 0.9);
			double force = RANDOM.random_real(0.1, 0.9);
			D.step(prob, force, GetInTriangle);
			if ((epoch&2047) == 0)
			{
				if (D.best_cost < best_overall_cost)
				{
					best_overall_cost = D.best_cost;
					best_overall_individual = D.best;
					std::cout << '(' << meta_epoch << '/' << meta_epochs << ')' << "[" << epoch << '/' << num_epochs << "] Best is now " << -best_overall_cost << std::endl;
				}
				else if (RANDOM.random_real(0.,1.) < 0.05)
				{
					D.InsertIndividual(random_points_in_triangle(k));
				}
			}
		}
		num_epochs *= 1.3;
		num_epochs = std::min(num_epochs,1000000);
	}
    

    std::cout << std::setprecision(10);
    double r = -best_overall_cost;
    std::cout << "Best achieved: " << r << "\n\n\n\n";

    std::cout << "P = line([" << T[0] << "," << T[1] << ", " << T[2] << ", "
              << T[0] << "],color=\"black\")\n";

    for (auto& x : best_overall_individual)
    {
        std::cout << "P += circle(" << x << ", " << r << ")\n";
    }
    std::cout << "P += text(\"r = " << r << "\",(0.2,0.8),color=\"red\")" << std::endl;
    std::cout << "P.show()" << std::endl;
    std::cout << std::endl;
    std::cout << "Total time taken: " << C.Peek() << "s\n";

    return 0;
}
