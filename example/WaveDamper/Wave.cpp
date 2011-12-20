#include "Wave.hpp"

#include <algorithm>
#include <cmath>

//---------------------------------------------------------------------------//
Wave::Wave(double x_min,
	   double x_max,
	   int num_x)
{
    // Create the grid.
    grid.resize(num_x);
    double x_size = (x_max - x_min) / (num_x - 1);

    std::vector<double>::iterator grid_iterator;
    int i = 0;

    for (grid_iterator = grid.begin();
	 grid_iterator != grid.end();
	 ++grid_iterator, ++i)
    {
	*grid_iterator = i*x_size + x_min;
    }

    // Set initial conditions.
    damping.resize(num_x);
    std::fill(damping.begin(), damping.end(), 0.0);
    f.resize(num_x);
    std::vector<double>::iterator f_iterator;
    for (f_iterator = f.begin(), grid_iterator = grid.begin();
	 f_iterator != f.end();
	 ++f_iterator, ++grid_iterator)
    {
	*f_iterator = cos( *grid_iterator );
    }
}

//---------------------------------------------------------------------------//
Wave::~Wave()
{ /* ... */ }

//---------------------------------------------------------------------------//
double Wave::solve()
{
    // Apply the dampened component.
    double l2_norm_residual = 0.0;
    double f_old = 0.0;
    std::vector<double>::iterator f_iterator;
    std::vector<double>::const_iterator damping_iterator;
    for (f_iterator = f.begin(), damping_iterator = damping.begin();
	 f_iterator != f.end();
	 ++f_iterator, ++damping_iterator)
    {
	f_old = f_iterator;
	*f_iterator -= *damping_iterator;
	l2_norm_residual += (*f_iterator - f_old)*(*f_iterator - f_old);
    }
    
    // Return the l2 norm of the local residual.
    return pow(l2_norm_residual, 0.5);
}

//---------------------------------------------------------------------------//
void Wave::output(int label)
{
    std::stringstream convert;
    convert << label;
    std::string filename = "time" + convert.str() + ".dat";
    std::ofstream output;
    output.open(&filename[0]);
    sor_output << "# f      time step\n";
    for (int i = 0; i < (int) f.size(); ++i)
    {
	output << f[i] << " " << grid[i] << "\n";
    }
    output.close();
}

//---------------------------------------------------------------------------//
