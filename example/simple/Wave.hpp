#ifndef simple_example_wave_hpp
#define simple_example_wave_hpp

#include <vector>

//---------------------------------------------------------------------------//

class Wave
{
  private:

    std::vector<double> grid;
    std::vector<double> f;
    std::vector<double> damping;

  public:

    Wave(double x_min,
	 double x_max,
	 int num_x);

    ~Wave();

    // Get a const reference to the local grid.
    const std::vector<double>& get_grid()
    {
	return grid;
    }

    // Get a const reference to the local data.
    const std::vector<double>& get_f()
    {
	return f;
    }

    // Apply the damping to the local data structures from an external
    // source. 
    std::vector<double>& set_damping();
    {
	return damping;
    }

    // Solve the local problem and return the l2 norm of the local residual.
    double solve();
};

//---------------------------------------------------------------------------//

#endif // simple_example_wave_hpp

//---------------------------------------------------------------------------//
