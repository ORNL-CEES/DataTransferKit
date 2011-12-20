#ifndef simple_example_damper_hpp
#define simple_example_damper_hpp

#include <vector>

//---------------------------------------------------------------------------//

class Damper
{
  private:

    std::vector<double> wave_data;
    std::vector<double> damping;
    std::vector<double> grid;

  public:

    Damper(double x_min,
	   double x_max,
	   int num_x);

    ~Damper();

    // Get a const reference to the local damping data.
    const std::vector<double>& get_damping()
    {
	return damping;
    }

    // Get a const reference to the local grid.
    const std::vector<double>& get_grid()
    {
	return grid;
    }

    // Get the wave data to apply damping to from an external source.
    std::vector<double>& set_wave_data()
    {
	return wave_data;
    }

    // Apply damping to the local problem.
    void solve();
};

//---------------------------------------------------------------------------//

#endif // simple_example_damper_hpp

//---------------------------------------------------------------------------//
