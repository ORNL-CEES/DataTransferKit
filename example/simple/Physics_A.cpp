#include <vector>
#include <cmath>

class Physics_A
{
  private:

    std::vector<double> grid;
    std::vector<double> f;

  public:

    Physics_A(double x_min,
	      double x_max,
	      int num_x)
    {
	// Create the grid.
	f.resize(num_x);
	grid.resize(num_x);
	double x_size = (x_max - x_min) / (num_x - 1);

	std::vector<double>::iterator grid_iterator;
	int i = 0;

	for (grid_iterator = grid.begin();
	     grid_iterator != grid.end();
	     ++grid_iterator, ++i)
	{
	    *grid_iterator = i*x_size;
	}
    }

    void solve()
    {
	
    }

    void output(int label)
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
};
