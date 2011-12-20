#ifndef simple_example_damper_hpp
#define simple_example_damper_hpp

#include <vector>

//---------------------------------------------------------------------------//

class Damper
{
  private:

    std::vector<double> wave_data;
    std::vector<double> damping;

  public:

    Damper();

    ~Damper();

    // Get a const reference to the local damping data.
    const std::vector<double>& get_damping()
    {
	return damping;
    }

    // Get the wave data to apply damping to from an external source.
    void set_wave_data( const std::vector<double>& external_wave )
    {
	wave = external_wave;
    }

    // Apply damping to the local problem.
    void solve();
};

//---------------------------------------------------------------------------//

#endif // simple_example_damper_hpp

//---------------------------------------------------------------------------//
