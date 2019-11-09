//Unsteady Aerodynamics Assignment - Method of Characteristics
//Simulation of Non-simple characteristic regions

#include<iostream>
#include<cmath>

struct characteristic_lines
{
    int index;
    double slope, u, a, J;
};

void main()
{
    void write_gnuplot(double x, double t, double u, double a);
    void compute_point(int i);
    
}
