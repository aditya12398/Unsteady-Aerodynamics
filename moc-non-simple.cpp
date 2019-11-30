//Unsteady Aerodynamics Assignment - Method of Characteristics
//Simulation of Non-simple characteristic regions

#include <iostream>
#include <cmath>
#include <string>

//Basic Thermodynamic parameters
double p3_p4 = 0.2, gamma = 1.4, R = 287.00;
double u3, T_ref = 300, a_ref;

struct characteristic_lines
{
    int index;
    double slope, u, T, a, p, J;
};
struct characteristic_lines *incident_line, *reflected_line;

void main()
{
    void write_gnuplot(std::string, double, double, double, double);
    void compute_point(int);
    void compute_line_data(double);
    
    int n;
    double u, a, T, p;

    std::string filename;
    std::cout << "Enter the number of characteristic lines needed: ";
    std::cin >> n;

    //Generating n characteristic lines
    incident_line = new characteristic_lines[n];
    reflected_line = new characteristic_lines[n];

    a_ref = sqrt(gamma * R * T_ref);
    u3 = (2 * a_ref / (gamma - 1)) * (1 - pow(p3_p4, (2 * gamma / (gamma - 1))));
    compute_line_data(n);
}

void compute_line_data(double n)
{
    double initial_slope, final_slope;
    initial_slope = - 1 / a_ref;
    final_slope = - 1 / (a_ref - u3);

    for(int i = 0; i < n; i++)
    {
        incident_line[i].slope = tan(atan(initial_slope) + i * (atan(final_slope) - atan(initial_slope))/(n - 1));
        incident_line[i].u = 2 / (gamma + 1) * (a_ref + incident_line[i].slope);
        incident_line[i].T = T_ref * (1 - 0.5 * (gamma -1) * pow((incident_line[i].u / a_ref), 2));
        incident_line[i].a = sqrt(gamma * R * incident_line[i].T);
        incident_line[i].J = incident_line[i].u - incident_line[i].a * 2 / (gamma - 1);
    }
}