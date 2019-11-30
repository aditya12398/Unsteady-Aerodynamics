//Unsteady Aerodynamics Assignment - Method of Characteristics
//Simulation of Non-simple characteristic regions

#include <iostream>
#include <cmath>
#include <string>

//Basic Thermodynamic parameters
double p3_p4 = 0.2, gamma = 1.4, R = 287.00;
double u3, T_ref = 300, a_ref, p_ref, rho_ref = 1.225;

struct characteristic_lines
{
    int index;
    double slope, u, T, a, p, rho, J;
};
struct characteristic_lines *incident_line, *reflected_line;

void main()
{
    void write_gnuplot(std::string, double, double, double, double); //Functions writes the data such that it can be viewed in GNUPLOT
    void compute_point(int);                                         //Function computes the points in x-t plane and changes the data of intersecting characteristic lines
    void compute_line_data(double);                                  //Function Calculates initial data for all incident characteristic lines

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
    p_ref = rho_ref * R * T_ref;
    compute_line_data(n);
}

void compute_line_data(double n)
{
    double initial_slope, final_slope;
    initial_slope = - a_ref;
    final_slope = -1 * (a_ref - u3);

    for (int i = 0; i < n; i++)
    {
        //Divide the lines at equal angles from initial to final slope
        incident_line[i].slope = tan(atan(initial_slope) + i * (atan(final_slope) - atan(initial_slope)) / (n - 1));

        //Relations can be found in "Modern Compressible flow by J.D. Anderson"
        incident_line[i].u = 2 / (gamma + 1) * (a_ref + incident_line[i].slope);
        incident_line[i].T = T_ref * pow((1 - 0.5 * (gamma - 1) * (incident_line[i].u / a_ref)), 2);
        incident_line[i].a = sqrt(gamma * R * incident_line[i].T);
        incident_line[i].p = p_ref * pow((1 - 0.5 * (gamma - 1) * (incident_line[i].u / a_ref)), (2 * gamma / (gamma - 1)));
        incident_line[i].rho = rho_ref * pow((1 - 0.5 * (gamma - 1) * (incident_line[i].u / a_ref)), (2 / (gamma - 1)));
        incident_line[i].J = incident_line[i].u - incident_line[i].a * 2 / (gamma - 1);
    }
}