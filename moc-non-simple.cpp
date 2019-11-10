//Unsteady Aerodynamics Assignment - Method of Characteristics
//Simulation of Non-simple characteristic regions

#include<iostream>
#include<cmath>
#include<string>

struct characteristic_lines
{
    int index;
    double slope, u, T, a, p, J;
};

void main()
{
    struct characteristic_lines *incident_line, *reflected_line;
    void write_gnuplot(double, double, double, double);
    void compute_point(int);
    int n;

    //Basic Thermodynamic parameters
    double p3_p4 = 0.2, gamma = 1.4, R = 287.00;
    double u3, T_ref = 300, a_ref;
    double u, a, T, p;
    std::string filename;
    std::cout << "Enter the number of characteristic lines needed: ";
    std::cin >> n;

    //Generating n characteristic lines
    incident_line=new characteristic_lines[n];
    reflected_line=new characteristic_lines[n];

    a_ref = sqrt(gamma * R * T_ref);
    u3 = (2 * a_ref / (gamma - 1)) * (1 - pow(p3_p4, (2 * gamma / (gamma - 1))));
    /*u = 2 / (gamma + 1) * (a_ref + x /t); //x and t yet to be defined. Syntax error inevitable until resolved; x = 0 refers to wall
    T = T_ref * (1 - 0.5 * (gamma -1) * pow((u / a_ref), 2));*/
}
