//Unsteady Aerodynamics Assignment - Method of Characteristics
//Simulation of Non-simple characteristic regions

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

//Basic Thermodynamic parameters
double p3_p4 = 0.4, gma = 1.4, R = 287.00;
double u3, T_ref = 300, a_ref, p_ref, rho_ref = 1.225;

struct characteristic_lines
{
    int index;
    double slope, u, T, a, p, rho, J;
    double x1, t1;
};
struct characteristic_lines *incident_line, *reflected_line;

int main()
{
    void compute_point(int);                                         //Function computes the points in x-t plane and changes the data of intersecting characteristic lines
    void compute_line_data(int);                                     //Function Calculates initial data for all incident characteristic lines
    void print_data(int);

    int n;
    double u, a, T, p;

    std::string filename;
    std::cout << "Enter the number of characteristic lines needed: ";
    std::cin >> n;

    //Generating n characteristic lines
    incident_line = new characteristic_lines[n];
    reflected_line = new characteristic_lines[n];

    a_ref = sqrt(gma * R * T_ref);
    u3 = (2 * a_ref / (gma - 1)) * (1 - pow(p3_p4, (0.5 * (gma - 1) / gma)));
    p_ref = rho_ref * R * T_ref;
    compute_line_data(n);
    compute_point(n);
    //print_data(n);
    return 0;
}

void compute_line_data(int n)
{
    double initial_slope, final_slope;
    initial_slope = -a_ref;
    final_slope = -1 * (a_ref - u3);

    for (int i = 0; i < n; i++)
    {
        //Divide the lines at equal angles from initial to final slope
        incident_line[i].slope = tan(atan(initial_slope) + i * (atan(final_slope) - atan(initial_slope)) / (n - 1));

        //Relations can be found in "Modern Compressible flow by J.D. Anderson"
        incident_line[i].u = 2 / (gma + 1) * (a_ref + incident_line[i].slope);
        incident_line[i].T = T_ref * pow((1 - 0.5 * (gma - 1) * (incident_line[i].u / a_ref)), 2);
        incident_line[i].a = sqrt(gma * R * incident_line[i].T);
        incident_line[i].p = p_ref * pow((1 - 0.5 * (gma - 1) * (incident_line[i].u / a_ref)), (2 * gma / (gma - 1)));
        incident_line[i].rho = rho_ref * pow((1 - 0.5 * (gma - 1) * (incident_line[i].u / a_ref)), (2 / (gma - 1)));
        incident_line[i].J = incident_line[i].u - incident_line[i].a * 2 / (gma - 1);

        //Writing x1 and t1 to make the equation of the line
        incident_line[i].x1 = incident_line[i].t1 = 0;
    }
}

void compute_point(int n)
{
    void write_gnuplot(std::string, double, double, double, double, double); //Functions writes the data such that it can be viewed in GNUPLOT

    double x1 ,t1;
    double temp1, temp2;
    double temp3, temp4;
    for (int i = 0; i < n; i++)
    {
        //First interaction with the wall. Only u,a, and J are calculated for the time being.
        reflected_line[i].J = -1 * incident_line[i].J;
        reflected_line[i].u = incident_line[i].u = 0;
        reflected_line[i].a = incident_line[i].a = (reflected_line[i].J - incident_line[i].J) * (gma - 1) / 4;
        reflected_line[i].slope = reflected_line[i].u + reflected_line[i].a;
        incident_line[i].x1 = reflected_line[i].x1 = -1;
        incident_line[i].t1 = reflected_line[i].t1 = (1 / incident_line[i].slope) * (-1 - incident_line[i].x1) + incident_line[i].t1;
        write_gnuplot("./Incident_wave.tsv", incident_line[i].x1, incident_line[i].t1, incident_line[i].u, incident_line[i].a, incident_line[i].J);
        write_gnuplot("./Reflected_wave.tsv", reflected_line[i].x1, reflected_line[i].t1, reflected_line[i].u, reflected_line[i].a, reflected_line[i].J);
        for (int j = i + 1; j < n; j++)
        {
            temp1 = incident_line[j].slope * reflected_line[i].x1 - incident_line[j].x1 * reflected_line[i].slope;
            temp2 = (incident_line[j].slope * reflected_line[i].slope) * (incident_line[j].t1 - reflected_line[i].t1);

            x1 = (temp1 + temp2) / (incident_line[j].slope - reflected_line[i].slope);

            temp3 = reflected_line[i].x1 - incident_line[j].x1;
            temp4 = (incident_line[j].slope * incident_line[j].t1 - reflected_line[i].slope * reflected_line[i].t1);

            t1 = (temp3 + temp4) / (incident_line[j].slope - reflected_line[i].slope);
            //Mid-plane interactions beyond the wall. Only u,a, and J are calculated for the time being.
            incident_line[j].u = reflected_line[i].u = (reflected_line[i].J + incident_line[j].J) * 0.5;
            incident_line[j].a = reflected_line[i].a = (reflected_line[i].J - incident_line[j].J) * (gma - 1) / 4;
            reflected_line[i].J = reflected_line[i].u + reflected_line[i].a * 2 / (gma - 1);
            incident_line[j].J = incident_line[j].u - incident_line[j].a * 2 / (gma - 1);
            incident_line[j].slope = incident_line[j].u - incident_line[j].a;
            reflected_line[i].slope = reflected_line[i].u + reflected_line[i].a;

            incident_line[j].x1 = reflected_line[i].x1 = x1;
            incident_line[j].t1 = reflected_line[i].t1 = t1;
            write_gnuplot("./Incident_wave.tsv", incident_line[j].x1, incident_line[j].t1, incident_line[j].u, incident_line[j].a, incident_line[j].J);
            write_gnuplot("./Reflected_wave.tsv", reflected_line[i].x1, reflected_line[i].t1, reflected_line[i].u, reflected_line[i].a, reflected_line[i].J);
        }
    }
}
void print_data(int n)
{
    std::cout << "Incident Wave\n";
    std::cout << "u\ta\tJ\n";
    for (int i = 0; i < n; i++)
    {
        std::cout << incident_line[i].u << "\t" << incident_line[i].a << "\t" << incident_line[i].J << std::endl;
        std::cout << reflected_line[i].u << "\t" << reflected_line[i].a << "\t" << reflected_line[i].J << std::endl;
    }
}

void write_gnuplot(std::string filename, double x, double t, double u, double a, double J)
{
    std::ofstream outfile;
    outfile.open(filename, std::fstream::app);
    if (u == 0)// && J > 0)
        outfile << std::endl;
    outfile << x << "\t" << t << "\t" << u << "\t" << a << "\t" << J << std::endl;
    //if (u == 0 && J < 0)
    //    outfile << std::endl;
}