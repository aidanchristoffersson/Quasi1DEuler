#include <iostream>
#include <cmath> // to get trig and power functions
#include <numbers> // to get pi
#include <vector>
#include <fstream>

// ideal gas assumption for fluid flow
const double gamma {1.4};
const double R {1716}; // units = ft lb / (slug R)

// constant static variables
const double Tt = 531.2; // unit = R
constexpr double Pt = 2117; // unit = lb / ft^2

// commonly used variables
const double cv {R / (gamma - 1)};
const double gamma_ratio {(gamma - 1) / (gamma + 1)};
const double gamma_ratio_sub {(gamma - 1) / (gamma)};
const double alpha_star {2 * gamma * gamma_ratio * cv * Tt};

// area parameters
const static double h = 0.15; // nozzle depth
const static double t1 = 0.8; // location
const static double t2 = 3;  // bump width

// nozzle length parameters
static constexpr double a {0};
static constexpr double b {1};
// I would have put these variables inside the classes for better organization/clarity, but I believe that they would be
// copied as I am using templates, so I thought this would be more efficient.

// function initialization
void Log(auto p_val) { std::cout << p_val << std::endl; } // used to print to console for debugging

// Struct initialization
struct InputParameters
{
    const double PressureRatio, M0, CFL, Epsilon;
};

// Class initialization:
template<size_t G>
class Grid
{
private: // variable parameters, only dependent on size

    double dx;
    size_t m_size;

    double m_nodes[(2 * G + 1)] {}; // initialize array of nodes, want G node points and G+1 boundary points
    double m_area[(2 * G + 1)] {};  // final index = 2 * G

public:
    // constructor to initialize object
    Grid()
        : dx((b - a) / G), m_size((2 * G) + 1) // cell size is dx, number of nodes is 2 * G + 1 w/ dist 0.5 * dx
    {
        auto areaCalc = [&] (double x)
        { return 1 - h * pow(sin(std::numbers::pi * pow(x, t1)), t2); }; //

        // loop to initialize cell values and boundaries as well as the area at all these points
        for (int i = 0; i < m_size; i++) // # of grid nodes = 2 * G + 1, index: i=0, 1, ... 2 * G
        {
            m_nodes[i] = a + ((dx * i) / 2); // even indices are boundaries, odd are cells, distance to boundary is dx/2

            m_area[i] = areaCalc(m_nodes[i]); // even indices are boundaries, odd are cells
        }
    }

    // Getters
    // Index to get specific area values
    double& operator [](size_t index) { return m_area[index]; } // use size_t instead of int to avoid negative
    const double& operator[](size_t index) const { return m_area[index]; } // works for constant data, doesn't modify data

    // method to return spacing (dx)
    double& cellSize() { return dx; }

    // method to extract the length of the grid
    [[nodiscard]] constexpr size_t& size() {return m_size;}

    // method to write x values to file
    void write(std::string& classification)
    {
        // write grid
        std::ofstream xVals("Data/Grid, " + classification + ".txt");

        if (!xVals.is_open())
        { Log("Error: Could not open file."); }
        else
        {
            for (size_t i = 0; i < size(); i++)
            {
                xVals << m_nodes[i] << " ";
            }
        }

        xVals.close();
    }
};

// Class to hold all variables and iterate
template<size_t L>
class Properties
{
private:
    // initialize data structure
    // add two ghost cells on either side of range
    double rho[L + 2];  // density
    double w2[L + 2];   // momentum (rho * u)
    double e[L + 2];    // energy

    double P[L + 2];    // pressure
    double c[L + 2];    // speed of sound

    double Mi;      // inlet mach number

private:
    struct PastIteration // needed so can both calculate residual and iterate within the same loop
    {
        double rho[2];  // density
        double w2[2];   // momentum (rho * u)
        double e[2];    // energy

        double P[2];    // pressure
        double c[2];    // speed of sound
    } prev;

    // is this a recommended way to store a set of vectors?
    struct fluxVec
    {
        double one[L + 1] {};
        double two[L + 1] {};
        double three[L + 1] {};
    } flux;

    const size_t indexL {L + 1}; // last index = size() - 1 = L + 2 - 1 = L + 1
    double dt; // time step

private: // vector used as number of iterations unknown
    std::vector<double> rhoResidual;

public: // methods
    // length of data (unused)
    constexpr size_t size() { return (L + 2); }

    void CalcC(const size_t& index)
    { c[index] = sqrt(gamma * P[index] / rho[index]); }

    // initial condition (Constructor)
    explicit Properties<L>(InputParameters& input)
    {
        // initialize constant static properties
        const double Pe = input.PressureRatio * Pt;
        const double Te = Tt * pow(input.PressureRatio, gamma_ratio_sub);

        // isentropic calculations
        auto isentropicMachP = [&] (double staticP, double M) constexpr
        { return staticP * pow((1 + ((gamma - 1) / 2) * pow(M, 2)), (-1 / gamma_ratio_sub)); };

        auto isentropicMachT = [&] (double staticT, double M) constexpr
        { return staticT * pow((1 + ((gamma - 1) / 2) * pow(M, 2)), (-1)); };

        // initializing calculations for density and energy
        auto rho_calc = [&] (size_t index, double T)  // better to pass already indexed or index
        { rho[index] = P[index] / (R * T); };

        auto e_calc = [&] (size_t index, double T)
        { e[index] =  rho[index] * ((cv * T) + (0.5 * pow((w2[index] / rho[index]), 2))); };

        // initialize dynamic temperature
        const double Ti {isentropicMachT(Tt, input.M0)};
        const double TL {isentropicMachT(Te, input.M0)};

        // set all initial values (except for final element)
        for (size_t i = 0; i < indexL; i++)
        {
                P[i] = isentropicMachP(Pt, input.M0);
                rho_calc(i, Ti);
                CalcC(i);
                w2[i] = (input.M0 * c[i]) * rho[i];
                e_calc(i, Ti);
        }

        // update last element
        P[indexL] = isentropicMachP(Pe, input.M0);
        rho_calc(indexL, TL);
        CalcC(indexL);
        w2[indexL] = (input.M0 * c[indexL]) * rho[indexL]; // rho*u = rho * (M * c)
        e_calc(indexL, TL);

        Mi = input.M0;

        // initialized all property values

        // record past properties at the neighbouring cells to the ghost cells (for use in BC iteration)
        prev = {{rho[1], rho[indexL - 1]}, // last index = size() - 1
                {w2[1], w2[indexL - 1]},
                {e[1], e[indexL - 1]},
                {P[1], P[indexL - 1]},
                {c[1], c[indexL - 1]} };

        // all variables initialized
    }

    // non-isentropic pressure calculation
    void CalcP(size_t& index)
    { P[index] = (gamma - 1) * (e[index] - (0.5 * (pow(w2[index], 2) / rho[index]))); }

    // max eigenvalue calculation
    double eigenCalc(unsigned long index)
    { return (w2[index] / rho[index]) +  c[index]; }

    // numerical flux calculation
    void fluxCalc(InputParameters& input, Grid<L>& grid) // function or subclass?
    {
        auto BoundaryEigenCalc = [&](size_t index)  // average of eigenvalues for surrounding cells
                { return 0.5 * (eigenCalc(index) + eigenCalc((index + 1))); }; // each index + 1

        // initialize variables
        double eigenMax;

        // kinetic energy calc (f2 term)
        auto kinE = [&] (size_t index) { return pow(w2[index], 2) / rho[index]; };

        // f3 term
        auto f3 = [&] (size_t index) { return (e[index] + P[index]) * (w2[index] / rho[index]); };

        // iterate over gridPts + 1 values to update all boundary terms
        for (size_t z = 0; z < indexL; ++z)  // L + 1 boundaries (size - 1) (cells from 0 -> L - 1)
        {
            eigenMax = BoundaryEigenCalc(z);
            flux.one[z] = 0.5 * grid[2 * z] * ((w2[z] + w2[z + 1]) -  // only area boundaries, so only even indices
                                               input.Epsilon * eigenMax * (rho[z + 1] - rho[z]));

            flux.two[z] = 0.5 * grid[2 * z] * ((kinE(z) + P[z]) + (kinE(z + 1) + P[z + 1]) -
                                               input.Epsilon * eigenMax * (w2[z + 1] - w2[z]));

            flux.three[z] = 0.5 * grid[2 * z] * ((f3(z) + f3(z + 1)) -
                                                 input.Epsilon * eigenMax * (e[z + 1] - e[z]));
        }
    }

    // method to update inlet condition
    void UpdateInlet(InputParameters& input, Grid<L>& grid)
    {
        // initialize all necessary variables: better to use lambdas or = f(x)? (See dPdU for example)
        double u0 = w2[0] / rho[0];
        double u1 {prev.w2[0] / prev.rho[0]};
        double dPdU { (Pt / gamma_ratio_sub) * pow(1 - gamma_ratio * (pow(u0, 2) / alpha_star),
                                                   (1 / (gamma - 1))) * (-2 * gamma_ratio * (u0 / alpha_star)) };

        double dt0 = input.CFL * grid.cellSize() / (u0 + c[0]);
        double eigen0 = 0.5 * ((u0 + u1) - (c[0] + prev.c[0])) * dt0 / grid.cellSize();

        double du = -eigen0 * ((prev.P[0] - P[0]) - rho[0] * c[0] * (u1 - u0)) / (dPdU - (rho[0] * c[0]));

        // isentropic calculations (based on u)
        auto isentropicT = [&](double u)
        { return Tt * (1 - gamma_ratio * pow(u, 2) / alpha_star); };

        // increment speed at inlet
        u0 += du;

        // update all properties based on above
        double T0 = isentropicT(u0);
        P[0] = Pt * pow((T0 / Tt), 1 / gamma_ratio_sub);
        rho[0] = P[0] / (R * T0);
        w2[0] = rho[0] * u0;
        e[0] = rho[0] * (cv * T0 + 0.5 * pow(u0, 2));
        c[0] = sqrt(gamma * P[0] / rho[0]);
        Mi = u0 / c[0];

        // save properties of neighbouring cell for next iteration
        prev.rho[0] = rho[1];
        prev.w2[0] = w2[1];
        prev.e[0] = e[1];
        prev.P[0] = P[1];
        prev.c[0] = c[1];
    }

    // method to calculate exit eigenvalues
    auto eigenEnd(InputParameters& input, Grid<L>& grid, double& uL, double& uLminus1) // checked
    {
        // time-step
        double dtL { input.CFL * grid.cellSize() / (uL + c[indexL]) };

        double ave_u = 0.5 * (uL + uLminus1);
        double ave_c = 0.5 * (c[indexL] + prev.c[1]);

        struct Lambda{
            double one;
            double two;
            double three;
        } l {};

        // could also put this in the above {} to initialize instantly but thought this was more readable
        l.one = ave_u * dtL / grid.cellSize();
        l.two = (ave_u + ave_c) * ( dtL / grid.cellSize());
        l.three = (ave_u - ave_c) * ( dtL / grid.cellSize());

        return l;
    }

    // method to update exit properties
    void UpdateExit(InputParameters& input, Grid<L>& grid) // checked
    {
        // variables (changing to u for legibility of code)
        double uL { w2[indexL] / rho[indexL] };
        double uLminus1 { prev.w2[1] / prev.rho[1] };

        // compute eigenvals (see above method)
        auto[l1, l2, l3] = eigenEnd(input, grid, uL, uLminus1);

        // calculate characteristic relations
        double r1 = -l1 * ((rho[indexL] - prev.rho[1]) - (P[indexL] - prev.P[1]) / pow(c[indexL], 2));
        double r2 = -l2 * ((P[indexL] - prev.P[1]) + (rho[indexL] * c[indexL] * (uL - uLminus1)));
        double r3 = -l3 * ((P[indexL] - prev.P[1]) - (rho[indexL] * c[indexL] * (uL - uLminus1)));

        // calculate mach number
        double ML { (uL + uLminus1) / (c[indexL] + prev.c[1]) };

        double dp; // change in pressure (should I put "double dp" in the if else below or is this better?)
        if (ML > 1)     // supersonic outlet
        {
            dp = 0.5 * (r2 + r3);
        }
        else { dp = 0; }    // subsonic outlet

        // calculate changes
        double drho = r1 + (dp / pow(c[indexL], 2));
        double du = (r2 - dp) / (rho[indexL] * c[indexL]);

        // update properties
        rho[indexL] += drho;
        uL += du;
        P[indexL] += dp;

        // update w2
        w2[indexL] = rho[indexL] * uL;

        double TL = P[indexL] / (R * rho[indexL]);
        e[indexL] = rho[indexL] * ((cv * TL) + (0.5 * pow(uL, 2)));
        c[indexL] = sqrt(gamma * P[indexL] / rho[indexL]);

        // save properties for next iteration
        prev.rho[1] = rho[indexL - 1]; // last index = size() - 1 = indexL
        prev.w2[1] = w2[indexL - 1];
        prev.e[1] = e[indexL - 1];
        prev.P[1] = P[indexL - 1];
        prev.c[1] = c[indexL - 1];
    }

    // method to iterate forward in time (Explicit Euler)
    double iterate(InputParameters& input, Grid<L>& grid)
    {
        fluxCalc(input, grid);  // must first iterate over all values to calculate flux (length = L + 1)
        // is there another faster way? Need flux[i] and flux[i + 1] value.

        // residual vector: in this case is a struct better or worse? could add a method that returns r1[end]
        double r1, r2, r3, residual = 0;

        for (size_t i = 1; i < indexL; ++i) // iterate from cell = 1 -> L - 1 (exclude ghost cells)
        {
            // calculate time step based on largest eigenvalues
            dt = (input.CFL * grid.cellSize()) / eigenCalc(i);  // eigenMax = u + c

            // calculate residual terms
            r1 = (flux.one[i] - flux.one[i - 1]) / grid.cellSize();  // i = 0, first boundary. i = L,

            r2 = (flux.two[i] - flux.two[i - 1]) / grid.cellSize() -
                    P[i] * (grid[2*i] - grid[2*(i - 1)]) / grid.cellSize();

            r3 = (flux.three[i] - flux.three[i - 1]) / grid.cellSize();

            residual += pow(r1, 2); // norm 2 residual calculation

            // Euler iteration
            rho[i] -= (dt / grid[(2 * i) - 1]) * r1;  // index = 1, 3, ... 2 * L - 1 (odd indices are cell centres)

            w2[i] -= (dt / grid[(2 * i) - 1]) * r2;

            e[i] -= (dt / grid[(2 * i) - 1]) * r3;

            // update dependent properties
            CalcP(i);
            CalcC(i);
        }
        // store residual term
        rhoResidual.emplace_back(sqrt(residual)); // append rho residual

        if (Mi < 1) // if subsonic, update inlet using u - c characteristic
        {
            UpdateInlet(input,grid);
            Mi = (w2[0] / rho[0]) / c[0]; // M = u / c
        }

        UpdateExit(input, grid); // update outlet

        return sqrt(residual); // return residual value to check for convergence
    }

    // method to save properties to file
    void write(std::string& classification)
    {
        // write properties:
        std::ofstream properties("Data/Properties, " + classification + ".txt");

        if (!properties.is_open())
        {
            Log("Error: Could not open file.");
        }
        else
        {
            for (size_t i = 1; i < indexL; i++) // not including ghost cells currently
            {
                double M = w2[i] / rho[i];

                properties << P[i] << " " << M << std::endl;
            }
        }
        properties.close();

        // write residuals
        std::string residName = "Data/Residuals, " + classification + ".txt";
        std::ofstream residuals(residName);

        if (!residuals.is_open())
        {
            Log("Error: Could not open Residual file.");
        }
        else {
            for (double &resid: rhoResidual) {
                residuals << resid << " ";
            }
        }
        residuals.close();

////      printing final residual value, used to verify convergence before graphing in Python
//        Log(rhoResidual[rhoResidual.size() - 1]);
    }
};

// function to iterate for each set of input parameters (should this be a class? I wasn't sure what that would add)
template<size_t S>
void iteration(InputParameters& input) // function or class? seems identical to me here
{
    // initialize grid parameters and properties
    // Legend: // dx = grid.cellSize() // x* = grid.xNodes(); //S* = grid.Area(); // Si = S[i] = grid[]
    Grid<S> grid;
    Properties<S> data{input};

    // initialize iteration condition
    size_t count = 0;
    double resid = 1;

    // number of iterations before stopping
    const size_t lim = pow(10, 5);

    // iteration loop
    while ((count < lim) && (resid > 5 * std::numeric_limits<double>::epsilon()))
    {
        resid = data.iterate(input, grid);
        count += 1;
    }

    // save to file: (could have as a function, would look cleaner but would have to pass grid, input, S
    // naming
    std::string GridPts = std::to_string(S);
    std::string PressureRatio = std::to_string(input.PressureRatio).erase(4);
    std::string cfl = std::to_string(input.CFL).erase(4);
    std::string epsilon = std::to_string(input.Epsilon).erase(4);

    // naming convention so different iterations aren't overwritten
    std::string classification = "Grid " + GridPts + " Pr " + PressureRatio + " CFL " + cfl + " eps " + epsilon;

    { // writing classification to document
        std::ofstream classifications;
        classifications.open("Data/Iterations.txt", std::ios_base::app);

        if (!classifications.is_open())
        { Log("Error: Could not open file."); }
        else { classifications << classification << " -"; }

        classifications.close();
    }
    // writing x values and data to file
    grid.write(classification);
    data.write(classification);
}

// insertion point for the code: starts here
int main() {

    // Clearing iteration list: (not needed, makes it easier on the Python end)
    std::ofstream classifications("Data/Iterations.txt");
    classifications.close();

    // Defining Pressure Study:
    const int GridPts1 {50};
    const double pRange[4] {0.76, 0.72, 0.68, 0.60};
    for (const double & p : pRange)
    {
        InputParameters InputP{p, 0.2, 0.1, 0.08};
        iteration<GridPts1>(InputP);
    }

    // Defining Grid Study
    InputParameters InputG{0.72, 0.2, 0.1, 0.1};
    iteration<25>(InputG);
    iteration<50>(InputG);
    iteration<100>(InputG);
    iteration<200>(InputG);

    std::cout << "Completed" << std::endl;
    return 0;
}
