//
// Created by jaimy on 02/03/2023.
//

#ifndef DOUBLEMACHREFLECTION_SOLVER_H
#define DOUBLEMACHREFLECTION_SOLVER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <limits>
#include <algorithm>

#include "CellCentre.h"
#include "Interface.h"

class Solver2D
{
public:
    void setDomainData (double lengthX, int numGridCellsX, double lengthY, int numGridCellsY,
                        double endTime, double maxCFL)
    {
        Lx = lengthX; Nx = numGridCellsX; dx = Lx/Nx;
        Ly = lengthY; Ny = numGridCellsY; dy = Ly/Ny;
        T  = endTime; CFL = maxCFL;
        N = Nx*Ny;
    }
    void doAnalysis (const std::string& filePath, int numberOutputFiles)
    {
        // Set grid data and initialize
        setGenericGridData();
        initializeCellCentredData();

        // Export initial data
        exportData(filePath,0);

        // Get time step size and set initial time
        getTimeStepSize(); t = 0.0;
        std::cout << "Initial time step size is set to: " << dt << std::endl;
        if ( T/dt < numberOutputFiles)
        {
            throw std::runtime_error("Currently more output files are requested than the expected number of time steps, "
                                     "please use a smaller maximum CFL-number or decrease the number of requested output files.");
        }

        std::cout << "Solver started" << std::endl;
        while (t < T)
        {
            // Inform user
            std::cout << "Solving for t = " << t << std::endl;

            // Do calculation
            getTimeStepSolution();

            // Update time step size and current time
            getTimeStepSize();
            t += dt;

            // Export data if required
            double I_dt = dt/T * numberOutputFiles;         // Index step size
            double I_d = t/T * numberOutputFiles;           // Index in double format
            int I_i = static_cast<int>(std::round(I_d));    // Index in int format
            double d_id = I_d - I_i;                        // Difference between double and int format
            if (std::abs(d_id) < std::abs(d_id - I_dt) && std::abs(d_id) < std::abs(d_id + I_dt) && I_i > 0)
            {   /* Only export data if:
                 * The difference is smaller than the exact previous difference
                 * The difference is smaller than the approximate next difference
                 * The integer index is larger than the initial condition index
                 */
                std::cout << "\nExporting time step: " << I_i << ", t = " << t << std::endl << std::endl;
                exportData(filePath, I_i);
            }
        }

        if (t >= T)
        {
            std::cout << "Solver finished successfully!" << std::endl;
        }
        else
        {
            std::cout << "Shit hit the fan." << std::endl;
        }
    }
private:
    ////////// Data Storage //////////
    //// Solution data
    std::vector<cellCentreData> gridCell;
    std::vector<interfaceData> interfaceX;
    std::vector<interfaceData> interfaceY;
    //// Domain constants
    int N {};                                           // Total number of grid cells
    double CFL {};                                      // Maximum CFL
    double Lx {};                                       // Domain size in x-direction
    int Nx {};                                          // Number of grid cells in x-direction
    double dx {};                                       // Grid size in x-direction
    double Ly {};                                       // Domain size in y-direction
    int Ny {};                                          // Number of grid cells in y-direction
    double dy {};                                       // Grid size in y-direction
    double t {};                                        // Current time
    double T {};                                        // Analysis end time
    double dt {};                                       // Time step size
    //// Constants
    const double PI {2*std::acos(0.0)};
    //// Boundary constants
    double shockAngle {60.0 * PI/180};                  // Angle of shock w.r.t. positive x-axis
    double wallLocation{1.0/6.0};                       // Starting location of the lower boundary wall
    //// Physics constants
    double S{20.0/std::sqrt(3.0)};                   // Shock velocity
    double gamma {1.4};                                 // Specific heat ratio
    std::array<double,4> preConditions                  // Conditions in front (right) of moving shock
            {gamma,0.0,0.0,1.0};                        // {rho, u, v, p}
    std::array<double,4> postConditions                 // Conditions in after (left of) the moving shock
            {8.0, 4.125*std::sqrt(3), -4.125, 116.5};// {rho, u, v, p}

    ////////// Functions //////////
    //// Domain data and pre-calculation functions
    // Store generic grid data
    void setGenericGridData()
    {
        // Initialize data storage
        gridCell = std::vector<cellCentreData> (N);
        interfaceX = std::vector<interfaceData> (N+Ny);
        interfaceY = std::vector<interfaceData> (N+Nx);

        for(int j{0}; j < Ny; j++)
        {
            for (int i{0}; i < Nx; i++)
            {
                gridCell[indexCc(i,j)].Position = {(i+0.5)*dx,(j+0.5)*dy,0.0};
            }
        }
    }

    // Index functions
    int indexCc (int i, int j) const    // Cell centre index (gridCell)
    {return j + i*Ny;}
    int indexIe (int i, int j) const    // Eastern interface index (interfaceX)
    {return j + (i+1)*Ny;}
    int indexIw (int i, int j) const    // Western interface index (interfaceX)
    {return j + i*Ny;}
    int indexIn (int i, int j) const    // Northern interface index (interfaceY)
    {return j + i*(Ny+1) + 1;}
    int indexIs (int i, int j) const    // Southern interface index (interfaceY)
    {return j + i*(Ny+1);}


    // Initialization of cell centred data {rho, rho*u, rho*v, rho*E}
    void initializeCellCentredData ()
    {
        for(int j{0}; j < Ny; j++)
        {   for (int i{0}; i < Nx; i++)
            {
                cellCentreData &CCD = gridCell[indexCc(i,j)];
                double x{(i+0.5)*dx};double y{(j+0.5)*dy};
                if ( x < wallLocation )
                {   // Conditions left of the shock
                    CCD.U = {postConditions[0],
                             postConditions[0]*postConditions[1],
                             postConditions[0]*postConditions[2],
                             getDensityEnergy(postConditions[3],postConditions[0],postConditions[1],postConditions[2])};
                }
                else if ( x > (wallLocation + Ly/std::tan(shockAngle)) )
                {   // Conditions right of the shock
                    CCD.U = {preConditions[0],
                             preConditions[0]*preConditions[1],
                             preConditions[0]*preConditions[2],
                             getDensityEnergy(preConditions[3],preConditions[0],preConditions[1],preConditions[2])};
                }
                else
                {   // Conditions within shock vicinity
                    if ( y >= (x - wallLocation)*std::tan(shockAngle))
                    {   // Conditions above the shock
                        CCD.U = {postConditions[0],
                                 postConditions[0]*postConditions[1],
                                 postConditions[0]*postConditions[2],
                                 getDensityEnergy(postConditions[3],postConditions[0],postConditions[1],postConditions[2])};
                    }
                    else
                    {   // Conditions below the shock
                        CCD.U = {preConditions[0],
                                 preConditions[0]*preConditions[1],
                                 preConditions[0]*preConditions[2],
                                 getDensityEnergy(preConditions[3],preConditions[0],preConditions[1],preConditions[2])};
                    }
                }
            }
        }
        std::cout << "Data initialized" << std::endl;
    }

    //// Variable conversion functions
    // Returns rho*E
    double getDensityEnergy (double &p, double &rho, double &u, double &v) const
    {
        return p/(gamma-1) + 0.5*rho*(u*u+v*v);
    }
    // Returns p
    double getPressure (std::array<double,4> &U) const
    {
        return (gamma-1)*(U[3] - 0.5*(U[1]*U[1] + U[2]*U[2])/U[0]);
    }
    // Returns c
    double getSpeedOfSound (std::array<double,4> &U) const
    {
        return std::sqrt( gamma*(gamma - 1) * (U[3]/U[0] - (U[1]*U[1] + U[2]*U[2])/(2*U[0]*U[0])) );
    }
    // Returns s
    double getEntropy (double &p, double &rho) const
    {
        return std::log( p / std::pow(rho, gamma) );
    }
    double getEntropy (std::array<double,4> &U)
    {
        return std::log(getPressure(U) / std::pow(U[0],gamma));
    }
    // Returns H
    double getEnthalpy (std::array<double,4> &U)
    {
        return (U[3] + getPressure(U))/U[0];
    }

    //// Interface flux
    // Determines all interface fluxes according to the referenced conservedVariables (U,K1,K2)
    void getInterfaceFlux(std::array<double,4> cellCentreData::*conservedVariables)
    {
        //// Transform the cell centred conserved variables into the primitive variables
        for (int n{0}; n < N; n++)
        {
            gridCell[n].getCellCentredPrimitiveVariables(gamma,conservedVariables);
        }
        //// Interface interpolation of primitive variables
        /// Boundary faces
        // X-direction
        for (int j{0}; j < Ny; j++)
        {   // West boundary
            // Pointers
            interfaceData &IFX = interfaceX[indexIw(0,j)];

            // Set primitive variables from boundary conditions
            IFX.ZL = {postConditions[3], postConditions[1], postConditions[2], getEntropy(postConditions[3], postConditions[0])};
            getMUSCL_WEST_BC(IFX, gridCell[indexCc(0,j)], gridCell[indexCc(1,j)]);
        }
        for (int j{0}; j < Ny; j++)
        {   // East boundary
            // Pointers
            interfaceData &IFX = interfaceX[indexIe(Nx-1,j)];

            // Set primitive variables from boundary conditions
            IFX.ZR = {preConditions[3], preConditions[1], preConditions[2], getEntropy(preConditions[3], preConditions[0])};
            getMUSCL_EAST_BC(IFX, gridCell[indexCc(Nx-2,j)], gridCell[indexCc(Nx-1,j)]);
        }
        // Y-direction
        for (int i{0}; i < Nx; i++)
        {   // south boundary
            // Pointers
            interfaceData &IFY = interfaceY[indexIs(i,0)];

            // Set primitive variables from boundary conditions
            if ( (i+0.5)*dx < wallLocation )
            {   // Pre wall
                IFY.ZL = {postConditions[3], postConditions[1], postConditions[2], getEntropy(postConditions[3], postConditions[0])};
            }
            else
            {   // Reflective slip wall (invert v)
                IFY.ZL = {gridCell[indexCc(i,0)].Z[0], gridCell[indexCc(i,0)].Z[1], -gridCell[indexCc(i,0)].Z[2], gridCell[indexCc(i,0)].Z[3]};
            }

            getMUSCL_WEST_BC(IFY, gridCell[indexCc(i,0)], gridCell[indexCc(i,0)]);

            if ( (i+0.5)*dx >= wallLocation )
            {   // Reflective slip wall (invert v)
                IFY.ZL = {IFY.ZR[0], IFY.ZR[1], -IFY.ZR[2], IFY.ZR[3]};
            }
        }
        for (int i{0}; i < Nx; i++)
        {   // North boundary
            // Pointers
            interfaceData &IFY = interfaceY[indexIn(i,Ny-1)];

            // Set primitive variables from boundary conditions
            if ( (i+0.5)*dx <= (wallLocation + Ly/std::tan(shockAngle)) + S * t )
            {   // Left of shock where shock position is determined using x_1 + S * T
                IFY.ZR = {postConditions[3], postConditions[1], postConditions[2], getEntropy(postConditions[3], postConditions[0])};
            }
            else
            {   // Right of shock
                IFY.ZR = {preConditions[3], preConditions[1], preConditions[2], getEntropy(preConditions[3], preConditions[0])};
            }
            getMUSCL_EAST_BC(IFY, gridCell[indexCc(i,Ny-2)], gridCell[indexCc(i,Ny-1)]);
        }
        /// Interior interfaces
        // X-direction
        for (int j{0}; j < Ny; j++)
        {   // West near boundary
            getMUSCL_WEST_NB(interfaceX[indexIe(0,j)],interfaceX[indexIw(0,j)],gridCell[indexCc(0,j)]
                    ,gridCell[indexCc(1,j)],gridCell[indexCc(2,j)]);
            // East near boundary
            getMUSCL_EAST_NB(interfaceX[indexIe(Nx-2,j)],gridCell[indexCc(Nx-3,j)],
                                 gridCell[indexCc(Nx-2,j)],gridCell[indexCc(Nx-1,j)],interfaceX[indexIe(Nx-1,j)]);
        }
        for (int i{1}; i < (Nx-2); i++)
        {   for (int j{0}; j < Ny; j++)
            {   // Interior
                getMUSCL(interfaceX[indexIe(i,j)],gridCell[indexCc(i-1,j)],
                         gridCell[indexCc(i,j)],gridCell[indexCc(i+1,j)],gridCell[indexCc(i+2,j)]);
            }
        }
        // Y-direction
        for (int i{0}; i < Nx; i++)
        {   // South near boundary
            getMUSCL_WEST_NB(interfaceY[indexIn(i,0)],interfaceY[indexIs(i,0)],gridCell[indexCc(i,0)]
                    ,gridCell[indexCc(i,1)],gridCell[indexCc(i,2)]);
            // North near boundary
            getMUSCL_EAST_NB(interfaceY[indexIn(i,Ny-2)],gridCell[indexCc(i,Ny-3)],
                                 gridCell[indexCc(i,Ny-2)],gridCell[indexCc(i,Ny-1)],interfaceY[indexIn(i,Ny-1)]);
        }
        for (int i{0}; i < Nx; i++)
        {   for (int j{1}; j < (Ny-2); j++)
            {   // Interpolate y-direction interface values
                getMUSCL(interfaceY[indexIn(i,j)],gridCell[indexCc(i,j-1)],
                         gridCell[indexCc(i,j)],gridCell[indexCc(i,j+1)],gridCell[indexCc(i,j+2)]);
            }
        }
        //// Get Roe's fluxes
        // X-direction
        for (int nx {0}; nx < (N + Ny); nx++)
        {
            // Pointers
            interfaceData &IFX = interfaceX[nx];

            // Transform interface primitive variables into conserved variables
            IFX.getInterfaceConservedVariables(gamma);

            // Get averaged state variables
            IFX.getAveragedStateVariables();

            // Get Roe's averaged Jacobian matrices
            IFX.getAveragedStateMatrixX(gamma);

            // Get Roe's averaged fluxes
            IFX.getRoeFluxX(gamma);
        }
        // Y-direction
        for (int ny {0}; ny < (N + Nx); ny++)
        {
            // Pointers
            interfaceData &IFY = interfaceY[ny];

            // Transform interface primitive variables into conserved variables
            IFY.getInterfaceConservedVariables(gamma);

            // Get averaged state variables
            IFY.getAveragedStateVariables();

            // Get Roe's averaged Jacobian matrices
            IFY.getAveragedStateMatrixY(gamma);

            // Get Roe's averaged fluxes
            IFY.getRoeFluxY(gamma);
        }
    }
    //// Interface interpolation Schemes
    // MC flux limiter to ensure TVD flux
    static double getLimiter (double r)
    {
        return std::max( 0.0, std::min((1.0+r)/2.0, std::min(2.0*r,2.0)) );
    }
    // MUSCL scheme for complete eastern/northern interface
    static void getMUSCL (interfaceData &e, const cellCentreData &W,const cellCentreData &P, const cellCentreData &E, const cellCentreData &EE)
    {
        for (int i{0}; i < 4; i++)
        {
            e.ZL[i] = P.Z[i] + 0.5*getLimiter((P.Z[i]-W.Z[i])/(E.Z[i]-P.Z[i]))*(E.Z[i]-P.Z[i]);
            e.ZR[i] = E.Z[i] - 0.5*getLimiter((E.Z[i]-P.Z[i])/(EE.Z[i]-E.Z[i]))*(EE.Z[i]-E.Z[i]);
        }
    }
    static void getMUSCL_WEST_NB (interfaceData &e, const interfaceData &w, const cellCentreData &P, const cellCentreData &E, const cellCentreData &EE)
    {
        for (int i{0}; i < 4; i++)
        {
            e.ZL[i] = P.Z[i] + 0.5*getLimiter((P.Z[i]-w.ZL[i])/(E.Z[i]-P.Z[i]))*(E.Z[i]-P.Z[i]);
            e.ZR[i] = E.Z[i] - 0.5*getLimiter((E.Z[i]-P.Z[i])/(EE.Z[i]-E.Z[i]))*(EE.Z[i]-E.Z[i]);
        }
    }
    static void getMUSCL_EAST_NB (interfaceData &e, const cellCentreData &W,const cellCentreData &P, const cellCentreData &E, const interfaceData &ee)
    {
        for (int i{0}; i < 4; i++)
        {
            e.ZL[i] = P.Z[i] + 0.5*getLimiter((P.Z[i]-W.Z[i])/(E.Z[i]-P.Z[i]))*(E.Z[i]-P.Z[i]);
            e.ZR[i] = E.Z[i] - 0.5*getLimiter((E.Z[i]-P.Z[i])/(ee.ZR[i]-E.Z[i]))*(ee.ZR[i]-E.Z[i]);
        }
    }
    static void getMUSCL_WEST_BC (interfaceData &w, const cellCentreData &P, const cellCentreData &E)
    {
        for (int i{0}; i < 4; i++)
        {
            w.ZR[i] = P.Z[i] - 0.5*getLimiter((P.Z[i]-w.ZL[i])/(E.Z[i]-P.Z[i]))*(E.Z[i]-P.Z[i]);
        }
    }
    static void getMUSCL_EAST_BC (interfaceData &e, const cellCentreData &W,const cellCentreData &P)
    {
        for (int i{0}; i < 4; i++)
        {
            e.ZL[i] = P.Z[i] + 0.5*getLimiter((P.Z[i]-W.Z[i])/(e.ZR[i]-P.Z[i]))*(e.ZR[i]-P.Z[i]);
        }
    }
    //// Return the numerical operator for the semi-discrete system
    std::array<double,4> getOperator (int &i, int &j)
    {
        std::array<double,4> numericalOperator =
                LinearAlgebra::VectorVectorAddition(
                LinearAlgebra::ScalarVectorMultiplication(1/dx, LinearAlgebra::VectorVectorSubtraction(interfaceX[indexIw(i,j)].F, interfaceX[indexIe(i,j)].F)),
                LinearAlgebra::ScalarVectorMultiplication(1/dy, LinearAlgebra::VectorVectorSubtraction(interfaceY[indexIs(i,j)].F, interfaceY[indexIn(i,j)].F))
        );
        return numericalOperator;
    }
    //// Determine step size using the max CFL
    void getTimeStepSize ()
    {
        std::vector<double> timeStepSize (N);
        for (int n{0}; n < N; n++)
        {
            // Current cell pointer
            cellCentreData &CCD = gridCell[n];

            // Get current cell max time step size
            timeStepSize[n] = CFL / ( (std::abs(CCD.U[1])/dx + std::abs(CCD.U[2])/dy)/CCD.U[0] + getSpeedOfSound(CCD.U)*(1/dx + 1/dy) );
            CCD.dt = timeStepSize[n];
        }
        // Get minimum step size
        dt = *std::min_element(timeStepSize.begin(), timeStepSize.end());
    }
    //// Explicit time integration step using SSP-RK3
    void getTimeStepSolution ()
    {
        // Get first kutta coefficient
        getInterfaceFlux(&cellCentreData::U);
        for (int i{0}; i < Nx; i++)
        {   for (int j{0}; j < Ny; j++)
            {
                gridCell[indexCc(i,j)].K1 = LinearAlgebra::VectorVectorAddition(gridCell[indexCc(i,j)].U,
                                                                                LinearAlgebra::ScalarVectorMultiplication(dt,getOperator(i,j)));
                gridCell[indexCc(i,j)].entropyCorrection = interfaceX[indexIw(i,j)].entropyCorrection + interfaceX[indexIe(i,j)].entropyCorrection
                                                         + interfaceY[indexIs(i,j)].entropyCorrection + interfaceY[indexIn(i,j)].entropyCorrection;
            }
        }
        // Get second kutta coefficient
        getInterfaceFlux(&cellCentreData::K1);
        for (int i{0}; i < Nx; i++)
        {   for (int j{0}; j < Ny; j++)
            {
                gridCell[indexCc(i,j)].K2 = LinearAlgebra::VectorVectorVectorAddition(LinearAlgebra::ScalarVectorMultiplication(3.0/4.0, gridCell[indexCc(i,j)].U),
                                                                                      LinearAlgebra::ScalarVectorMultiplication(1.0/4.0, gridCell[indexCc(i,j)].K1),
                                                                                      LinearAlgebra::ScalarVectorMultiplication(dt/4.0,getOperator(i,j)));
            }
        }
        // Get ze final solution
        getInterfaceFlux(&cellCentreData::K2);
        for (int i{0}; i < Nx; i++)
        {   for (int j{0}; j < Ny; j++)
            {
                gridCell[indexCc(i,j)].U = LinearAlgebra::VectorVectorVectorAddition(LinearAlgebra::ScalarVectorMultiplication(1.0/3.0, gridCell[indexCc(i,j)].U),
                                                                                      LinearAlgebra::ScalarVectorMultiplication(2.0/3.0, gridCell[indexCc(i,j)].K2),
                                                                                      LinearAlgebra::ScalarVectorMultiplication(dt*2.0/3.0,getOperator(i,j)));
            }
        }
    }

    //// Export data
    void exportData (const std::string& filePath, int timeStepIndex)
    {
        std::ofstream DataFile;
        DataFile.open(filePath + "DoubleMachReflection_Nx" + std::to_string(Nx) +
                      "_Ny" + std::to_string(Ny) +
                      "_CFL" + std::to_string(CFL) +
                      "_nt" + std::to_string(timeStepIndex) + ".vts");
        DataFile << "<VTKFile type='StructuredGrid' version='1.0' byte_order='LittleEndian'>\n";
        DataFile << "<StructuredGrid WholeExtent='0 " << Nx-1 << " 0 " << Ny-1 << " 0 0'>\n";
        DataFile << "<Piece Extent='0 " << Nx-1 << " 0 " << Ny-1 << " 0 0'>\n";

        DataFile << "<Points>\n";
        DataFile << "<DataArray type='Float64' Name='Points' NumberOfComponents='3'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            cellCentreData &CCD = gridCell[n];
            DataFile << CCD.Position[0] << " " << CCD.Position[1] << " " << CCD.Position[2] << "\n";
        }
        DataFile << "</DataArray>\n";
        DataFile << "</Points>\n";
        DataFile << "<PointData>\n";

        // Export velocity
        DataFile << "<DataArray type='Float64' Name='velocity' NumberOfComponents='2' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            cellCentreData &CCD = gridCell[n];
            DataFile << CCD.U[1]/CCD.U[0] << " " << CCD.U[2]/CCD.U[0] << "\n";
        }
        DataFile << "</DataArray>\n";

        // Export density
        DataFile << "<DataArray type='Float64' Name='density' NumberOfComponents='1' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            DataFile << gridCell[n].U[0] << "\n";
        }
        DataFile << "</DataArray>\n";

        // Export pressure
        DataFile << "<DataArray type='Float64' Name='pressure' NumberOfComponents='1' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            DataFile << getPressure(gridCell[n].U) << "\n";
        }
        DataFile << "</DataArray>\n";

        // Export speed of sound
        DataFile << "<DataArray type='Float64' Name='speed of sound' NumberOfComponents='1' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            DataFile << getSpeedOfSound(gridCell[n].U) << "\n";
        }
        DataFile << "</DataArray>\n";

        // Export total energy
        DataFile << "<DataArray type='Float64' Name='total energy' NumberOfComponents='1' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            cellCentreData &CCD = gridCell[n];
            DataFile << CCD.U[3]/CCD.U[0] << "\n";
        }
        DataFile << "</DataArray>\n";

        // Export enthalpy
        DataFile << "<DataArray type='Float64' Name='enthalpy' NumberOfComponents='1' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            DataFile << getEnthalpy(gridCell[n].U) << "\n";
        }
        DataFile << "</DataArray>\n";

        // Export entropy
        DataFile << "<DataArray type='Float64' Name='entropy' NumberOfComponents='1' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            DataFile << getEntropy(gridCell[n].U) << "\n";
        }
        DataFile << "</DataArray>\n";

        // Export time step size
        DataFile << "<DataArray type='Float64' Name='dt' NumberOfComponents='1' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            DataFile << gridCell[n].dt << "\n";
        }
        DataFile << "</DataArray>\n";

        // Export entropy correction
        DataFile << "<DataArray type='Float64' Name='entropyCorrection' NumberOfComponents='1' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            DataFile << gridCell[n].entropyCorrection << "\n";
        }
        DataFile << "</DataArray>\n";

        // Export eigenvalue_x
        DataFile << "<DataArray type='Float64' Name='eigenvalue_x' NumberOfComponents='2' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            cellCentreData &CCD = gridCell[n];
            double c {getSpeedOfSound(gridCell[n].U)};
            DataFile << CCD.U[1]/CCD.U[0] + c << " " << CCD.U[1]/CCD.U[0] - c << "\n";
        }
        DataFile << "</DataArray>\n";

        // Export eigenvalue_y
        DataFile << "<DataArray type='Float64' Name='eigenvalue_y' NumberOfComponents='2' format='ascii'>\n";
        for (int n{0}; n < N; n++)
        {   // Export grid positions
            cellCentreData &CCD = gridCell[n];
            double c {getSpeedOfSound(gridCell[n].U)};
            DataFile << CCD.U[2]/CCD.U[0] + c << " " << CCD.U[1]/CCD.U[0] - c << "\n";
        }
        DataFile << "</DataArray>\n";

        DataFile << "</PointData>\n";
        DataFile << "</Piece>\n";
        DataFile << "</StructuredGrid>\n";
        DataFile << "</VTKFile>\n";
        DataFile.close();
    }
};

#endif //DOUBLEMACHREFLECTION_SOLVER_H
