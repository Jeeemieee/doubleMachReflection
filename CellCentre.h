//
// Data storage for the cell centres
//

#ifndef DOUBLEMACHREFLECTION_CELLCENTRE_H
#define DOUBLEMACHREFLECTION_CELLCENTRE_H

#include <array>
#include <cmath>

struct cellCentreData
{
    // Conserved variable
    std::array<double,4> U{};
    // Primitive variable
    std::array<double,4> Z{};
    // Time integration results
    std::array<double,4> K1{};
    std::array<double,4> K2{};
    // Generic data
    std::array<double,3> Position {};
    double dt {};
    int entropyCorrection {};

    // Conserved to primitive variable transformation
    void getCellCentredPrimitiveVariables (const double &gamma, std::array<double,4> cellCentreData::*conservedVariables)
    {
        // Get pointer to appropriate conserved variable array (U, K1, or K2)
        const std::array<double,4>& U_ptr = this->*conservedVariables;

        // Set pressure from energy, u-mom, v-mom, and mass
        Z[0] = (gamma-1)*(U_ptr[3] - 0.5*(U_ptr[1]*U_ptr[1] + U_ptr[2]*U_ptr[2])/U_ptr[0]);
        // Set u velocity from u-mom and mass
        Z[1] = U_ptr[1]/U_ptr[0];
        // Set v velocity from u-mom and mass
        Z[2] = U_ptr[2]/U_ptr[0];
        // Set entropy from pressure and mass
        Z[3] = std::log( Z[0] / std::pow(U_ptr[0],gamma));
    }
};

#endif //DOUBLEMACHREFLECTION_CELLCENTRE_H
