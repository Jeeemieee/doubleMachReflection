//
// Created by jaimy on 12 Mar 2023.
//

#ifndef DOUBLEMACHREFLECTION_INTERFACE_H
#define DOUBLEMACHREFLECTION_INTERFACE_H

#include <array>
#include <cmath>

#include "LinearAlgebra.h"

class interfaceData
{
public:
    // Interface interpolated conserved variables (U,K1,K2)
    std::array<double,4> VL;
    std::array<double,4> VR;
    std::array<double,4> V;
    // Interface interpolated primitive variables
    std::array<double,4> ZL;
    std::array<double,4> ZR;
    // Roe's matrix
    std::array<std::array<double,4>,4> A;
    // Eigen values
    std::array<double,3> L;
    // Interface flux
    std::array<double,4> F;
    // Entropy correction
    bool entropyCorrection;

    // Primitive to conserved variable transformation
    void getInterfaceConservedVariables (const double &gamma)
    {
        getConservedVariables(VL,ZL,gamma);
        getConservedVariables(VR,ZR,gamma);
    }

    // Roe's averaged state conserved variables
    void getAveragedStateVariables()
    {
        for (int i{0}; i < 4; i++)
        {
            V[i] = 0.5*(VL[i] + VR[i]);
        }
    }

    // Roe's averaged state matrix for x-direction
    void getAveragedStateMatrixX(const double &gamma)
    {
        // Set input variables
        double u {V[1]/V[0]};
        double uu = u*u;
        double v {V[2]/V[0]};
        double vv = v*v;
        double u_sum = uu + vv;
        double cc { gamma*(gamma - 1) * (V[3]/V[0] - u_sum/2) };
        double c = std::sqrt(cc);
        double H = cc / (gamma-1) + (u_sum)/2;

        // Update eigenvalues according to entropy fix
        //L = {std::abs(u), std::abs(u-c), std::abs(u+c)};
        getFixedEigenValuesX(u, c, gamma);

        // Get reduced variables
        double H_red = -2*H + u_sum;
        double L_diff = (L[1]) - (L[2]);
        double L_sum = (L[1]) + (L[2]);
        double L_red = -2*(L[0]) + L_sum;

        // Update linearized Jacobian matrix
        A = {{
                     {L_diff*u/(2*c) + (L[0]) - L_red*u_sum/(2*H_red), -L_diff/(2*c) + L_red*u/H_red, L_red*v/H_red, -L_red/H_red},
                     {L_diff*uu/(2*c) - L_red*u/2 + L_diff*c*u_sum/(2*H_red) - L_red*u_sum*u/(2*H_red), -L_diff*u/(2*c) + L_sum/2 - L_diff*c*u/H_red + L_red*uu/H_red, v*(-L_diff*c + L_red*u)/H_red, (L_diff*c - L_red*u)/H_red},
                     {v*(H_red*L_diff*u-L_red*c*u_sum)/(2*H_red*c), -L_diff*v/(2*c) + L_red*u*v/H_red, (L[0]) + L_red*vv/H_red, -L_red*v/H_red},
                     {H*L_diff*u/(2*c) - L_red*(uu-vv)/4 + L_diff*c*u_sum*u/(2*H_red) - L_red*u_sum*u_sum/(4*H_red), -H*L_diff/(2*c) - L_diff*c*uu/H_red + L_red*u_sum*u/(2*H_red), v*(H*L_red - L_diff*c*u)/H_red, (-H*L_sum + L_diff*c*u + u_sum*(L[0]))/(H_red)}
             }};
    }

    // Roe's averaged state matrix for y-direction
    void getAveragedStateMatrixY(const double &gamma)
    {
        // Set input variables
        double u {V[1]/V[0]};
        double uu = u*u;
        double v {V[2]/V[0]};
        double vv = v*v;
        double u_sum = uu + vv;
        double cc { gamma*(gamma - 1) * (V[3]/V[0] - u_sum/2) };
        double c = std::sqrt(cc);
        double H = cc / (gamma-1) + (u_sum)/2;

        // Update eigenvalues according to entropy fix
        //L = {std::abs(v), std::abs(v-c), std::abs(v+c)};
        getFixedEigenValuesY(v, c, gamma);

        // Get reduced variables
        double H_red = -2*H + u_sum;
        double L_diff = (L[1]) - (L[2]);
        double L_sum = (L[1]) + (L[2]);
        double L_red = -2*(L[0]) + L_sum;

        // Update linearized Jacobian matrix
        A = {{
                     {L_diff*v/(2*c) + (L[0]) - L_red*u_sum/(2*H_red), L_red*u/H_red, -L_diff/(2*c) + L_red*v/H_red, -L_red/H_red},
                     {u*(H_red*L_diff*v-L_red*c*u_sum)/(2*H_red*c), (L[0]) + L_red*uu/H_red, -L_diff*u/(2*c) + L_red*u*v/H_red, -L_red*u/H_red},
                     {L_diff*vv/(2*c) - L_red*v/2 + L_diff*c*u_sum/(2*H_red) - L_red*u_sum*v/(2*H_red), u*(-L_diff*c + L_red*v)/H_red, -L_diff*v/(2*c) + L_sum/2 - L_diff*c*v/H_red + L_red*vv/H_red, (L_diff*c - L_red*v)/H_red},
                     {H*L_diff*v/(2*c) + L_red*(uu-vv)/4 + L_diff*c*u_sum*v/(2*H_red) - L_red*u_sum*u_sum/(4*H_red), u*(H*L_red - L_diff*c*v)/H_red, -H*L_diff/(2*c) - L_diff*c*vv/H_red + L_red*u_sum*v/(2*H_red), (-H*L_sum + L_diff*c*v + u_sum*(L[0]))/(H_red)}
            }};
    }

    // Roe flux for x-direction
    void getRoeFluxX (const double &gamma)
    {
        F = LinearAlgebra::VectorVectorSubtraction(
                LinearAlgebra::ScalarVectorMultiplication(0.5, LinearAlgebra::VectorVectorAddition(getInterfaceFluxX(gamma,VR),getInterfaceFluxX(gamma,VL))),
                LinearAlgebra::ScalarVectorMultiplication(0.5, LinearAlgebra::MatrixVectorMultiplication(A,LinearAlgebra::VectorVectorSubtraction(VR,VL)))
        );
    }

    // Roe flux for y-direction
    void getRoeFluxY (const double &gamma)
    {
        F = LinearAlgebra::VectorVectorSubtraction(
                LinearAlgebra::ScalarVectorMultiplication(0.5, LinearAlgebra::VectorVectorAddition(getInterfaceFluxY(gamma,VR),getInterfaceFluxY(gamma,VL))),
                LinearAlgebra::ScalarVectorMultiplication(0.5, LinearAlgebra::MatrixVectorMultiplication(A,LinearAlgebra::VectorVectorSubtraction(VR,VL)))
        );
    }

private:
    static void getConservedVariables(std::array<double,4> &V, const std::array<double,4> &Z, const double &gamma)
    {
        // Get density from the entropy relation
        V[0] = std::exp((std::log(Z[0]) - Z[3])/ gamma );
        // Get u-mom
        V[1] = V[0]*Z[1];
        // Get V-mom
        V[2] = V[0]*Z[2];
        // Get energy from the pressure relation
        V[3] = Z[0]/(gamma-1) + 0.5*V[0]*(Z[1]*Z[1] + Z[2]*Z[2]);
    }
    // Eigen values using constant entropy fix
    // Eigen values using constant entropy fix
    // x-direction
    void getFixedEigenValuesX (double &u, double &c, const double &gamma)
    {
        // Set initial eigenvalues
        L = {(u), (u-c), (u+c)};
        // Get left and right state eigenvalues
        double cL = getSpeedOfSound(VL, gamma);
        std::array<double,3> LL {ZL[1], ZL[1] - cL, ZL[1] + cL};
        double cR = getSpeedOfSound(VR, gamma);
        std::array<double,3> LR {ZR[1], ZR[1] - cR, ZL[1] + cR};

        for (int k{0}; k < 3; k++)
        {
            double delta = std::max( 0.0, std::max( L[k] - LL[k],  LR[k] - L[k]) );//0.2*(std::abs(u) + c);//
            if ( std::abs(L[k]) < delta )
            {
                L[k] = (L[k]*L[k]/delta + delta)/2.0;
                entropyCorrection = true;
            }
            else
            {
                L[k] = std::abs(L[k]);
                entropyCorrection = false;
            }
        }
    }
    // y-direction
    void getFixedEigenValuesY (double &v, double &c, const double &gamma)
    {
        // Set initial eigenvalues
        L = {(v), (v-c), (v+c)};
        // Get left and right state eigenvalues
        double cL = getSpeedOfSound(VL, gamma);
        std::array<double,3> LL {ZL[2], ZL[2] - cL, ZL[2] + cL};
        double cR = getSpeedOfSound(VR, gamma);
        std::array<double,3> LR {ZR[2], ZR[2] - cR, ZL[2] + cR};

        for (int k{0}; k < 3; k++)
        {
            double delta = std::max( 0.0, std::max( L[k] - LL[k],  LR[k] - L[k]) );//0.2*(std::abs(v) + c);//
            if ( std::abs(L[k]) < delta )
            {
                L[k] = (L[k]*L[k]/delta + delta)/2.0;
                entropyCorrection = true;
            }
            else
            {
                L[k] = std::abs(L[k]);
                entropyCorrection = false;
            }
        }
    }
    static double getSpeedOfSound (std::array<double,4> &U, const double &gamma)
    {
        return std::sqrt( gamma*(gamma - 1) * (U[3]/U[0] - (U[1]*U[1] + U[2]*U[2])/(2*U[0]*U[0])) );
    }

    // Interface fluxes
    // x-direction
    static std::array<double,4> getInterfaceFluxX(const double &gamma, std::array<double,4> &U)
    {
        double p {(gamma-1)*(U[3] - 0.5*(U[1]*U[1] + U[2]*U[2])/U[0])};
        std::array<double,4> F {U[1], U[1]*U[1]/U[0]+p, U[1]*U[2]/U[0], (U[3]+p)*U[1]/U[0]};
        return F;
    }
    // y-direction
    static std::array<double,4> getInterfaceFluxY(const double &gamma, std::array<double,4> &U)
    {
        double p {(gamma-1)*(U[3] - 0.5*(U[1]*U[1] + U[2]*U[2])/U[0])};
        std::array<double,4> G {U[2], U[1]*U[2]/U[0], U[2]*U[2]/U[0]+p, (U[3]+p)*U[2]/U[0]};
        return G;
    }
};


#endif //DOUBLEMACHREFLECTION_INTERFACE_H
