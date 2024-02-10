//
// Linear algebra namespace including for easy linear algebra manipulations
//

#ifndef DOUBLEMACHREFLECTION_LINEARALGEBRA_H
#define DOUBLEMACHREFLECTION_LINEARALGEBRA_H

#include <array>
#include <cmath>

namespace LinearAlgebra
{
    // Returns scalar * vector
    std::array<double,4> ScalarVectorMultiplication (double scalar, std::array<double,4> vector)
    {
        std::array<double,4> result = vector;
        for (int i{0}; i < 4; i++)
        {
            result[i] *= scalar;
        }
        return result;
    }
    // Returns vectorA + vectorB
    std::array<double,4> VectorVectorAddition (std::array<double,4> vectorA, std::array<double,4> vectorB)
    {
        std::array<double,4> result = vectorA;
        for (int i{0}; i < 4; i++)
        {
            result[i] += vectorB[i];
        }
        return result;
    }
    // Returns vectorA + vectorB + vectorC
    std::array<double,4> VectorVectorVectorAddition (std::array<double,4> vectorA,std::array<double,4> vectorB,std::array<double,4> vectorC)
    {
        std::array<double,4> result = vectorA;
        for (int i{0}; i < 4; i++)
        {
            result[i] += vectorB[i] + vectorC[i];
        }
        return result;
    }
    // Returns vectorA - vectorB
    std::array<double,4> VectorVectorSubtraction (std::array<double,4> vectorA, std::array<double,4> vectorB)
    {
        std::array<double,4> result = vectorA;
        for (int i{0}; i < 4; i++)
        {
            result[i] -= vectorB[i];
        }
        return result;
    }
    // Returns matrix * vector
    std::array<double,4> MatrixVectorMultiplication (std::array<std::array<double,4>,4> &matrix, std::array<double,4> vector)
    {
        std::array<double,4> result = {};
        for (int i{0}; i < 4; i++)
        {   for (int j{0}; j < 4; j++)
            {
                result[i] += matrix[i][j]*vector[j];
            }
        }
        return result;
    }
}

#endif //DOUBLEMACHREFLECTION_LINEARALGEBRA_H
