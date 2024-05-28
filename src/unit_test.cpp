///
/// \file unit_test.cpp
/// \brief Unit test
/// \author Lionel GENEVE
/// \date 15/03/2022
/// \version 1.0
///

// ---------- Headers ----------
// STD
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <utility>
#include <cmath>

// EIGEN
#include <Eigen/Dense>

// Project
#include "FilterUtils.hpp"
#include "Measurement.hpp"
#include "State3d.hpp"

using namespace ekf;

int main(int argc, char** argv) {
    VectorX z1 = VectorX::Zero(3);
    MatrixX R1 = MatrixX::Zero(3,3);
    int id1 = 123;
    std::string name1 = "dummy_sensor";
    double timestamp1 = 0.12345;
    double mahaThresh1 = 1.0;

    Measurement m1(z1, R1, id1, name1, timestamp1, mahaThresh1);

    Measurement m2(m1);

    VectorX z3 = VectorX::Zero(6);
    MatrixX R3 = MatrixX::Zero(6, 6);
    double ts3 = 2.3456;
    Measurement m3(z3, R3, 456, "accelero", ts3, 2.0);

    std::string m1_str = m1.toString();
    std::cout << m1_str << std::endl;

    std::cout << "M1: " << m1 << std::endl;
    std::cout << "M2: " << m2 << std::endl;

    const int nStates = State3d::N_STATES;
    const int nErrorStates = State3d::N_ERROR_STATES;

    VectorX x0 = VectorX::Zero(nStates);
    MatrixX P0 = MatrixX::Zero(nErrorStates, nErrorStates);
    double ts0 = 0.111;
    State3d s0(x0, P0, ts0);
    std::string s0_str = s0.toString();
    std::cout << s0_str << std::endl;

    State3d s0bis(s0);

    State3d s0ter = s0;

    std::cout << "s0: " << s0 << std::endl;
    std::cout << "s0bis: " << s0bis << std::endl;
    std::cout << "s0ter: " << s0ter << std::endl;
    
    s0.addMeasurement(m1);
    s0.addMeasurement(m3);

    s0.checkStateVector();
    s0.checkCovarianceMatrix();

    std::cout << "s0: " << s0 << std::endl;
    std::cout << "S0 nMeas=" << s0.nMeasurements() << std::endl;

    VectorX x1 = VectorX::Zero(nStates);
    MatrixX P1 = MatrixX::Zero(nErrorStates, nErrorStates);
    double ts1 = 0.222;
    State3d s1(x1, P1, ts1);

    s1 > s0 ? std::cout << "s1 > s0\n" : std::cout << "s1 <= s0\n";
    
    return 0;
}