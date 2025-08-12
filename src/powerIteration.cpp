#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <Eigen/Dense>


Eigen::VectorXd randomVector(int size, double epsilon){
    
    Eigen::VectorXd vector = Eigen::VectorXd(size);

    std::normal_distribution<double> normalDistribution(0, 1);
    std::mt19937 generator(123456);

    for (int i = 0; i < size; i++){
        vector(i) = normalDistribution(generator);
    }

    double norm = vector.norm();
    if (norm < epsilon){
        vector = Eigen::VectorXd::Ones(size);
        return vector;
    }

    return vector / norm;
}


Eigen::MatrixXd powerIterationChannel(Eigen::MatrixXd A, int k, int maxIterations, double epsilon){

    /*
    A: input channel matrix (R/G/B)
    k: k-biggest singular values
    epsilon: precision for computations
    */

    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(A.rows(), k);
    Eigen::VectorXd S = Eigen::VectorXd::Zero(k);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(A.cols(), k);

    Eigen::MatrixXd A_copy = A;
    Eigen::MatrixXd A_final = Eigen::MatrixXd::Zero(A.rows(), A.cols());

    for (int idx = 0; idx < k; idx++){

        std::cout << "Now at singular value number " << idx + 1 << std::endl;

        //v -> the right singular vector
        Eigen::VectorXd v = randomVector(A.cols(), epsilon);

        //u -> the left singular vector
        Eigen::VectorXd u;
        for (int iteration = 0; iteration < maxIterations; iteration++){

            u = A_copy * v;
            double uNorm = u.norm();

            if (uNorm < epsilon){
                break;
            }
            u /= uNorm;

            Eigen::VectorXd v_new = A_copy.transpose() * u;
            double vNorm = v_new.norm();

            if (vNorm < epsilon){
                break;
            }
            v_new /= vNorm;
            
            double normDiff = (v_new - v).norm();
            if (normDiff < epsilon){
                v = v_new;
                break;
            }

            v = v_new;
        }

        u = A_copy * v;
        double singularValue = u.norm();
        u /= singularValue;
        std::cout << "singular value = " << singularValue << std::endl;
        

        U.col(idx) = u; // or U.block(0, idx, U.rows(), 1) = u; 
        S(idx) = singularValue;
        V.col(idx) = v;


        // deflation step: remove rank-1 component
        A_copy -= singularValue * u * v.transpose();
        // reconstruction
        A_final += singularValue * u * v.transpose();

        //std::cout << "u = :" << std::endl << u << std::endl;
        //std::cout << "v = :" << std::endl << v << std::endl;
    }

    //std::cout << "A_final " << std::endl << A_final << std::endl;
    return A_final;
}

std::vector<Eigen::MatrixXd> compressImage(std::vector<Eigen::MatrixXd> channels, int k, int maxIterations, double epsilon){

    std::vector<Eigen::MatrixXd> finalMatrices(3);

    for (int i = 0; i < finalMatrices.size(); i++){
        Eigen::MatrixXd newChannel = powerIterationChannel(channels[i], k, maxIterations, epsilon);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}


