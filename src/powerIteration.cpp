#include "include/power_iteration.hpp"
#include "include/utils.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

/**
 * @brief Applies the power iteration algorithm with deflation on a given channel
 *
 * @param A Input channel from matrix (R/G/B)
 * @param k Keep only the k-biggest singular values
 * @param maxIterations Max iterations for the algorithm
 * @param epsilon Floating point precision for computations
 * @return The truncated channel as an Eigen::MatrixXd.
 *
 */
Eigen::MatrixXd SVD::powerIteration(Eigen::MatrixXd A, int k, int maxIterations, double epsilon){

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


/// Applies the powerIteration on all 3 channels from the image
std::vector<Eigen::MatrixXd> SVD::apply(){

    std::vector<Eigen::MatrixXd> finalMatrices(3);

    for (int i = 0; i < finalMatrices.size(); i++){
        Eigen::MatrixXd newChannel = powerIteration(channels[i], k, maxIterations, epsilon);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}


/**
 * @brief Applies the SVD decomposition on smaller submatrices from a given channel
 *
 * @param A Input channel from matrix (R/G/B)
 * @param patchSize The dimensions of the smaller submatrices (patchSize X patchSize)
 * @param stride The step size then traversing the channel
 * @param k Keep only the k-biggest singular values
 * @param maxIterations Max iterations for the algorithm
 * @param epsilon Floating point precision for computations
 * @return The truncated channel as an Eigen::MatrixXd.
 *
 */
Eigen::MatrixXd SVD::applyPatches(Eigen::MatrixXd A, int patchSize, int stride, int k, int maxIterations, double epsilon){

    Eigen::MatrixXd A_final = Eigen::MatrixXd::Zero(A.rows(), A.cols());
    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(A.rows(), A.cols());

    for (int row = 0; row < A.rows(); row += stride){
        int rowStart = 0;
        if (row + patchSize > A.rows()){
            rowStart = A.rows() - patchSize;
        } else {
            rowStart = row;
        }

        for (int col = 0; col < A.cols(); col += stride){
            int colStart = 0;
            if (col + patchSize > A.cols()){
                colStart = A.cols() - patchSize;
            } else {
                colStart = col;
            }

            Eigen::MatrixXd patch = A.block(rowStart, colStart, patchSize, patchSize);
            Eigen::MatrixXd newPatch = powerIteration(patch, k, maxIterations, epsilon);

            // there might be multiple patches for a given pixel
            A_final.block(rowStart, colStart, patchSize, patchSize) += newPatch;
            cumulateWeights(weights, rowStart, colStart, patchSize);

        }
    }

    return averagePixels(A_final, weights);
}



/// Applies the patch SVD on all 3 channels from the image
std::vector<Eigen::MatrixXd> SVD::apply(int patchSize, int stride){

    std::vector<Eigen::MatrixXd> finalMatrices(3);

    for (int i = 0; i < finalMatrices.size(); i++){
        Eigen::MatrixXd newChannel = applyPatches(channels[i], patchSize, stride, k, maxIterations, epsilon);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}



