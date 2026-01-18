#include "include/power_iteration.hpp"
#include "include/utils.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
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
std::vector<Eigen::MatrixXd> SVD::powerIteration(Eigen::MatrixXd A, std::vector<int> k, int maxIterations, double epsilon){

    std::sort(k.begin(), k.end()); // ensure that k's are in order

    int biggestK = k[k.size() - 1];

    std::vector<Eigen::MatrixXd> U = zeroMatrix(A.rows(), k);
    std::vector<Eigen::VectorXd> S = zeroVector(k);
    std::vector<Eigen::MatrixXd> V = zeroMatrix(A.cols(), k);

    Eigen::MatrixXd A_copy = A;
    std::vector<Eigen::MatrixXd> A_final = zeroMatrix(A.rows(), A.cols(), k.size());

    for (int idx = 0; idx < biggestK; idx++){

        //std::cout << "Now at singular value number " << idx + 1;

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
        //std::cout << "; singular value = " << singularValue << std::endl;
        
        for (int updateIdx = 0; updateIdx < k.size(); updateIdx++){
            if (k[updateIdx] > idx) {
                U[updateIdx].col(idx) = u; // or U.block(0, idx, U.rows(), 1) = u; 
                S[updateIdx](idx) = singularValue;
                V[updateIdx].col(idx) = v;
            }
        }

        // deflation step: remove rank-1 component
        A_copy -= singularValue * u * v.transpose();

        // reconstruction for each k
        for (int updateIdx = 0; updateIdx < k.size(); updateIdx++){
            if (k[updateIdx] > idx){
                A_final[updateIdx] += singularValue * u * v.transpose();
            }
        }

        //std::cout << "u = :" << std::endl << u << std::endl;
        //std::cout << "v = :" << std::endl << v << std::endl;
    }

    //std::cout << "A_final " << std::endl << A_final << std::endl;
    return A_final; // vec k1, vec k2, vec k3
}


/// Applies the powerIteration on all 3 channels from the image
std::vector<std::vector<Eigen::MatrixXd>> SVD::apply(){

    std::vector<std::vector<Eigen::MatrixXd>> finalMatrices(3, std::vector<Eigen::MatrixXd>(kVals.size()));

    /// TODO: modificat aici
    // tine cont de faptul ca powerIteration returneaza: // vec k1, vec k2, vec k3

    for (int i = 0; i < finalMatrices.size(); i++){
        std::cout << "Now at channel " << i + 1 << std::endl;

        /// aici logica e gresita, deoarece daca ai un singur k / sau doi,  kreconstruction.size = 1, respectiv 2
        // si nu apuci sa muti matricea in 'finalMatrices' 
        std::vector<Eigen::MatrixXd> kReconstructions = powerIteration(channels[i], kVals, maxIterations, epsilon);
        std::cout << " kReconstructions.size() = " << kReconstructions.size() << std::endl;
        for (int j = 0; j < kReconstructions.size(); j++){
            finalMatrices[i][j] = kReconstructions[j];
        }
    }

    std::cout << "apply in svd ok " << std::endl; 
    return finalMatrices;
}


void SVD::cumulateWeights(Eigen::MatrixXd &weights, int rowStart, int colStart, int patchSize){

    for (int i = 0; i < patchSize; i++){
        for (int j = 0; j < patchSize; j++){
            weights(rowStart + i, colStart + j) += 1.0;
        }
    }
}

std::vector<Eigen::MatrixXd> SVD::averagePixels(std::vector<Eigen::MatrixXd> A, std::vector<Eigen::MatrixXd> weights){

    assert (A.size() == weights.size());

    for (int idx = 0; idx < A.size(); idx++){

        Eigen::MatrixXd innerA = A[idx];
        Eigen::MatrixXd innerW = weights[idx];
        
        assert (innerA.rows() == innerW.rows() && innerA.cols() == innerW.cols());

        for (int i = 0; i < innerA.rows(); i++){
            for (int j = 0; j < innerA.cols(); j++){
                if (innerW(i, j) > 0){
                    innerA(i, j) /= innerW(i, j);
                }
            }
        }
    }

    return A;
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
 * @return The truncated channel as a vector of Eigen::MatrixXd.
 *
 */
std::vector<Eigen::MatrixXd> SVD::applyPatches(Eigen::MatrixXd A, int patchSize, int stride, std::vector<int> k, int maxIterations, double epsilon){

    std::vector<Eigen::MatrixXd> A_final = zeroMatrix(A.rows(), A.cols(), k.size());
    std::vector<Eigen::MatrixXd> weights = zeroMatrix(A.rows(), A.cols(), k.size());

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
            std::vector<Eigen::MatrixXd> newPatches = powerIteration(patch, k, maxIterations, epsilon);

            // there might be multiple patches for a given pixel

            for (int updateIdx = 0; updateIdx < k.size(); updateIdx++){
                A_final[updateIdx].block(rowStart, colStart, patchSize, patchSize) += newPatches[updateIdx];
                cumulateWeights(weights[updateIdx], rowStart, colStart, patchSize);
            }

        }
    }

    return averagePixels(A_final, weights);
}


/// Applies the patch SVD on all 3 channels from the image
std::vector<std::vector<Eigen::MatrixXd>> SVD::apply(int patchSize, int stride){

    std::vector<std::vector<Eigen::MatrixXd>> finalMatrices(3);

    for (int i = 0; i < kVals.size(); i++){
        std::vector<Eigen::MatrixXd> currentReconstruction(3);        
        finalMatrices[i] = applyPatches(channels[i], patchSize, stride, kVals, maxIterations, epsilon);
    }

    return  finalMatrices;
}



