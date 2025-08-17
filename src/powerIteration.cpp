#include "powerIteration.hpp"
#include "utils.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>
#include <Eigen/Dense>


static Eigen::MatrixXd powerIteration(Eigen::MatrixXd A, int k, int maxIterations, double epsilon){

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

std::vector<Eigen::MatrixXd> approximateImage(std::vector<Eigen::MatrixXd> channels, int k, int maxIterations, double epsilon){

    std::vector<Eigen::MatrixXd> finalMatrices(3);

    for (int i = 0; i < finalMatrices.size(); i++){
        Eigen::MatrixXd newChannel = powerIteration(channels[i], k, maxIterations, epsilon);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}

static Eigen::MatrixXd denoiseWithPatches(Eigen::MatrixXd A, int patchSize, int stride, int k, int maxIterations, double epsilon){

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

std::vector<Eigen::MatrixXd> denoiseImage(std::vector<Eigen::MatrixXd> channels, int patchSize, int stride, int k, int maxIterations, double epsilon){

    std::vector<Eigen::MatrixXd> finalMatrices(3);

    for (int i = 0; i < finalMatrices.size(); i++){
        Eigen::MatrixXd newChannel = denoiseWithPatches(channels[i], patchSize, stride, k, maxIterations, epsilon);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}

static Eigen::MatrixXd medianFilter(Eigen::MatrixXd A, int windowSize){

    int originalHeight = A.rows();
    int originalWidth = A.cols();

    int padSize = windowSize / 2;

    Eigen::MatrixXd paddedMatrix = padMatrix(A, padSize);
    Eigen::MatrixXd A_final = A;

    for (int row = 0; row < A.rows(); row++){
        for (int col = 0; col < A.cols(); col++){
            
            std::cout << "Now at (row, col) = " << row << ", " << col << std::endl;
            std::vector<double> window;
            // initialize window
            for (int i = 0; i < padSize; i++){
                for (int j = 0; j < padSize; j++){
                    window.push_back(A_final(row + i, col + j));
                }
            }

            std::sort(window.begin(), window.end());
            double medianValue = window[window.size() / 2];

            A_final(row,col) = medianValue;
        }
    }

    return A_final;
}

std::vector<Eigen::MatrixXd> applyMedianFilter(std::vector<Eigen::MatrixXd> channels, int windowSize){

    std::vector<Eigen::MatrixXd> finalMatrices(3);

    for (int i = 0; i < finalMatrices.size(); i++){
        std::cout << "Now at index = " << i << std::endl;
        Eigen::MatrixXd newChannel = medianFilter(channels[i], windowSize);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}



