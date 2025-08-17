# image_denoising

### image compression

In the singular value decomposition (SVD), a matrix A can be decomposed using 3 matrices, `A = U * Σ * Vᵀ`, where U and V are orthogonal, meaning `U * Uᵀ = I`, and Σ is a diagonal matrix, whose non-zero values are the eigenvalues of A.

Power iteration is an iterative algorithm for finding the biggest eigenvalue from a matrix, starting from a random, but normalized vector `v`. If the rank of the matrix is deflated at each iteration, then it can be used to approximate all the eigenvalues.

Aditionally, the left and right singular vectors can also be approximated in the process, making matrix reconstruction possible.

`A * v1 = U * Σ * Vᵀ * v1`. v1 is the first column of V, and since V is orthogonal,  `Vᵀ * v1 = e1`, where e1 is the first column from the identity matrix. Σ * e1 will be just the first singular value, σ1​.

The equation is therefore reduced to `A * v1 = U * σ1​ * e1` (with σ1​ being just a scaler), so `A * v1 = u1 * σ1​`. Take the norm of both sides and the result is `||A * v1|| = σ1` (the norm of a singular vector is always normalized). Reintroduce in the equation above and `u1 = (A * v1) / ||A * v1||`.

Similar steps are taken to deduce the update formula for v, with the subtle change that U needs to be on the right side. We can accomplish this by taking the transposde of A `Aᵀ = V * Σᵀ * Uᵀ`. It then follows that `v1 = (Aᵀ * u1) / σ1`.


#### If we use the k biggest eigenvalues and singular vectors from the input, we esentially compress the image.

k=5
<img src="images/compressed_k5.png" alt="k5" width=400/> 

k=25
<img src="images/compressed_k25.png" alt="k25" width=400/>

k=50
<img src="images/compressed_k50.png" alt="k50" width=400/>

k=75
<img src="images/compressed_k75.png" alt="k75" width=400/>


