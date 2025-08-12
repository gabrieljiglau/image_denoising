# image_denoising

### image compression

Power iteration is an iterative algorithm for finding the biggest eigenvalue from a matrix. If the rank of the matrix is deflated at each iteration, then it can be used to approximate all the eigenvalues.

Aditionally, the left and right singular vectors can be approximated too in the process, making matrix reconstruction possible.

#### If only we use the k biggest eigenvalues and singular values,compression of the input is possible.

original
<img src="images/original.png" alt="original" width=400/>

k=25
<img src="images/compressed_k25.png" alt="original" width=400/>

k=50
<img src="images/compressed_k50.png" alt="original" width=400/>

k=75
<img src="images/compressed_k75.png" alt="original" width=400/>

k=100
<img src="images/compressed_k100.png" alt="original" width=400/>

