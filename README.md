### system-solver

#### It can:
- Solve equation AX=B using Gauss methods.
- Find inverse matrix using Gauss.
- Make DLUP decomposition, where D - a rows permutations matrix, P - columns permutations matrix(https://en.wikipedia.org/wiki/LU_decomposition#LDU_decomposition). Famous LU decomposition is a particular case of this method, when D = P = Identical Matrix (https://en.wikipedia.org/wiki/LU_decomposition).
- Solve AX=B and find inversed matrix using DLUP decomposition
- make LDL^T decomposition for square symmetric matrices(https://en.wikipedia.org/wiki/Cholesky_decomposition).
- make QR decomposition of matrix(https://en.wikipedia.org/wiki/QR_decomposition).

Some methods have some optimizations, for example, Gauss method with choosing next row/column.

The repository was originally designed for passing lab work at university, so it has a lot of unneeded files.
```c++
// @TODO refactor everything
```
