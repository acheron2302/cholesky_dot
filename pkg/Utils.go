package cholesky

import (
	"gonum.org/v1/gonum/mat"
)

func createUMatrix(t, n int) *mat.Dense {
	I := createIdentityMatrix(n)

	for i := 0; i < t+1; i++ {
		if i == t {
			arr := make([]float64, n)
			arr[0] = 1
			I.SetRow(i, arr)
			continue
		}
		arr := make([]float64, n)
		arr[i+1] = 1
		I.SetRow(i, arr)
	}

	return I
}

func createUInvMatrix(t, n int) *mat.Dense {
	I := createIdentityMatrix(n)

	for i := 0; i < t+1; i++ {
		if i == t {
			arr := make([]float64, n)
			arr[0] = 1
			I.SetRow(i, arr)
			continue
		}
		arr := make([]float64, n)
		arr[i+1] = 1
		I.SetRow(i, arr)
	}

	temp := I.T()
	I.CloneFrom(temp)
	return I
}

func createSMatrix(companion mat.Matrix, r, c int) (S *mat.Dense) {
	n, _ := companion.Dims()
	S = createIdentityMatrix(n)

	arr := make([]float64, n)
	arr[r-1] = 1
	for i := c; i < n; i++ {
		arr[i] = companion.At(r, i)
	}
	S.SetRow(r-1, arr)

	return
}

func createSinvMatrix(companion mat.Matrix, r, c int) (S_inv *mat.Dense) {
	n, _ := companion.Dims()
	S_inv = createIdentityMatrix(n)

	arr := make([]float64, n)
	arr[r-1] = 1
	for i := c; i < n; i++ {
		arr[i] = -companion.At(r, i)
	}
	S_inv.SetRow(r-1, arr)

	return
}

func createIdentityMatrix(n int) *mat.Dense {
	arr := make([]float64, n*n)

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				arr[i*n+j] = 1
				continue
			}
			arr[i*n+j] = 0
		}
	}

	matrix := mat.NewDense(n, n, arr)
	return matrix
}

func createExchangeIndentityMatrix(r1, r2, dims int) *mat.Dense {
	mat := createIdentityMatrix(dims)

	mat.Set(r1, r1, 0)
	mat.Set(r2, r2, 0)
	mat.Set(r2, r1, 1)
	mat.Set(r1, r2, 1)

	return mat
}
