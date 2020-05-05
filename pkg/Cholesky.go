package cholesky

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
	"math"
)

func cholUpper(triMatrix *mat.TriDense, sym mat.Symmetric) error {
	for i := 0; i < sym.Symmetric(); i++ {
		for j := 0; j < i; j++ {
			Uji := 1 / triMatrix.At(j, j)
			s := 0.0

			for k := 0; k < j; k++ {
				s += triMatrix.At(k, i) * triMatrix.At(k, j)
			}

			Uji *= sym.At(i, j) - s
			triMatrix.SetTri(j, i, Uji)
		}

		s := 0.0
		for k := 0; k < i; k++ {
			s += math.Pow(triMatrix.At(k, i), 2)
		}
		Uii := math.Sqrt(sym.At(i, i) - s)
		if Uii == 0 || math.IsNaN(Uii) {
			return fmt.Errorf("The matrix is not positive definite")
		}

		triMatrix.SetTri(i, i, Uii)
	}

	return nil
}

func cholLower(triMatrix *mat.TriDense, sym mat.Symmetric) error {
	for i := 0; i < sym.Symmetric(); i++ {
		for j := 0; j < i; j++ {
			Lij := 1 / triMatrix.At(j, j)
			s := 0.0

			for k := 0; k < j; k++ {
				s += triMatrix.At(i, k) * triMatrix.At(j, k)
			}

			Lij *= sym.At(i, j) - s
			triMatrix.SetTri(i, j, Lij)
		}

		s := 0.0

		for k := 0; k < i; k++ {
			s += math.Pow(triMatrix.At(i, k), 2)
		}

		Lii := math.Sqrt(sym.At(i, i) - s)

		if Lii == 0 || math.IsNaN(Lii) {
			return fmt.Errorf("The matrix is not positive definite")
		}

		triMatrix.SetTri(i, i, Lii)
	}

	return nil
}

func Chol(sym mat.Symmetric, kind mat.TriKind) (*mat.TriDense, error) {
	// early stop
	if sym.At(0, 0) <= 0 {
		return nil, fmt.Errorf("a_00 must be real and positive")
	}

	if sym.Symmetric() == 0 {
		return &mat.TriDense{}, nil
	}

	triMatrix := mat.NewTriDense(sym.Symmetric(), kind, nil)

	// if it is upper matrix
	if kind == mat.Upper {
		if err := cholUpper(triMatrix, sym); err != nil {
			return nil, err
		}

		return triMatrix, nil
	}

	// if it is lower matrix
	if err := cholLower(triMatrix, sym); err != nil {
		return nil, err
	}

	return triMatrix, nil
}

func SolveChol(b mat.Vector, tri mat.Triangular) *mat.VecDense {
	resultVec := mat.NewVecDense(b.Len(), nil)

	// forward
	lower := tri
	resultVec.SetVec(0, b.AtVec(0)/lower.At(0, 0))
	for i := 1; i < b.Len(); i++ {
		s := 0.0
		for j := 0; j < i; j++ {
			s += lower.At(i, j) * resultVec.AtVec(j)
		}
		resultVec.SetVec(i, (b.AtVec(i)-s)/lower.At(i, i))
	}

	// backward
	n := b.Len()
	upper := tri.TTri()
	resultVec.SetVec(n-1, resultVec.AtVec(n-1)/upper.At(n-1, n-1))
	for i := n - 2; i >= 0; i-- {
		s := 0.0
		for j := i + 1; j < n; j++ {
			s += upper.At(i, j) * resultVec.AtVec(j)
		}
		resultVec.SetVec(i, (resultVec.AtVec(i)-s)/upper.At(i, i))
	}

	return resultVec
}

func NewSym(n int, data []float64) (*mat.SymDense, error) {
	matrix := mat.NewDense(n, n, data)

	if mat.Equal(matrix, matrix.T()) {
		return mat.NewSymDense(n, data), nil
	}
	return nil, fmt.Errorf("The matrix is not Symmetric")
}
