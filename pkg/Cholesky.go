package cholesky

import (
	"fmt"
	"math/cmplx"

	"gonum.org/v1/gonum/mat"
)

func cholUpper(triMatrix *CTriDense, sym mat.Symmetric) error {
	for i := 0; i < sym.Symmetric(); i++ {
		for j := 0; j < i; j++ {
			Uji := 1 / triMatrix.At(j, j)
			s := 0 + 0i

			for k := 0; k < j; k++ {
				s += triMatrix.At(k, i) * cmplx.Conj(triMatrix.At(k, j))
			}

			Uji *= complex(sym.At(i, j), 0) - s
			triMatrix.SetTri(j, i, Uji)
		}

		s := 0 + 0i
		for k := 0; k < i; k++ {
			s += triMatrix.At(k, i) * cmplx.Conj(triMatrix.At(k, i))
		}
		Uii := cmplx.Sqrt(complex(sym.At(i, i), 0) - s)
		if Uii == 0 {
			return fmt.Errorf("The matrix is not positive definite")
		}

		triMatrix.SetTri(i, i, Uii)
	}

	return nil
}

func cholLower(triMatrix *CTriDense, sym mat.Symmetric) error {
	for i := 0; i < sym.Symmetric(); i++ {
		for j := 0; j < i; j++ {
			Lij := 1 / triMatrix.At(j, j)
			s := 0 + 0i

			for k := 0; k < j; k++ {
				s += triMatrix.At(i, k) * cmplx.Conj(triMatrix.At(j, k))
			}

			Lij *= complex(sym.At(i, j), 0) - s
			triMatrix.SetTri(i, j, Lij)
		}

		s := 0 + 0i
		for k := 0; k < i; k++ {
			s += triMatrix.At(i, k) * cmplx.Conj(triMatrix.At(i, k))
		}
		Lii := cmplx.Sqrt(complex(sym.At(i, i), 0) - s)
		if Lii == 0 {
			return fmt.Errorf("The matrix is not positive definite")
		}

		triMatrix.SetTri(i, i, Lii)
	}

	return nil
}

func Chol(sym mat.Symmetric, kind mat.TriKind) (*CTriDense, error) {
	// early stop
	if sym.At(0, 0) <= 0 {
		return nil, fmt.Errorf("a_00 must be real and positive")
	}

	if sym.Symmetric() == 0 {
		return &CTriDense{}, nil
	}

	triMatrix := newCTriDense(sym.Symmetric(), kind)

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

	fmt.Printf("The triangle matrix is:\n% v\n", Formatted(triMatrix))

	return triMatrix, nil
}

func SolveChol(b mat.Vector, tri *CTriDense) *CVectorDense {
	resultVec := newCVector(b.Len(), nil)

	// forward
	lower := tri
	resultVec.SetVec(0, complex(b.AtVec(0), 0)/lower.At(0, 0))
	for i := 1; i < b.Len(); i++ {
		s := 0 + 0i
		for j := 0; j < i; j++ {
			s += lower.At(i, j) * resultVec.AtVec(j)
		}
		resultVec.SetVec(i, (complex(b.AtVec(i), 0)-s)/lower.At(i, i))
	}

	// backward
	n := b.Len()
	upper := tri.H()
	resultVec.SetVec(n-1, resultVec.AtVec(n-1)/upper.At(n-1, n-1))
	for i := n - 2; i >= 0; i-- {
		s := 0 + 0i
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

func ChangeToSymMatrix(matrix *mat.Dense, vec *mat.VecDense) (*mat.SymDense, *mat.VecDense) {
	if !mat.Equal(matrix, matrix.T()) {
		vec.MulVec(matrix.T(), vec)
		matrix.Mul(matrix.T(), matrix)
	}

	arr := []float64{}
	for i := 0; i < vec.Len(); i++ {
		for j := 0; j < vec.Len(); j++ {
			arr = append(arr, matrix.At(i, j))
		}
	}
	resultMat := mat.NewSymDense(vec.Len(), arr)
	fmt.Printf("The matrix is: \n% v\n", mat.Formatted(matrix))
	fmt.Printf("The vector is: \n% v\n", mat.Formatted(vec))

	return resultMat, vec
}

func SolveGeneral(vec *mat.VecDense, matrix *mat.Dense) (*mat.VecDense, error) {
	sym, vec := ChangeToSymMatrix(matrix, vec)
	tri, err := Chol(sym, mat.Lower)
	if err != nil {
		return nil, err
	}
	complexVector := SolveChol(vec, tri)
	result := mat.NewVecDense(complexVector.Len(), nil)
	for i := 0; i < complexVector.Len(); i++ {
		result.SetVec(i, real(complexVector.AtVec(i)))
	}

	return result, nil
}
