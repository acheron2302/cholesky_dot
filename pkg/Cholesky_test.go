package cholesky

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"gonum.org/v1/gonum/mat"
)

func TestCholesky(t *testing.T) {
	matrix1 := mat.NewSymDense(3, []float64{
		3, 4, 3,
		4, 8, 6,
		3, 6, 9,
	})

	matrix2 := mat.NewSymDense(3, []float64{
		25, 10, -5,
		15, 18, 0,
		-5, 0, 11,
	})

	LowerMatrix, err := Chol(matrix1, mat.Upper)
	if err != nil {
		t.Fatal(err)
	}

	LowerMatrix, err = Chol(matrix2, mat.Lower)
	if err != nil {
		t.Fatal(err)
	}
	_ = LowerMatrix
}

func TestSolveCholesky(t *testing.T) {
	A := mat.NewSymDense(3, []float64{
		4, 10, 8,
		10, 26, 26,
		8, 26, 61,
	})

	b := mat.NewVecDense(3, []float64{
		44, 128, 214,
	})

	tri, err := Chol(A, mat.Lower)
	if err != nil {
		t.Error(err)
	}
	result := SolveChol(b, tri)
	expected := mat.NewVecDense(3, []float64{
		-8, 6, 2,
	})
	assert.True(t, mat.Equal(result, expected), "The result is different from expeceted")
}
