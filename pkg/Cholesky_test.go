package cholesky

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
	"gonum.org/v1/gonum/mat"
)

func TestCholesky(t *testing.T) {
	/*
	 *matrix1 := mat.NewSymDense(3, []float64{
	 *    3, 4, 3,
	 *    4, 8, 6,
	 *    3, 6, 9,
	 *})
	 */

	matrix2 := mat.NewSymDense(3, []float64{
		25, 451, -251,
		451, 18, 98,
		-251, 98, 11,
	})

	LowerMatrix, err := Chol(matrix2, mat.Lower)
	if err != nil {
		t.Fatal(err)
	}

	fmt.Printf("a = \n% v\n", Formatted(LowerMatrix))
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
	fmt.Printf("The result x is: \n% v\n", Formatted(result))
}

func TestChangeToSymMatrix(t *testing.T) {
	A := mat.NewDense(3, 3, []float64{
		4, 10, 19,
		71, 26, 26,
		28, 36, 61,
	})

	b := mat.NewVecDense(3, []float64{
		44, 128, 214,
	})

	result1, result2 := ChangeToSymMatrix(A, b)

	expected1 := mat.NewSymDense(3, []float64{
		5841, 2894, 3630,
		2894, 2072, 3062,
		3630, 3062, 4758,
	})

	expected2 := mat.NewVecDense(3, []float64{
		15256,
		11472,
		17218,
	})

	assert.Equal(t, expected1, result1)
	assert.Equal(t, expected2, result2)
}

func TestSolveGeneral(t *testing.T) {
	A := mat.NewDense(3, 3, []float64{
		4, 10, 19,
		71, 26, 26,
		28, 36, 61,
	})

	b := mat.NewVecDense(3, []float64{
		44, 128, 214,
	})

	result, err := SolveGeneral(b, A)
	if err != nil {
		t.Fatal(err)
	}

	expected := newCVector(3, []complex128{
		-9.219662058369668,
		54.55760368662457,
		-24.45775729646123,
	})

	assert.Equal(t, expected, result)
}
