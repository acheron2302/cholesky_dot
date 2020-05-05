package main

import (
	"fmt"

	"cholesky/pkg"

	"gonum.org/v1/gonum/mat"
)

func main() {
	var dim int

	fmt.Printf("Insert the dim of the square matrix: ")
	if _, err := fmt.Scanf("%d", &dim); err != nil {
		fmt.Println("The dim must be integer")
		return
	}

	data := make([]float64, dim*dim)
	for i := 0; i < dim; i++ {
		fmt.Printf("Insert row number %d: ", i+1)
		for j := 0; j < dim; j++ {
			fmt.Scanf("%f", &data[i*dim+j])
		}
	}
	sym, err := cholesky.NewSym(dim, data)
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Printf("The matrix is:\n% v\n", mat.Formatted(sym, mat.Squeeze()))

	triMat, err := cholesky.Chol(sym, mat.Lower)
	if err != nil {
		fmt.Println(err)
		return
	}

	fmt.Printf("The lower triangle matrix is:\n% v\n", mat.Formatted(triMat, mat.Squeeze()))

	/*
	 *inp := make([]float64, dim)
	 *fmt.Printf("Input the vector: ")
	 *for i := 0; i < dim; i++ {
	 *    fmt.Scanf("%f", &inp[i])
	 *}
	 *vec := mat.NewVecDense(dim, inp)
	 *result := cholesky.SolveChol(vec, triMat)
	 *fmt.Printf("The result x is: \n%v\n", mat.Formatted(result))
	 */
}
