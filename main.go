package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	cholesky "cholesky/pkg"
	"gonum.org/v1/gonum/mat"
)

func main() {
	matrix, err := inputMatrix("./TestData/matrixTest.txt")
	if err != nil {
		fmt.Printf("Error when try to read matrix: %v\n", err)
	}

	vec, err := inputVec("./TestData/vectorTest.txt")
	if err != nil {
		fmt.Printf("Error when try to read vec: %v\n", err)
	}

	result, err := cholesky.SolveGeneral(vec, matrix)
	if err != nil {
		fmt.Printf("Error when try to solve cholesky: %v\n", err)
	}

	fmt.Printf("The result vector x: \n% v\n", mat.Formatted(result))
}

func inputMatrix(file string) (*mat.Dense, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var arr []float64
	var n int
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		arrText := strings.Split(scanner.Text(), " ")
		n = len(arrText)
		for _, text := range arrText {
			value, err := strconv.ParseFloat(text, 64)
			if err != nil {
				return nil, err
			}
			arr = append(arr, value)
		}
	}

	return mat.NewDense(n, n, arr), nil
}

func inputVec(file string) (*mat.VecDense, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var arr []float64
	var n int
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		value, err := strconv.ParseFloat(scanner.Text(), 64)
		if err != nil {
			return nil, err
		}
		arr = append(arr, value)
		n++
	}

	return mat.NewVecDense(n, arr), nil
}
