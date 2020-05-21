package cholesky

import (
	"gonum.org/v1/gonum/mat"
)

type CTriDense struct {
	Matrix *mat.CDense
	Kind   mat.TriKind
	Len    int
}

func newCTriDense(n int, kind mat.TriKind) *CTriDense {
	return &CTriDense{
		Matrix: mat.NewCDense(n, n, nil),
		Kind:   kind,
		Len:    n,
	}
}

func (tri *CTriDense) SetTri(i, j int, value complex128) {
	if tri.Kind == mat.Lower {
		tri.Matrix.Set(i, j, value)
	}
}

func (tri *CTriDense) Dims() (r, c int) {
	return tri.Len, tri.Len
}

func (tri *CTriDense) At(i, j int) complex128 {
	return tri.Matrix.At(i, j)
}

func (tri *CTriDense) H() mat.CMatrix {
	if tri.Kind == mat.Upper {
		for i := 0; i < tri.Len; i++ {
			for j := i; j < tri.Len; j++ {
				if i == j {
					continue
				}

				tri.Matrix.Set(j, i, tri.Matrix.At(i, j))
				tri.Matrix.Set(i, j, 0)
			}
		}

		tri.Kind = mat.Lower
		return tri
	}

	for i := 0; i < tri.Len; i++ {
		for j := 0; j < i; j++ {
			if i == j {
				continue
			}

			tri.Matrix.Set(j, i, tri.Matrix.At(i, j))
			tri.Matrix.Set(i, j, 0)
		}
	}
	tri.Kind = mat.Upper
	return tri
}

type CVector interface {
	mat.CMatrix
	AtVec(int) complex128
	Len() int
}

type CVectorDense struct {
	Matrix *mat.CDense
}

func newCVector(len int, values []complex128) *CVectorDense {
	return &CVectorDense{
		Matrix: mat.NewCDense(len, 1, values),
	}
}

func (vec *CVectorDense) SetVec(i int, value complex128) {
	vec.Matrix.Set(i, 0, value)
}

func (vec *CVectorDense) AtVec(i int) complex128 {
	return vec.Matrix.At(i, 0)
}

func (vec *CVectorDense) At(i, j int) complex128 {
	return vec.Matrix.At(i, j)
}

func (vec *CVectorDense) Len() int {
	l, _ := vec.Matrix.Dims()
	return l
}

func (vec *CVectorDense) Dims() (r, c int) {
	r, c = vec.Matrix.Dims()
	return
}

func (vec *CVectorDense) H() mat.CMatrix {
	result := mat.NewCDense(1, vec.Len(), nil)
	for i := 0; i < vec.Len(); i++ {
		result.Set(0, i, vec.AtVec(i))
	}
	return result
}
