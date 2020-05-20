package cholesky

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

// Formatted returns a fmt.Formatter for the matrix m using the given options.
func Formatted(m mat.CMatrix) fmt.Formatter {
	f := formatter{
		matrix: m,
		dot:    '.',
	}

	return f
}

type formatter struct {
	matrix mat.CMatrix
	dot    byte

	format func(m mat.CMatrix, dot byte, fs fmt.State, c rune)
}

func (f formatter) Format(fs fmt.State, c rune) {
	if f.format == nil {
		f.format = format
	}
	f.format(f.matrix, f.dot, fs, c)
}

func format(m mat.CMatrix, dot byte, fs fmt.State, c rune) {
	rows, cols := m.Dims()

	skipZero := fs.Flag(' ')

	for i := 0; i < rows; i++ {
		var el string
		switch {
		case rows == 1:
			fmt.Fprint(fs, "[")
			el = "]"
		case i == 0:
			fmt.Fprint(fs, "⎡")
			el = "  \t⎤\n"
		case i < rows-1:
			fmt.Fprint(fs, "⎢")
			el = "  \t⎥\n"
		default:
			fmt.Fprint(fs, "⎣")
			el = "  \t⎦"
		}

		for j := 0; j < cols; j++ {
			v := m.At(i, j)
			if v == 0 && skipZero {
				fmt.Fprint(fs, ".     ")
			} else {
				fmt.Fprintf(fs, "%.4v", v)
			}

			if j < cols-1 {
				fs.Write([]byte("   \t"))
			}
		}

		fmt.Fprint(fs, el)
	}
}
