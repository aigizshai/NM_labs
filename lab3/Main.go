package main

import (
	"fmt"
	"math"
)

type Matr [][]float64
type Vector []float64

func Regul(n int, a Matr, b Vector, eps float64) Vector {
	a1 := makeMatrix(n, n)
	a2 := makeMatrix(n, n)
	b1 := makeVector(n)
	b2 := makeVector(n)
	x := makeVector(n)
	x0 := makeVector(n)

	for i := 0; i < n; i++ {
		for k := 0; k < n; k++ {
			s := 0.0
			for j := 0; j < n; j++ {
				s += a[j][i] * a[j][k]
			}
			a1[i][k] = s
		}
	}

	for i := 0; i < n; i++ {
		s := 0.0
		for j := 0; j < n; j++ {
			s += a[j][i] * b[j]
		}
		b1[i] = s
	}

	alpha := 0.0
	k := 0
	maxDelta := 0.0

	for {
		alpha += 1e-8
		k++
		copyMatrix(a2, a1)
		for i := 0; i < n; i++ {
			a2[i][i] += alpha
		}
		for i := 0; i < n; i++ {
			b2[i] = b1[i] + alpha*x0[i]
		}

		x = solveSimq(n, a2, b2)
		copyVector(x0, x)

		// Проверка на сходимость
		b2 = solveSimq(n, a2, b2)
		maxDelta = math.Abs(b2[0] - x[0])
		for i := 1; i < n; i++ {
			delta := math.Abs(b2[i] - x[i])
			if delta > maxDelta {
				maxDelta = delta
			}
		}

		if maxDelta < eps {
			break
		}
	}

	return x
}

func makeMatrix(rows, cols int) Matr {
	matrix := make(Matr, rows)
	for i := range matrix {
		matrix[i] = make([]float64, cols)
	}
	return matrix
}

func makeVector(size int) Vector {
	return make([]float64, size)
}

func copyMatrix(dest, src Matr) {
	for i := range src {
		for j := range src[i] {
			dest[i][j] = src[i][j]
		}
	}
}

func copyVector(dest, src Vector) {
	for i := range src {
		dest[i] = src[i]
	}
}

// Решение СЛАУ методом Гаусса
func solveSimq(n int, A Matr, B Vector) Vector {
	for i := 0; i < n; i++ {
		A[i] = append(A[i], B[i])
	}

	for k := 1; k < n; k++ {
		for j := k; j < n; j++ {
			m := A[j][k-1] / A[k-1][k-1]
			for i := 0; i < n+1; i++ {
				A[j][i] = A[j][i] - m*A[k-1][i]
			}
		}
	}

	X := make([]float64, n)

	for i := n - 1; i >= 0; i-- {
		X[i] = A[i][n] / A[i][i]
		for c := n - 1; c > i; c-- {
			X[i] = X[i] - A[i][c]*X[c]/A[i][i]
		}
	}

	return X
}

func Vrash(n int, a Matr, b Vector) Vector {
	// Создаем расширенную матрицу b
	aExt := make(Matr, n)
	for i := 0; i < n; i++ {
		aExt[i] = make([]float64, n+1)
		copy(aExt[i], a[i])
		aExt[i][n] = b[i]
	}

	for i := 0; i < n-1; i++ {
		for k := i + 1; k < n; k++ {
			var m, l float64
			if aExt[i][i] != 0.0 || aExt[k][i] != 0.0 {
				m = math.Sqrt(aExt[i][i]*aExt[i][i] + aExt[k][i]*aExt[k][i])
				l = -aExt[k][i] / m
				m = aExt[i][i] / m
			} else {
				m = 1.0
				l = 0.0
			}

			for j := i; j < n+1; j++ {
				r := m*aExt[i][j] - l*aExt[k][j]
				aExt[k][j] = l*aExt[i][j] + m*aExt[k][j]
				aExt[i][j] = r
			}
		}
	}

	x := make(Vector, n)
	for i := n - 1; i >= 0; i-- {
		m := 0.0
		for k := 0; k < n-i-1; k++ {
			m += aExt[i][n-k-1] * x[n-k-1]
		}
		x[i] = (aExt[i][n] - m) / aExt[i][i]
	}

	return x
}

func main() {
	n := 2
	eps := 0.0001
	a := Matr{
		{1.03, 0.991},
		{0.991, 0.943},
	}
	b := Vector{2.58, 2.48}

	fmt.Println("Исходная матрица")
	for i := 0; i < n; i++ {
		fmt.Println(a[i])
	}
	fmt.Println("Столбец B")
	fmt.Println(b[0])
	fmt.Println(b[1], "\n")

	x := Regul(n, a, b, eps)
	fmt.Println("Решение методом регуляризации:")
	fmt.Printf("x1=%.4f\n", x[0])
	fmt.Printf("x2=%.4f\n", x[1])

	fmt.Println("Решение методом вращения:")
	x = Vrash(n, a, b)
	fmt.Printf("x1=%.4f\n", x[0])
	fmt.Printf("x2=%.4f\n", x[1])
}
