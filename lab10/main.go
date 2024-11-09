package main

import (
	"fmt"
)

type Matrix [][]float64
type Vector []float64

// Функция для умножения двух квадратных матриц размера n x n
func mulMat(n int, a, b Matrix) Matrix {
	c := make(Matrix, n)
	for i := range c {
		c[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				c[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	return c
}

// Функция для умножения матрицы на вектор
func matVecMul(n int, a Matrix, b Vector) Vector {
	c := make(Vector, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			c[i] += a[i][j] * b[j]
		}
	}
	return c
}

// Функция для умножения матрицы на скаляр
func mulConst(n int, h float64, a Matrix) Matrix {
	b := make(Matrix, n)
	for i := range a {
		b[i] = make([]float64, n)
		for j := range a[i] {
			b[i][j] = a[i][j] * h
		}
	}
	return b
}

// Метод Ракитского
func metRak(n int, h *float64, A Matrix, Y *Vector, Q *Matrix, ind, nts int) {
	if ind == 0 { // Вычисление матричной экспоненты
		// Вычисление нормы матрицы A
		s := 0.0
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				s += A[i][j] * A[i][j]
			}
		}
		// s = math.Sqrt(s)
		// *h = 0.1 / s

		// Инициализация Q как нулевой матрицы
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				(*Q)[i][j] = 0
			}
		}

		// Умножение матрицы A на шаг h
		A = mulConst(n, *h, A)

		// Вычисление экспоненты матрицы с помощью разложения в ряд
		for k := 5; k >= 1; k-- {
			*Q = mulConst(n, 1/float64(k), *Q)
			for i := 0; i < n; i++ {
				(*Q)[i][i] += 1
			}
			*Q = mulMat(n, *Q, A)
		}
		for i := 0; i < n; i++ {
			(*Q)[i][i] += 1
		}

		// Построение матрицы экспоненты exp(A*h)
		for i := 0; i < n; i++ {
			y0 := make(Vector, n)
			y0[i] = 1
			for j := 0; j < int(0.1 / *h); j++ {
				y0 = matVecMul(n, *Q, y0)
			}
			for j := 0; j < n; j++ {
				(*Q)[j][i] = y0[j]
			}
		}

	} else if nts > 0 { // Вычисление решения на текущем шаге
		*Y = matVecMul(n, *Q, *Y)
	}
}

func main() {
	// Исходные данные
	n := 4
	a := Matrix{
		{233.70, 359.13, 2461.63, 235.38},
		{13.424, 14.239, 1.7391, 12.554},
		{105.79, 167.88, 124.13, 106.22},
		{232.77, 357.39, 244.89, 235.33},
	}
	y := Vector{-0.8462, -0.5328, -0.4322, 0.9018}
	q := make(Matrix, n)
	for i := range q {
		q[i] = make([]float64, n)
	}

	h := 0.1
	x := 0.0

	// Вычисление матричной экспоненты
	metRak(n, &h, a, &y, &q, 0, 0)

	fmt.Println(" x\tY1\tY2\tY3\tY4")
	fmt.Printf("x=%.1f\t%.5e\t%.5e\t%.5e\t%.5e\n", x, y[0], y[1], y[2], y[3])
	// Численное решение задачи методом Ракитского
	for nts := 1; nts < 21; nts++ {
		metRak(n, &h, a, &y, &q, 1, nts)
		x += h
		fmt.Printf("x=%.1f\t%.5e\t%.5e\t%.5e\t%.5e\n", x, y[0], y[1], y[2], y[3])
	}
}

//добавить метод Гира
