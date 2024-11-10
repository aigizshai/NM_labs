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

// Метод Рунге-Кутты четвертого порядка для инициализации первых 4 шагов
func rungeKutta(a, b float64, n int, x float64, y Vector, A Matrix) {
	h := (b - a) / float64(n)
	k := make([]Vector, 4)
	for i := range k {
		k[i] = make(Vector, len(y))
	}

	for j := 0; j < len(y); j++ {
		// Вычисляем промежуточные значения k1, k2, k3 и k4
		k1 := matVecMul(len(A), A, y)
		for i := range k1 {
			k[0][i] = h * k1[i]
		}

		yTemp := make(Vector, len(y))
		for i := range y {
			yTemp[i] = y[i] + 0.5*k[0][i]
		}
		k2 := matVecMul(len(A), A, yTemp)
		for i := range k2 {
			k[1][i] = h * k2[i]
		}

		for i := range y {
			yTemp[i] = y[i] + 0.5*k[1][i]
		}
		k3 := matVecMul(len(A), A, yTemp)
		for i := range k3 {
			k[2][i] = h * k3[i]
		}

		for i := range y {
			yTemp[i] = y[i] + k[2][i]
		}
		k4 := matVecMul(len(A), A, yTemp)
		for i := range k4 {
			k[3][i] = h * k4[i]
		}

		// Обновляем y по формуле Рунге-Кутты
		for i := range y {
			y[i] += (k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i]) / 6
		}
	}
}

// Метод Гира для решения системы ОДУ
func gear(a, b float64, n int, x float64, y Vector, A Matrix) {
	h := (b - a) / float64(n)
	fun := make([][]float64, len(y))
	for i := range fun {
		fun[i] = make([]float64, 4)
	}

	// Инициализация с помощью метода Рунге-Кутты
	for i := 0; i < 4; i++ {
		rungeKutta(a, b, n, x, y, A)
		x += h
		for j := 0; j < len(y); j++ {
			fun[j][i] = y[j]
		}
	}

	// Основной цикл метода Гира
	for i := 0; i < n; i++ {
		// Вычисляем коэффициенты c1 и c2 для метода Гира
		c := make(Vector, len(y))
		for j := range y {
			c[j] = (-48*fun[j][3] + 36*fun[j][2] - 16*fun[j][1] + 3*fun[j][0]) / 1.2
		}

		// Решаем систему линейных уравнений для y
		denom := (25/1.2+11)*(25/1.2+11) - 9*9
		y[0] = (-9*c[1] - (25/1.2+11)*c[0]) / denom
		y[1] = (9*y[0] - c[1]) / (25/1.2 + 11)

		// Обновляем прошлые значения для метода Гира
		for j := 0; j < len(y); j++ {
			copy(fun[j][:3], fun[j][1:])
			fun[j][3] = y[j]
		}

		// Вывод значений
		fmt.Printf("x=%.1f \t y1=%.4e \t y2=%.4e\t y3=%.4e\t y4=%.4e\n", x, y[0], y[1], y[2], y[3])
		x += h
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
	x = 0.0
	a = Matrix{
		{233.70, 359.13, 2461.63, 235.38},
		{13.424, 14.239, 1.7391, 12.554},
		{105.79, 167.88, 124.13, 106.22},
		{232.77, 357.39, 244.89, 235.33},
	}
	y = Vector{-0.8462, -0.5328, -0.4322, 0.9018}
	fmt.Println("Решение системы методом Гира:")
	gear(0.0, 2.0, 20, x, y, a)
}

//добавить метод Гира, проверить и исправить
