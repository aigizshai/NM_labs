package newton

import (
	"fmt"
	"math"
)

const (
	n   = 2
	eps = 1e-11
)

type (
	Matrix [n][n]float64
	Vector [n]float64
)

// Задание функций системы уравнений
func Func(x Vector, i int) float64 {
	switch i {
	case 1:
		return math.Cos(x[0]+0.5) + x[1] - 0.8
	case 2:
		return math.Sin(x[1]) - 2*x[0] - 1.6
	default:
		return 0.0
	}
}

// Задание элементов матрицы Якоби
func Jacobian(x Vector, i, j int) float64 {
	switch i {
	case 1:
		switch j {
		case 1:
			return -math.Sin(x[0] + 0.5)
		case 2:
			return 1.0
		}
	case 2:
		switch j {
		case 1:
			return -2.0
		case 2:
			return math.Cos(x[1])
		}
	}
	return 0.0
}

// Вывод вектора
func PrintVector(vector Vector) {
	for i := 0; i < n; i++ {
		fmt.Printf("x%d = %10.5f\n", i+1, vector[i])
	}
}

// Решение системы линейных уравнений методом Гаусса
func SolveLinearSystem(a Matrix, b Vector) Vector {
	var x Vector
	augmented := make([][]float64, n)

	for i := 0; i < n; i++ {
		augmented[i] = make([]float64, n+1)
		for j := 0; j < n; j++ {
			augmented[i][j] = a[i][j]
		}
		augmented[i][n] = b[i]
	}

	for i := 0; i < n; i++ {
		maxElem := math.Abs(augmented[i][i])
		maxRow := i

		for k := i + 1; k < n; k++ {
			if math.Abs(augmented[k][i]) > maxElem {
				maxElem = math.Abs(augmented[k][i])
				maxRow = k
			}
		}

		for k := i; k < n+1; k++ {
			tmp := augmented[maxRow][k]
			augmented[maxRow][k] = augmented[i][k]
			augmented[i][k] = tmp
		}

		for k := i + 1; k < n; k++ {
			c := -augmented[k][i] / augmented[i][i]
			for j := i; j < n+1; j++ {
				if i == j {
					augmented[k][j] = 0
				} else {
					augmented[k][j] += c * augmented[i][j]
				}
			}
		}
	}

	for i := n - 1; i >= 0; i-- {
		x[i] = augmented[i][n] / augmented[i][i]
		for k := i - 1; k >= 0; k-- {
			augmented[k][n] -= augmented[k][i] * x[i]
		}
	}

	return x
}

func Newton() (Vector, int) {
	var (
		a        Matrix
		x, f, dx Vector
		iter     int
		max      float64
	)

	// Начальное приближение
	x[0] = 0.0
	x[1] = 0.0
	iter = 0

	for {
		// PrintVector(x)
		// fmt.Printf("Итерация № %d\n", iter)
		// fmt.Println("=================")

		// Подсчет количества итераций

		// Вычисление значений элементов матрицы Якоби
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				a[i][j] = Jacobian(x, i+1, j+1)
			}
		}

		// Вычисление правой части СЛАУ
		for i := 0; i < n; i++ {
			f[i] = -Func(x, i+1)
		}

		// Решение системы для нахождения поправок
		dx = SolveLinearSystem(a, f)

		// Вычисление максимального изменения
		max = math.Abs(dx[0])
		for i := 1; i < n; i++ {
			if math.Abs(dx[i]) > max {
				max = math.Abs(dx[i])
			}
		}

		// Обновление приближения
		for i := 0; i < n; i++ {
			x[i] += dx[i]
		}

		// Условие завершения
		if max < eps {
			return x, iter
		}
		iter++
	}
}
