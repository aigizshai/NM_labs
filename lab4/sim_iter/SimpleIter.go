package sim_iter

import (
	"fmt"
	"math"
)

const (
	n   = 2
	eps = 1e-4
)

type (
	Matrix [n][n]float64
	Vector [n]float64
)

// Вычисление нормы матрицы
func Norma(a Matrix, n int) float64 {
	res := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			res += a[i][j] * a[i][j]
		}
	}
	return math.Sqrt(res)
}

// Задание функции
func Func(x Vector, i int) float64 {
	switch i {
	case 1:
		return (math.Sin(x[1]) - 1.6) / 2
	case 2:
		return 0.8 - math.Cos(x[0]+0.5)
	default:
		return 0.0
	}
}

// Построение матрицы Якоби
func MatrJacobi(x Vector, i, j int) float64 {
	switch i {
	case 1:
		switch j {
		case 1:
			return math.Sin(x[0] + 0.5)
		case 2:
			return 0.0
		}
	case 2:
		switch j {
		case 1:
			return 0.0
		case 2:
			return 0.5 * math.Cos(x[1])
		}
	}
	return 0.0
}

// Вывод матрицы
func PrintMatrix(mat Matrix) {
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			fmt.Printf("%10.5f ", mat[i][j])
		}
		fmt.Println()
	}
}

// Вывод вектора
func PrintVector(vector Vector) {
	for j := 0; j < n; j++ {
		fmt.Printf("x%d = %10.5f\n", j+1, vector[j])
	}
}

func SimpleIter() (Vector, float64, int) {
	//х1,х2, норма, колво итераций
	var (
		a     Matrix
		x, x0 Vector
		iter  int
		max   float64
	)

	// Начальные приближения
	x0[0] = 0
	x0[1] = 0
	iter = 0

	for {
		// Построение матрицы Якоби
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				a[i][j] = MatrJacobi(x0, i+1, j+1)
			}
		}

		// Вывод вектора x0
		// PrintVector(x0)

		// Вычисление нормы матрицы A
		norma := Norma(a, n)
		// fmt.Printf("Норма = %10.5f\n", norma)

		// Подсчет количества итераций
		// fmt.Printf("Номер итерации - %d\n", iter)
		// fmt.Println("=================")

		// Нахождение нового приближения функции
		for i := 0; i < n; i++ {
			x[i] = Func(x0, i+1)
		}

		// Нахождение максимального изменения
		max = math.Abs(x[0] - x0[0])
		for i := 1; i < n; i++ {
			diff := math.Abs(x[i] - x0[i])
			if diff > max {
				max = diff
			}
		}

		// Обновляем x0
		x0 = x

		// Условие остановки
		if max < eps || iter > 20 {
			return x0, norma, iter

		}
		iter++
	}

}
