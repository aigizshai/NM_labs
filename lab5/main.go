package main

import (
	"fmt"
	"math"
)

// Определяем типы для матрицы и вектора
type Matrix [][]float64
type Vector []float64

// Main - метод запуска Фадеева и вывода результатов
func main() {
	// Пример ввода матрицы (замените на реальное значение)
	matrix := Matrix{
		{1, -1, -1, 2},
		{2, 3, 0, -4},
		{1, 1, -2, -2},
		{1, 1, 0, -1},
	}

	characteristicPolynom := Fadeev(matrix)
	fmt.Println("Characteristic Polynomial Coefficients:", characteristicPolynom)

	eigenValues := Solve(characteristicPolynom)
	fmt.Println("Eigenvalues:", eigenValues)

	for _, eigenvalue := range eigenValues {
		eigenvector := FindEigenvector(matrix, eigenvalue)
		fmt.Printf("Eigenvector for eigenvalue %.4f: %v\n", eigenvalue, eigenvector)
	}
}

// Fadeev - Метод Фадеева для нахождения характеристического многочлена
func Fadeev(matrix Matrix) Vector {
	n := len(matrix)
	S := make(Vector, n)
	P := make(Vector, n+1)
	Matr1 := CopyMatrix(matrix)
	Matr2 := CopyMatrix(matrix)

	for i := 0; i < n; i++ {
		S[i] = trace(Matr2) / float64(i+1)
		for j := 0; j < n; j++ {
			Matr1[j][j] = Matr2[j][j] - S[i]
		}
		Matr2 = MatrixMultiplication(matrix, Matr1)
	}

	for i := 0; i < n; i++ {
		P[i+1] = -S[i]
	}
	P[0] = 1
	return P
}

// MatrixMultiplication - Перемножение двух матриц
func MatrixMultiplication(a, b Matrix) Matrix {
	n := len(a)
	c := make(Matrix, n)
	for i := 0; i < n; i++ {
		c[i] = make(Vector, n)
		for j := 0; j < n; j++ {
			sum := 0.0
			for k := 0; k < n; k++ {
				sum += a[i][k] * b[k][j]
			}
			c[i][j] = sum
		}
	}
	return c
}

// trace - След матрицы (сумма диагональных элементов)
func trace(matrix Matrix) float64 {
	sum := 0.0
	for i := 0; i < len(matrix); i++ {
		sum += matrix[i][i]
	}
	return sum
}

// CopyMatrix - Копирование матрицы
func CopyMatrix(matrix Matrix) Matrix {
	n := len(matrix)
	copy := make(Matrix, n)
	for i := 0; i < n; i++ {
		copy[i] = make(Vector, n)
		for j := 0; j < n; j++ {
			copy[i][j] = matrix[i][j]
		}
	}
	return copy
}

func Solve(poly Vector) Vector {
	n := len(poly) - 1
	results := make(Vector, n)

	for i := 0; i < n; i++ {
		x0 := 0.1 // Первое начальное приближение
		x1 := 0.2 // Второе начальное приближение (для метода секущей)
		iterCount := 0
		maxIterations := 10000

		var xn float64

		for {
			fx0 := Substitute(poly, x0)
			fx1 := Substitute(poly, x1)

			// Проверка на деление на слишком малую разницу
			if math.Abs(fx1-fx0) < 1e-8 {
				fmt.Println("Разница между fx1 и fx0 слишком мала, останавливаем вычисления.")
				break
			}

			// Метод секущей
			xn = x1 - fx1*(x1-x0)/(fx1-fx0)

			// Сдвигаем значения для следующей итерации
			x0, x1 = x1, xn

			// Проверка на сходимость или превышение количества итераций
			if math.Abs(x1-x0) < 0.00001 || iterCount > maxIterations {
				break
			}
			iterCount++
		}

		if iterCount >= maxIterations {
			fmt.Println("Достигнуто максимальное количество итераций.")
		}

		// Обновляем многочлен делением на (x - xn)
		poly = PolynomialDivision(poly, xn)
		results[i] = xn
	}
	return results
}

// Substitute - Подстановка значения в многочлен
func Substitute(poly Vector, x float64) float64 {
	result := 0.0
	for i := len(poly) - 1; i >= 0; i-- {
		result = result*x + poly[i]
	}
	return result
}

// Derivative - Производная многочлена
func Derivative(poly Vector) Vector {
	n := len(poly) - 1
	derivative := make(Vector, n)
	for i := 0; i < n; i++ {
		derivative[i] = poly[i] * float64(n-i)
	}
	return derivative
}

// PolynomialDivision - Деление многочлена на корень
func PolynomialDivision(poly Vector, root float64) Vector {
	n := len(poly)
	result := make(Vector, n-1)
	result[0] = poly[0]
	for i := 1; i < n-1; i++ {
		result[i] = poly[i] + result[i-1]*root
	}
	return result
}

// FindEigenvector - Нахождение собственных векторов
func FindEigenvector(matrix Matrix, eigenvalue float64) Vector {
	n := len(matrix)
	Matr1 := CopyMatrix(matrix)
	for i := 0; i < n; i++ {
		Matr1[i][i] -= eigenvalue
	}
	return GaussianElimination(Matr1)
}

// GaussianElimination - Метод Гаусса для нахождения решений системы уравнений
func GaussianElimination(matrix Matrix) Vector {
	n := len(matrix)
	x := make(Vector, n)

	// Initialize x with zeroes
	for i := 0; i < n; i++ {
		x[i] = 0
	}

	// Perform Gaussian elimination to row echelon form
	for i := 0; i < n; i++ {
		// Find the pivot element in the current column
		maxRow := i
		for k := i + 1; k < n; k++ {
			if math.Abs(matrix[k][i]) > math.Abs(matrix[maxRow][i]) {
				maxRow = k
			}
		}

		// Swap the current row with the maxRow if needed
		matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]

		// Check for a zero pivot (if zero, continue to the next column)
		if math.Abs(matrix[i][i]) < 1e-8 {
			continue
		}

		// Normalize the current row
		for k := i + 1; k < n; k++ {
			factor := matrix[k][i] / matrix[i][i]
			for j := i; j < n; j++ {
				matrix[k][j] -= factor * matrix[i][j]
			}
		}
	}

	// Check for free variables in row echelon form
	for i := n - 1; i >= 0; i-- {
		if math.Abs(matrix[i][i]) < 1e-8 {
			// Free variable, we can assign it any non-zero value for a non-trivial solution
			x[i] = 1 // Assign 1 or any non-zero value to suggest a non-trivial solution
		} else {
			// Back-substitution
			sum := 0.0
			for j := i + 1; j < n; j++ {
				sum += matrix[i][j] * x[j]
			}
			x[i] = -sum / matrix[i][i]
		}
	}

	return x
}
