package main

import (
	"fmt"
)

// lagrangeInterpolation выполняет интерполяцию с расширением границ.
func lagrange(xPoints, yPoints []float64, x float64) float64 {

	//расширяем границы
	xExtended, yExtended := extendBounds(xPoints, yPoints)
	n := len(xExtended)
	result := 0.0
	for i := 0; i < n; i++ {

		L := 1.0
		for j := 0; j < n; j++ {
			if i != j {
				L *= (x - xExtended[j]) / (xExtended[i] - xExtended[j])
			}
		}

		result += L * yExtended[i]

	}

	return result
}

// добавляет точки за пределы интервала
func extendBounds(xPoints, yPoints []float64) ([]float64, []float64) {
	n := len(xPoints)

	dxLeft := xPoints[1] - xPoints[0] // Шаг между первыми точками
	xMin1 := xPoints[0] - dxLeft      // Первая левая точка
	xMin2 := xPoints[0] - 2*dxLeft    // Вторая левая точка
	yMin1 := linearExtrapolation(xPoints[0], xPoints[1], yPoints[0], yPoints[1], xMin1)
	yMin2 := linearExtrapolation(xPoints[0], xPoints[1], yPoints[0], yPoints[1], xMin2)

	dxRight := xPoints[n-1] - xPoints[n-2] // Шаг между последними точками
	xMax1 := xPoints[n-1] + dxRight        // Первая правая точка
	xMax2 := xPoints[n-1] + 2*dxRight      // Вторая правая точка
	yMax1 := linearExtrapolation(xPoints[n-2], xPoints[n-1], yPoints[n-2], yPoints[n-1], xMax1)
	yMax2 := linearExtrapolation(xPoints[n-2], xPoints[n-1], yPoints[n-2], yPoints[n-1], xMax2)

	// Добавляем новые точки к массивам
	xExtended := append([]float64{xMin2, xMin1}, xPoints...)
	xExtended = append(xExtended, xMax1, xMax2)

	yExtended := append([]float64{yMin2, yMin1}, yPoints...)
	yExtended = append(yExtended, yMax1, yMax2)

	return xExtended, yExtended
}

// выполняет линейную экстраполяцию между двумя точками
func linearExtrapolation(x1, x2, y1, y2, x float64) float64 {
	return y1 + (y2-y1)*(x-x1)/(x2-x1)
}

func Progon(a, b, c, d []float64, n int) []float64 {
	u := make([]float64, n+1)
	v := make([]float64, n+1)
	x := make([]float64, n)

	v[1] = 0
	u[1] = 0

	for i := 1; i < n; i++ {
		z := 1.0 / (b[i] - a[i]*v[i])
		v[i+1] = c[i] * z
		u[i+1] = (a[i]*u[i] - d[i]) * z
	}

	x[n-1] = u[n]
	for i := n - 2; i >= 0; i-- {
		x[i] = v[i+1]*x[i+1] + u[i+1]
	}

	return x
}

// Функция интерполяции с помощью кубического сплайна
func Spline(n int, x, y []float64, xx float64) float64 {
	a := make([]float64, n)
	b := make([]float64, n)
	c := make([]float64, n)
	d := make([]float64, n)
	z := make([]float64, n)

	// Построение коэффициентов трехдиагональной матрицы
	a[1] = 0
	b[1] = -2
	c[1] = 1
	d[1] = 3 * (y[1] - y[0]) / (x[1] - x[0])

	for j := 1; j < n-1; j++ {
		hj := x[j+1] - x[j]
		hj1 := x[j] - x[j-1]
		am := hj1 / (hj1 + hj)
		al := 1 - am
		a[j] = al
		b[j] = -2
		c[j] = am
		d[j] = 3 * (am*(y[j+1]-y[j])/hj + al*(y[j]-y[j-1])/hj1)
	}

	a[n-1] = 1
	b[n-1] = -2
	c[n-1] = 0
	d[n-1] = 3 * (y[n-1] - y[n-2]) / (x[n-1] - x[n-2])

	// Решение системы методом прогонки
	z = Progon(a, b, c, d, n)

	// Если xx равен последнему узлу, возвращаем значение функции в последнем узле
	if xx == x[n-1] {
		return y[n-1]
	}

	// Вычисление значения функции в точке xx
	var s float64
	for j := 1; j < n; j++ {
		if x[j] > xx {
			hj := x[j] - x[j-1]
			t := (xx - x[j-1]) / hj
			t1 := 1 - t
			t2 := t * t
			t12 := t1 * t1
			s = y[j-1]*t12*(1+2*t) + y[j]*t2*(3-2*t) + z[j-1]*hj*t*t12 - z[j]*hj*t2*t1
			break
		}
	}
	return s
}

// Функция для интерполяции методом Ньютона
func NewtonInterpolation(n int, a, b []float64, x float64) float64 {
	// Создаем массив для разделенных разностей
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
	}

	// Заполняем первую строку разделенными разностями
	for i := 1; i < n; i++ {
		m[0][i-1] = b[i] - b[i-1]
	}

	// Вычисляем разделенные разности для каждого порядка
	for i := 1; i < n-1; i++ {
		f := factorial(i + 1)
		for j := 1; j < n-i; j++ {
			m[i][j-1] = (m[i-1][j] - m[i-1][j-1]) / float64(f)
		}
	}

	// Определяем, где начнется интерполяция (с начала или с конца массива)
	var s float64
	if (a[n-1]-x)-(x-a[0]) < 0 {
		// Интерполяция от конца массива
		s = b[n-1]
		for i := 1; i < n; i++ {
			p := 1.0
			for j := 0; j < i; j++ {
				p *= (x - a[n-1-j])
			}
			s += p * m[i-1][n-i-1]
		}
	} else {
		// Интерполяция от начала массива
		s = b[0]
		for i := 1; i < n; i++ {
			p := 1.0
			for j := 0; j < i; j++ {
				p *= (x - a[j])
			}
			s += p * m[i-1][0]
		}
	}

	return s
}

// Вспомогательная функция для вычисления факториала
func factorial(n int) int {
	if n == 0 {
		return 1
	}
	f := 1
	for i := 2; i <= n; i++ {
		f *= i
	}
	return f
}

// Функция для аппроксимации методом наименьших квадратов
func LeastSquares(n int, x, y []float64, q float64) (float64, float64, float64, float64) {
	// Матрица коэффициентов для системы уравнений
	a := make([][]float64, 3)
	for i := range a {
		a[i] = make([]float64, 3)
	}

	// Вектор правых частей
	b := make([]float64, 3)

	// Вычисление элементов матрицы `a` и вектора `b`
	for i := 0; i < n; i++ {
		a[0][0] += 1
		a[0][1] += x[i]
		a[0][2] += x[i] * x[i]
		a[1][0] += x[i]
		a[1][1] += x[i] * x[i]
		a[1][2] += x[i] * x[i] * x[i]
		a[2][0] += x[i] * x[i]
		a[2][1] += x[i] * x[i] * x[i]
		a[2][2] += x[i] * x[i] * x[i] * x[i]

		b[0] += y[i]
		b[1] += y[i] * x[i]
		b[2] += y[i] * x[i] * x[i]
	}

	// Решаем систему уравнений методом Гаусса
	coefficients := gaussSolve(a, b)

	// Возвращаем коэффициенты c0, c1, c2 и значение полинома в точке q
	c0, c1, c2 := coefficients[0], coefficients[1], coefficients[2]
	result := c0 + c1*q + c2*q*q
	return result, c0, c1, c2
}

// Функция для решения системы линейных уравнений методом Гаусса
func gaussSolve(a [][]float64, b []float64) []float64 {
	n := len(b)

	// Прямой ход метода Гаусса
	for i := 0; i < n; i++ {
		// Нормализация строки
		for j := i + 1; j < n; j++ {
			ratio := a[j][i] / a[i][i]
			for k := i; k < n; k++ {
				a[j][k] -= ratio * a[i][k]
			}
			b[j] -= ratio * b[i]
		}
	}

	// Обратный ход метода Гаусса
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		x[i] = b[i]
		for j := i + 1; j < n; j++ {
			x[i] -= a[i][j] * x[j]
		}
		x[i] /= a[i][i]
	}

	return x
}

func main() {
	x := []float64{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0}
	y := []float64{-1.92, -1.60, -1.57, -1.41, -1.36, -0.97, -0.59, -0.71, -0.15, 0.01, 0.22, 0.63, 1.07, 1.42, 1.68, 2.49, 2.57, 3.09, 3.40, 4.0}
	//y := []float64{1.99, 2.03, 2.20, 2.39, 2.19, 2.61, 2.35, 2.60, 2.55, 2.49, 2.50, 2.52, 2.44, 2.35, 2.26, 2.19, 2.24, 2.34, 1.96, 2.19}
	dx := 0.05

	q := x[0]

	fmt.Println("Интерполяция методом Лагранжа")
	for i := 0; i < 39; i++ {
		resultLagr := lagrange(x, y, q)
		resultSpline := Spline(len(x), x, y, q)
		resultNewton := NewtonInterpolation(len(x), x, y, q)
		resultLeastSqr, _, _, _ := LeastSquares(len(x), x, y, q)

		fmt.Printf("x= %3.2f y = %6.3f, y = %6.3f, y = %6.3f, y = %6.3f\n", q, resultLagr, resultSpline, resultNewton, resultLeastSqr)
		//		fmt.Println("q=", q)
		q += dx

	}
	_, c0, c1, c2 := LeastSquares(len(x), x, y, q)

	fmt.Printf("Коэффициенты: c0 = %.5f, c1 = %.5f, c2 = %.5f\n", c0, c1, c2)

	// n := len(x)

	// fmt.Println("Интерполяция функций с помощью кубического сплайна:")
	// xx := 1.0
	// for i := 0; i < 39; i++ {
	// 	s := Spline(n, x, y, xx)
	// 	fmt.Printf("x = %.2f, y = %.3f\n", xx, s)
	// 	xx += dx
	// }

	// q = x[0]

	// fmt.Println("Интерполяция методом Ньютона")
	// for i := 0; i < 39; i++ {
	// 	result := NewtonInterpolation(len(x), x, y, q)
	// 	fmt.Printf("x = %.2f y = %.3f\n", q, result)
	// 	q += dx

	// }

	// fmt.Println("Аппроксимация методом наименьших квадратов")
	// q = x[0]
	// result, c0, c1, c2 := LeastSquares(len(x), x, y, q)

	// fmt.Printf("Коэффициенты: c0 = %.5f, c1 = %.5f, c2 = %.5f\n", c0, c1, c2)
	// fmt.Printf("x = %.2f y = %.3f\n", q, result)
	// q += dx
	// for i := 0; i < 38; i++ {
	// 	result, c0, c1, c2 = LeastSquares(len(x), x, y, q)
	// 	fmt.Printf("x = %.2f y = %.3f, c0=%.4f,c1=%.4f,c2=%.4f\n", q, result, c0, c1, c2)
	// 	q += dx

	// }

	//сдедлать вывод в несколько колонок и поменять нач значения
}
