package main

import "fmt"

// Функция вычисления значения интерполяционного многочлена Лагранжа
func Lagrange(n int, x, y []float64, q float64) float64 {
	L := 0.0

	for i := 0; i < n; i++ {
		s := 1.0
		for j := 0; j < n; j++ {
			if j != i {
				s *= (q - x[j]) / (x[i] - x[j])
			}
		}
		L += y[i] * s
	}

	return L
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
	// Пример использования функции Lagrange
	x := []float64{1.0, 2.0, 3.0, 4, 5, 6, 7, 8, 9, 10}
	y := []float64{2.05, 1.94, 1.92, 1.87, 1.77, 1.88, 1.71, 1.60, 1.56, 1.40}
	dx := 0.5

	q := x[0]

	fmt.Println("Интерполяция методом Лагранжа")
	for i := 0; i < 19; i++ {
		result := Lagrange(len(x), x, y, q)
		fmt.Printf("x= %.2f y=%.3f\n", q, result)
		q += dx

	}

	n := len(x)

	fmt.Println("Интерполяция функций с помощью кубического сплайна:")
	xx := 1.0
	for i := 0; i < 19; i++ {
		s := Spline(n, x, y, xx)
		fmt.Printf("x = %.2f, y = %.3f\n", xx, s)
		xx += 0.5
	}

	q = x[0]

	fmt.Println("Интерполяция методом Ньютона")
	for i := 0; i < 19; i++ {
		result := NewtonInterpolation(len(x), x, y, q)
		fmt.Printf("x= %.2f y=%.3f\n", q, result)
		q += dx

	}

	fmt.Println("Аппроксимация методом наименьших квадратов")
	q = x[0]
	result, c0, c1, c2 := LeastSquares(n, x, y, q)

	fmt.Printf("Коэффициенты: c0 = %.5f, c1 = %.5f, c2 = %.5f\n", c0, c1, c2)
	fmt.Printf("x = %.2f: %.3f\n", q, result)
	q += dx
	for i := 0; i < 18; i++ {
		result, c0, c1, c2 = LeastSquares(n, x, y, q)
		fmt.Printf("x = %.2f: %.3f\n", q, result)
		q += dx

	}

	//сдедлать вывод в несколько колонок и поменять нач значения
}
