package main

import (
	"fmt"
)

// интерполяция Лагранжа
func lagrange(xPoints, yPoints []float64, x float64) float64 {

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

// добавляем значения к массиву за пределы
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

// линейная экстраполяция
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

// интерполяция сплайном
func Spline(n int, x, y []float64, xx float64) float64 {
	a := make([]float64, n)
	b := make([]float64, n)
	c := make([]float64, n)
	d := make([]float64, n)
	z := make([]float64, n)

	//коэффициенты трехдиагональной матрицы
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

	z = Progon(a, b, c, d, n)

	if xx == x[n-1] {
		return y[n-1]
	}

	//значение функции в точке xx
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

// метод Ньютона
func NewtonInterpolation(x []float64, y []float64, x0 float64) float64 {
	n := len(x)

	dividedDifferences := make([][]float64, n)
	for i := range dividedDifferences {
		dividedDifferences[i] = make([]float64, n)
		dividedDifferences[i][0] = y[i]
	}

	//разделенные разности
	for j := 1; j < n; j++ {
		for i := 0; i < n-j; i++ {
			dividedDifferences[i][j] = (dividedDifferences[i+1][j-1] - dividedDifferences[i][j-1]) / (x[i+j] - x[i])
		}
	}

	//значение в точке x0
	result := dividedDifferences[0][0]
	term := 1.0
	for i := 1; i < n; i++ {
		term *= (x0 - x[i-1])
		result += dividedDifferences[0][i] * term
	}

	return result
}

// Функция для аппроксимации методом наименьших квадратов
func LeastSquares(n int, x, y []float64, q float64) (float64, float64, float64, float64) {
	//матрица коэффициентов
	a := make([][]float64, 3)
	for i := range a {
		a[i] = make([]float64, 3)
	}

	//вектор правых частей
	b := make([]float64, 3)

	//вычисление элементов
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

	//метод Гаусса
	coefficients := gaussSolve(a, b)

	c0, c1, c2 := coefficients[0], coefficients[1], coefficients[2]
	result := c0 + c1*q + c2*q*q
	return result, c0, c1, c2
}

// метод Гаусса
func gaussSolve(a [][]float64, b []float64) []float64 {
	n := len(b)

	//прямой ход
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			ratio := a[j][i] / a[i][i]
			for k := i; k < n; k++ {
				a[j][k] -= ratio * a[i][k]
			}
			b[j] -= ratio * b[i]
		}
	}

	//обратный ход
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
	//x := []float64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
	//y := []float64{2.05, 1.94, 1.92, 1.87, 1.77, 1.88, 1.71, 1.60, 1.56, 1.40}
	x := []float64{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0}
	y := []float64{-1.92, -1.60, -1.57, -1.41, -1.36, -0.97, -0.59, -0.71, -0.15, 0.01, 0.22, 0.63, 1.07, 1.42, 1.68, 2.49, 2.57, 3.09, 3.40, 4.0}
	//y := []float64{1.99, 2.03, 2.20, 2.39, 2.19, 2.61, 2.35, 2.60, 2.55, 2.49, 2.50, 2.52, 2.44, 2.35, 2.26, 2.19, 2.24, 2.34, 1.96, 2.19}
	dx := 0.05

	q := x[0]
	c := 0
	var resultSpline float64
	fmt.Printf("Знач х  |Исходные данные|Лагранж\t\t|Сплайн\t\t| Ньютон\t|МНК\n")
	for i := 0; i < 39; i++ {

		resultLagr := lagrange(x, y, q)
		if i == 38 {
			resultSpline = Spline(len(x), x, y, x[len(x)-1])
		} else {
			resultSpline = Spline(len(x), x, y, q)
		}
		//resultNewton := NewtonInterpolation(len(x), x, y, q)
		resultNewton := NewtonInterpolation(x, y, q)
		resultLeastSqr, _, _, _ := LeastSquares(len(x), x, y, q)

		//fmt.Printf("x= %3.2f, y = %6.3f, y = %6.3f, y = %6.3f, y = %6.5f\n", q, resultLagr, resultSpline, resultNewton, resultLeastSqr)

		fmt.Printf("x= %3.2f | y =%8.2f | y = %8.3f | y = %8.3f | y = %8.3f | y = %8.3f\n", q, y[c], resultLagr, resultSpline, resultNewton, resultLeastSqr)

		if i%2 == 1 {
			c++
		}
		q += dx

	}
	_, c0, c1, c2 := LeastSquares(len(x), x, y, q)
	fmt.Println()
	fmt.Printf("Коэффициенты: c0 = %.5f, c1 = %.5f, c2 = %.5f\n", c0, c1, c2)

}
