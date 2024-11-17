package main

import (
	"fmt"
	"math"
)

type CubicSpline struct {
	a, b, c, d []float64 //коэффициенты
	x          []float64 //узлы
}

// кубический сплайн по заданным узлам x и значениям y
func buildCubicSpline(x, y []float64) CubicSpline {
	n := len(x) - 1
	h := make([]float64, n)
	for i := 0; i < n; i++ {
		h[i] = x[i+1] - x[i]
	}

	alpha := make([]float64, n)
	for i := 1; i < n; i++ {
		alpha[i] = (3/h[i])*(y[i+1]-y[i]) - (3/h[i-1])*(y[i]-y[i-1])
	}

	l := make([]float64, n+1)
	mu := make([]float64, n+1)
	z := make([]float64, n+1)
	l[0], mu[0], z[0] = 1, 0, 0

	for i := 1; i < n; i++ {
		l[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*mu[i-1]
		mu[i] = h[i] / l[i]
		z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i]
	}

	l[n], z[n] = 1, 0

	c := make([]float64, n+1)
	b := make([]float64, n)
	d := make([]float64, n)
	a := make([]float64, n)

	for j := n - 1; j >= 0; j-- {
		c[j] = z[j] - mu[j]*c[j+1]
		b[j] = (y[j+1]-y[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3
		d[j] = (c[j+1] - c[j]) / (3 * h[j])
		a[j] = y[j]
	}

	return CubicSpline{a, b, c, d, x}
}

// первая производная
func (s CubicSpline) firstDerivative(x0 float64) float64 {
	i := findSegment(s.x, x0)
	dx := x0 - s.x[i]
	return s.b[i] + 2*s.c[i]*dx + 3*s.d[i]*dx*dx
}

// вторая производная
func (s CubicSpline) secondDerivative(x0 float64) float64 {
	i := findSegment(s.x, x0)
	dx := x0 - s.x[i]
	return 2*s.c[i] + 6*s.d[i]*dx
}

func findSegment(x []float64, x0 float64) int {
	for i := 0; i < len(x)-1; i++ {
		if x0 >= x[i] && x0 <= x[i+1] {
			return i
		}
	}
	panic("x0 вне диапазона интервала")
}

func main() {
	var a, b float64
	var n int
	var x0 float64

	fmt.Print("Введите границы интервала (a и b): ")
	fmt.Scan(&a, &b)
	fmt.Print("Введите количество точек (n): ")
	fmt.Scan(&n)

	if n < 2 {
		fmt.Println("Количество точек должно быть не менее 2")
		return
	}

	fmt.Print("Введите точку для вычисления производных (x0): ")
	fmt.Scan(&x0)

	if x0 < a || x0 > b {
		fmt.Println("x0 должна быть внутри интервала [a, b]")
		return
	}

	x := make([]float64, n)
	y := make([]float64, n)

	// Генерация узлов и значений функции
	for i := 0; i < n; i++ {
		x[i] = a + float64(i)*(b-a)/float64(n-1)
		y[i] = 1.0 / x[i] * math.Exp(x[i]) //задание функци
	}

	spline := buildCubicSpline(x, y)

	first := spline.firstDerivative(x0)
	second := spline.secondDerivative(x0)

	fmt.Printf("Первая производная в точке x0 = %.2f: %.3f\n", x0, first)
	fmt.Printf("Вторая производная в точке x0 = %.2f: %.3f\n", x0, second)
}
