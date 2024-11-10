package main

import (
	"fmt"
	"math"
)

// Определение системы дифференциальных уравнений

func derivatives(x, y1, y2 float64) (float64, float64) {
	dy1 := math.Cos(x + y2)
	dy2 := math.Sin(x - y2)
	return dy1, dy2
}

func f1(x, y1, y2 float64) float64 {
	return math.Cos(x + y2)
}

func f2(x, y1, y2 float64) float64 {
	return math.Sin(x - y2)
}

// Метод Рунге-Кутта 4-го порядка для начальных значений
func rungeKutta4(a, b, h float64, y1Init, y2Init float64) ([]float64, []float64, []float64) {
	//h := (b - a) / float64(n) // Шаг
	n := int((b - a) / h)
	// Массивы для хранения значений x, y1 и y2 на каждом шаге
	x := make([]float64, n+1)
	y1 := make([]float64, n+1)
	y2 := make([]float64, n+1)

	// Начальные условия
	x[0] = a
	y1[0] = y1Init
	y2[0] = y2Init

	// Основной цикл метода Рунге-Кутта
	for i := 0; i < n; i++ {
		k1_y1, k1_y2 := derivatives(x[i], y1[i], y2[i])
		k2_y1, k2_y2 := derivatives(x[i]+h/2, y1[i]+h/2*k1_y1, y2[i]+h/2*k1_y2)
		k3_y1, k3_y2 := derivatives(x[i]+h/2, y1[i]+h/2*k2_y1, y2[i]+h/2*k2_y2)
		k4_y1, k4_y2 := derivatives(x[i]+h, y1[i]+h*k3_y1, y2[i]+h*k3_y2)

		// Обновляем значения y1 и y2
		y1[i+1] = y1[i] + h/6*(k1_y1+2*k2_y1+2*k3_y1+k4_y1)
		y2[i+1] = y2[i] + h/6*(k1_y2+2*k2_y2+2*k3_y2+k4_y2)
		x[i+1] = x[i] + h // Обновляем значение x
	}

	return x, y1, y2
}

// Метод Адамса четвертого порядка
func adams4(a, b, h float64, y1Init, y2Init float64) ([]float64, []float64, []float64) {
	n := int((b - a) / h)

	x := make([]float64, n+1)
	y1 := make([]float64, n+1)
	y2 := make([]float64, n+1)

	// Инициализация начальных условий
	x[0] = a
	y1[0] = y1Init
	y2[0] = y2Init

	// Используем метод Рунге-Кутта для первых четырех значений
	x, y1, y2 = rungeKutta4(a, b, h, y1Init, y2Init)
	//rungeKutta4(a, h, y1, y2, n)

	for i := 4; i <= n; i++ {
		y1[i] = y1[i-1] + h/24*(55*f1(x[i-1], y1[i-1], y2[i-1])-
			59*f1(x[i-2], y1[i-2], y2[i-2])+
			37*f1(x[i-3], y1[i-3], y2[i-3])-
			9*f1(x[i-4], y1[i-4], y2[i-4]))

		y2[i] = y2[i-1] + h/24*(55*f2(x[i-1], y1[i-1], y2[i-1])-
			59*f2(x[i-2], y1[i-2], y2[i-2])+
			37*f2(x[i-3], y1[i-3], y2[i-3])-
			9*f2(x[i-4], y1[i-4], y2[i-4]))

		x[i] = x[i-1] + h
	}

	return x, y1, y2
}

func main() {
	a := 0.0       // начальная точка
	b := 4.0       // конечная точка
	h := 0.1       // шаг
	y1Init := 0.7  // начальное значение y1
	y2Init := -0.5 // начальное значение y2

	// метод Адамса
	x, y1, y2 := adams4(a, b, h, y1Init, y2Init)

	fmt.Println("Метод Адамса")
	for i := range x {
		fmt.Printf("x=%4.1f\ty1=%10.5f\ty2=%10.5f\n", x[i], y1[i], y2[i])
	}
}
