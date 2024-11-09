package main

import (
	"fmt"
	"math"
)

// Функция для вычисления производных dy1/dx и dy2/dx
func derivatives(x, y1, y2 float64) (float64, float64) {
	dy1 := math.Cos(x + y2)
	dy2 := math.Sin(x - y2)
	return dy1, dy2
}

// Метод Эйлера
func eulerMethod(a, b, h float64, y1Init, y2Init float64) ([]float64, []float64, []float64) {
	// Количество шагов
	n := int((b - a) / h)

	// Массивы для хранения значений x, y1 и y2 на каждом шаге
	x := make([]float64, n+1)
	y1 := make([]float64, n+1)
	y2 := make([]float64, n+1)

	// Начальные условия
	x[0] = a
	y1[0] = y1Init
	y2[0] = y2Init

	// Основной цикл метода Эйлера
	for i := 0; i < n; i++ {
		// Вычисляем производные
		dy1, dy2 := derivatives(x[i], y1[i], y2[i])

		// Вычисляем значения на следующем шаге
		y1[i+1] = y1[i] + h*dy1
		y2[i+1] = y2[i] + h*dy2
		x[i+1] = x[i] + h
	}

	return x, y1, y2
}

// Метод Эйлера-Коши
func eulerCauchyMethod(a, b, h float64, y1Init, y2Init float64) ([]float64, []float64, []float64) {
	// Количество шагов
	n := int((b - a) / h)

	// Массивы для хранения значений x, y1 и y2 на каждом шаге
	x := make([]float64, n+1)
	y1 := make([]float64, n+1)
	y2 := make([]float64, n+1)

	// Начальные условия
	x[0] = a
	y1[0] = y1Init
	y2[0] = y2Init

	// Основной цикл метода Эйлера-Коши
	for i := 0; i < n; i++ {
		// Вычисляем производные в начальной точке
		dy1, dy2 := derivatives(x[i], y1[i], y2[i])

		// Предсказание (метод Эйлера)
		y1Pred := y1[i] + h*dy1
		y2Pred := y2[i] + h*dy2
		xPred := x[i] + h

		// Вычисляем производные в предсказанной точке
		dy1Pred, dy2Pred := derivatives(xPred, y1Pred, y2Pred)

		// Корректировка
		y1[i+1] = y1[i] + (h/2)*(dy1+dy1Pred)
		y2[i+1] = y2[i] + (h/2)*(dy2+dy2Pred)
		x[i+1] = xPred
	}

	return x, y1, y2
}

// Метод Рунге-Кутта 4 порядка
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

func main() {
	// Параметры задачи
	a := 0.0       // начальная точка
	b := 4.0       // конечная точка
	h := 0.1       // шаг
	y1Init := 0.7  // начальное значение y1
	y2Init := -0.5 // начальное значение y2

	// Решение задачи методом Эйлера
	x, y1, y2 := eulerMethod(a, b, h, y1Init, y2Init)

	// Вывод результатов
	for i := range x {
		fmt.Printf("x=%.1f\ty1=%.4f\ty2=%.4f\n", x[i], y1[i], y2[i])
	}

	y1Init = 0.7  // начальное значение y1
	y2Init = -0.5 // начальное значение y2
	fmt.Println("Метод Эйлера-Коши")
	// Решение задачи методом Эйлера-Коши
	x, y1, y2 = eulerCauchyMethod(a, b, h, y1Init, y2Init)

	// Вывод результатов
	for i := range x {
		fmt.Printf("x=%.1f\ty1=%.4f\ty2=%.4f\n", x[i], y1[i], y2[i])
	}

	fmt.Println("Метод Рунге-Кутта")
	// Решение задачи методом Рунге-Кутта 4-го порядка
	x, y1, y2 = rungeKutta4(a, b, h, y1Init, y2Init)

	// Вывод результатов
	for i := range x {
		fmt.Printf("x=%.1f\ty1=%.5f\ty2=%.5f\n", x[i], y1[i], y2[i])
	}

	//сделать красивый вывод
}
