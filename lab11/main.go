package main

import (
	"fmt"
	"math"
)

type Array []float64

// Spline3 вычисляет коэффициенты сплайна
func Spline3(N int, X, Y Array, S0, SN float64, A, B, C, D *Array) {
	// Временные массивы
	F := make(Array, N+1)
	h2 := X[1] - X[0]
	h3 := X[2] - X[1]

	// Заполняем массивы коэффициентов нулями
	*A = make(Array, N+1)
	*B = make(Array, N+1)
	*C = make(Array, N+1)
	*D = make(Array, N+1)

	// Заполнение значений для A и F на начальных точках
	(*A)[1] = 2 * (h2 + h3) / h3
	F[1] = (6/h3)*((Y[2]-Y[1])/h3-(Y[1]-Y[0])/h2) - (h2 * S0 / h3)

	// Заполнение значений для A, B и F для остальных точек
	for i := 3; i <= N-1; i++ {
		h2 = X[i-1] - X[i-2]
		h3 = X[i] - X[i-1]
		(*A)[i-2] = 2 * (h2 + h3) / h3
		(*B)[i-2] = h2 / h3
		F[i-2] = (6 / h3) * ((Y[i]-Y[i-1])/h3 - (Y[i-1]-Y[i-2])/h2)
	}

	// Заполнение значений для последних точек
	h2 = X[N-1] - X[N-2]
	h3 = X[N] - X[N-1]
	p := 2 * (h2 + h3)
	(*B)[1] = h2 / p
	F[N-2] = (6/p)*((Y[N]-Y[N-1])/h3-(Y[N-1]-Y[N-2])/h2) - (h3 * SN / p)

	// Прямой проход для вычисления d и c
	DTemp := make(Array, N+1)
	CTemp := make(Array, N+1)
	DTemp[1] = 1 / (*A)[1]
	CTemp[1] = F[1]

	for i := 2; i <= N-3; i++ {
		DTemp[i] = 1 / ((*A)[i] - (*B)[i]*DTemp[i-1])
		CTemp[i] = F[i] - (*B)[i]*DTemp[i-1]*CTemp[i-1]
	}
	DTemp[N-2] = (F[N-2] - (*B)[1]*DTemp[N-3]*CTemp[N-3]) / (1 - (*B)[1]*DTemp[N-3])

	// Обратный проход для вычисления значений d
	for i := N - 3; i >= 1; i-- {
		DTemp[i] *= (CTemp[i] - DTemp[i+1])
	}

	// Установка значений c
	(*C)[0] = S0
	(*C)[N] = SN
	for i := 1; i <= N-2; i++ {
		(*C)[i] = DTemp[i-1]
	}

	// Вычисление окончательных значений A, B, D для каждого отрезка
	for i := 1; i <= N; i++ {
		h2 = X[i] - X[i-1]
		(*D)[i] = ((*C)[i] - (*C)[i-1]) / h2
		(*B)[i] = h2*(*C)[i]/2 - h2*h2*(*D)[i]/6 + (Y[i]-Y[i-1])/h2
		(*A)[i] = Y[i]
	}
}

// DifSpline находит первую производную сплайна в точке x
func DifSpline(x float64, X Array, A, B, C, D Array) float64 {
	// Находим нужный интервал [X[i-1], X[i]]
	i := 1
	for i < len(X) && x > X[i] {
		i++
	}
	i-- // Переходим на индекс предыдущего узла

	// Длина интервала
	h := x - X[i]

	// Вычисление производной сплайна
	return B[i] + h*C[i] + h*h*D[i]/2
}

func main() {
	// Задаем функцию e^(x^2) и интерполируемые точки
	N := 5
	X := Array{0.0, 0.25, 0.5, 0.75, 1.0}
	Y := make(Array, len(X))
	for i := range X {
		Y[i] = math.Exp(X[i] * X[i])
	}

	// Задаем начальные и конечные производные
	S0 := 2 * X[0] * math.Exp(X[0]*X[0])       // Производная e^(x^2) при x = X[0]
	SN := 2 * X[N-1] * math.Exp(X[N-1]*X[N-1]) // Производная e^(x^2) при x = X[N-1]

	// Коэффициенты сплайна
	var A, B, C, D Array

	// Вычисляем коэффициенты сплайна
	Spline3(N, X, Y, S0, SN, &A, &B, &C, &D)

	// Находим производную функции e^(x^2) в точке x = 0.55
	x := 0.55
	derivative := DifSpline(x, X, A, B, C, D)

	fmt.Printf("Производная функции e^(x^2) в точке x=%.2f: %.6f\n", x, derivative)
}
