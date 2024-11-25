package main

import (
	"fmt"
	"math"
	"math/cmplx"
)

// выполняет быстрое преобразование Фурье
func FFT(x []complex128) []complex128 {
	n := len(x)
	if n <= 1 {
		return x
	}

	//разделение на четные и нечетные элементы
	even := make([]complex128, n/2)
	odd := make([]complex128, n/2)
	for i := 0; i < n/2; i++ {
		even[i] = x[i*2]
		odd[i] = x[i*2+1]
	}

	//рекурсивные вызовы
	Feven := FFT(even)
	Fodd := FFT(odd)

	//объединение
	combined := make([]complex128, n)
	for k := 0; k < n/2; k++ {
		t := cmplx.Exp(complex(0, -2*math.Pi*float64(k)/float64(n))) * Fodd[k]
		combined[k] = Feven[k] + t
		combined[k+n/2] = Feven[k] - t
	}
	return combined
}

// создает дискретные значения функции f(x)
func generateSamples(f func(float64) float64, a, b float64, n int) ([]float64, []complex128) {
	samples := make([]complex128, n)
	xValues := make([]float64, n) //сохраняем координаты x
	step := (b - a) / float64(n)
	for i := 0; i < n; i++ {
		x := a + float64(i)*step
		xValues[i] = x
		samples[i] = complex(f(x), 0) //реальная часть f(x),мнимая — 0
	}
	return xValues, samples
}

// вычисляет интеграл для каждого значения омега
func computeIntegralForOmega(f func(float64) float64, a, b float64, n int) []complex128 {
	step := (b - a) / float64(n) // шаг по x
	integrals := make([]complex128, n)

	for k := 0; k < n; k++ {
		omega := 2 * math.Pi * float64(k) / (b - a)
		sum := complex(0, 0)

		//расчет интеграла
		for i := 0; i < n; i++ {
			x := a + float64(i)*step
			term := complex(f(x), 0) * cmplx.Exp(1i*complex(omega*x, 0))
			sum += term
		}

		integrals[k] = sum * complex(step, 0) // умножаем на шаг дискретизации
	}
	return integrals
}

func main() {
	// границы интервала
	var a, b float64
	a = -math.Pi
	b = math.Pi

	const n = 512 //количество точек
	fmt.Printf("Количество точек n: %d\n", n)

	//определение функции f(x)
	f := func(x float64) float64 {
		return math.Exp(-10 * (x - math.Pi) * (x - math.Pi))
		//return math.Sin(x*x) + x
	}

	//дискретные значения функции
	_, samples := generateSamples(f, a, b, n)

	fftResult := FFT(samples)

	//вывод спектра
	fmt.Println("Быстрое преобразование Фурье:")
	for k := 0; k < 10; k++ {
		omega := 2 * math.Pi * float64(k) / (b - a) //частота
		fmt.Printf("w_k = %.4f, F(w_k) = %.4f + %.4fi\n", omega, real(fftResult[k]), imag(fftResult[k]))
	}

	// нтегралы для каждого w_k
	integrals := computeIntegralForOmega(f, a, b, n)

	fmt.Println("\nЗначение интегралов w_k:")
	for k := 0; k < n; k++ {
		omega := 2 * math.Pi * float64(k) / (b - a)
		fmt.Printf("w_k = %5.0f : %10.4f + %10.4fi\n", omega, real(integrals[k]), imag(integrals[k]))
	}
}
