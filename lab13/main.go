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
// func computeIntegralForOmega(f func(float64) float64, a, b float64, n int) []complex128 {
// 	n = 16384
// 	step := (b - a) / float64(n) // шаг по x
// 	integrals := make([]complex128, n)

// 	for k := 0; k < n; k++ {
// 		omega := float64(k)
// 		sum := complex(0, 0)

// 		//расчет интеграла
// 		for i := 0; i < n; i++ {
// 			x := a + float64(i)*step
// 			term := complex(f(x), 0) * cmplx.Exp(1i*complex(omega*x, 0))
// 			sum += term
// 		}

//			integrals[k] = sum * complex(step, 0) // умножаем на шаг дискретизации
//		}
//		return integrals
//	}
// func computeIntegralForOmega(f func(float64) float64, a, b float64, n int) []complex128 {

// 	n = 1048576 * 16
// 	step := (b - a) / float64(n) // шаг по x
// 	fmt.Println("step=", step)
// 	integrals := make([]complex128, n)

// 	for k := 0; k < 512; k++ {

// 		sum := complex(0, 0)
// 		omega := float64(k)
// 		// Расчет интеграла с использованием формулы Симпсона
// 		for i := 0; i <= n; i++ {
// 			x := a + float64(i)*step
// 			weight := 1.0

// 			// Определяем вес для Симпсона
// 			if i == 0 || i == n {
// 				weight = 1.0
// 			} else if i%2 == 0 {
// 				weight = 2.0
// 			} else {
// 				weight = 4.0
// 			}

// 			term := complex(f(x), 0) * cmplx.Exp(1i*complex(omega*x, 0)) * complex(weight, 0)
// 			sum += term
// 		}

// 		integrals[k] = sum * complex(step/3.0, 0) // Умножаем на шаг по формуле Симпсона
// 	}
// 	return integrals
// }

func computeIntegral(f func(float64) float64, a, b float64, k int, epsilon float64, results chan<- IntegralResult) {
	var result complex128
	n := 4096 // Начальное количество точек (должно быть четным)
	step := (b - a) / float64(n)
	omega := float64(k)

	for {
		sum := complex(0, 0)
		for i := 0; i <= n; i++ {
			x := a + float64(i)*step
			weight := 1.0

			// Определяем вес для Симпсона
			if i == 0 || i == n {
				weight = 1.0
			} else if i%2 == 0 {
				weight = 2.0
			} else {
				weight = 4.0
			}

			term := complex(f(x), 0) * cmplx.Exp(1i*complex(omega*x, 0)) * complex(weight, 0)
			sum += term
		}

		newResult := sum * complex(step/3.0, 0) // Умножаем на шаг по формуле Симпсона

		// Проверка сходимости
		if cmplx.Abs(newResult-result) < epsilon {
			result = newResult
			break
		}

		// Увеличиваем количество точек
		result = newResult
		n *= 2
		step = (b - a) / float64(n)
	}

	results <- IntegralResult{k: k, integral: result}
}

type IntegralResult struct {
	k        int
	integral complex128
}

func main() {
	// границы интервала
	var a, b float64
	a = -math.Pi
	b = math.Pi
	const epsilon = 1e-15
	const n = 512 //количество точек
	fmt.Printf("Количество точек n: %d\n", n)

	//определение функции f(x)
	f := func(x float64) float64 {
		return math.Exp(-10 * (x - math.Pi) * (x - math.Pi))
		//return math.Sin(x*x) + x
		//return math.Exp(x * x)
	}

	//дискретные значения функции
	_, samples := generateSamples(f, a, b, n)

	fftResult := FFT(samples)

	//вывод спектра
	fmt.Println("Быстрое преобразование Фурье:")
	for k := 0; k < 10; k++ {
		omega := float64(k) //частота
		fmt.Printf("w_k = %.4f, F(w_k) = %.4f + %.4fi\n", omega, real(fftResult[k]), imag(fftResult[k]))
	}

	results := make(chan IntegralResult, n)

	// Запуск горутин
	for k := 0; k < n; k++ {
		go computeIntegral(f, a, b, k, epsilon, results)
	}

	// Сбор результатов
	integrals := make([]complex128, n)
	for i := 0; i < n; i++ {
		result := <-results
		integrals[result.k] = result.integral
	}

	// нтегралы для каждого w_k
	//integrals := computeIntegralForOmega(f, a, b, n)

	fmt.Println("\nЗначение интегралов w_k:")
	for k := 0; k < n; k++ {
		omega := float64(k)
		fmt.Printf("w_k = %5.0f : %12.5f + %12.5fi\n", omega, real(integrals[k]), imag(integrals[k]))
	}
}
