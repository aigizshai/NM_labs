package main

import (
	"fmt"
	"math"
	"math/cmplx"
)

// Функция для выполнения быстрого преобразования Фурье (FFT)
func FFT(x []complex128, inverse bool) ([]complex128, []complex128) {
	n := len(x)
	if n <= 1 {
		return x, nil
	}

	// Разбиение на четные и нечетные элементы
	even, _ := FFT(x[0:n/2], inverse)
	odd, _ := FFT(x[n/2:], inverse)

	// Инициализация результата
	result := make([]complex128, n)
	angle := 2 * math.Pi / float64(n)
	if inverse {
		angle = -angle
	}

	omega := make([]complex128, n)
	// Преобразование Фурье
	for k := 0; k < n/2; k++ {
		// Расчет комплексных коэффициентов
		omega[k] = cmplx.Exp(complex(0, float64(k)*angle)) // e^(2πik/n)
		result[k] = even[k] + omega[k]*odd[k]
		result[k+n/2] = even[k] - omega[k]*odd[k]
	}

	return result, omega
}

// Функция для выполнения обратного преобразования Фурье (IFFT)
func IFFT(x []complex128) ([]complex128, []complex128) {
	// В обратном преобразовании используем инвертированное значение (inverse = true)
	return FFT(x, true)
}

// Функция для демонстрации работы
func main() {
	// Пример использования с функцией f(x) = e^(-x^2)
	n := 512
	x := make([]complex128, n)

	// Инициализация значений функции f(x)
	for i := 0; i < n; i++ {
		x[i] = complex(math.Exp(-float64(i)*float64(i)), 0)
	}

	// Выполнение FFT
	//fmt.Println("FFT:")
	i := 1
	fftResult, omega := FFT(x, false)
	for _, val := range fftResult {
		fmt.Printf("n=%d omega=%.5f val=%.5f\n", i, omega[i-1], val)
		i++
	}

	// Выполнение IFFT (обратное преобразование)
	fmt.Println("\nIFFT:")
	// ifftResult, _ := IFFT(fftResult)
	// i = 1
	//
	//	for _, val := range ifftResult {
	//		fmt.Printf("n=%d val=%.5f\n", i, val)
	//		i++
	//	}
}
