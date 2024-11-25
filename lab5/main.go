package main

import (
	"fmt"
	"math"
)

func Leverrier(n int, a [][]float64) []float64 {
	p := make([]float64, n)
	p[0] = Trace(a)
	b
}

// Matrix type as a 2D slice
type Matrix [][]float64

// CreateIdentity creates an identity matrix of size n x n
func CreateIdentity(n int) Matrix {
	I := make(Matrix, n)
	for i := range I {
		I[i] = make([]float64, n)
		I[i][i] = 1.0
	}
	return I
}

// Multiply multiplies two matrices A and B
func Multiply(A, B Matrix) Matrix {
	n := len(A)
	m := len(A[0])
	p := len(B[0])
	result := make(Matrix, n)
	for i := range result {
		result[i] = make([]float64, p)
	}

	for i := 0; i < n; i++ {
		for j := 0; j < p; j++ {
			for k := 0; k < m; k++ {
				result[i][j] += A[i][k] * B[k][j]
			}
		}
	}
	return result
}

// ScaleAndAdd scales matrix A by alpha and adds it to matrix B: result = alpha * A + B
func ScaleAndAdd(A Matrix, alpha float64, B Matrix) Matrix {
	n := len(A)
	m := len(A[0])
	result := make(Matrix, n)
	for i := range result {
		result[i] = make([]float64, m)
		for j := 0; j < m; j++ {
			result[i][j] = alpha*A[i][j] + B[i][j]
		}
	}
	return result
}

// Trace computes the trace of a square matrix A
func Trace(A Matrix) float64 {
	trace := 0.0
	for i := 0; i < len(A); i++ {
		trace += A[i][i]
	}
	return trace
}

// CopyMatrix creates a deep copy of a matrix
func CopyMatrix(A Matrix) Matrix {
	n := len(A)
	m := len(A[0])
	B := make(Matrix, n)
	for i := 0; i < n; i++ {
		B[i] = make([]float64, m)
		copy(B[i], A[i])
	}
	return B
}

// LeVerrier computes the coefficients of the characteristic polynomial of a square matrix
func LeVerrier(A Matrix) []float64 {
	n := len(A)
	I := CreateIdentity(n)
	Mk := CopyMatrix(A)
	coeffs := make([]float64, n+1)
	coeffs[n] = 1.0 // Leading coefficient is 1

	for k := 1; k <= n; k++ {
		trace := Trace(Mk)
		coeffs[n-k] = (-1.0 / float64(k)) * trace

		// Update Mk = A * (Mk + coeffs[n-k] * I)
		temp := ScaleAndAdd(I, coeffs[n-k], Mk)
		Mk = Multiply(A, temp)
	}

	return coeffs
}

func EvaluatePolynomial(coeffs []float64, x float64) float64 {
	result := 0.0
	n := len(coeffs) - 1
	for i := 0; i <= n; i++ {
		result += coeffs[i] * math.Pow(x, float64(n-i))
	}
	return result
}

// EvaluatePolynomialDerivative computes the value of the derivative P'(x) at x
func EvaluatePolynomialDerivative(coeffs []float64, x float64) float64 {
	result := 0.0
	n := len(coeffs) - 1
	for i := 0; i < n; i++ {
		result += float64(n-i) * coeffs[i] * math.Pow(x, float64(n-i-1))
	}
	return result
}

// NewtonRaphsonRootFinder finds a single root of the polynomial using the Newton-Raphson method
func NewtonRaphsonRootFinder(coeffs []float64, initialGuess float64, tolerance float64, maxIterations int) float64 {
	x := initialGuess
	for i := 0; i < maxIterations; i++ {
		Px := EvaluatePolynomial(coeffs, x)
		Pdx := EvaluatePolynomialDerivative(coeffs, x)

		if math.Abs(Pdx) < tolerance {
			// Avoid division by zero or near-zero derivative
			break
		}

		xNew := x - Px/Pdx
		if math.Abs(xNew-x) < tolerance {
			return xNew
		}
		x = xNew
	}
	return x
}

// FindRootsStub computes all roots of the polynomial using Newton-Raphson
func FindRootsStub(coeffs []float64) []float64 {
	n := len(coeffs) - 1
	roots := make([]float64, 0, n)
	tolerance := 1e-6
	maxIterations := 100

	// Use evenly spaced initial guesses around 0
	initialGuesses := make([]float64, n)
	for i := 0; i < n; i++ {
		initialGuesses[i] = float64(i) - float64(n)/2.0
	}

	for _, guess := range initialGuesses {
		root := NewtonRaphsonRootFinder(coeffs, guess, tolerance, maxIterations)

		// Check if the root is already found (to handle duplicates)
		isUnique := true
		for _, r := range roots {
			if math.Abs(r-root) < tolerance {
				isUnique = false
				break
			}
		}

		if isUnique {
			roots = append(roots, root)
		}
	}

	return roots
}

func main() {
	// Example matrix
	A := Matrix{
		{1, -1, -1, 2},
		{2, 3, 0, -4},
		{1, 1, -2, -2},
		{1, 1, 0, -1},
	}

	// Step 1: Compute characteristic polynomial
	coeffs := LeVerrier(A)
	fmt.Println("Characteristic polynomial coefficients:", coeffs)

	// Step 2: Compute eigenvalues (roots of the polynomial)
	eigenvalues := FindRootsStub(coeffs)
	fmt.Println("Eigenvalues:", eigenvalues)

	// Note: Eigenvector computation would follow using the eigenvalues
}
