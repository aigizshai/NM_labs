package main

import (
	"fmt"
	n "lab4/Newton"
	si "lab4/sim_iter"
)

func main() {
	x, norma, iter := si.SimpleIter()
	fmt.Println("Вычисление методом простых итераций")
	fmt.Printf("Корни: %.4f\n", x)
	fmt.Printf("Норма: %.5f\n", norma)
	fmt.Println("Количество итераций", iter)

	fmt.Println("")
	x_newton, iter_newton := n.Newton()
	fmt.Println("Вычисление методом Ньютона")
	fmt.Printf("Корни: %.4f\n", x_newton)
	fmt.Println("Количество итераций", iter_newton)

}
