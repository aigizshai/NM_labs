import numpy as np

# данные для примера
xi = np.array([0.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0])
yi = np.array([-1.92, -1.60, -1.57, -1.41, -1.36, -0.97, -0.59, -0.71, -0.15, 0.01, 0.22, 0.63, 1.07, 1.42, 1.68, 2.49, 2.57, 3.09, 3.40, 4.0])


def get_coefficients(_pl: int, _xi: np.ndarray):
    '''
    Определение коэффициентов для множителей базисных полиномов l_i
    :param _pl: индекс базисного полинома
    :param _xi: массив значений x
    :return:
    '''
    n = int(_xi.shape[0])
    coefficients = np.empty((n, 2), dtype=float)
    for i in range(n):
        if i == _pl:
            coefficients[i][0] = float('inf')
            coefficients[i][1] = float('inf')
        else:
            coefficients[i][0] = 1 / (_xi[_pl] - _xi[i])
            coefficients[i][1] = -_xi[i] / (_xi[_pl] - _xi[i])
    filtered_coefficients = np.empty((n - 1, 2), dtype=float)
    j = 0
    for i in range(n):
        if coefficients[i][0] != float('inf'):
            # изменение последовательности, степень увеличивается
            filtered_coefficients[j][0] = coefficients[i][1]
            filtered_coefficients[j][1] = coefficients[i][0]
            j += 1
    return filtered_coefficients


def get_polynomial_l(_xi: np.ndarray):
    '''
    Определение базисных полиномов
    :param _xi: массив значений x
    :return:
    '''
    n = int(_xi.shape[0])
    pli = np.zeros((n, n), dtype=float)
    for pl in range(n):
        coefficients = get_coefficients(pl, _xi)
        for i in range(1, n - 1):  # проходим по массиву coefficients
            if i == 1:  # на второй итерации занимаются 0-2 степени
                pli[pl][0] = coefficients[i - 1][0] * coefficients[i][0]
                pli[pl][1] = coefficients[i - 1][1] * coefficients[i][0] + coefficients[i][1] * coefficients[i - 1][0]
                pli[pl][2] = coefficients[i - 1][1] * coefficients[i][1]
            else:
                clone_pli = np.zeros(n, dtype=float)
                for val in range(n):
                    clone_pli[val] = pli[pl][val]
                zeros_pli = np.zeros(n, dtype=float)
                for j in range(n-1):  # проходим по строке pl массива pli
                    product_1 = clone_pli[j] * coefficients[i][0]
                    product_2 = clone_pli[j] * coefficients[i][1]
                    zeros_pli[j] += product_1
                    zeros_pli[j+1] += product_2
                for val in range(n):
                    pli[pl][val] = zeros_pli[val]
    return pli

def get_polynomial(_xi: np.ndarray, _yi: np.ndarray):
    '''
    Определение интерполяционного многочлена Лагранжа
    :param _xi: массив значений x
    :param _yi: массив значений y
    :return:
    '''
    n = int(_xi.shape[0])
    polynomial_l = get_polynomial_l(_xi)
    for i in range(n):
        for j in range(n):
            polynomial_l[i][j] *= _yi[i]
    L = np.zeros(n, dtype=float)
    for i in range(n):
        for j in range(n):
            L[i] += polynomial_l[j][i]
    return L


print(get_polynomial(xi,yi))