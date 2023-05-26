# PROYECTO FINAL ANALISIS NUMÃ‰RICO
# Gabriela Linares y Henry Velandia

import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
import io
import base64


def replaceDiagonal(matrix, replacementList):
    for i in range(len(replacementList)):
        matrix[i][i + 1] = replacementList[i]


def generar_matriz_cuadrada(j, betas):
    matriz = np.zeros((j, j))
    np.fill_diagonal(matriz, betas)
    listareemplazo = []

    for i in range(j - 1):
        listareemplazo.append(1)

    replaceDiagonal(matriz, listareemplazo)
    np.fill_diagonal(matriz[1:], 1)
    matriz[-1, -2] = 2

    return matriz


def jacobi(A, B, error):
    n = len(B)
    x = np.zeros(n)
    error_norm = float('inf')

    while error_norm > error:
        x_new = np.zeros(n)

        for i in range(n):
            sum_term = np.dot(A[i, :i], x[:i]) + np.dot(A[i, i+1:], x[i+1:])
            x_new[i] = (B[i] - sum_term) / A[i, i]

        error_norm = np.linalg.norm(x_new - x, np.inf)
        x = x_new

    return x


def sol_estacionaria(longitud, num_nodos_espacio, mconst, temp_contorno):
    delta_x = longitud / (num_nodos_espacio - 1)
    betas = -2 - mconst**2 * delta_x**2
    A = generar_matriz_cuadrada(num_nodos_espacio, betas)
    B = np.zeros((num_nodos_espacio, 1))
    B[0] = -temp_contorno
    x = jacobi(A, B, 0.0001)
    return x


def mostrar_grafica(x, longitud):
    plt.plot(np.linspace(0, longitud, len(x)), x)
    plt.xlabel('Posición')
    plt.ylabel('Temperatura')
    plt.title('Posición - Temperatura')

    img_bytes = io.BytesIO()
    plt.savefig(img_bytes, format='png')
    plt.close()

    img_bytes.seek(0)
    encoded_img = base64.b64encode(img_bytes.getvalue()).decode('utf-8')

    return encoded_img


layout_metodo = [
    [sg.Text("Seleccione el metodo de solucion:")],
    [sg.Radio("Estado estacionario", "metodo", key="metodo_1", default=True)],
    [sg.Radio("Estado transitorio", "metodo", key="metodo_2")],
    [sg.Button("Siguiente")],
]

ventana_metodo = sg.Window("Seleccion de metodo", layout_metodo, finalize=True)

metodo = None

while True:
    evento_metodo, valores_metodo = ventana_metodo.read()

    if evento_metodo == sg.WINDOW_CLOSED:
        break

    if evento_metodo == "Siguiente":
        if valores_metodo["metodo_1"]:
            metodo = 1
        elif valores_metodo["metodo_2"]:
            metodo = 2
        break

ventana_metodo.close()

if metodo == 1:
    layout_inputs = [
        [
            sg.Text("Longitud de la superficie de difusion:"),
            sg.InputText(key="longitud"),
        ],
        [sg.Text("Numero de nodos (espacio):"), sg.InputText(key="num_nodos_espacio")],
        [sg.Text("Constante de material:"), sg.InputText(key="mconst")],
        [sg.Text("Temperatura del contorno:"), sg.InputText(key="temp_contorno")],
    ]
elif metodo == 2:
    layout_inputs = [
        [
            sg.Text("Longitud de la superficie de difusion:"),
            sg.InputText(key="longitud"),
        ],
        [sg.Text("Tiempo total de simulacion:"), sg.InputText(key="t_total")],
        [sg.Text("Numero de nodos (espacio):"), sg.InputText(key="num_nodos_espacio")],
        [sg.Text("Numero de nodos (tiempo):"), sg.InputText(key="num_nodos_tiempo")],
        [sg.Text("Constante de material:"), sg.InputText(key="mconst")],
        [sg.Text("Temperatura del contorno:"), sg.InputText(key="temp_contorno")],
        [sg.Text("Temperatura inicial:"), sg.InputText(key="temp_inicial")],
    ]

layout_inputs.append([sg.Button("Resolver ecuacion")])
layout_inputs.append([sg.Button("Volver")])

ventana_inputs = sg.Window("Ingreso de inputs", layout_inputs, finalize=True)

while True:
    evento_inputs, valores_inputs = ventana_inputs.read()

    if evento_inputs == sg.WINDOW_CLOSED:
        break

    if evento_inputs == "Volver":
        ventana_inputs.close()
        ventana_metodo.un_hide()
        break

    if evento_inputs == "Resolver ecuacion":
        if metodo == 1:
            longitud = float(valores_inputs["longitud"])
            num_nodos_espacio = int(valores_inputs["num_nodos_espacio"])
            mconst = float(valores_inputs["mconst"])
            temp_contorno = float(valores_inputs["temp_contorno"])
            x = sol_estacionaria(longitud, num_nodos_espacio, mconst, temp_contorno)
            encoded_img = mostrar_grafica(x, longitud)
            layout_grafica = [[sg.Image(data=encoded_img)]]
            ventana_grafica = sg.Window("Gráfica", layout_grafica, finalize=True)

            while True:
                evento_grafica, valores_grafica = ventana_grafica.read()
                if evento_grafica == sg.WINDOW_CLOSED:
                    break
                elif metodo == 2:
                    pass

ventana_inputs.close()


# DATOS A USAR

# L = 1
# Tt = 2
# Nodos espacio = 10
# Nodos t = 5
# alpha = 0.01
# Temp contorno = 90
# Temp inicial = 20

# L/h = nodos
