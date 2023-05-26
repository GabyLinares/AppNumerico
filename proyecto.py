# PROYECTO FINAL ANALISIS NUMÃ‰RICO
# Gabriela Linares y Henry Velandia

import tkinter as tk
from tkinter import Tk, Frame, Label, Entry, Radiobutton, Button
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np

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

def mostrar_grafica(x, longitud, window):
    fig = plt.figure()
    plt.plot(np.linspace(0, longitud, len(x)), x)
    plt.xlabel('Posición')
    plt.ylabel('Temperatura')
    plt.title('Posición - Temperatura')

    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().pack()
    plt.show()

def mostrar_ventana_metodo():
    window_metodo = Tk()
    window_metodo.title("Seleccion de metodo")

    metodo = tk.IntVar()
    metodo.set(1)

    metodo_frame = Frame(window_metodo)
    metodo_frame.pack(pady=10)

    radio_1 = Radiobutton(metodo_frame, text="Estado estacionario", variable=metodo, value=1)
    radio_1.pack()

    radio_2 = Radiobutton(metodo_frame, text="Estado transitorio", variable=metodo, value=2)
    radio_2.pack()

    button_siguiente = Button(window_metodo, text="Siguiente", command=lambda: mostrar_ventana_inputs(metodo.get(), window_metodo))
    button_siguiente.pack(pady=10)

    window_metodo.mainloop()

def mostrar_ventana_inputs(metodo, window_metodo):
    window_metodo.destroy()

    window_inputs = Tk()
    window_inputs.title("Ingreso de inputs")

    inputs_frame = Frame(window_inputs)
    inputs_frame.pack(pady=10)

    if metodo == 1:
        longitud_label = Label(inputs_frame, text="Longitud de la superficie de difusion:")
        longitud_label.pack()

        longitud_entry = Entry(inputs_frame)
        longitud_entry.pack()

        num_nodos_espacio_label = Label(inputs_frame, text="Numero de nodos (espacio):")
        num_nodos_espacio_label.pack()

        num_nodos_espacio_entry = Entry(inputs_frame)
        num_nodos_espacio_entry.pack()

        mconst_label = Label(inputs_frame, text="Constante de material:")
        mconst_label.pack()

        mconst_entry = Entry(inputs_frame)
        mconst_entry.pack()

        temp_contorno_label = Label(inputs_frame, text="Temperatura del contorno:")
        temp_contorno_label.pack()

        temp_contorno_entry = Entry(inputs_frame)
        temp_contorno_entry.pack()

        button_resolver = Button(window_inputs, text="Resolver ecuacion", command=lambda: resolver_ecuacion_estacionaria(float(longitud_entry.get()), int(num_nodos_espacio_entry.get()), float(mconst_entry.get()), float(temp_contorno_entry.get()), window_inputs))
        button_resolver.pack(pady=10)
    elif metodo == 2:
        longitud_label = Label(inputs_frame, text="Longitud de la superficie de difusion:")
        longitud_label.pack()

        longitud_entry = Entry(inputs_frame)
        longitud_entry.pack()

        t_total_label = Label(inputs_frame, text="Tiempo total de simulacion:")
        t_total_label.pack()

        t_total_entry = Entry(inputs_frame)
        t_total_entry.pack()

        num_nodos_espacio_label = Label(inputs_frame, text="Numero de nodos (espacio):")
        num_nodos_espacio_label.pack()

        num_nodos_espacio_entry = Entry(inputs_frame)
        num_nodos_espacio_entry.pack()

        num_nodos_tiempo_label = Label(inputs_frame, text="Numero de nodos (tiempo):")
        num_nodos_tiempo_label.pack()

        num_nodos_tiempo_entry = Entry(inputs_frame)
        num_nodos_tiempo_entry.pack()

        mconst_label = Label(inputs_frame, text="Constante de material:")
        mconst_label.pack()

        mconst_entry = Entry(inputs_frame)
        mconst_entry.pack()

        temp_contorno_label = Label(inputs_frame, text="Temperatura del contorno:")
        temp_contorno_label.pack()

        temp_contorno_entry = Entry(inputs_frame)
        temp_contorno_entry.pack()

        temp_inicial_label = Label(inputs_frame, text="Temperatura inicial:")
        temp_inicial_label.pack()

        temp_inicial_entry = Entry(inputs_frame)
        temp_inicial_entry.pack()

        button_resolver = Button(window_inputs, text="Resolver ecuacion", command=lambda: resolver_ecuacion_transitoria(float(longitud_entry.get()), float(t_total_entry.get()), int(num_nodos_espacio_entry.get()), int(num_nodos_tiempo_entry.get()), float(mconst_entry.get()), float(temp_contorno_entry.get()), float(temp_inicial_entry.get()), window_inputs))
        button_resolver.pack(pady=10)

    button_volver = Button(window_inputs, text="Volver", command=lambda: volver(window_metodo, window_inputs))
    button_volver.pack()

    window_inputs.mainloop()

def resolver_ecuacion_estacionaria(longitud, num_nodos_espacio, mconst, temp_contorno, window_inputs):
    x = sol_estacionaria(longitud, num_nodos_espacio, mconst, temp_contorno)
    mostrar_grafica(x, longitud, window_inputs)

def resolver_ecuacion_transitoria(longitud, t_total, num_nodos_espacio, num_nodos_tiempo, mconst, temp_contorno, temp_inicial, window_inputs):
    # Implementar la resolución de la ecuación para el estado transitorio
    pass

def volver(window_metodo, window_inputs):
    window_inputs.destroy()
    mostrar_ventana_metodo()

mostrar_ventana_metodo()

# DATOS A USAR

# L = 1
# Tt = 2
# Nodos espacio = 10
# Nodos t = 5
# alpha = 0.01
# Temp contorno = 90
# Temp inicial = 20

# L/h = nodos
