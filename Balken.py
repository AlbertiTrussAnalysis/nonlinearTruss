import copy
from tkinter import messagebox
import numpy as np
import math
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk

class Knotenn:

    def __init__(self, x, y, ID, h=0, v=0, Fx=0, Fy=0):
        self.x = x      #aktuelle-Koordinaten
        self.y = y      #aktuelle-Koordinaten
        self.ID = ID
        self.h = h
        self.v = v
        self.Fx = Fx
        self.Fy = Fy
        self.x0 = x     #Ursprüngliche Koordinaten
        self.y0 = y     #Ursprüngliche Koordinaten
        self.x_linear = x
        self.y_linear = y
        self.dof = [2 * self.ID, 2 * self.ID + 1]

    def auf_Startposition_linear(self):
        self.x_linear, self.y_linear = self.x0, self.y0

    def auf_Startposition(self):
        self.x, self.y = self.x0, self.y0

class Stab:

    def __init__(self, knoten1, knoten2, e_modul, flaechenquerschnitt, ID):
        self.ID = ID
        self.knoten1 = knoten1
        self.knoten2 = knoten2
        self.e_modul = e_modul
        self.flaechenquerschnitt = flaechenquerschnitt
        self.dof = [2 * self.knoten1.ID, 2 * self.knoten1.ID + 1,
                    2 * self.knoten2.ID, 2 * self.knoten2.ID + 1]
        self.L = math.sqrt((-self.knoten1.x0 + self.knoten2.x0) ** 2.0 + (-self.knoten1.y0 + self.knoten2.y0) ** 2.0)


    def F_int(self):
        u_1, v_1, u_2, v_2 = self.knoten1.x - self.knoten1.x0, self.knoten1.y - self.knoten1.y0, self.knoten2.x - self.knoten2.x0, self.knoten2.y - self.knoten2.y0
        x1_0, y1_0, x2_0, y2_0 = self.knoten1.x0, self.knoten1.y0, self.knoten2.x0, self.knoten2.y0


        res1 = self.e_modul * self.flaechenquerschnitt * self.L * ((-0.5*(-x1_0 + x2_0)**2 - 0.5*(-y1_0 + y2_0)**2 + 0.5*(-u_1 + u_2 - x1_0 + x2_0)**2 + 0.5*(-v_1 + v_2 - y1_0 + y2_0)**2)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)) * \
               np.array(
                   [[(1.0 * u_1 - 1.0 * u_2 + 1.0 * x1_0 - 1.0 * x2_0) / ((-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2)],
                    [(1.0 * v_1 - 1.0 * v_2 + 1.0 * y1_0 - 1.0 * y2_0) / ((-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2)],
                    [(-1.0 * u_1 + 1.0 * u_2 - 1.0 * x1_0 + 1.0 * x2_0) / ((-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2)],
                    [(-1.0 * v_1 + 1.0 * v_2 - 1.0 * y1_0 + 1.0 * y2_0) / ((-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2)]])

        return res1

    def Elementsteifigkeitsmatrix(self):
        u_1, v_1, u_2, v_2 = self.knoten1.x - self.knoten1.x0, self.knoten1.y - self.knoten1.y0, self.knoten2.x - self.knoten2.x0, self.knoten2.y - self.knoten2.y0
        x1_0, y1_0, x2_0, y2_0 = self.knoten1.x0, self.knoten1.y0, self.knoten2.x0, self.knoten2.y0


        res2 = self.e_modul * self.flaechenquerschnitt * self.L * \
               np.array([[(1.0*u_1 - 1.0*u_2 + 1.0*x1_0 - 1.0*x2_0)**2/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2 + 1.0*(-0.5*(-x1_0 + x2_0)**2 - 0.5*(-y1_0 + y2_0)**2 + 0.5*(-u_1 + u_2 - x1_0 + x2_0)**2 + 0.5*(-v_1 + v_2 - y1_0 + y2_0)**2)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (1.0*u_1 - 1.0*u_2 + 1.0*x1_0 - 1.0*x2_0)*(1.0*v_1 - 1.0*v_2 + 1.0*y1_0 - 1.0*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (-1.0*u_1 + 1.0*u_2 - 1.0*x1_0 + 1.0*x2_0)*(1.0*u_1 - 1.0*u_2 + 1.0*x1_0 - 1.0*x2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2 - 1.0*(-0.5*(-x1_0 + x2_0)**2 - 0.5*(-y1_0 + y2_0)**2 + 0.5*(-u_1 + u_2 - x1_0 + x2_0)**2 + 0.5*(-v_1 + v_2 - y1_0 + y2_0)**2)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (1.0*u_1 - 1.0*u_2 + 1.0*x1_0 - 1.0*x2_0)*(-1.0*v_1 + 1.0*v_2 - 1.0*y1_0 + 1.0*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2], [(1.0*u_1 - 1.0*u_2 + 1.0*x1_0 - 1.0*x2_0)*(1.0*v_1 - 1.0*v_2 + 1.0*y1_0 - 1.0*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (1.0*v_1 - 1.0*v_2 + 1.0*y1_0 - 1.0*y2_0)**2/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2 + 1.0*(-0.5*(-x1_0 + x2_0)**2 - 0.5*(-y1_0 + y2_0)**2 + 0.5*(-u_1 + u_2 - x1_0 + x2_0)**2 + 0.5*(-v_1 + v_2 - y1_0 + y2_0)**2)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (-1.0*u_1 + 1.0*u_2 - 1.0*x1_0 + 1.0*x2_0)*(1.0*v_1 - 1.0*v_2 + 1.0*y1_0 - 1.0*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (-1.0*v_1 + 1.0*v_2 - 1.0*y1_0 + 1.0*y2_0)*(1.0*v_1 - 1.0*v_2 + 1.0*y1_0 - 1.0*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2 - 1.0*(-0.5*(-x1_0 + x2_0)**2 - 0.5*(-y1_0 + y2_0)**2 + 0.5*(-u_1 + u_2 - x1_0 + x2_0)**2 + 0.5*(-v_1 + v_2 - y1_0 + y2_0)**2)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2], [(-1.0*u_1 + 1.0*u_2 - 1.0*x1_0 + 1.0*x2_0)*(1.0*u_1 - 1.0*u_2 + 1.0*x1_0 - 1.0*x2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2 - 1.0*(-0.5*(-x1_0 + x2_0)**2 - 0.5*(-y1_0 + y2_0)**2 + 0.5*(-u_1 + u_2 - x1_0 + x2_0)**2 + 0.5*(-v_1 + v_2 - y1_0 + y2_0)**2)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (-1.0*u_1 + 1.0*u_2 - 1.0*x1_0 + 1.0*x2_0)*(1.0*v_1 - 1.0*v_2 + 1.0*y1_0 - 1.0*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (-1.0*u_1 + 1.0*u_2 - 1.0*x1_0 + 1.0*x2_0)**2/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2 + 1.0*(-0.5*(-x1_0 + x2_0)**2 - 0.5*(-y1_0 + y2_0)**2 + 0.5*(-u_1 + u_2 - x1_0 + x2_0)**2 + 0.5*(-v_1 + v_2 - y1_0 + y2_0)**2)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (-1.0*u_1 + 1.0*u_2 - 1.0*x1_0 + 1.0*x2_0)*(-1.0*v_1 + 1.0*v_2 - 1.0*y1_0 + 1.0*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2], [(1.0*u_1 - 1.0*u_2 + 1.0*x1_0 - 1.0*x2_0)*(-1.0*v_1 + 1.0*v_2 - 1.0*y1_0 + 1.0*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (-1.0*v_1 + 1.0*v_2 - 1.0*y1_0 + 1.0*y2_0)*(1.0*v_1 - 1.0*v_2 + 1.0*y1_0 - 1.0*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2 - 1.0*(-0.5*(-x1_0 + x2_0)**2 - 0.5*(-y1_0 + y2_0)**2 + 0.5*(-u_1 + u_2 - x1_0 + x2_0)**2 + 0.5*(-v_1 + v_2 - y1_0 + y2_0)**2)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (-1.0*u_1 + 1.0*u_2 - 1.0*x1_0 + 1.0*x2_0)*(-1.0*v_1 + 1.0*v_2 - 1.0*y1_0 + 1.0*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2, (-1.0*v_1 + 1.0*v_2 - 1.0*y1_0 + 1.0*y2_0)**2/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2 + 1.0*(-0.5*(-x1_0 + x2_0)**2 - 0.5*(-y1_0 + y2_0)**2 + 0.5*(-u_1 + u_2 - x1_0 + x2_0)**2 + 0.5*(-v_1 + v_2 - y1_0 + y2_0)**2)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2)**2]])

        return res2

    def lineare_Elementsteifigkeitsmatrix(self):
        theta = np.arctan2(self.knoten2.y0 - self.knoten1.y0, self.knoten2.x0 - self.knoten1.x0)
        c = np.cos(theta)
        s = np.sin(theta)
        dr = (self.flaechenquerschnitt * self.e_modul / self.L) * np.array([[c ** 2, c * s, -c ** 2, -c * s],
                                                                                 [c * s, s ** 2, -s * c, -s ** 2],
                                                                                [-c ** 2, -s * c, c ** 2, s * c],
                                                                                 [-s * c, -s ** 2, s * c, s ** 2]])
        return dr

class Fachwerk:

    def __init__(self, Stabliste, Knotenliste):
        self.Stabliste = Stabliste
        self.Knotenliste = Knotenliste
        self.Norm_history = []

    def Knoten_linear_zuruecksetzten(self):
        for Knoten in self.Knotenliste:
            Knoten.auf_Startposition_linear()

    def Knoten_zuruecksetzten(self):
        for Knoten in self.Knotenliste:
            Knoten.auf_Startposition()

    def Knoten_aktualisieren(self, u):
        for Knoten in self.Knotenliste:
            dof = Knoten.dof
            Knoten.x += u[dof[0]]
            Knoten.y += u[dof[1]]

    def Knoten_aktualisieren_linear(self, u):
        for Knoten in self.Knotenliste:
            dof = Knoten.dof
            Knoten.x_linear += u[dof[0]]
            Knoten.y_linear += u[dof[1]]

    def Residuum(self):
        r_ges = np.zeros(2 * len(self.Knotenliste))

        k = 0
        for Knoten in self.Knotenliste:
            r_ges[k] = r_ges[k] - Knoten.Fx
            r_ges[k + 1] = r_ges[k + 1] - Knoten.Fy
            k = k + 2

        for Stab in self.Stabliste:
            r = Stab.F_int()
            for i, dof in enumerate(Stab.dof):
                r_ges[dof] += r[i]
            if Stab.knoten1.h == 1:
                r_ges[Stab.dof[0]] = 0
            if Stab.knoten1.v == 1:
                r_ges[Stab.dof[1]] = 0
            if Stab.knoten2.h == 1:
                r_ges[Stab.dof[2]] = 0
            if Stab.knoten2.v == 1:
                r_ges[Stab.dof[3]] = 0

        return r_ges

    def globale_Steifigkeitsmatrix(self):
        dr_ges = np.zeros((2 * len(self.Knotenliste), 2 * len(self.Knotenliste)))

        # Assemblierung der globalen Steifigkeitsmatrix
        for Stab in self.Stabliste:
            Steifigkeitsmatrix = Stab.Elementsteifigkeitsmatrix()
            for i, dof1 in enumerate(Stab.dof):
                for j, dof2 in enumerate(Stab.dof):
                    dr_ges[dof1, dof2] += Steifigkeitsmatrix[i, j]

        # Anwendung der Randbedingungen
        for Stab in self.Stabliste:
            if Stab.knoten1.h == 1:
                dr_ges[Stab.dof[0], :] = 0
                dr_ges[:, Stab.dof[0]] = 0
                dr_ges[Stab.dof[0], Stab.dof[0]] = 1
            if Stab.knoten1.v == 1:
                dr_ges[Stab.dof[1], :] = 0
                dr_ges[:, Stab.dof[1]] = 0
                dr_ges[Stab.dof[1], Stab.dof[1]] = 1
            if Stab.knoten2.h == 1:
                dr_ges[Stab.dof[2], :] = 0
                dr_ges[:, Stab.dof[2]] = 0
                dr_ges[Stab.dof[2], Stab.dof[2]] = 1
            if Stab.knoten2.v == 1:
                dr_ges[Stab.dof[3], :] = 0
                dr_ges[:, Stab.dof[3]] = 0
                dr_ges[Stab.dof[3], Stab.dof[3]] = 1

        return dr_ges

    def F_bestimmen(self):
        return (self.Residuum()[3])/2


    def linear_F(self):
        F_ges = np.zeros(2 * len(self.Knotenliste))
        for Knoten in self.Knotenliste:
            F_ges[Knoten.ID * 2] += Knoten.Fx
            F_ges[Knoten.ID * 2 + 1] += Knoten.Fy
            if Knoten.h == 1:
                F_ges[Knoten.ID * 2] = 0
            if Knoten.v == 1:
                F_ges[Knoten.ID * 2 + 1] = 0
        return F_ges

    def linear_globale_Steifigkeitsmatrix(self):
        dr_ges = np.zeros((2 * len(self.Knotenliste), 2 * len(self.Knotenliste)))

        # Assemblierung der globalen Steifigkeitsmatrix
        for Stab in self.Stabliste:
            Steifigkeitsmatrix = Stab.lineare_Elementsteifigkeitsmatrix()
            for i, dof1 in enumerate(Stab.dof):
                for j, dof2 in enumerate(Stab.dof):
                    dr_ges[dof1, dof2] += Steifigkeitsmatrix[i, j]

        # Anwendung der Randbedingungen
        for Stab in self.Stabliste:
            if Stab.knoten1.h == 1:
                dr_ges[Stab.dof[0], :] = 0
                dr_ges[:, Stab.dof[0]] = 0
                dr_ges[Stab.dof[0], Stab.dof[0]] = 1
            if Stab.knoten1.v == 1:
                dr_ges[Stab.dof[1], :] = 0
                dr_ges[:, Stab.dof[1]] = 0
                dr_ges[Stab.dof[1], Stab.dof[1]] = 1
            if Stab.knoten2.h == 1:
                dr_ges[Stab.dof[2], :] = 0
                dr_ges[:, Stab.dof[2]] = 0
                dr_ges[Stab.dof[2], Stab.dof[2]] = 1
            if Stab.knoten2.v == 1:
                dr_ges[Stab.dof[3], :] = 0
                dr_ges[:, Stab.dof[3]] = 0
                dr_ges[Stab.dof[3], Stab.dof[3]] = 1
        return dr_ges

    def loesen(self):

        self.Norm_history=[]

        dr = self.linear_globale_Steifigkeitsmatrix()
        F = self.linear_F()

        try:
            u_linear = np.linalg.solve(dr, F)
        except np.linalg.LinAlgError:
            messagebox.showerror("Fehler", "Das System ist kinematisch!")
            return

        self.Knoten_aktualisieren_linear(u_linear)

        while True:
            r = self.Residuum()
            dr = self.globale_Steifigkeitsmatrix()
            u = np.linalg.solve(dr, -r)

            self.Knoten_aktualisieren(u)

            res = np.linalg.norm(r)

            self.Norm_history.append(res)

            if abs(res) < 1e-3:
                break

    def plot_norm(self):

        iterations = np.arange(len(self.Norm_history))

        # Filtere die Daten, sodass nur Werte kleiner als 1 berücksichtigt werden
        mask = np.array(self.Norm_history) < 1
        filtered_iterations = iterations[mask]
        filtered_log_errors = np.log(self.Norm_history)[mask]

        # Berechne die lineare Regression nur für die gefilterten Werte
        if len(filtered_iterations) > 0 and len(filtered_log_errors) > 0:
            slope, intercept = np.polyfit(filtered_iterations, filtered_log_errors, 1)
            plt.plot(self.Norm_history, marker='o')
            plt.yscale('log')
            plt.title('Konvergenzverlauf')
            plt.xlabel(r'$k$')
            plt.ylabel(r'$ \|\| \mathbf{r}^k \|\|_2 $')

            # Zeichne die Gerade nur für die gefilterten Werte
            plt.plot(filtered_iterations, np.exp(intercept + slope * filtered_iterations), 'r--',
                     label=f'Steigung für r < 1  = {slope:.2f}')

            plt.legend()
            plt.grid(True)
            plt.show()
        else:
            print("Warning: filtered_iterations or filtered_log_errors is empty. Cannot fit a line.")

K1 = Knotenn(0, 4, ID=0, v=1, h=1)
K2 = Knotenn(3, 0, ID=1)
K3 = Knotenn(6, 4, ID=2, v=1, h=1)
Knotenliste = [K1, K2, K3]
S1 = Stab(K1, K2, 20, 4, 0)
S2 = Stab(K2, K3, 20, 4, 1)
Stabliste = [S1, S2]
test = Fachwerk(Stabliste, Knotenliste)

Fy_werte = np.linspace(0, 20, 400)
Fy_ergebnisse = []
u1_werte = []
u2_werte = []
u3_werte = np.linspace(0, 10, 200)

for Fyy in Fy_werte:
    K2.Fy = 2 * Fyy
    test.Knoten_linear_zuruecksetzten()
    test.loesen()
    u1_werte.append(K2.y)
    u2_werte.append(K2.y_linear)
    if abs(Fyy - 4.96240506329114) < 1e-6:
        test.plot_norm()


test.Knoten_zuruecksetzten()
K2.Fy = 0

for u in u3_werte:
    K2.y = u
    Fy_ergebnisse.append(test.F_bestimmen())
    test.Knoten_zuruecksetzten()

u_values = np.linspace(0, 12, 400)
F_equation_values = (8 * (u_values - 8) * (u_values - 4) * u_values) / 25

plt.figure(figsize=(10, 10))
plt.plot(u_values, F_equation_values, label="Analytische Lösung", color='yellow')
plt.plot(u1_werte, Fy_werte, label="Lösung der Eigenentwicklung (Nichtlinear, Kraftgesteuert)", linestyle='--', color='blue')
plt.plot(u2_werte, Fy_werte, label="Lineare Berechnung", linestyle='--', color='C2')
plt.plot(u3_werte, Fy_ergebnisse, label="Lösung der Eigenentwicklung (Nichtlinear, Verschiebungsgesteuert)", linestyle='--', color='red')
plt.title("Kraft-Verschiebungsverlauf")
plt.xlim(0, 10)
plt.ylim(-10, 20)
plt.xlabel('u')
plt.ylabel('F')
plt.legend()
plt.grid(True)
plt.show()
