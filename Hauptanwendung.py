import copy
from tkinter import messagebox
import numpy as np
import math
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk

class Knotenn:

    def __init__(self, x, y, ID, h=0, v=0, Fx=0, Fy=0):
        self.x = x     #Ziel-Koordinaten
        self.y = y     #Ziel-Koordinaten
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

    def auf_Startposition(self):
        self.x_linear, self.y_linear = self.x0, self.y0

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

    def plot_structure(self):

        for Knoten in self.Knotenliste:
            plt.plot(Knoten.x, Knoten.y, "go")

        for s in self.Stabliste:
            plt.plot([s.knoten1.x, s.knoten2.x], [s.knoten1.y, s.knoten2.y], "green")

        a = 8
        h = a * np.sqrt(3) / 2

        for Knoten_instanz in self.Knotenliste:
            auflager_punkt = np.array([Knoten_instanz.x, Knoten_instanz.y])
            if Knoten_instanz.h == 1:
                punkt1 = auflager_punkt
                punkt2 = auflager_punkt + np.array([-h, -a / 2])
                punkt3 = auflager_punkt + np.array([-h, a / 2])
                plt.plot([punkt1[0], punkt2[0]], [punkt1[1], punkt2[1]], "green")
                plt.plot([punkt1[0], punkt3[0]], [punkt1[1], punkt3[1]], "green")
                plt.plot([punkt2[0], punkt3[0]], [punkt2[1], punkt3[1]], "green")
            if Knoten_instanz.v == 1:
                punkt1 = auflager_punkt
                punkt2 = auflager_punkt + np.array([-a / 2, -h])
                punkt3 = auflager_punkt + np.array([a / 2, -h])
                plt.plot([punkt1[0], punkt2[0]], [punkt1[1], punkt2[1]], "green")
                plt.plot([punkt1[0], punkt3[0]], [punkt1[1], punkt3[1]], "green")
                plt.plot([punkt2[0], punkt3[0]], [punkt2[1], punkt3[1]], "green")

        plt.show()


    def loesen(self):

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


stabListe = []
KnotenListe = []

ersterKnoten = None
zweiterKnoten = None
ID_counter = 0
stabID_counter = 0
letzterStab = None


class Fenster:

    def __init__(self, Fenster):
        self.Fenster = Fenster
        Fenster.title("geometrisch Nichtlineare Fachwerke")

        self.Oberfläche = tk.Canvas(Fenster, width=1400, height=800)
        self.Oberfläche.pack()

        self.Oberfläche.bind("<Button-1>", self.Knoten_hinzufügen)
        self.Oberfläche.bind("<Button-2>", self.Knoten_bearbeiten)
        self.Oberfläche.bind("<Button-3>", self.Stab_hinzufuegen)

        # self.button2 = tk.Button(Fenster, text="Visualisieren", command=self.plot_results)
        # self.button2.pack(side=tk.BOTTOM, anchor=tk.S)

        # self.button = tk.Button(Fenster, text="Neustart", command=self.restart_program)
        # self.button.pack(side=tk.BOTTOM, anchor=tk.SE)

        self.plot_norm_var = tk.BooleanVar(value=False)
        self.plot_geo_var = tk.BooleanVar(value=False)
        self.button = tk.Button(Fenster, text="Berechnen", command=self.calculate)
        self.checkbox = tk.Checkbutton(Fenster, text="Konvergenz anzeigen", variable=self.plot_norm_var)
        self.checkbox.pack(side=tk.BOTTOM, anchor=tk.S)
        self.button.pack(side=tk.BOTTOM, anchor=tk.S)
        self.checkbox2 = tk.Checkbutton(Fenster, text="Geometrieverlauf anzeigen", variable=self.plot_geo_var)
        self.checkbox2.pack(side=tk.BOTTOM, anchor=tk.S)

        raster_größe = 40

        for x in range(0, 1400, raster_größe):
            self.Oberfläche.create_line(x, 0, x, 800, fill="gray")

        for y in range(0, 800, raster_größe):
            self.Oberfläche.create_line(0, y, 1400, y, fill="gray")

        pfeil_y = self.Oberfläche.create_line(20, 760, 20, 720, arrow=tk.LAST, fill="black")
        pfeil_x = self.Oberfläche.create_line(20, 760, 60, 760, arrow=tk.LAST, fill="black")
        self.Oberfläche.create_text(10, 760, text="y", fill="black")
        self.Oberfläche.create_text(60, 770, text="x", fill="black")


        self.Oberfläche.create_text(180, 90, text="Linke Maustaste erstellt Knoten", font=("Arial", 12, "bold"),
                                    fill="black")
        self.Oberfläche.create_text(180, 110, text="Mittlere Maustaste bearbeitet Knoten", font=("Arial", 12, "bold"),
                                    fill="black")
        self.Oberfläche.create_text(180, 130, text="Rechte Maustaste auf zwei Knoten erstellt Stab",
                                    font=("Arial", 12, "bold"), fill="black")

    def Knoten_hinzufügen(self, Event):
        x = Event.x
        y = Event.y

        # Pop-up-Fenster für die Eingabe der Koordinaten
        coordFenster = tk.Toplevel()
        coordFenster.title("Knoten-Koordinaten eingeben")

        x_Schriftzug = tk.Label(coordFenster, text="X-Koordinate:")
        x_Schriftzug.grid(row=0, column=0)
        x_Eingabe = tk.Entry(coordFenster)
        x_Eingabe.insert(0, str(x))  # Anfangswert ist der visuell ausgewählte Punkt
        x_Eingabe.grid(row=0, column=1)

        y_Schriftzug = tk.Label(coordFenster, text="Y-Koordinate:")
        y_Schriftzug.grid(row=1, column=0)
        y_Eingabe = tk.Entry(coordFenster)
        y_Eingabe.insert(0, str(800 - y))  # Anfangswert ist der visuell ausgewählte Punkt
        y_Eingabe.grid(row=1, column=1)

        def hinzufuegen():
            # Benutzen Sie die eingegebenen Koordinaten, um einen Knoten hinzuzufügen
            x_knoten = int(x_Eingabe.get())
            y_knoten = 800 - int(y_Eingabe.get())

            Punkt = self.Oberfläche.create_oval(x_knoten - 5, y_knoten - 5, x_knoten + 5, y_knoten + 5, outline="black",
                                                width=1, fill="white")
            global ID_counter
            knoten = Knotenn(x_knoten, y_knoten, ID_counter)
            knoten.oval = Punkt
            ID_counter = ID_counter + 1
            KnotenListe.append(knoten)
            coordFenster.destroy()

        Hinzufüge_Knopf = tk.Button(coordFenster, text="Hinzufügen", command=hinzufuegen)
        Hinzufüge_Knopf.grid(row=2, column=0, columnspan=2)
        Hinzufüge_Knopf.bind("<Return>", hinzufuegen)
        Hinzufüge_Knopf.focus_set()

    def Knoten_abfrage(self):
        return [(knoten.x, 800 - knoten.y, knoten.ID, knoten.h, knoten.v, knoten.Fx, knoten.Fy) for knoten in
                KnotenListe]

    def Stab_abfrage(self):
        return [(Stab.knoten1.ID, Stab.knoten2.ID, Stab.e_modul, Stab.flaechenquerschnitt, Stab.ID) for Stab in
                stabListe]

    def Knoten_bearbeiten(self, Event):
        x = Event.x
        y = Event.y

        KnotenID = None
        geringste_Entfernung = float("inf")

        for knoten in KnotenListe:
            gefundene_Entfernung = math.sqrt((x - knoten.x) ** 2 + (y - knoten.y) ** 2)
            if gefundene_Entfernung < geringste_Entfernung:
                geringste_Entfernung = gefundene_Entfernung
                KnotenID = knoten.ID

        Extrafenster = tk.Toplevel()
        Extrafenster.title("Knoten bearbeiten")

        Betrag_Schriftzug = tk.Label(Extrafenster, text="Betrag in [kN]:")
        Betrag_Schriftzug.grid(row=0, column=0)
        Betrag_Eingabe = tk.Entry(Extrafenster)
        Betrag_Eingabe.grid(row=0, column=1)

        Richtung_Schriftzug = tk.Label(Extrafenster, text="Richtung in [°]:")
        Richtung_Schriftzug.grid(row=1, column=0)
        Richtung_Eingabe = tk.Entry(Extrafenster)
        Richtung_Eingabe.grid(row=1, column=1)

        Auflager_Schriftzug = tk.Label(Extrafenster, text="Auflager:")
        Auflager_Schriftzug.grid(row=2, column=0)
        horizontal_var = tk.IntVar()
        vertikal_var = tk.IntVar()
        horizontal_Auswahl = ttk.Checkbutton(Extrafenster, text="horizontal", variable=horizontal_var)
        horizontal_Auswahl.grid(row=2, column=1)
        vertikal_Auswahl = ttk.Checkbutton(Extrafenster, text="vertikal", variable=vertikal_var)
        vertikal_Auswahl.grid(row=2, column=2)

        bearbeiteter_Knoten = None
        for k in KnotenListe:
            if k.ID == KnotenID:
                bearbeiteter_Knoten = k
                break

        def hinzufügen():
            Betrag = float(Betrag_Eingabe.get() or 0)
            Richtung = math.radians(float(Richtung_Eingabe.get() or 0))

            if bearbeiteter_Knoten:
                self.Oberfläche.delete(bearbeiteter_Knoten.oval)

                bearbeiteter_Knoten.oval = self.Oberfläche.create_oval(bearbeiteter_Knoten.x - 5,
                                                                       bearbeiteter_Knoten.y - 5,
                                                                       bearbeiteter_Knoten.x + 5,
                                                                       bearbeiteter_Knoten.y + 5,
                                                                       outline="black", width=1, fill="white")
                Pfeil = self.Oberfläche.create_line(bearbeiteter_Knoten.x, bearbeiteter_Knoten.y,
                                                    bearbeiteter_Knoten.x + Betrag * math.cos(Richtung),
                                                    bearbeiteter_Knoten.y - Betrag * math.sin(Richtung),
                                                    arrow=tk.LAST,
                                                    fill="red")
                bearbeiteter_Knoten.Fy = Betrag * math.sin(Richtung)
                bearbeiteter_Knoten.Fx = Betrag * math.cos(Richtung)

            if vertikal_var.get() == 1:
                Auflager = self.Oberfläche.create_polygon(
                    [bearbeiteter_Knoten.x - 5, bearbeiteter_Knoten.y + 10, bearbeiteter_Knoten.x,
                     bearbeiteter_Knoten.y, bearbeiteter_Knoten.x + 5, bearbeiteter_Knoten.y + 10], fill='white',
                    outline='black', width=2.1)
                bearbeiteter_Knoten.v = 1

            if horizontal_var.get() == 1:
                Auflager = self.Oberfläche.create_polygon(
                    [bearbeiteter_Knoten.x - 10, bearbeiteter_Knoten.y - 5, bearbeiteter_Knoten.x,
                     bearbeiteter_Knoten.y, bearbeiteter_Knoten.x - 10, bearbeiteter_Knoten.y + 5], fill='white',
                    outline='black', width=2.1)
                bearbeiteter_Knoten.h = 1

            Extrafenster.destroy()

        Hinzufüge_Knopf = tk.Button(Extrafenster, text="Hinzufügen", command=hinzufügen)
        Hinzufüge_Knopf.grid(row=5, column=0, columnspan=4)

    def Stab_hinzufuegen(self, Event):
        global ID_counter
        global stabID_counter
        global ersterKnoten
        global zweiterKnoten
        global letzterStab

        x = Event.x
        y = Event.y

        KnotenID = None
        geringste_Entfernung = float("inf")

        for knoten in KnotenListe:
            gefundene_Entfernung = math.sqrt((x - knoten.x) ** 2 + (y - knoten.y) ** 2)
            if gefundene_Entfernung < geringste_Entfernung:
                geringste_Entfernung = gefundene_Entfernung
                KnotenID = knoten.ID

        if KnotenID == None:
            return

        if ersterKnoten == None:
            ersterKnoten = KnotenID
            return
        elif zweiterKnoten == None:
            zweiterKnoten = KnotenID
            stabFenster = tk.Toplevel()
            stabFenster.title("Stab hinzufügen")
            e_modul_Schriftzug = tk.Label(stabFenster, text="E-Modul:")
            e_modul_Schriftzug.grid(row=0, column=0)
            e_modul_Eingabe = tk.Entry(stabFenster)
            e_modul_Eingabe.insert(0, letzterStab.e_modul if letzterStab else "4")
            e_modul_Eingabe.grid(row=0, column=1)

            flaechenquerschnitt_Schriftzug = tk.Label(stabFenster, text="Flächenquerschnitt:")
            flaechenquerschnitt_Schriftzug.grid(row=1, column=0)
            flaechenquerschnitt_Eingabe = tk.Entry(stabFenster)
            flaechenquerschnitt_Eingabe.insert(0, letzterStab.flaechenquerschnitt if letzterStab else "44")
            flaechenquerschnitt_Eingabe.grid(row=1, column=1)

            def hinzufuegen():
                global ersterKnoten
                global zweiterKnoten
                global stabID_counter
                global letzterStab

                if ersterKnoten is None or zweiterKnoten is None:
                    messagebox.showerror("Fehler", "Bitte wählen Sie zwei Knoten aus, bevor Sie einen Stab hinzufügen.")
                    return

                # Finden Sie die Knoten mit den entsprechenden IDs
                knoten1 = next((k for k in KnotenListe if k.ID == ersterKnoten), None)
                knoten2 = next((k for k in KnotenListe if k.ID == zweiterKnoten), None)

                # Überprüfen Sie, ob die Knoten gefunden wurden
                if knoten1 is None or knoten2 is None:
                    messagebox.showerror("Fehler",
                                         "Einer oder beide der ausgewählten Knoten konnten nicht gefunden werden.")
                    return

                e_modul = float(e_modul_Eingabe.get() or 0)
                flaechenquerschnitt = float(flaechenquerschnitt_Eingabe.get() or 0)

                stab = Stab(KnotenListe[ersterKnoten], KnotenListe[zweiterKnoten], e_modul, flaechenquerschnitt,
                            stabID_counter)
                stabID_counter += 1
                stabListe.append(stab)
                letzterStab = stab
                self.Oberfläche.create_line(KnotenListe[ersterKnoten].x, KnotenListe[ersterKnoten].y,
                                            KnotenListe[zweiterKnoten].x,
                                            KnotenListe[zweiterKnoten].y,
                                            width=2.5)
                ersterKnoten = None
                zweiterKnoten = None

                stabFenster.destroy()

            Hinzufüge_Knopf = tk.Button(stabFenster, text="Hinzufügen", command=hinzufuegen)
            Hinzufüge_Knopf.grid(row=2, column=0, columnspan=2)
            Hinzufüge_Knopf.bind("<Return>", hinzufuegen)
            Hinzufüge_Knopf.focus_set()

    def calculate(self):

        knoten_data = GUI.Knoten_abfrage()
        stab_data = GUI.Stab_abfrage()
        Knotenliste = []

        Fx_informationen = []
        Fy_informationen = []

        for knoten in knoten_data:
            x, y, ID, h, v, Fx, Fy = knoten
            Knotenliste.append(Knotenn(x, y, ID, h, v))

            Fx_informationen.append([np.linspace(0, Fx, 100),ID])
            Fy_informationen.append([np.linspace(0, Fy, 100),ID])

        Stabliste = []

        for stab in stab_data:
            knoten1_ID, knoten2_ID, e_modul, flaechenquerschnitt, ID = stab

            knoten1 = next((k for k in Knotenliste if k.ID == knoten1_ID), None)
            knoten2 = next((k for k in Knotenliste if k.ID == knoten2_ID), None)

            if knoten1 is not None and knoten2 is not None:
                Stabliste.append(Stab(knoten1, knoten2, e_modul, flaechenquerschnitt, ID))

        fachwerk = Fachwerk(Stabliste, Knotenliste)

        Knotenliste_neu = copy.deepcopy(Knotenliste)
        Stabliste_neu = copy.deepcopy(Stabliste)



        if self.plot_geo_var.get():

            for i in range(100):

                for Knoten in fachwerk.Knotenliste:
                # Searching for the corresponding Fx and Fy information based on Knoten.ID
                    matching_Fx_info = next((info for info in Fx_informationen if info[1] == Knoten.ID), None)
                    matching_Fy_info = next((info for info in Fy_informationen if info[1] == Knoten.ID), None)

                    if matching_Fx_info:
                        Knoten.Fx = matching_Fx_info[0][i]

                    if matching_Fy_info:
                        Knoten.Fy = matching_Fy_info[0][i]

                fachwerk.Knoten_zuruecksetzten()
                fachwerk.Norm_history = []
                fachwerk.loesen()

                if i == 0 or i == 20 or i == 40 or i == 60 or i == 80:
                    self.plot_results(Knotenliste_neu, Stabliste_neu, Knotenliste, Stabliste)

                if i == 33:
                    if self.plot_norm_var.get():
                        fachwerk.plot_norm()



        else:

            for i in range(100):

                for Knoten in fachwerk.Knotenliste:
                    # Searching for the corresponding Fx and Fy information based on Knoten.ID
                    matching_Fx_info = next((info for info in Fx_informationen if info[1] == Knoten.ID), None)
                    matching_Fy_info = next((info for info in Fy_informationen if info[1] == Knoten.ID), None)

                    if matching_Fx_info:
                        Knoten.Fx = matching_Fx_info[0][i]

                    if matching_Fy_info:
                        Knoten.Fy = matching_Fy_info[0][i]

                fachwerk.Knoten_zuruecksetzten()
                fachwerk.Norm_history = []
                fachwerk.loesen()

                if i == 33:
                    if self.plot_norm_var.get():
                        fachwerk.plot_norm()

        self.plot_results(Knotenliste_neu, Stabliste_neu, Knotenliste, Stabliste)

    def plot_results(self, Knotenliste, Stabliste, Knotenliste_neu, Stabliste_neu):

        self.plot_structure(Knotenliste, Stabliste, 'ko', 'k-', Auflagerfarbe="black")

        self.plot_structure(Knotenliste_neu, Stabliste_neu, 'go', 'g-', Auflagerfarbe="green")

        self.plot_structure2(Knotenliste_neu, Stabliste_neu, 'bo', 'b-', Auflagerfarbe="blue")

        plt.title(
            'Fachwerk vor der Belastung (schwarz) und nach der Belastung nicht linear berechnet (grün) und linear berechnet (blau)')
        plt.show()

    def plot_structure(self, Knotenliste, Stabliste, node_color, stab_color, Auflagerfarbe):
        for k in Knotenliste:
            plt.plot(k.x, k.y, node_color)

        for s in Stabliste:
            plt.plot([s.knoten1.x, s.knoten2.x], [s.knoten1.y, s.knoten2.y], stab_color)

        a = 8
        h = a * np.sqrt(3) / 2

        for Knoten in Knotenliste:
            auflager_punkt = np.array([Knoten.x, Knoten.y])
            if Knoten.h == 1:
                punkt1 = auflager_punkt
                punkt2 = auflager_punkt + np.array([-h, -a / 2])
                punkt3 = auflager_punkt + np.array([-h, a / 2])
                plt.plot([punkt1[0], punkt2[0]], [punkt1[1], punkt2[1]], Auflagerfarbe)
                plt.plot([punkt1[0], punkt3[0]], [punkt1[1], punkt3[1]], Auflagerfarbe)
                plt.plot([punkt2[0], punkt3[0]], [punkt2[1], punkt3[1]], Auflagerfarbe)
            if Knoten.v == 1:
                punkt1 = auflager_punkt
                punkt2 = auflager_punkt + np.array([-a / 2, -h])
                punkt3 = auflager_punkt + np.array([a / 2, -h])
                plt.plot([punkt1[0], punkt2[0]], [punkt1[1], punkt2[1]], Auflagerfarbe)
                plt.plot([punkt1[0], punkt3[0]], [punkt1[1], punkt3[1]], Auflagerfarbe)
                plt.plot([punkt2[0], punkt3[0]], [punkt2[1], punkt3[1]], Auflagerfarbe)

    def plot_structure2(self, Knotenliste, Stabliste, node_color, stab_color, Auflagerfarbe):

        for k in Knotenliste:
            plt.plot(k.x_linear, k.y_linear, node_color)

        for s in Stabliste:
            plt.plot([s.knoten1.x_linear, s.knoten2.x_linear], [s.knoten1.y_linear, s.knoten2.y_linear], stab_color)

        a = 8
        h = a * np.sqrt(3) / 2

        for Knoten in Knotenliste:
            auflager_punkt = np.array([Knoten.x_linear, Knoten.y_linear])
            if Knoten.h == 1:
                punkt1 = auflager_punkt
                punkt2 = auflager_punkt + np.array([-h, -a / 2])
                punkt3 = auflager_punkt + np.array([-h, a / 2])
                plt.plot([punkt1[0], punkt2[0]], [punkt1[1], punkt2[1]], Auflagerfarbe)
                plt.plot([punkt1[0], punkt3[0]], [punkt1[1], punkt3[1]], Auflagerfarbe)
                plt.plot([punkt2[0], punkt3[0]], [punkt2[1], punkt3[1]], Auflagerfarbe)
            if Knoten.v == 1:
                punkt1 = auflager_punkt
                punkt2 = auflager_punkt + np.array([-a / 2, -h])
                punkt3 = auflager_punkt + np.array([a / 2, -h])
                plt.plot([punkt1[0], punkt2[0]], [punkt1[1], punkt2[1]], Auflagerfarbe)
                plt.plot([punkt1[0], punkt3[0]], [punkt1[1], punkt3[1]], Auflagerfarbe)
                plt.plot([punkt2[0], punkt3[0]], [punkt2[1], punkt3[1]], Auflagerfarbe)

testFenster = tk.Tk()
GUI = Fenster(testFenster)
testFenster.mainloop()


