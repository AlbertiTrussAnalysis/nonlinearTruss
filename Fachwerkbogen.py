import numpy as np
import math
import matplotlib.pyplot as plt


class Knoten:

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
        self.dof = [2 * self.ID, 2 * self.ID + 1]

    def auf_Startposition(self):
        self.x, self.y = self.x0, self.y0

class Stab:

    def __init__(self, knoten1, knoten2, e_modul, flaechenquerschnitt):
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
        self.counter = 0

    def Knoten_zuruecksetzten(self):
        for Knoten in self.Knotenliste:
            Knoten.auf_Startposition()

    def Knoten_aktualisieren(self, u):
        for Knoten in self.Knotenliste:
            dof = Knoten.dof
            Knoten.x += u[dof[0]]
            Knoten.y += u[dof[1]]

    def Residuum(self):
        r_ges = np.zeros(2 * len(self.Knotenliste))

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

        k = 0
        for Knoten in self.Knotenliste:
            r_ges[k] = r_ges[k] - Knoten.Fx
            r_ges[k + 1] = r_ges[k + 1] - Knoten.Fy
            k = k + 2


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
            if Knoten.y0 == 50:
                plt.plot(Knoten.x, Knoten.y, "ro")

        for s in self.Stabliste:
            plt.plot([s.knoten1.x, s.knoten2.x], [s.knoten1.y, s.knoten2.y], "green")

        a = 0.5
        h = a * np.sqrt(3) / 2

        for Knoten_instanz in self.Knotenliste:
            auflager_punkt = np.array([Knoten_instanz.x, Knoten_instanz.y])
            if Knoten_instanz.h == 1:
                punkt1 = auflager_punkt
                punkt2 = auflager_punkt + np.array([-h, -a / 2])
                punkt3 = auflager_punkt + np.array([-h, a / 2])
                plt.plot([punkt1[0], punkt2[0]], [punkt1[1], punkt2[1]], "black")
                plt.plot([punkt1[0], punkt3[0]], [punkt1[1], punkt3[1]], "black")
                plt.plot([punkt2[0], punkt3[0]], [punkt2[1], punkt3[1]], "black")
            if Knoten_instanz.v == 1:
                punkt1 = auflager_punkt
                punkt2 = auflager_punkt + np.array([-a / 2, -h])
                punkt3 = auflager_punkt + np.array([a / 2, -h])
                plt.plot([punkt1[0], punkt2[0]], [punkt1[1], punkt2[1]], "black")
                plt.plot([punkt1[0], punkt3[0]], [punkt1[1], punkt3[1]], "black")
                plt.plot([punkt2[0], punkt3[0]], [punkt2[1], punkt3[1]], "black")

        plt.show()

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





    def loesen(self):

        self.Norm_history=[]



        while True:
            r = self.Residuum()
            dr = self.globale_Steifigkeitsmatrix()
            u = np.linalg.solve(dr, -r)

            #self.plot_structure()
            self.Knoten_aktualisieren(u)

            res = np.linalg.norm(r)
            self.Norm_history.append(res)



            if res < 1e-4:
                break


    def linear_loesen(self):
        dr = self.linear_globale_Steifigkeitsmatrix()
        F = self.linear_F()
        u_linear = np.linalg.solve(dr, F)
        self.Knoten_aktualisieren(u_linear)






################################################################Instanzierung#########################################

with open("Koordinaten-Fachwerkbogen.txt", "r") as file:
    lines = file.readlines()

# Initialisierung der Variablen
reading_nodes = False
reading_elements = False

# Knoten-Instanzen erneut erstellen
knoten_list = []
stab_list = []

# Fortlaufende ID für Knoten und Stäbe zurücksetzen
knoten_id = 0
stab_id = 0

for line in lines:
    line = line.strip()
    if line == "Begin Nodes":
        reading_nodes = True
        continue
    elif line == "End Nodes":
        reading_nodes = False
        continue
    elif "Begin Elements" in line:
        reading_elements = True
        continue
    elif line == "End Elements":
        reading_elements = False
        continue
    if reading_nodes:
        data = line.split()
        x, y = float(data[1]), float(data[2])
        knoten = Knoten(x, y, knoten_id)
        if abs(y - 33.9) < 0.2:
            knoten.h = 1
            knoten.v = 1
        knoten_list.append(knoten)
        knoten_id += 1
    elif reading_elements:
        data = line.split()
        # Überprüfe, ob genügend Werte in der Zeile vorhanden sind
        if len(data) < 4:
            continue
        k1, k2 = map(int, data[2:])
        # Indizes korrigieren
        stab = Stab(knoten_list[k1 - 1], knoten_list[k2 - 1], flaechenquerschnitt=0.01, e_modul=5E9)
        stab_list.append(stab)
        stab_id += 1



##########################Berechnung################################################################

instanz = Fachwerk(stab_list, knoten_list)

a = 1

if a == 2:

    for k in instanz.Knotenliste:
        if k.y0 == 50:
            k.h = 1
            k.v = 1

    for k in instanz.Knotenliste:
        if k.y0 == 50:
            k.y = k.y0 - 0.8

    instanz.loesen()
    instanz.plot_structure()


else:
    N1_values = np.linspace(0, 39, 300)
    N2_values = []
    Fy2_values = np.linspace(0, 1e7, 50)
    Fy_values = []

    for N1 in N1_values:

        for k in instanz.Knotenliste:
            if k.y0 == 50:
                k.h = 1
                k.v = 1

        for k in instanz.Knotenliste:
            if k.y0 == 50:
                k.y = k.y0 - N1

        instanz.loesen()

        for k in instanz.Knotenliste:
            if k.y0 == 50:
                k.h = 0
                k.v = 0

        Fy_values.append(-instanz.Residuum()[43])

        for k in instanz.Knotenliste:
            if k.y0 == 50:
                k.h = 1
                k.v = 1


    instanz.Knoten_zuruecksetzten()
    for k in instanz.Knotenliste:
        if k.y0 == 50:
            k.h = 0
            k.v = 0

    for F in Fy2_values:

        for k in instanz.Knotenliste:
            if k.y0 == 50:
                k.Fy = -F

        instanz.loesen()

        for k in instanz.Knotenliste:
            if k.y0 == 50:
                N2_values.append(-k.y + k.y0)




    plt.figure(figsize=(10, 6))
    #plt.plot(x, load * 1e7, '-', label="Referenzlösung (Kraftgesteuert)", color='gray')
    #plt.plot(x, dis*1e7,  '-', label="Referenzlösung (Verschiebungsgesteuert)", color='lightskyblue')
    plt.plot(N2_values, Fy2_values, '--', label="Eigenentwicklung (Kraftgesteuert)", color='black')
    plt.plot(N1_values, Fy_values, '--', label="Eigenentwicklung (Verschiebungsgesteuert)", color='blue')
    plt.title("Kraft-Verschiebungsverlauf des Scheitelpunktes")
    plt.xlabel('Verschiebung in [m]')
    plt.ylabel('vertikale Kraft am Knoten in [N]')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

