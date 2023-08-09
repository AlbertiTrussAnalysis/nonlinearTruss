import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import numpy as np
import math
import matplotlib.pyplot as plt
import copy


class Knoten:
    def __init__(self, x, y, ID, h=0, v=0, Fx=0, Fy=0, vV=0, uV=0):
        self.x = x
        self.y = y
        self.ID = ID
        self.h = h
        self.v = v
        self.Fx = Fx
        self.Fy = Fy
        self.oval = None
        self.vV = vV
        self.uV = uV

class Stab_linear:
    def __init__(self, knoten1, knoten2, e_modul, flaechenquerschnitt, ID):
        self.knoten1 = knoten1
        self.knoten2 = knoten2
        self.dof = [2 * self.knoten1.ID, 2 * self.knoten1.ID + 1,
                    2 * self.knoten2.ID, 2 * self.knoten2.ID + 1]
        self.e_modul = e_modul
        self.flaechenquerschnitt = flaechenquerschnitt
        self.ID = ID
        self.L = math.sqrt((-self.knoten1.x + self.knoten2.x) ** 2.0 + (-self.knoten1.y + self.knoten2.y) ** 2.0)
        theta = np.arctan2(self.knoten2.y - self.knoten1.y, self.knoten2.x - self.knoten1.x)
        c = np.cos(theta)
        s = np.sin(theta)
        self.dr = self.flaechenquerschnitt * self.e_modul / self.L * np.array([[c ** 2, c * s, -c ** 2, -c *s],
                                                                                [c* s, s ** 2, -s * c, -s **2],
                                                                                [-c **2, -s *c, c ** 2, s * c],
                                                                                [-s *c , -s **2, s * c, s **2 ]])

class Fachwerk_linear:

    def __init__(self, Stabliste, Knotenliste):
        self.Knotenliste = Knotenliste
        self.Stabliste = Stabliste

    def F(self):
        F_ges = np.zeros(2 * len(self.Knotenliste))
        for Knoten in self.Knotenliste:
            F_ges[Knoten.ID * 2] += Knoten.Fx
            F_ges[Knoten.ID * 2 + 1] += Knoten.Fy
        return F_ges

    def dr(self):
        dr_ges = np.zeros((2 * len(self.Knotenliste), 2 * len(self.Knotenliste)))
        for Stab in self.Stabliste:
            for i, dof1 in enumerate(Stab.dof):
                for j, dof2 in enumerate(Stab.dof):
                    dr_ges[dof1, dof2] = dr_ges[dof1, dof2] + Stab.dr[i, j]
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



class Stab:

    def __init__(self, knoten1, knoten2, e_modul, flaechenquerschnitt, ID):
        self.knoten1 = knoten1
        self.knoten2 = knoten2
        self.dof = [2 * self.knoten1.ID, 2 * self.knoten1.ID + 1,
                    2 * self.knoten2.ID, 2 * self.knoten2.ID + 1]
        self.e_modul = e_modul
        self.flaechenquerschnitt = flaechenquerschnitt
        self.ID = ID
        self.L = math.sqrt((-self.knoten1.x + self.knoten2.x) ** 2.0 + (-self.knoten1.y + self.knoten2.y) ** 2.0)
        self.update()


    def update(self):
            self.e = (-0.5 * (-self.knoten1.x + self.knoten2.x) ** 2 - 0.5 * (
                    -self.knoten1.y + self.knoten2.y) ** 2 + 0.5 * (
                              -self.knoten1.uV - self.knoten1.x + self.knoten2.uV + self.knoten2.x) ** 2 + 0.5 * (
                              -self.knoten1.vV - self.knoten1.y + self.knoten2.vV + self.knoten2.y) ** 2) / (
                             (-self.knoten1.x + self.knoten2.x) ** 2 + (-self.knoten1.y + self.knoten2.y) ** 2)
            de_du1 = (1.0 * self.knoten1.uV + 1.0 * self.knoten1.x - 1.0 * self.knoten2.uV - 1.0 * self.knoten2.x) / (
                    (-self.knoten1.x + self.knoten2.x) ** 2 + (-self.knoten1.y + self.knoten2.y) ** 2)
            de_dv1 = (1.0 * self.knoten1.vV + 1.0 * self.knoten1.y - 1.0 * self.knoten2.vV - 1.0 * self.knoten2.y) / (
                    (-self.knoten1.x + self.knoten2.x) ** 2 + (-self.knoten1.y + self.knoten2.y) ** 2)
            de_du2 = (-1.0 * self.knoten1.uV - 1.0 * self.knoten1.x + 1.0 * self.knoten2.uV + 1.0 * self.knoten2.x) / (
                    (-self.knoten1.x + self.knoten2.x) ** 2 + (-self.knoten1.y + self.knoten2.y) ** 2)
            de_dv2 = (-1.0 * self.knoten1.vV - 1.0 * self.knoten1.y + 1.0 * self.knoten2.vV + 1.0 * self.knoten2.y) / (
                    (-self.knoten1.x + self.knoten2.x) ** 2 + (-self.knoten1.y + self.knoten2.y) ** 2)
            d2e_d2x = 1.0 / ((-self.knoten1.x + self.knoten2.x) ** 2 + (-self.knoten1.y + self.knoten2.y) ** 2)
            d2e_du1du2 = -1.0 / ((-self.knoten1.x + self.knoten2.x) ** 2 + (-self.knoten1.y + self.knoten2.y) ** 2)
            self.r = self.e * self.e_modul * self.flaechenquerschnitt * self.L * np.array(
                [de_du1, de_dv1, de_du2, de_dv2]) - np.array(
                [self.knoten1.Fx, self.knoten1.Fy, self.knoten2.Fx, self.knoten2.Fy])
            self.dr = self.e_modul * self.flaechenquerschnitt * self.L * \
                      np.array([[de_du1 * de_du1 + self.e * d2e_d2x, de_dv1 * de_du1,
                                 de_du2 * de_du1 + self.e * d2e_du1du2, de_dv2 * de_du1],
                                [de_du1 * de_dv1, de_dv1 * de_dv1 + self.e * d2e_d2x,
                                 de_du2 * de_dv1, de_dv2 * de_dv1 + self.e * d2e_du1du2],
                                [de_du1 * de_du2 + self.e * d2e_du1du2, de_dv1 * de_du2,
                                 de_du2 * de_du2 + self.e * d2e_d2x, de_dv2 * de_du2],
                                [de_du1 * de_dv2, de_dv1 * de_dv2 + self.e * d2e_du1du2,
                                 de_du2 * de_dv2, de_dv2 * de_dv2 + self.e * d2e_d2x]])


class Fachwerk:

    def __init__(self, Stabliste, Knotenliste):
        self.Knotenliste = Knotenliste
        self.Stabliste = Stabliste

    def r(self):
        r_ges = np.zeros(2 * len(self.Knotenliste))
        for Stab in self.Stabliste:
            for i, dof in enumerate(Stab.dof):
                r_ges[dof] += Stab.r[i]
            if Stab.knoten1.h == 1:
                r_ges[Stab.dof[0]] = 0
            if Stab.knoten1.v == 1:
                r_ges[Stab.dof[1]] = 0
            if Stab.knoten2.h == 1:
                r_ges[Stab.dof[2]] = 0
            if Stab.knoten2.v == 1:
                r_ges[Stab.dof[3]] = 0
        return r_ges

    def dr(self):
        dr_ges = np.zeros((2 * len(self.Knotenliste), 2 * len(self.Knotenliste)))
        for Stab in self.Stabliste:
            for i, dof1 in enumerate(Stab.dof):
                for j, dof2 in enumerate(Stab.dof):
                    dr_ges[dof1, dof2] = dr_ges[dof1, dof2] + Stab.dr[i, j]
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

class NewtonRaphson:
    def __init__(self, Fachwerk, e=1e-3):
        self.Fachwerk = Fachwerk
        self.e = e
        self.norm_history = []

    def newton_raphson(self):
        u = np.zeros(2 * len(self.Fachwerk.Knotenliste))
        while True:
            ru = self.Fachwerk.r()
            dru = self.Fachwerk.dr()
            norm_ru = np.linalg.norm(ru)
            self.norm_history.append(norm_ru)
            if norm_ru < self.e:
                return u
            du = np.linalg.solve(dru, ru)
            u = u - du
            for Knoten in self.Fachwerk.Knotenliste:
                Knoten.uV -= du[Knoten.ID * 2]
                Knoten.vV -= du[Knoten.ID * 2 + 1]
            for Stab in self.Fachwerk.Stabliste:
                Stab.update()

    def plot_norm(self):
        plt.figure(figsize=(10, 6))
        plt.plot(self.norm_history, marker='o')
        plt.title('Die Zweite Euklidische Norm von r(u) in Abhängigkeit der Iteration')
        plt.xlabel('Iterationen')
        plt.ylabel('Die zweite Euklidische Norm von r(u)')
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
        Fenster.title("geometrisch nicht lineare Fachwerke")

        self.Oberfläche = tk.Canvas(Fenster, width=1400, height=800)
        self.Oberfläche.pack()

        self.Oberfläche.bind("<Button-1>", self.Knoten_hinzufügen)
        self.Oberfläche.bind("<Button-2>", self.Knoten_bearbeiten)
        self.Oberfläche.bind("<Button-3>", self.Stab_hinzufuegen)

        #self.button2 = tk.Button(Fenster, text="Visualisieren", command=self.plot_results)
        #self.button2.pack(side=tk.BOTTOM, anchor=tk.S)


        #self.button = tk.Button(Fenster, text="Neustart", command=self.restart_program)
        #self.button.pack(side=tk.BOTTOM, anchor=tk.SE)

        self.button = tk.Button(Fenster, text="Berechnen", command=self.calculate)
        self.button.pack(side=tk.BOTTOM, anchor=tk.S)


        raster_größe = 40

        for x in range(0, 1400, raster_größe):
            self.Oberfläche.create_line(x, 0, x, 800, fill="gray")

        for y in range(0, 800, raster_größe):
            self.Oberfläche.create_line(0, y, 1400, y, fill="gray")

        pfeil_y = self.Oberfläche.create_line(20, 760, 20, 720, arrow=tk.LAST, fill="black")
        pfeil_x = self.Oberfläche.create_line(20, 760, 60, 760, arrow=tk.LAST, fill="black")
        self.Oberfläche.create_text(10, 760, text="y", fill="black")
        self.Oberfläche.create_text(60, 770, text="x", fill="black")

        self.Oberfläche.create_text(700, 785, text="Bitte beachten sie, dass die Stäbe erst nach den Knoten erstellt werden dürfen!", font=("Arial", 12, "bold"), fill="black")
        self.Oberfläche.create_text(180, 90, text="Linke Maustaste erstellt Knoten", font=("Arial", 12, "bold"), fill="black")
        self.Oberfläche.create_text(180, 110, text="Mittlere Maustaste bearbeitet Knoten", font=("Arial", 12, "bold"), fill="black")
        self.Oberfläche.create_text(180, 130, text="Rechte Maustaste auf zwei Knoten erstellt Stab", font=("Arial", 12, "bold"), fill="black")


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
            knoten = Knoten(x_knoten, y_knoten, ID_counter)
            knoten.oval = Punkt
            ID_counter = ID_counter + 1
            KnotenListe.append(knoten)
            coordFenster.destroy()

        Hinzufüge_Knopf = tk.Button(coordFenster, text="Hinzufügen", command=hinzufuegen)
        Hinzufüge_Knopf.grid(row=2, column=0, columnspan=2)
        Hinzufüge_Knopf.bind("<Return>", hinzufuegen)
        Hinzufüge_Knopf.focus_set()

    def Knoten_abfrage(self):
        return [(knoten.x, 800 - knoten.y, knoten.ID, knoten.h, knoten.v, knoten.Fx, knoten.Fy) for knoten in KnotenListe]

    def Stab_abfrage(self):
        return [(Stab.knoten1.ID, Stab.knoten2.ID, Stab.e_modul, Stab.flaechenquerschnitt, Stab.ID) for Stab in stabListe]

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
        #Neigung_Schriftzug = tk.Label(Extrafenster, text="Neigung in [°]:")
        #Neigung_Schriftzug.grid(row=3, column=1)
        #Neigung_Eingabe = tk.Entry(Extrafenster)
        #Neigung_Eingabe.grid(row=3, column=2)

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
            flaechenquerschnitt_Eingabe.insert(0,letzterStab.flaechenquerschnitt if letzterStab else "44")
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
##############################nonlinear################################################################################
        global Knoten

        knoten_data = GUI.Knoten_abfrage()
        stab_data = GUI.Stab_abfrage()
        Knotenliste = []

        for knoten in knoten_data:
            x, y, ID, h, v, Fx, Fy = knoten
            Knotenliste.append(Knoten(x, y, ID, h, v, Fx, Fy))

        Stabliste = []

        for stab in stab_data:
            knoten1_ID, knoten2_ID, e_modul, flaechenquerschnitt, ID = stab

            knoten1 = next((k for k in Knotenliste if k.ID == knoten1_ID), None)
            knoten2 = next((k for k in Knotenliste if k.ID == knoten2_ID), None)

            if knoten1 is not None and knoten2 is not None:
                Stabliste.append(Stab(knoten1, knoten2, e_modul, flaechenquerschnitt, ID))

        fachwerk = Fachwerk(Stabliste, Knotenliste)
        N1 = NewtonRaphson(fachwerk).newton_raphson()

        Knotenliste_neu = copy.deepcopy(Knotenliste)

        for Knoten in Knotenliste_neu:
            Knoten.x += Knoten.uV
            Knoten.y += Knoten.vV


        #diese Stabliste ist notwendig, weil die neuen Koordinaten der Stabentpunkte aktualisiert werden muss
        #damit die visualiesierung des stabes als linie richtig ist. und nich nur die knoten neu visualiert werden
        Stabliste_neu = []

        for stab in Stabliste:
            knoten1_ID, knoten2_ID, e_modul, flaechenquerschnitt, ID = stab.knoten1.ID, stab.knoten2.ID, stab.e_modul, stab.flaechenquerschnitt, stab.ID

            knoten1 = next((k for k in Knotenliste_neu if k.ID == knoten1_ID), None)
            knoten2 = next((k for k in Knotenliste_neu if k.ID == knoten2_ID), None)

            if knoten1 is not None and knoten2 is not None:
                Stabliste_neu.append(Stab(knoten1, knoten2, e_modul, flaechenquerschnitt, ID))

#################################linear#######################################################################

        Knotenliste_neu_linear = copy.deepcopy(Knotenliste)

        fachwerk_linear = Fachwerk_linear(Stabliste, Knotenliste)
        dr = fachwerk_linear.dr()
        F = fachwerk_linear.F()
        u_linear = np.linalg.solve(dr, F)

        for Knoten in Knotenliste_neu_linear:
            Knoten.x = Knoten.x + u_linear[Knoten.ID * 2]
            Knoten.y = Knoten.y + u_linear[Knoten.ID * 2 + 1]


        Stabliste_neu_linear = []

        for stab in Stabliste:
            knoten1_ID, knoten2_ID, e_modul, flaechenquerschnitt, ID = stab.knoten1.ID, stab.knoten2.ID, stab.e_modul, stab.flaechenquerschnitt, stab.ID

            knoten1 = next((k for k in Knotenliste_neu_linear if k.ID == knoten1_ID), None)
            knoten2 = next((k for k in Knotenliste_neu_linear if k.ID == knoten2_ID), None)

            if knoten1 is not None and knoten2 is not None:
                Stabliste_neu_linear.append(Stab(knoten1, knoten2, e_modul, flaechenquerschnitt, ID))


######################################################################################################################

        self.plot_results(Knotenliste, Stabliste, Knotenliste_neu, Stabliste_neu, Knotenliste_neu_linear, Stabliste_neu_linear)



    def plot_results(self, Knotenliste, Stabliste, Knotenliste_neu, Stabliste_neu, Knotenliste_neu_linear, Stabliste_neu_linear):

        self.plot_structure(Knotenliste, Stabliste, 'ko', 'k-', Auflagerfarbe="black")

        self.plot_structure(Knotenliste_neu, Stabliste_neu, 'go', 'g-', Auflagerfarbe="green")

        self.plot_structure(Knotenliste_neu_linear, Stabliste_neu_linear, 'bo', 'b-', Auflagerfarbe="green")

        plt.title('Fachwerk vor der Belastung (schwarz) und nach der Belastung nicht linear berechnet (grün) und linear berechnet (blau)')
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

testFenster = tk.Tk()
GUI = Fenster(testFenster)
testFenster.mainloop()

#print(GUI.Knoten_abfrage())
#print(GUI.Stab_abfrage())
