import copy
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import patches


class Knoten:

    def __init__(self, x, y, z, ID, h=0, v=0, w=0, Fx=0, Fy=0, Fz=0):
        self.x = x
        self.y = y
        self.z = z
        self.ID = ID
        self.h = h
        self.v = v
        self.w = w
        self.Fx = Fx
        self.Fy = Fy
        self.Fz = Fz
        self.x0 = x
        self.y0 = y
        self.z0 = z
        self.dof = [3 * self.ID, 3 * self.ID + 1, 3 * self.ID + 2]

    def auf_Startposition(self):
        self.x, self.y, self.z = self.x0, self.y0, self.z0

class Stab:

    def __init__(self, knoten1, knoten2, e_modul, flaechenquerschnitt, ID):
        self.knoten1 = knoten1
        self.knoten2 = knoten2
        self.e_modul = e_modul
        self.ID = ID
        self.flaechenquerschnitt = flaechenquerschnitt
        self.dof = [3 * self.knoten1.ID, 3 * self.knoten1.ID + 1, 3 * self.knoten1.ID + 2,
                    3 * self.knoten2.ID, 3 * self.knoten2.ID + 1, 3 * self.knoten2.ID + 2]
        self.L = math.sqrt((self.knoten1.x0 - self.knoten2.x0) ** 2 + (self.knoten1.y0 - self.knoten2.y0) ** 2 + (self.knoten1.z0 - self.knoten2.z0) ** 2)


    def F_int(self):
        u_1, v_1, w_1, u_2, v_2, w_2 = (self.knoten1.x - self.knoten1.x0,
                               self.knoten1.y - self.knoten1.y0,
                               self.knoten1.z - self.knoten1.z0,
                               self.knoten2.x - self.knoten2.x0,
                               self.knoten2.y - self.knoten2.y0,
                               self.knoten2.z - self.knoten2.z0)

        x1_0, y1_0, z1_0, x2_0, y2_0, z2_0 = (self.knoten1.x0,
                                     self.knoten1.y0,
                                     self.knoten1.z0,
                                     self.knoten2.x0,
                                     self.knoten2.y0,
                                     self.knoten2.z0)

        #print("y0_Wert links:    ", y1_0)
        #print("y0_Wert rechts:    ", y2_0)


        res1 = self.e_modul * self.flaechenquerschnitt * self.L * (0.5*(-(-x1_0 + x2_0)**2 - (-y1_0 + y2_0)**2 - (-z1_0 + z2_0)**2 + (-u_1 + u_2 - x1_0 + x2_0)**2 + (-v_1 + v_2 - y1_0 + y2_0)**2 + (-w_1 + w_2 - z1_0 + z2_0)**2)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2 + (-z1_0 + z2_0)**2)) * \
               np.array([[0.5*(2*u_1 - 2*u_2 + 2*x1_0 - 2*x2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2 + (-z1_0 + z2_0)**2)], [0.5*(2*v_1 - 2*v_2 + 2*y1_0 - 2*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2 + (-z1_0 + z2_0)**2)], [0.5*(2*w_1 - 2*w_2 + 2*z1_0 - 2*z2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2 + (-z1_0 + z2_0)**2)], [0.5*(-2*u_1 + 2*u_2 - 2*x1_0 + 2*x2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2 + (-z1_0 + z2_0)**2)], [0.5*(-2*v_1 + 2*v_2 - 2*y1_0 + 2*y2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2 + (-z1_0 + z2_0)**2)], [0.5*(-2*w_1 + 2*w_2 - 2*z1_0 + 2*z2_0)/((-x1_0 + x2_0)**2 + (-y1_0 + y2_0)**2 + (-z1_0 + z2_0)**2)]])
        #print("F_int:   ")
        #print(res1)
        #print("          ")
        #print("          ")
        return res1

    def Elementsteifigkeitsmatrix(self):
        u_1, v_1, w_1, u_2, v_2, w_2 = (self.knoten1.x - self.knoten1.x0,
                                        self.knoten1.y - self.knoten1.y0,
                                        self.knoten1.z - self.knoten1.z0,
                                        self.knoten2.x - self.knoten2.x0,
                                        self.knoten2.y - self.knoten2.y0,
                                        self.knoten2.z - self.knoten2.z0)

        x1_0, y1_0, z1_0, x2_0, y2_0, z2_0 = (self.knoten1.x0,
                                              self.knoten1.y0,
                                              self.knoten1.z0,
                                              self.knoten2.x0,
                                              self.knoten2.y0,
                                              self.knoten2.z0)

        #print("Ursprüngliche Geometrie:                  aktuelle Geometrie:")
        #print("x0_Wert Knoten", self.knoten1.ID, "ist", self.knoten1.x0, "              x_Wert Knoten", self.knoten1.ID,
              #"ist", self.knoten1.x)
        #print("y0_Wert Knoten", self.knoten1.ID, "ist", self.knoten1.y0, "              y_Wert Knoten", self.knoten1.ID,
              #"ist", self.knoten1.y)
        #print("x0_Wert Knoten", self.knoten2.ID, "ist", self.knoten2.x0, "              x_Wert Knoten", self.knoten2.ID,
              #"ist", self.knoten2.x)
        #print("y0_Wert Knoten", self.knoten2.ID, "ist", self.knoten2.y0, "              y_Wert Knoten", self.knoten2.ID,
              #"ist", self.knoten2.y)


        res2 = self.e_modul * self.flaechenquerschnitt * self.L * \
               np.array([[0.25 * (2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) ** 2 / (
                           (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 + 0.5 * (
                                    -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                        -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                    (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                        0.25 * (2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) * (2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) / (
                                    (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                        0.25 * (2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) * (2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) / (
                                    (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                        0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) * (
                                    2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) / (
                                    (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 - 0.5 * (
                                    -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                        -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                    (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                        0.25 * (2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) * (
                                    -2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) / (
                                    (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                        0.25 * (2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) * (
                                    -2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) / (
                                    (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2], [
                           0.25 * (2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) * (
                                       2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) ** 2 / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 + 0.5 * (
                                       -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                           -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                   -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) * (
                                       2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) * (
                                       2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) * (
                                       2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 - 0.5 * (
                                       -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                           -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                   -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) * (
                                       -2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2], [
                           0.25 * (2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) * (
                                       2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) * (
                                       2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) ** 2 / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 + 0.5 * (
                                       -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                           -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                   -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) * (
                                       2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) * (
                                       2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) * (
                                       2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 - 0.5 * (
                                       -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                           -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                   -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2], [
                           0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) * (
                                       2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 - 0.5 * (
                                       -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                           -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                   -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) * (
                                       2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) * (
                                       2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) ** 2 / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 + 0.5 * (
                                       -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                           -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                   -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) * (
                                       -2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) * (
                                       -2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2], [
                           0.25 * (2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) * (
                                       -2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) * (
                                       2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 - 0.5 * (
                                       -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                           -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                   -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) * (
                                       2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) * (
                                       -2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) ** 2 / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 + 0.5 * (
                                       -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                           -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                   -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) * (
                                       -2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2], [
                           0.25 * (2 * u_1 - 2 * u_2 + 2 * x1_0 - 2 * x2_0) * (
                                       -2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (2 * v_1 - 2 * v_2 + 2 * y1_0 - 2 * y2_0) * (
                                       -2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) * (
                                       2 * w_1 - 2 * w_2 + 2 * z1_0 - 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 - 0.5 * (
                                       -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                           -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                   -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * u_1 + 2 * u_2 - 2 * x1_0 + 2 * x2_0) * (
                                       -2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * v_1 + 2 * v_2 - 2 * y1_0 + 2 * y2_0) * (
                                       -2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2,
                           0.25 * (-2 * w_1 + 2 * w_2 - 2 * z1_0 + 2 * z2_0) ** 2 / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2 + 0.5 * (
                                       -(-x1_0 + x2_0) ** 2 - (-y1_0 + y2_0) ** 2 - (-z1_0 + z2_0) ** 2 + (
                                           -u_1 + u_2 - x1_0 + x2_0) ** 2 + (-v_1 + v_2 - y1_0 + y2_0) ** 2 + (
                                                   -w_1 + w_2 - z1_0 + z2_0) ** 2) / (
                                       (-x1_0 + x2_0) ** 2 + (-y1_0 + y2_0) ** 2 + (-z1_0 + z2_0) ** 2) ** 2]])
        #print("Elementsteifigkeitsmatrix:   ")
        #print(res2)
        #print("          ")
        #print("          ")
        return res2

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
            Knoten.z += u[dof[2]]

    def Residuum(self):
        r_ges = np.zeros(3 * len(self.Knotenliste))

        for Stab in self.Stabliste:
            r = Stab.F_int()
            for i, dof in enumerate(Stab.dof):
                r_ges[dof] += r[i]
            if Stab.knoten1.h == 1:
                r_ges[Stab.dof[0]] = 0
            if Stab.knoten1.v == 1:
                r_ges[Stab.dof[1]] = 0
            if Stab.knoten1.w == 1:  # Neue z-Richtungsbedingung
                r_ges[Stab.dof[2]] = 0
            if Stab.knoten2.h == 1:
                r_ges[Stab.dof[3]] = 0
            if Stab.knoten2.v == 1:
                r_ges[Stab.dof[4]] = 0
            if Stab.knoten2.w == 1:  # Neue z-Richtungsbedingung
                r_ges[Stab.dof[5]] = 0

        k = 0
        for Knoten in self.Knotenliste:
            r_ges[k] = r_ges[k] - Knoten.Fx
            r_ges[k + 1] = r_ges[k + 1] - Knoten.Fy
            r_ges[k + 2] = r_ges[k + 2] - Knoten.Fz  # Neue z-Richtung
            k = k + 3

        #print("Residualvektor:  ", r_ges)
        #print("          ")
        #print("          ")
        #print("          ")

        return r_ges

    def globale_Steifigkeitsmatrix(self):
        dr_ges = np.zeros((3 * len(self.Knotenliste), 3 * len(self.Knotenliste)))

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
            if Stab.knoten1.w == 1:  # Neue z-Richtungsbedingung
                dr_ges[Stab.dof[2], :] = 0
                dr_ges[:, Stab.dof[2]] = 0
                dr_ges[Stab.dof[2], Stab.dof[2]] = 1
            if Stab.knoten2.h == 1:
                dr_ges[Stab.dof[3], :] = 0
                dr_ges[:, Stab.dof[3]] = 0
                dr_ges[Stab.dof[3], Stab.dof[3]] = 1
            if Stab.knoten2.v == 1:
                dr_ges[Stab.dof[4], :] = 0
                dr_ges[:, Stab.dof[4]] = 0
                dr_ges[Stab.dof[4], Stab.dof[4]] = 1
            if Stab.knoten2.w == 1:  # Neue z-Richtungsbedingung
                dr_ges[Stab.dof[5], :] = 0
                dr_ges[:, Stab.dof[5]] = 0
                dr_ges[Stab.dof[5], Stab.dof[5]] = 1

        #print("globale Steifikgeitsmatrix:  ")
        #print(dr_ges)
        #print("          ")
        #print("          ")


        return dr_ges

    def loesen(self):

        self.Norm_history=[]



        while True:
            r = self.Residuum()
            dr = self.globale_Steifigkeitsmatrix()
            u = np.linalg.solve(dr, -r)




            self.Knoten_aktualisieren(u)

            res = np.linalg.norm(r)

            self.Norm_history.append(res)



            if abs(res) < 1e-4:
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

#####################################################################################################


# Lists to store instances
nodes = []
elements = []

with open("koordinaten3D.txt", 'r') as file:
    lines = file.readlines()
    mode = None

    for line in lines:
        line = line.strip()
        if not line: continue

        if "Begin Nodes" in line:
            mode = "Nodes"
            continue
        elif "End Nodes" in line:
            mode = None
        elif "Begin Elements" in line:
            mode = "Elements"
            continue
        elif "End Elements" in line:
            mode = None

        if mode == "Nodes":
            parts = line.split()
            node = Knoten(float(parts[1]), float(parts[2]), float(parts[3]), int(parts[0]))
            nodes.append(node)
        elif mode == "Elements":
            parts = line.split()
            stab = Stab(nodes[int(parts[1])], nodes[int(parts[2])], flaechenquerschnitt=99, e_modul=68947572800, ID=int(parts[0]))
            elements.append(stab)



for k in nodes:
    if k.ID==0 or k.ID==1 or k.ID==2 or k.ID==6 or k.ID==11 or k.ID==15 or k.ID==18 or k.ID==17 or k.ID==16 or k.ID==12 or k.ID==7 or k.ID==3:
        k.h=1
        k.v=1
        k.w=1

for s in elements:
    if s.ID == 9:
        s.flaechenquerschnitt = 0.001091030076
    if s.ID == 10:
        s.flaechenquerschnitt = 0.001091030076
    if s.ID == 15:
        s.flaechenquerschnitt = 0.001091030076
    if s.ID == 20:
        s.flaechenquerschnitt = 0.001091030076
    if s.ID == 19:
        s.flaechenquerschnitt = 0.001091030076
    if s.ID == 14:
        s.flaechenquerschnitt = 0.001091030076
    if s.ID == 8:
        s.flaechenquerschnitt = 0.0008926434
    if s.ID == 5:
        s.flaechenquerschnitt = 0.0008926434
    if s.ID == 11:
        s.flaechenquerschnitt = 0.0008926434
    if s.ID == 18:
        s.flaechenquerschnitt = 0.0008926434
    if s.ID == 24:
        s.flaechenquerschnitt = 0.0008926434
    if s.ID == 21:
        s.flaechenquerschnitt = 0.0008926434
    if s.ID == 13:
        s.flaechenquerschnitt = 0.0001586448
    if s.ID == 0:
        s.flaechenquerschnitt = 0.0001586448
    if s.ID == 3:
        s.flaechenquerschnitt = 0.0001586448
    if s.ID == 16:
        s.flaechenquerschnitt = 0.0001586448
    if s.ID == 29:
        s.flaechenquerschnitt = 0.0001586448
    if s.ID == 26:
        s.flaechenquerschnitt = 0.0001586448
    if s.ID == 7:
        s.flaechenquerschnitt = 0.000064516 #1
    if s.ID == 4:
        s.flaechenquerschnitt = 0.000064516 #2
    if s.ID == 1:
        s.flaechenquerschnitt = 0.000064516 #3
    if s.ID == 2:
        s.flaechenquerschnitt = 0.000064516 #4
    if s.ID == 6:
        s.flaechenquerschnitt = 0.000064516 #5
    if s.ID == 12:
        s.flaechenquerschnitt = 0.000064516 #6
    if s.ID == 22:
        s.flaechenquerschnitt = 0.000064516 #7
    if s.ID == 25:
        s.flaechenquerschnitt = 0.000064516 #8
    if s.ID == 28:
        s.flaechenquerschnitt = 0.000064516 #9
    if s.ID == 27:
        s.flaechenquerschnitt = 0.000064516 #10
    if s.ID == 23:
        s.flaechenquerschnitt = 0.000064516 #11
    if s.ID == 17:
        s.flaechenquerschnitt = 0.000064516 #12


instanz = Fachwerk(elements, nodes)


###############################PDF-visualisierung######################

def visualize_structure(Fachwerk1, Fachwerk2):
    quadrat_groesse = 0.5
    quadrat_groesse_z = 0.1

    fig = plt.figure(figsize=(15, 15))

    # 3D plot for the structure
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    ax1.set_title("3D Ansicht")

    # Plotting for Fachwerk1 object
    for stab in Fachwerk1.Stabliste:
        ax1.plot([stab.knoten1.x, stab.knoten2.x],
                 [stab.knoten1.y, stab.knoten2.y],
                 [stab.knoten1.z, stab.knoten2.z], color="green")
    for Knoten in Fachwerk1.Knotenliste:
        ax1.text(Knoten.x, Knoten.y, Knoten.z, f"{Knoten.ID}", fontsize=10, color="blue")

    # Plotting for Fachwerk2 object
    for stab in Fachwerk2.Stabliste:
        ax1.plot([stab.knoten1.x, stab.knoten2.x],
                 [stab.knoten1.y, stab.knoten2.y],
                 [stab.knoten1.z, stab.knoten2.z], color="black", linestyle='--')
    #for Knoten in Fachwerk2.Knotenliste:
        #ax1.text(Knoten.x, Knoten.y, Knoten.z, f"{Knoten.ID}", fontsize=10, color="blue")

    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')

    # XY projection
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.set_title("XY Projektion")
    ax2.title.set_position([.5, 1.1])

    # Plotting for Fachwerk1 object in XY projection
    for Stab in Fachwerk1.Stabliste:
        if Stab.flaechenquerschnitt == 0.000064516:
            ax2.plot([Stab.knoten1.x, Stab.knoten2.x],
                     [Stab.knoten1.y, Stab.knoten2.y], color="green")
        else:
            ax2.plot([Stab.knoten1.x, Stab.knoten2.x],
                     [Stab.knoten1.y, Stab.knoten2.y], color="red")

    for Knoten in Fachwerk1.Knotenliste:
        ax2.text(Knoten.x, Knoten.y, f"{Knoten.ID}", fontsize=10, color="blue")

        if Knoten.h == 1 and Knoten.v == 1 and Knoten.w == 1:
            x_ecke = Knoten.x - quadrat_groesse / 2
            y_ecke = Knoten.y - quadrat_groesse / 2
            rect = patches.Rectangle((x_ecke, y_ecke), quadrat_groesse, quadrat_groesse, color="black")
            ax2.add_patch(rect)

    # Plotting for Fachwerk2 object in XY projection
    for Stab in Fachwerk2.Stabliste:
        ax2.plot([Stab.knoten1.x, Stab.knoten2.x],
                 [Stab.knoten1.y, Stab.knoten2.y], color="green")
    #for Knoten in Fachwerk2.Knotenliste:
        #ax2.text(Knoten.x, Knoten.y, f"{Knoten.ID}", fontsize=10, color="blue")

    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')

    # XZ projection
    ax3 = fig.add_subplot(2, 2, 3)
    ax3.set_title("XZ Projektion")
    ax3.title.set_position([.5, 1.1])

    # Plotting for Fachwerk1 object in XZ projection
    for Stab in Fachwerk1.Stabliste:
        ax3.plot([Stab.knoten1.x, Stab.knoten2.x],
                 [Stab.knoten1.z, Stab.knoten2.z], color="green")
    for Knoten in Fachwerk1.Knotenliste:
        ax3.text(Knoten.x, Knoten.z, f"{Knoten.ID}", fontsize=10, color="blue")

        if Knoten.h == 1 and Knoten.v == 1 and Knoten.w == 1:
            x_ecke = Knoten.x - quadrat_groesse / 2
            z_ecke = Knoten.z - quadrat_groesse_z / 2  # Anpassung für die Z-Richtung
            rect = patches.Rectangle((x_ecke, z_ecke), quadrat_groesse, quadrat_groesse_z,
                                     color="black")  # Anpassung für die Z-Richtung
            ax3.add_patch(rect)

    # Plotting for Fachwerk2 object in XZ projection
    for Stab in Fachwerk2.Stabliste:
        ax3.plot([Stab.knoten1.x, Stab.knoten2.x],
                 [Stab.knoten1.z, Stab.knoten2.z], color="black", linestyle='--')
    #for Knoten in Fachwerk2.Knotenliste:
        #ax3.text(Knoten.x, Knoten.z, f"{Knoten.ID}", fontsize=10, color="blue")

    ax3.set_xlabel('X')
    ax3.set_ylabel('Z')

    # YZ projection
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.set_title("YZ Projektion")
    ax4.title.set_position([.5, 1.1])

    # Plotting for Fachwerk1 object in YZ projection
    for Stab in Fachwerk1.Stabliste:
        ax4.plot([Stab.knoten1.y, Stab.knoten2.y],
                 [Stab.knoten1.z, Stab.knoten2.z], color="green")
    for Knoten in Fachwerk1.Knotenliste:
        ax4.text(Knoten.y, Knoten.z, f"{Knoten.ID}", fontsize=10, color="blue")

        if Knoten.h == 1 and Knoten.v == 1 and Knoten.w == 1:
            y_ecke = Knoten.y - quadrat_groesse / 2
            z_ecke = Knoten.z - quadrat_groesse_z / 2  # Anpassung für die Z-Richtung
            rect = patches.Rectangle((y_ecke, z_ecke), quadrat_groesse, quadrat_groesse_z,
                                     color="black")  # Anpassung für die Z-Richtung
            ax4.add_patch(rect)

    # Plotting for Fachwerk2 object in YZ projection
    for Stab in Fachwerk2.Stabliste:
        ax4.plot([Stab.knoten1.y, Stab.knoten2.y],
                 [Stab.knoten1.z, Stab.knoten2.z], color="black", linestyle='--')
    #for Knoten in Fachwerk2.Knotenliste:
        #ax4.text(Knoten.y, Knoten.z, f"{Knoten.ID}", fontsize=10, color="blue")

    ax4.set_xlabel('Y')
    ax4.set_ylabel('Z')

    # Display the plot
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4, wspace=0.3, bottom=0.15)
    plt.show()


#hier einstellen ob verlauf oder plot!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
a = 1

if a == 1:
    N1_values = np.linspace(0, 0.3, 300)
    N2_values = []
    Fy2_values = np.linspace(0, 9000, 500)
    Fy_values = []

    for N1 in N1_values:

        for k in instanz.Knotenliste:
            if k.ID == 9:
                k.h = 1
                k.v = 1
                k.w = 1


        for k in instanz.Knotenliste:
            if k.ID == 9:
                k.z = k.z0 - N1

        instanz.loesen()


        for k in instanz.Knotenliste:
            if k.ID == 9:
                k.h = 0
                k.v = 0
                k.w = 0

        Fy_values.append(-instanz.Residuum()[29])

        for k in instanz.Knotenliste:
            if k.ID == 9:
                k.h = 1
                k.v = 1
                k.w = 1

    instanz.Knoten_zuruecksetzten()
    for k in instanz.Knotenliste:
        if k.ID == 9:
            k.h = 0
            k.v = 0
            k.w = 0

    for F in Fy2_values:
        for k in instanz.Knotenliste:
            if k.ID == 9:
                k.Fz = -F
        instanz.loesen()
        N2_values.append(abs(nodes[9].z-nodes[9].z0))



        # Daten mit Matplotlib plotten
    plt.figure(figsize=(10, 6))
    #plt.plot(y_values, x_values, '-', label="Referenzlösung", color='grey')
    plt.plot(N1_values, Fy_values, '--', label="Eigenentwicklung (Verschiebungsgesteuert)", color='C0')
    plt.plot(N2_values, Fy2_values, '--', label="Eigenentwicklung (Kraftgesteuert)", color='C1')
    plt.title("Kraft-Verschiebungsverlauf des höchsten Punktes")
    plt.xlabel('Verschiebung in [m]')
    plt.ylabel('vertikale Kraft am Knoten in [N]')
    plt.xlim(0, 0.275)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if a == 2:
    #8896.44
    Fy_values = np.linspace(0, 9000, 300)
    N1_values = []

    i = 0

    for Fyy in Fy_values:
        for Knoten in nodes:
            if Knoten.ID == 9:
                Knoten.Fz = - Fyy
        instanz.loesen()
        i = i + 1
        if i == 50:
            instanz.plot_norm()
        N1_values.append(abs(nodes[9].z-nodes[9].z0))




    # Daten mit Matplotlib plotten
    plt.figure(figsize=(10, 6))
    #plt.plot(y_values, x_values, '-', label="Referenzlösung", color='C0')
    plt.plot(N1_values, Fy_values, '--', label="Eigenentwicklung", color='C1')
    plt.title("Kraft-Verschiebungsverlauf des Scheitelpunktes")
    plt.xlabel('Verschiebung in [m]')
    plt.ylabel('vertikale Kraft am Knoten in [N]')
    plt.xlim(0, 0.275)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()



else:

    Knotenliste = copy.deepcopy(nodes)
    Stabliste = copy.deepcopy(elements)

    Instanz2 = Fachwerk(Stabliste, Knotenliste)

    Fy_values = np.linspace(14000, 15000, 1)

    for Fyy in Fy_values:
        for Knoten in nodes:
            if Knoten.ID == 9:
                Knoten.Fz = - Fyy
        instanz.loesen()
        visualize_structure(instanz, Instanz2)
        #instanz.Knoten_zuruecksetzten()
