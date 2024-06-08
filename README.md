# Geometrisch Nichtlineare Fachwerkstäbe

Die hier bereitgestellten Codes wurden im Rahmen meiner Bachelorarbeit an der Technischen Universität München am Lehrstuhl für Statik entwickelt.

Das "Hauptanwendung" ermöglicht die lineare und nichtlineare Analyse von beliebigen 2D-Fachwerken. Hierbei kann sowohl die Konvergenz, als auch Veränderungen in der Geometrie im Zuge der Laststeigerung betrachtet werden.

Mit "Balken" kann ein Dreigelenktragwerk sowohl linear als auch nicht-linear berechnet werden, und zwar kraft- und verschiebungsgesteuert. Die erhaltenen Ergebnisse werden im Anschluss mit der analytischen Lösung abgeglichen.

Folgende Codes benötigen es, dass die äquivalenten Koordinatendateien sich im gleichen Ordner wie die Pyhton Dateien befinden. 

"Fachwerkbogen" ist ein weiterer Benchmark, bei welcher die Verschiebungsgesteuerte Rechenmethode aufwendiger ist als bei dem Dreigelenktragwerk.

In "Raumfachwerk" können die beiden Berechnungsmethoden an einem dreidimensionalen Tragwerk nachvollzogen werden.

Die Gültigkeit der Ergebnisse wurden in der Arbeit nachgewiesen.

## Notwendige Python Packete

Es muss 'numpy' und 'matplotlib' installiert werden.

- Im Jupyter Notebook:
```
import sys
!{sys.executable} -m pip install numpy
!{sys.executable} -m pip install matplotlib
```

- In der Python Konsole:
```
pip install numpy matplotlib
```

## Erreichbarkeit
Bei Fragen zur Bachelorarbeit oder zu den veröffentlichen Codes stehe ich gerne zur Verfügung.

alberti.maxhof@gmx.de



