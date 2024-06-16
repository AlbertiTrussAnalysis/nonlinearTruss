# Geometrisch Nichtlineare Fachwerkstäbe

Die hier zur Verfügung gestellten Codes wurden im Rahmen meiner Bachelorarbeit an der Technischen Universität München am Lehrstuhl für Statik entwickelt.

Die "Hauptanwendung" ermöglicht die lineare und nichtlineare Analyse beliebiger 2D-Fachwerke. Dabei kann sowohl die Konvergenz als auch die Veränderung der Geometrie bei Laststeigerung untersucht werden.

Mit "Balken" kann ein Dreigelenksystem sowohl linear als auch nichtlinear kraft- und weggesteuert berechnet werden. Die erhaltenen Ergebnisse werden dann mit der analytischen Lösung verglichen.

Die folgenden Codes setzen voraus, dass sich die äquivalenten Koordinatendateien im gleichen Verzeichnis wie die Pyhton-Dateien befinden. 

Der "Fachwerkbogen" ist ein weiterer Benchmark, bei dem die verschiebungsgesteuerte Berechnungsmethode komplexer ist als beim Dreigelenksystem.

In "Raumfachwerk" können beide Berechnungsmethoden an einem dreidimensionalen Tragwerk nachvollzogen werden.

Die Gültigkeit der Ergebnisse wurde in der Arbeit nachgewiesen.

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



