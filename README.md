# Berechnungsmodell der Flaschenrakete

Nulldimensionales Kammermodell der ventilierten Deflagration einer PET-Flasche
mit Alkohol-Luft-Gemisch, mit anschliessender Schub- und Flugrechnung.

Das Modell bestimmt den optimalen Duesendurchmesser und die
sicherheitstechnische Druckobergrenze und validiert diese gegen die Messwerte der
eigenen Erprobung.

## Aufbau

| Abschnitt | Inhalt |
|-----------|--------|
| 1 | Stoffdaten und Parameter (alle Eingabewerte an einer Stelle) |
| 2 | Stoechiometrie und Dosierung |
| 3 | Isochore Druckobergrenze (geschlossene Flasche) |
| 4 | Instationaeres Kammermodell (Zeitschrittintegration) |
| 5 | Schub, Impuls und Brennschlussgeschwindigkeit |
| 6 | Flugmechanik mit Luftwiderstand |
| 7 | Duesenstudie und Kalibrierung |
| 8 | Ausgabe |

## Ausfuehrung

```bash
python3 berechnungsmodell_flaschenrakete.py
```

Abhaengigkeiten: nur die Python-Standardbibliothek.

Die Duesenstudie wird zusaetzlich als `duesenstudie.csv` neben dem Skript abgelegt.

## Autor
Sikken


## Beispielberechnung
==============================================================================
BERECHNUNGSMODELL FLASCHENRAKETE
==============================================================================

1  STOECHIOMETRIE UND DOSIERUNG
------------------------------------------------------------------------------
  Luftbedarf                             21.48 mol/mol
  Stoechiometrischer Volumenanteil        4.45 Vol.-%
  Stoechiometrisches Luftverhaeltnis     10.35 kg/kg
  Gemischheizwert                        2.683 MJ/kg
  Brennstoffmasse je Flasche            0.1112 g
  Brennstoffvolumen je Flasche          0.1414 ml

2  ISOCHORE DRUCKOBERGRENZE (geschlossene Flasche)
------------------------------------------------------------------------------
  Stoffmengenverhaeltnis n2/n1          1.0667 -
  Isochorer Grenzdruck                    9.25 bar
  Mindestberstdruck (Literatur)          14.00 bar
  Sicherheitsfaktor                       1.51 -
  Geforderter Sicherheitsfaktor           2.00 -
  Bewertung: unzulaessig, Duesenbohrung zwingend erforderlich

3  FLUGMECHANISCHE KENNGROESSEN
------------------------------------------------------------------------------
  Stirnflaeche                           60.82 cm^2
  Grenzgeschwindigkeit                   11.97 m/s

4  DUESENSTUDIE (chi = 2.00)
------------------------------------------------------------------------------
      d    p_max    F_max    Impuls     v_bo    Hoehe  Sicherh.
   [mm]    [bar]      [N]     [N s]    [m/s]      [m]       [-]
------------------------------------------------------------------------------
      4     6.09     7.22     0.551     15.2     6.99       2.3
      5     5.29     9.54     0.603     17.3     8.21       2.6
      6     4.44    11.06     0.614     17.9     8.56       3.2
      7     3.59    11.42     0.592     17.4     8.30       3.9
      8     2.78    10.42     0.546     16.1     7.56       5.0
     10     1.78     7.63     0.413     12.1     5.16       7.8
     12     1.41     6.30     0.284      8.1     2.78       9.9
     16     1.16     4.48     0.151      4.0     0.76      12.1

  Optimum bei d = 6 mm, Hoehe 8.56 m, Spitzendruck 4.44 bar

5  KALIBRIERUNG DES FLAMMENFALTUNGSFAKTORS
------------------------------------------------------------------------------
     chi   d_opt Modell     d_krit     Hoehe
     [-]           [mm]       [mm]       [m]
------------------------------------------------------------------------------
    1.50              5       8.31      7.26
    1.75              6       8.98      7.99
    2.00              6       9.60      8.56
    2.50              7      10.73      9.36
    3.00              7      11.76     10.06
 Ergebnisse zusaetzlich gespeichert in ------------------------
 ==============================================================================
 

