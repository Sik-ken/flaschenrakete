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
