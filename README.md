# Berechnungsmodell Flaschenrakete

Nulldimensionales Kammermodell der ventilierten Deflagration einer PET-Flasche
mit Alkohol/Luft-Gemisch, mit anschließender Schub- und Flugrechnung.

**Projekt:** Funkzündmodul für Methanolrakete und Wasserstoffballon
**Datei:** `berechnungsmodell_flaschenrakete.py` (466 Zeilen)
**Zweck:** Bestimmung des optimalen Düsendurchmessers und der sicherheits­technischen
Druckobergrenze; Validierung gegen die Messwerte der eigenen Erprobung.

---

## 1. Ausführung

```bash
python3 berechnungsmodell_flaschenrakete.py
```

Getestet mit Python 3.12. **Keine externen Abhängigkeiten** — es werden nur
`math`, `csv` und `os` aus der Standardbibliothek verwendet. Das Skript ist
damit ohne Installation von Fremdpaketen lauffähig und auch in mehreren Jahren
noch reproduzierbar.

Laufzeit: wenige Sekunden.

### Ausgaben

1. **Konsole** — fünf Abschnitte:
   Stöchiometrie und Dosierung, isochore Druckobergrenze, flugmechanische
   Kenngrößen, Düsenstudie, Kalibrierung des Flammenfaltungsfaktors
2. **CSV-Datei** `duesenstudie_<brennstoff>.csv` im Skriptverzeichnis,
   Semikolon-getrennt (Excel-kompatibel), Spalten:
   `d_mm; p_max_bar; F_max_N; impuls_Ns; v_bo; hoehe_m; sicherheitsfaktor`

Der Dateiname enthält den Brennstoffnamen, damit sich die Läufe für
Isopropanol und Methanol nicht gegenseitig überschreiben.

---

## 2. Brennstoff wechseln

Im Parameterblock am Skriptanfang stehen zwei Brennstoffblöcke untereinander.
Genau einer ist aktiv, der andere ist auskommentiert. Zum Wechseln den aktiven
Block auskommentieren und den anderen einkommentieren — sonst ist nichts zu
ändern, alle abgeleiteten Größen berechnen sich neu.

| Parameter | Isopropanol C₃H₈O | Methanol CH₄O | Herkunft |
|---|---|---|---|
| `M_B` molare Masse | 60,10·10⁻³ kg/mol | 32,04·10⁻³ kg/mol | Literatur |
| `RHO_B` Dichte flüssig (20 °C) | 786 kg/m³ | 792 kg/m³ | Literatur |
| `H_U` unterer Heizwert | 30,45 MJ/kg | 19,90 MJ/kg | Literatur |
| `N_O2` Sauerstoffbedarf | 4,5 mol/mol | 1,5 mol/mol | ν = c + h/4 − o/2 |
| `N_PROD` gasförmige Produkte | 7,0 (3 CO₂ + 4 H₂O) | 3,0 (1 CO₂ + 2 H₂O) | Reaktionsgleichung |
| `S_L` lam. Brenngeschwindigkeit | 0,40 m/s | 0,45 m/s | Literatur |
| `SIGMA` Expansionsverhältnis | 8,10 | 7,80 | berechnet |
| `T_AD` isochore Flammentemperatur | 2595,77 K | 2570,36 K | externes Tool |
| `P_ISO_EXT` isochorer Grenzdruck | 9,02 bar | 9,18 bar | externes Tool |

> **Wichtig zu `T_AD`:** Es muss die **isochore** adiabate Flammentemperatur
> sein, nicht die isobare. In der Beziehung p₂/p₁ = (n₂/n₁)·(T₂/T₁) steht
> zwingend die isochore Temperatur; der isobare Wert liegt rund 15 % niedriger
> und würde den Grenzdruck systematisch unterschätzen. Das Skript gibt zur
> Kontrolle den aus der eigenen Energiebilanz berechneten Wert daneben aus —
> weicht der externe Wert um mehr als etwa 10 % nach unten ab, ist es der
> isobare.

### Herkunft der externen Werte: NASA CEA (CEARUN)

`T_AD` und `P_ISO_EXT` stammen aus **NASA CEA** (*Chemical Equilibrium with
Applications*), bedient über die Weboberfläche **CEARUN**:
<https://cearun.grc.nasa.gov/>

CEA berechnet die Gleichgewichtszusammensetzung einer Verbrennung
**inklusive Dissoziation** (bei rund 2600 K zerfällt ein merklicher Anteil
der Produkte in CO, OH, H, O, NO …). Die einfache Molbilanz im Skript
unterstellt dagegen eingefrorene Produkte (3 CO₂ + 4 H₂O + N₂). Deshalb
weichen `p_modell` und `p_extern` in Abschnitt 3 voneinander ab.

Für den Grenzdruck der **geschlossenen** Flasche wird der Problemtyp
**`uv`** gewählt (*assigned internal energy and volume*) — die isochore
adiabate Verbrennung bei konstantem Volumen. Der Kammerdruck ist dabei
Ergebnis der Rechnung, nicht Vorgabe.

**Eingaben** (für beide Brennstoffe identisch bis auf den Fuel):

| Größe | Wert |
|---|---|
| Fuel | `CH3OH` (Methanol) bzw. `C3H8O,2propanol` (Isopropanol) |
| Oxidant | `Air` |
| Äquivalenzverhältnis φ | 1,0 (stöchiometrisch) |
| Anfangstemperatur | 298,15 K |
| Anfangsdruck | 1,013 bar |

**Ergebnisse** (Gleichgewicht mit Dissoziation, Problemtyp `uv`):

| Brennstoff | O/F | %FUEL | p isochor → `P_ISO_EXT` | T isochor → `T_AD` |
|---|---|---|---|---|
| Methanol CH₃OH | 6,473 | 13,381 | 9,1842 bar → `9.18e5` | 2570,36 K → `2570.36` |
| Isopropanol C₃H₈O | 10,354 | 8,807 | 9,0204 bar → `9.02e5` | 2595,77 K → `2595.77` |

Auszug der CEA-Ausgabe — Methanol:

```
 FUEL        CH3OH                1.0000000   -203418.971    298.150
 OXIDANT     Air                  1.0000000     -2604.501    298.150

 O/F=    6.47313  %FUEL= 13.381268  R,EQ.RATIO= 1.000000  PHI,EQ.RATIO= 1.000000

 THERMODYNAMIC PROPERTIES

 P, BAR            9.1842
 T, K            2570.36
```

Isopropanol:

```
 FUEL        C3H8O,2propanol      1.0000000   -275178.971    298.150
 OXIDANT     Air                  1.0000000     -2604.501    298.150

 O/F=   10.35416  %FUEL=  8.807342  R,EQ.RATIO= 1.000000  PHI,EQ.RATIO= 1.000000

 THERMODYNAMIC PROPERTIES

 P, BAR            9.0204
 T, K            2595.77
```

---

## 3. Aufbau des Skripts

| Abschnitt | Inhalt | Kernfunktionen |
|---|---|---|
| 1 | Stoffdaten und Parameter (alle Eingabewerte an einer Stelle) | — |
| 2 | Stöchiometrie und Dosierung | `stoechiometrie()` |
| 3 | Isochore Druckobergrenze | `isochore_flammentemperatur()`, `grenzdruck()` |
| 4 | Instationäres Kammermodell | `ausflussfunktion()`, `kammermodell()` |
| 5 | Brennschlussgeschwindigkeit | `brennschluss()` |
| 6 | Flugmechanik mit Luftwiderstand | `grenzgeschwindigkeit()`, `steighoehe()` |
| 7 | Düsenstudie und Kalibrierung | `kritischer_durchmesser()`, `duesenstudie()`, `kalibrierung()` |
| 8 | Ausgabe | `main()` |

Jeder Eingabewert ist mit seiner Herkunft gekennzeichnet:

| Kürzel | Bedeutung |
|---|---|
| `[S]` | Stoffwert aus der Literatur |
| `[G]` | eigene Messung am Bauteil |
| `[A]` | Annahme |
| `[K]` | an der eigenen Messung kalibriert |

---

## 4. Kalibrierung

Das Modell hat **zwei freie Parameter**:

| Parameter | Wert | Wirkung |
|---|---|---|
| `CHI` Flammenfaltungsfaktor | 2,25 | Brenndauer → **Lage** des Optimums |
| `H_WALL` Wandwärmeübergang | 100 W/(m²K) | Wärmeverlust → **Höhenmaßstab** |

Beide wurden gemeinsam an der eigenen Erprobung kalibriert (1,0-L-Flasche,
Isopropanol, 20 °C):

| Messgröße | Messung | Modell |
|---|---|---|
| optimaler Düsendurchmesser | 6 mm | 6 mm |
| Steighöhe | ≥ 8 m | 8,21 m |

Dass **zwei unabhängige Messgrößen mit zwei Parametern** getroffen werden,
und zwar nicht beliebig (χ verschiebt fast nur das Optimum, H_WALL fast nur
den Maßstab), ist der eigentliche Validierungsnachweis des Modells.

**Bei Änderung von `L_FL` muss `H_WALL` nachgeführt werden**, damit das
kalibrierte Produkt H_WALL · A_wall erhalten bleibt.

### Zeitschrittkonvergenz

| `DT` | Optimum | Höhe | p_max |
|---|---|---|---|
| 1·10⁻⁶ s | 6 mm | 8,21 m | 3,88 bar |
| 5·10⁻⁶ s | 6 mm | 8,21 m | 3,88 bar |
| **1·10⁻⁵ s** (eingestellt) | 6 mm | 8,21 m | 3,88 bar |
| 5·10⁻⁵ s | 6 mm | 8,21 m | 3,88 bar |

Über zwei Größenordnungen identisch bis zur zweiten Nachkommastelle. Der
gewählte Zeitschritt ist konvergiert.



## 5. Änderungsprotokoll gegenüber dem Ausgangsskript

### Version 2 — Korrekturen und Erweiterungen

#### Fehlerkorrekturen

| # | Stelle | Problem | Korrektur |
|---|---|---|---|
| 1 | `T_AD` |  Wert korrigiert und auf die **isochore** Flammentemperatur 2595,77 K aus externer Gleichgewichtsrechnung umgestellt |
| 2 | Frischgasdichte | `rho_u = p / (R_S * T_0)` hielt T_u = T_0 fest, entspricht **isothermer** Verdichtung; bei 4,4 bar 35 % zu hohe Frischgasdichte | isentrope Verdichtung `rho_u = rho_0 * (p/P_0)**(1/KAPPA_U)` mit eigenem Isentropenexponenten für das kalte Frischgas |
| 3 | `KAPPA` | für das Frischgas wurde derselbe Exponent 1,25 wie für das heiße Abgas verwendet | `KAPPA_U = 1.38` als eigener Parameter ergänzt |
| 4 | `A_wall` | aus der äquivalenten Länge `V_FL/A_quer` = 164 mm berechnet; reale Innenhöhe 270 mm, Fläche dadurch um Faktor 1,51 zu klein | `L_FL = 0.270 m` als gemessener Eingabewert; `H_WALL` entsprechend von 150 auf 100 W/(m²K) nachgeführt |
| 5 | `m_frisch` | Massenanteil wurde mit `m` gebildet, das drei Zeilen vorher bereits auf `m_neu` überschrieben worden war | `m_alt` vor dem Überschreiben gesichert |
| 6 | `rate_aus` | nicht auf die verfügbare Masse begrenzt | `min(rate_aus, m/DT)` ergänzt |
| 7 | Protokollierung | `int(t/DT) % 200` — durch Gleitkomma-Akkumulation von `t` konnten Stichproben übersprungen werden | eigener Schrittzähler `schritt` |
| 8 | `kritischer_durchmesser` | `T_b = 2400.0` hart kodiert, lief beim Brennstoffwechsel gegen `T_AD` auseinander | `T_b=None` koppelt auf `T_AD` |

**Auswirkung von Korrektur 2 und 4:** Spitzendruck bei 6 mm sinkt von 4,44 auf
3,88 bar. Die **Optimumslage bleibt bei 6 mm** — die Kalibrierung ist gegen
beide Korrekturen robust.

#### Erweiterungen

| # | Erweiterung |
|---|---|
| 9 | **Brennstoffauswahl.** Schaltbare Blöcke für Isopropanol und Methanol; `BRENNSTOFF` erscheint in Kopfzeile, Düsenstudie und CSV-Dateinamen |
| 10 | **Zwei Wege zum Grenzdruck.** `grenzdruck()` gibt den Wert aus der Molbilanz mit `T_AD` **und** den extern ermittelten `P_ISO_EXT` aus, dazu die prozentuale Abweichung. Für die Auslegung wird konservativ der größere Wert verwendet |
| 11 | **Plausibilitätsprüfung.** `isochore_flammentemperatur()` berechnet T_ad,V aus derselben Energiebilanz wie das Kammermodell und stellt sie dem externen Wert gegenüber |
| 12 | `SIGMA` und `S_L` in den Brennstoffblock verschoben (sind Stoffgrößen, keine Modellparameter) |
| 13 | `DT` von 1·10⁻⁶ auf 1·10⁻⁵ s (konvergiert, rund zehnfache Rechengeschwindigkeit) |
| 14 | Nachkalibrierung `CHI` 2,00 → 2,25 nach den Korrekturen 2 und 4 |
| 15 | CSV-Dateiname enthält den Brennstoff |
