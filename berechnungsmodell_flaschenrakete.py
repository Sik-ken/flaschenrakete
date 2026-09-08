"""
Berechnungsmodell der Flaschenrakete
====================================

Nulldimensionales Kammermodell der ventilierten Deflagration einer PET-Flasche
mit Alkohol-Luft-Gemisch, mit anschliessender Schub- und Flugrechnung.

Projekt : Funkzuendmodul fuer Methanolrakete und Wasserstoffballon
Autor   : Lukas Sikken
Zweck   : Bestimmung des optimalen Duesendurchmessers und der sicherheits-
          technischen Druckobergrenze; Validierung gegen die Messwerte der
          eigenen Erprobung.

Aufbau des Skripts
------------------
  Abschnitt 1  Stoffdaten und Parameter (alle Eingabewerte an einer Stelle)
  Abschnitt 2  Stoechiometrie und Dosierung
  Abschnitt 3  Isochore Druckobergrenze (geschlossene Flasche)
  Abschnitt 4  Instationaeres Kammermodell (Zeitschrittintegration)
  Abschnitt 5  Schub, Impuls und Brennschlussgeschwindigkeit
  Abschnitt 6  Flugmechanik mit Luftwiderstand
  Abschnitt 7  Duesenstudie und Kalibrierung
  Abschnitt 8  Ausgabe

Ausfuehrung:  python3 modell.py
Abhaengigkeit: nur die Standardbibliothek.
"""

import math
import csv
import os

# Verzeichnis dieses Skripts (unabhaengig vom Arbeitsverzeichnis)
SKRIPT_VERZEICHNIS = os.path.dirname(os.path.abspath(__file__))

# =====================================================================
# 1  STOFFDATEN UND PARAMETER
# =====================================================================
# Jeder Wert ist mit seiner Herkunft gekennzeichnet:
#   [S] Stoffwert aus der Literatur
#   [G] eigene Messung am Bauteil
#   [A] Annahme
#   [K] an der eigenen Messung kalibriert

# --- Brennstoff: Isopropanol C3H8O -----------------------------------
M_B      = 60.10e-3    # kg/mol   Molmasse Brennstoff                    [S]
RHO_B    = 786.0       # kg/m^3   Dichte fluessig bei 20 C               [S]
H_U      = 30.45e6     # J/kg     unterer Heizwert                       [S]
N_O2     = 4.5         # mol/mol  Sauerstoffbedarf je mol Brennstoff     [S]
N_PROD   = 7.0         # mol/mol  gasfoermige Produkte (3 CO2 + 4 H2O)   [S]
T_AD     = 2509.0      # K        adiabate Verbrennungstemperatur        [S]

# --- Luft und Umgebung ------------------------------------------------
M_L      = 28.96e-3    # kg/mol   Molmasse Luft                          [S]
Y_O2     = 0.2095      # -        Sauerstoffanteil der Luft              [S]
R_UNIV   = 8.314       # J/(mol K) universelle Gaskonstante              [S]
R_S      = 287.0       # J/(kg K) spez. Gaskonstante des Gemisches       [A]
KAPPA    = 1.25        # -        Isentropenexponent heisses Gas         [A]
CV       = R_S / (KAPPA - 1.0)   # J/(kg K)  berechnet
CP       = CV + R_S              # J/(kg K)  berechnet
T_0      = 293.0       # K        Umgebungstemperatur                    [G]
P_0      = 1.013e5     # Pa       Umgebungsdruck                         [G]
RHO_L    = 1.20        # kg/m^3   Luftdichte                             [S]
G        = 9.81        # m/s^2    Erdbeschleunigung                      [S]

# --- Modellparameter --------------------------------------------------
S_L      = 0.40        # m/s      laminare Brenngeschwindigkeit          [S]
CHI      = 2.00        # -        Flammenfaltungsfaktor                  [K]
SIGMA    = 7.50        # -        Expansionsverhaeltnis                  [S]
C_D_HOLE = 0.80        # -        Durchflusszahl der Bohrung             [A]
H_WALL   = 150.0       # W/(m^2 K) Wandwaermeuebergang                   [K]
T_WALL   = 293.0       # K        Wandtemperatur                         [A]
ETA_C    = 0.95        # -        Verbrennungswirkungsgrad               [A]

# --- Flasche und Flug -------------------------------------------------
V_FL     = 1.0e-3      # m^3      Flaschenvolumen                        [G]
D_FL     = 0.088       # m        Flaschendurchmesser                    [G]
M_FL     = 0.032       # kg       Leermasse der Flasche                  [G]
C_D_AIR  = 0.60        # -        Widerstandsbeiwert                     [A]
P_BURST  = 14.0e5      # Pa       Mindestberstdruck                      [S] Riedel et al. 2025
S_MIN    = 2.0         # -        geforderter Sicherheitsfaktor          [A]

# --- Numerik ----------------------------------------------------------
DT       = 1.0e-6      # s        Zeitschritt
T_MAX    = 1.0         # s        Abbruchzeit


# =====================================================================
# 2  STOECHIOMETRIE UND DOSIERUNG
# =====================================================================
def stoechiometrie():
    """Stoechiometrische Kenngroessen des Brennstoff-Luft-Gemisches."""
    n_L   = N_O2 / Y_O2                     # mol Luft je mol Brennstoff
    phi   = 1.0 / (1.0 + n_L)               # Volumenanteil Brennstoff
    L_st  = n_L * M_L / M_B                 # kg Luft je kg Brennstoff
    H_G   = H_U / (1.0 + L_st)              # Gemischheizwert J/kg
    V_m   = R_UNIV * T_0 / P_0              # molares Volumen bei 20 C, m^3/mol
    m_B   = phi / V_m * M_B                 # kg Brennstoff je m^3 Gemisch
    return dict(n_L=n_L, phi=phi, L_st=L_st, H_G=H_G, V_m=V_m,
                m_B_pro_m3=m_B,
                m_B_flasche=m_B * V_FL,
                V_B_flasche=m_B * V_FL / RHO_B)


# =====================================================================
# 3  ISOCHORE DRUCKOBERGRENZE
# =====================================================================
def grenzdruck(st):
    """Druck bei vollstaendiger Verbrennung in der geschlossenen Flasche."""
    n1    = 1.0 + st["n_L"]
    n2    = N_PROD + st["n_L"] * (1.0 - Y_O2)
    p2    = P_0 * (n2 / n1) * (T_AD / T_0)
    return dict(n1=n1, n2=n2, verhaeltnis=n2 / n1, p_isochor=p2,
                sicherheitsfaktor=P_BURST / p2,
                p_zulaessig=P_BURST / S_MIN)


# =====================================================================
# 4  INSTATIONAERES KAMMERMODELL
# =====================================================================
# Bilanzraum: Flascheninnenraum, konstantes Volumen, ideal durchmischt.
#
#   Massenbilanz    dm/dt  = -m_aus
#   Energiebilanz   d(m cv T)/dt = m_verbr * H_G * eta_c
#                                  - m_aus * cp * T
#                                  - h_W * A_W * (T - T_Wand)
#   Zustandsgl.     p = m * R_s * T / V
#
# Verbrennungsrate: ebene Flammenfront ueber den Flaschenquerschnitt,
# die Vergroesserung der Frontflaeche durch Turbulenz ist im Faltungs-
# faktor CHI enthalten.
#
#   m_verbr = rho_u * CHI * S_L * A_quer
#
# Ausstroemung durch die Bohrung, unterschieden nach ueberkritischem und
# unterkritischem Druckverhaeltnis.

def ausflussfunktion(p_innen, p_aussen):
    """Ausflussfunktion Gamma nach Saint-Venant und Wantzel."""
    pr_krit = ((KAPPA + 1.0) / 2.0) ** (KAPPA / (KAPPA - 1.0))
    if p_innen / p_aussen >= pr_krit:                       # ueberkritisch
        return math.sqrt(KAPPA) * (2.0 / (KAPPA + 1.0)) ** \
               ((KAPPA + 1.0) / (2.0 * (KAPPA - 1.0))), True
    pr = p_aussen / p_innen                                 # unterkritisch
    term = pr ** (2.0 / KAPPA) - pr ** ((KAPPA + 1.0) / KAPPA)
    return math.sqrt(2.0 * KAPPA / (KAPPA - 1.0) * max(term, 0.0)), False


def kammermodell(d_bohrung, chi=CHI, verlauf=False):
    """Zeitschrittintegration fuer einen Bohrungsdurchmesser in Metern."""
    st = stoechiometrie()
    H_G = st["H_G"]

    A_quer  = math.pi / 4.0 * D_FL ** 2          # Flaschenquerschnitt
    L_FL    = V_FL / A_quer                      # aequivalente Laenge
    A_hole  = math.pi / 4.0 * d_bohrung ** 2     # Bohrungsflaeche
    A_wall  = math.pi * D_FL * L_FL + 2 * A_quer # Waermeuebertragende Flaeche

    rho_0   = P_0 / (R_S * T_0)                  # Anfangsdichte des Gemisches
    m       = rho_0 * V_FL                       # Masse in der Flasche
    T       = T_0
    p       = P_0
    m_frisch = m                                 # noch unverbrannte Masse
    m_verbrannt_gesamt = 0.0

    t = 0.0
    impuls = 0.0
    p_max = P_0
    F_max = 0.0
    verlaufsdaten = []

    while t < T_MAX:
        rho_u = p / (R_S * T_0)                  # Dichte des unverbrannten Gases
        rate_verbr = rho_u * chi * S_L * A_quer  # kg/s
        rate_verbr = min(rate_verbr, m_frisch / DT if m_frisch > 0 else 0.0)

        gamma, kritisch = ausflussfunktion(p, P_0)
        rate_aus = C_D_HOLE * A_hole * p * gamma / math.sqrt(R_S * T)

        q_wand = H_WALL * A_wall * (T - T_WALL)

        dU = (rate_verbr * H_G * ETA_C
              - rate_aus * CP * T
              - q_wand) * DT
        dm = -rate_aus * DT

        U_neu = m * CV * T + dU
        m_neu = m + dm
        if m_neu <= 1e-9:
            break
        T = U_neu / (m_neu * CV)
        T = max(T, T_0)
        m = m_neu
        p = m * R_S * T / V_FL
        p = max(p, P_0)

        m_frisch = max(m_frisch - rate_verbr * DT + dm * (m_frisch / m if m > 0 else 0), 0.0)
        m_verbrannt_gesamt += rate_verbr * DT

        # Schub aus Impulssatz. Bei ueberkritischem Druckverhaeltnis liegt
        # im Muendungsquerschnitt Schallgeschwindigkeit an.
        if kritisch:
            T_hals = T * 2.0 / (KAPPA + 1.0)
            v_aus = math.sqrt(KAPPA * R_S * T_hals)
            p_hals = p * (2.0 / (KAPPA + 1.0)) ** (KAPPA / (KAPPA - 1.0))
        else:
            v_aus = math.sqrt(max(2.0 * CP * T * (1.0 - (P_0 / p) **
                                                  ((KAPPA - 1.0) / KAPPA)), 0.0))
            p_hals = P_0
        F = rate_aus * v_aus + (p_hals - P_0) * A_hole

        impuls += F * DT
        p_max = max(p_max, p)
        F_max = max(F_max, F)
        if verlauf and int(t / DT) % 200 == 0:
            verlaufsdaten.append((t * 1e3, p / 1e5, F))

        t += DT
        if m_frisch <= 1e-9 and p <= P_0 * 1.001:
            break

    return dict(d=d_bohrung, p_max=p_max, F_max=F_max, impuls=impuls,
                t_ende=t, verlauf=verlaufsdaten)


# =====================================================================
# 5  BRENNSCHLUSSGESCHWINDIGKEIT
# =====================================================================
def brennschluss(impuls, t_brenn):
    """Geschwindigkeit am Brennschluss aus Impuls abzueglich Schwerkraft."""
    m_ges = M_FL
    v = impuls / m_ges - G * t_brenn
    return max(v, 0.0)


# =====================================================================
# 6  FLUGMECHANIK MIT LUFTWIDERSTAND
# =====================================================================
def grenzgeschwindigkeit():
    A_stirn = math.pi / 4.0 * D_FL ** 2
    return math.sqrt(2.0 * M_FL * G / (RHO_L * C_D_AIR * A_stirn))


def steighoehe(v_bo):
    """Steighoehe bei quadratischem Luftwiderstand."""
    v_t = grenzgeschwindigkeit()
    if v_bo <= 0.0:
        return 0.0
    return (v_t ** 2 / (2.0 * G)) * math.log(1.0 + v_bo ** 2 / v_t ** 2)


def steighoehe_ohne_widerstand(v_bo):
    return v_bo ** 2 / (2.0 * G)


# =====================================================================
# 7  DUESENSTUDIE UND KALIBRIERUNG
# =====================================================================
def kritischer_durchmesser(chi=CHI, T_b=2400.0):
    """Geschlossene Abschaetzung des kritischen Bohrungsdurchmessers."""
    gamma_krit, _ = ausflussfunktion(10.0 * P_0, P_0)
    nenner = C_D_HOLE * gamma_krit * math.sqrt(R_S * T_b)
    zaehler = chi * S_L * (SIGMA - 1.0)
    return D_FL * math.sqrt(zaehler / nenner)


def duesenstudie(durchmesser_mm, chi=CHI):
    ergebnisse = []
    for d_mm in durchmesser_mm:
        r = kammermodell(d_mm / 1000.0, chi=chi)
        v_bo = brennschluss(r["impuls"], r["t_ende"])
        ergebnisse.append(dict(
            d_mm=d_mm,
            p_max_bar=r["p_max"] / 1e5,
            F_max_N=r["F_max"],
            impuls_Ns=r["impuls"],
            v_bo=v_bo,
            hoehe_m=steighoehe(v_bo),
            sicherheitsfaktor=P_BURST / r["p_max"],
        ))
    return ergebnisse


def kalibrierung(chi_werte, durchmesser_mm):
    """Optimaler Durchmesser je Flammenfaltungsfaktor."""
    out = []
    for chi in chi_werte:
        erg = duesenstudie(durchmesser_mm, chi=chi)
        best = max(erg, key=lambda e: e["hoehe_m"])
        out.append(dict(chi=chi, d_opt_mm=best["d_mm"],
                        hoehe_m=best["hoehe_m"],
                        d_krit_mm=kritischer_durchmesser(chi) * 1000.0))
    return out


# =====================================================================
# 8  AUSGABE
# =====================================================================
def linie(z="-", n=78):
    print(z * n)


def main():
    st = stoechiometrie()
    gd = grenzdruck(st)

    linie("=")
    print("BERECHNUNGSMODELL FLASCHENRAKETE")
    linie("=")

    print("\n1  STOECHIOMETRIE UND DOSIERUNG")
    linie()
    print(f"  Luftbedarf                        {st['n_L']:10.2f} mol/mol")
    print(f"  Stoechiometrischer Volumenanteil  {st['phi']*100:10.2f} Vol.-%")
    print(f"  Stoechiometrisches Luftverhaeltnis{st['L_st']:10.2f} kg/kg")
    print(f"  Gemischheizwert                   {st['H_G']/1e6:10.3f} MJ/kg")
    print(f"  Brennstoffmasse je Flasche        {st['m_B_flasche']*1e3:10.4f} g")
    print(f"  Brennstoffvolumen je Flasche      {st['V_B_flasche']*1e6:10.4f} ml")

    print("\n2  ISOCHORE DRUCKOBERGRENZE (geschlossene Flasche)")
    linie()
    print(f"  Stoffmengenverhaeltnis n2/n1      {gd['verhaeltnis']:10.4f} -")
    print(f"  Isochorer Grenzdruck              {gd['p_isochor']/1e5:10.2f} bar")
    print(f"  Mindestberstdruck (Literatur)     {P_BURST/1e5:10.2f} bar")
    print(f"  Sicherheitsfaktor                 {gd['sicherheitsfaktor']:10.2f} -")
    print(f"  Geforderter Sicherheitsfaktor     {S_MIN:10.2f} -")
    print("  Bewertung: " + ("zulaessig" if gd['sicherheitsfaktor'] >= S_MIN
          else "unzulaessig, Duesenbohrung zwingend erforderlich"))

    print("\n3  FLUGMECHANISCHE KENNGROESSEN")
    linie()
    print(f"  Stirnflaeche                      {math.pi/4*D_FL**2*1e4:10.2f} cm^2")
    print(f"  Grenzgeschwindigkeit              {grenzgeschwindigkeit():10.2f} m/s")

    print("\n4  DUESENSTUDIE (chi = %.2f)" % CHI)
    linie()
    kopf = (f"  {'d':>5} {'p_max':>8} {'F_max':>8} {'Impuls':>9} "
            f"{'v_bo':>8} {'Hoehe':>8} {'Sicherh.':>9}")
    einh = (f"  {'[mm]':>5} {'[bar]':>8} {'[N]':>8} {'[N s]':>9} "
            f"{'[m/s]':>8} {'[m]':>8} {'[-]':>9}")
    print(kopf); print(einh); linie()
    durchmesser = [4, 5, 6, 7, 8, 10, 12, 16]
    erg = duesenstudie(durchmesser)
    for e in erg:
        print(f"  {e['d_mm']:5.0f} {e['p_max_bar']:8.2f} {e['F_max_N']:8.2f} "
              f"{e['impuls_Ns']:9.3f} {e['v_bo']:8.1f} {e['hoehe_m']:8.2f} "
              f"{e['sicherheitsfaktor']:9.1f}")
    best = max(erg, key=lambda e: e["hoehe_m"])
    print(f"\n  Optimum bei d = {best['d_mm']:.0f} mm, Hoehe {best['hoehe_m']:.2f} m, "
          f"Spitzendruck {best['p_max_bar']:.2f} bar")

    print("\n5  KALIBRIERUNG DES FLAMMENFALTUNGSFAKTORS")
    linie()
    print(f"  {'chi':>6} {'d_opt Modell':>14} {'d_krit':>10} {'Hoehe':>9}")
    print(f"  {'[-]':>6} {'[mm]':>14} {'[mm]':>10} {'[m]':>9}")
    linie()
    for k in kalibrierung([1.5, 1.75, 2.0, 2.5, 3.0], durchmesser):
        print(f"  {k['chi']:6.2f} {k['d_opt_mm']:14.0f} "
              f"{k['d_krit_mm']:10.2f} {k['hoehe_m']:9.2f}")


    csv_pfad = os.path.join(SKRIPT_VERZEICHNIS, "duesenstudie.csv")
    with open(csv_pfad, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(erg[0].keys()), delimiter=";")
        w.writeheader()
        for e in erg:
            w.writerow(e)
    print(f"\n  Ergebnisse zusaetzlich gespeichert in {csv_pfad}")
    linie("=")


if __name__ == "__main__":
    main()
