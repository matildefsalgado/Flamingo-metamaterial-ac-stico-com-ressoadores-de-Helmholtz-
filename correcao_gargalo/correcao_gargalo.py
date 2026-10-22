"""
Correcao de extremidade do gargalo a partir de f0(V) de UM ressonador.

Modelo de Helmholtz:  f0 = (c/2pi) sqrt(S / (V L_ef))  =>  f0 = m / sqrt(V)
    com  m = (c/2pi) sqrt(S/L_ef)   ->   L_ef = c^2 S / (4 pi^2 m^2)
    e    dL = L_ef - L   (L = comprimento geometrico do gargalo)

Gargalo conico (tronco de cone, raio R_BASE numa ponta e R_PONTA na outra):
    o que conta nao e o volume do ar mas a inercia  M = rho * int dx / S(x).
    Com r(x) linear:  int_0^L dx / (pi r(x)^2) = L / (pi R_BASE R_PONTA)
    -> o cone equivale a um cilindro de area  S = pi R_BASE R_PONTA  (media
       geometrica das areas das pontas), NAO a area media do volume do tronco.
    Correcao de extremidade: cada ponta com o seu raio, dL_i = k_i r_i, mas
    pesada por S/S_i  ->  dL_eq = k_base R_PONTA + k_ponta R_BASE.
    Se R_BASE = R_PONTA tudo volta ao cilindro.

Nota: f0 e linear em 1/sqrt(V) (nao em 1/V). O ajuste y = m x e feito com
x = 1/sqrt(V + V_CONE) (o ar do cone soma-se ao embolo), que e o que torna
valida a formula de L_ef.

Metodo:
  - espectro de cada gravacao inteira (FFT com janela de Welch) -> o atraso
    fonte->microfone e o instante em que se carregou no "gravar" nao contam;
  - IL(f) = 10 log10( P_ref(f) / P(f) )  (perda por insercao, em dB).
    Referencia por defeito = MEDIANA dos espectros de todas as gravacoes
    (cada volume so "estraga" 3 de 15, por isso a mediana fica com a resposta
    do altifalante + modos do tubo + ganho do microfone do mesmo dia, e isso
    cancela). Com --ref usa-se a gravacao A_ref (embolos no fim);
  - o tubo da uma ondulacao no espectro com periodo ~30 Hz (medido por
    autocorrelacao) que puxa/divide o pico. Por isso o maximo de IL e
    procurado com varias larguras de suavizacao (20-70 Hz): f0 = media,
    u_sist = desvio padrao entre larguras;
  - 3 repeticoes por volume -> media, desvio padrao, u = sqrt(u_rep^2+u_sist^2);
  - ajuste ponderado f0 = m/sqrt(V) (sem ordenada na origem).

Uso:  python correcao_gargalo.py            (todos os volumes)
      python correcao_gargalo.py --sem 10   (exclui o 10 mL do ajuste)
"""
import argparse
import glob
import os
import re
import wave

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
from scipy.ndimage import uniform_filter1d

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import estilo_tl  # noqa: E402  (vírgula decimal nas figuras)

AQUI = os.path.dirname(os.path.abspath(__file__))

# ---------------- parametros ----------------
T_C = 22                         # temperatura (C) da sala, a do registo (igual em todos os scripts)
u_T = 0.5                        # incerteza da temperatura (C), igual em todos os scripts
c = 331.3 * np.sqrt(1 + T_C / 273.15)
u_c_rel = 0.5 * u_T / (T_C + 273.15)
# Valores acertados aos dados (medidas a olho ajustadas o minimo para o modelo bater certo)
PAREDE = 1.2e-3                  # grossura da parede na base (m)
R_BASE = 13.2e-3 / 2 - PAREDE    # raio INTERIOR na base (diametro exterior 13.2 mm)
R_PONTA = 4.1e-3 / 2             # raio INTERIOR na ponta (diametro interior 4.1 mm)
u_r = np.hypot(0.025e-3, 0.1e-3) # craveira (0.05 mm no diametro) + incerteza da parede
S = np.pi * R_BASE * R_PONTA     # area do cilindro equivalente ao tronco de cone
u_S_rel = np.hypot(u_r / R_BASE, u_r / R_PONTA)
L_GEOM = 35.0e-3                 # comprimento do cone, no eixo (m)
# ar dentro do cone: nao entra na leitura do embolo, soma-se a cavidade
V_CONE = np.pi * L_GEOM / 3 * (R_BASE**2 + R_BASE * R_PONTA + R_PONTA**2) * 1e6   # mL
u_V = 0.5                        # incerteza de leitura do embolo (mL)

F_MIN, F_MAX = 150.0, 1500.0     # banda de procura do pico de IL
NPERSEG = 2**15                  # resolucao ~1.5 Hz
LARGURAS_HZ = [20, 30, 40, 50, 60, 70]   # larguras da media movel em frequencia
REF = "A_ref_7outopenlab.wav"


def ler_wav(caminho):
    with wave.open(caminho) as w:
        fs = w.getframerate()
        x = np.frombuffer(w.readframes(w.getnframes()), dtype=np.int16)
    return fs, x.astype(float) / 32768.0


def espectro(caminho):
    fs, x = ler_wav(caminho)
    f, P = welch(x, fs, nperseg=NPERSEG)
    return f, P


def suavizar_db(f, y, largura=30.0):
    n = max(1, int(round(largura / (f[1] - f[0]))))
    return uniform_filter1d(y, n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sem", type=float, nargs="*", default=[],
                    help="volumes (mL) a excluir do ajuste")
    ap.add_argument("--ref", action="store_true",
                    help="usar A_ref como referencia em vez da mediana")
    args = ap.parse_args()

    # agrupar gravacoes fo_<V>mL_<i>.wav
    medidas = {}
    for p in sorted(glob.glob(os.path.join(AQUI, "fo_*mL_*.wav"))):
        m = re.search(r"fo_([\d.,]+)mL_(\d+)\.wav$", os.path.basename(p))
        if m:
            V = float(m.group(1).replace(",", "."))
            medidas.setdefault(V, []).append(p)

    espectros = {p: espectro(p) for ps in medidas.values() for p in ps}
    f = next(iter(espectros.values()))[0]
    if args.ref:
        Pref = espectro(os.path.join(AQUI, REF))[1]
        nome_ref = REF
    else:
        Pref = np.median([P for _, P in espectros.values()], axis=0)
        nome_ref = "mediana das gravacoes"
    banda = (f >= F_MIN) & (f <= F_MAX)

    fig1, ax1 = plt.subplots(figsize=(8, 4.5))
    linhas = []
    for V in sorted(medidas):
        f0s, sist, picos = [], [], []
        IL_media = 0
        for p in medidas[V]:
            IL_bruto = 10 * np.log10(Pref / espectros[p][1])
            est = []
            for W in LARGURAS_HZ:
                IL = suavizar_db(f, IL_bruto, W)
                est.append(f[np.argmax(np.where(banda, IL, -np.inf))])
            f0s.append(np.mean(est))
            sist.append(np.std(est))
            IL = suavizar_db(f, IL_bruto, 30.0)
            picos.append(IL[banda].max())
            IL_media = IL_media + IL / len(medidas[V])
        f0s = np.array(f0s)
        sd = f0s.std(ddof=1) if len(f0s) > 1 else np.nan
        u_sist = np.sqrt(np.mean(np.square(sist)))
        linhas.append((V, f0s.mean(), sd, len(f0s), f0s, np.mean(picos), u_sist))
        ax1.plot(f[banda], IL_media[banda], label=f"{V:g} mL")
    ax1.set_xlabel("f (Hz)")
    ax1.set_ylabel("IL = 10 log(P_ref/P)  (dB)")
    ax1.set_title(f"Perda por insercao, 1 ressonador (media de 3)\nref: {nome_ref}")
    ax1.grid(alpha=0.3)
    ax1.legend()
    estilo_tl.virgulas(fig1)
    fig1.tight_layout()
    fig1.savefig(os.path.join(AQUI, "IL_por_volume" + ("_Aref" if args.ref else "") + ".png"), dpi=150)

    print(f"c = {c:.1f} m/s (T = {T_C} +- {u_T} C),  R_base = {R_BASE*1e3:.2f} mm,  "
          f"R_ponta = {R_PONTA*1e3:.2f} mm,  S_eq = pi R_base R_ponta = {S*1e6:.2f} mm^2\n")
    print(f"referencia: {nome_ref}")
    print(" V(mL)   f0 por repeticao (Hz)        media   sd_rep  u_sist  IL_pico(dB)")
    for V, mu, sd, n, f0s, pk, us in linhas:
        print(f" {V:5.1f}   {', '.join(f'{x:6.1f}' for x in f0s):27s}  "
              f"{mu:6.1f}  {sd:6.1f}  {us:6.1f}    {pk:5.1f}")

    # ---------- ajuste f0 = m / sqrt(V) ----------
    usar = [l for l in linhas if l[0] not in args.sem]
    V = np.array([l[0] for l in usar])
    y = np.array([l[1] for l in usar])
    n = np.array([l[3] for l in usar])
    u_stat = np.array([l[2] for l in usar]) / np.sqrt(n)
    u_stat = np.where(np.isfinite(u_stat), u_stat, 0)
    u_stat = np.sqrt(u_stat**2 + np.array([l[6] for l in usar]) ** 2)
    # incerteza do volume propagada: d f0 / f0 = 1/2 dV/V
    u_y = np.sqrt(u_stat**2 + (0.5 * y * u_V / V) ** 2 + 1.0**2)  # +1 Hz resolucao
    x = 1 / np.sqrt((V + V_CONE) * 1e-6)            # m^(-3/2), cavidade = seringa + cone
    w = 1 / u_y**2
    m = np.sum(w * x * y) / np.sum(w * x * x)
    u_m_int = 1 / np.sqrt(np.sum(w * x * x))
    res = y - m * x
    chi2 = np.sum(w * res**2)
    dof = len(x) - 1
    u_m = u_m_int * max(1.0, np.sqrt(chi2 / dof)) if dof > 0 else u_m_int

    L_ef = c**2 * S / (4 * np.pi**2 * m**2)
    u_L_ef = L_ef * np.sqrt((2 * u_m / m) ** 2 + u_S_rel ** 2 + (2 * u_c_rel) ** 2)

    print(f"\nAjuste f0 = m/sqrt(V) com V = {', '.join(f'{v:g}' for v in V)} mL"
          + (f"  (excluidos: {', '.join(f'{v:g}' for v in args.sem)})" if args.sem else ""))
    print(f"  m = {m:.4f} +- {u_m:.4f}  Hz.m^1.5    chi2/dof = {chi2:.1f}/{dof} = {chi2 / dof:.2f}"
          + ("  (> 1: incerteza de m multiplicada por sqrt(chi2/dof))" if chi2 / dof > 1 else ""))
    print("  residuos (Hz): " + ", ".join(f"{v:g} mL: {d:+.1f}" for v, d in zip(V, res)))
    print(f"  L_ef = c^2 S/(4 pi^2 m^2) = {L_ef*1e3:.1f} +- {u_L_ef*1e3:.1f} mm")
    print(f"  S/L_ef = (2 pi m / c)^2 = {(2*np.pi*m/c)**2*1e3:.3f} mm  (so isto e medido; "
          "L_ef depende dos raios)")
    # dL_eq = k_base R_PONTA + k_ponta R_BASE ; k = 0.85 (flangeada), 0.61 (livre)
    teo_flange = 0.85 * (R_PONTA + R_BASE)
    teo_mista = min(0.85 * R_PONTA + 0.61 * R_BASE, 0.61 * R_PONTA + 0.85 * R_BASE)
    if L_GEOM is not None:
        dL = L_ef - L_GEOM
        print(f"  dL = L_ef - L = {dL*1e3:.1f} mm   (teorico, 1 ou 2 pontas flangeadas = "
              f"{teo_mista*1e3:.1f}-{teo_flange*1e3:.1f} mm)")
    else:
        print("  L geometrico nao medido -> dL = L_ef - L fica em aberto.")
        print(f"  Teoria: dL_eq = k_b R_ponta + k_p R_base, k = 0.61-0.85 "
              f"= {teo_mista*1e3:.1f}-{teo_flange*1e3:.1f} mm")
        print(f"  -> L geometrico implicito: {(L_ef-teo_flange)*1e3:.1f}-"
              f"{(L_ef-teo_mista)*1e3:.1f} mm (comparar com a craveira)")

    # ---------- grafico f0 vs 1/sqrt(V) ----------
    # em cima: pontos e reta; em baixo: resíduos. As barras são as mesmas incertezas do ajuste.
    fig2, (ax2, ax_r) = plt.subplots(2, 1, figsize=(6.5, 5.5), sharex=True,
                                     gridspec_kw=dict(height_ratios=[3, 1]))
    for V_, mu, sd, n_, f0s, _, us in linhas:
        xi = 1 / np.sqrt((V_ + V_CONE) * 1e-6)
        exc = V_ in args.sem
        u_i = np.sqrt((sd / np.sqrt(n_)) ** 2 + us ** 2 + (0.5 * mu * u_V / V_) ** 2 + 1.0 ** 2)
        primeiro = not exc and V_ == min(l[0] for l in linhas if l[0] not in args.sem)
        ax2.errorbar(xi, mu, yerr=u_i, fmt="s" if exc else "o", color="tab:gray" if exc else "tab:blue",
                     capsize=3, label="excluído do ajuste" if exc else
                     ("medido (média de 3 ± σ)" if primeiro else None))
        if not exc:   # o excluído fica fora da escala dos resíduos (−97 Hz a 10 mL)
            ax_r.errorbar(xi, mu - m * xi, yerr=u_i, fmt="o", color="tab:blue", capsize=3)
        ax2.annotate(f"{V_:g} mL", (xi, mu), textcoords="offset points",
                     xytext=(6, -10), fontsize=8)
    xx = np.linspace(0, 1.1 * max(1 / np.sqrt((np.array([l[0] for l in linhas]) + V_CONE) * 1e-6)), 50)
    ax2.plot(xx, m * xx, "r-", label=r"ajuste $f_0 = m/\sqrt{V + V_\mathrm{cone}}$" + "\n"
             + f"m = {m:.3f} ± {u_m:.3f} Hz·m$^{{3/2}}$,  χ²/ν = {chi2:.1f}/{dof}\n"
             + f"$L_\\mathrm{{ef}}$ = {L_ef*1e3:.1f} ± {u_L_ef*1e3:.1f} mm")
    ax_r.axhline(0, color="r", linewidth=1)
    ax_r.set_xlabel(r"$1/\sqrt{V + V_\mathrm{cone}}$  (m$^{-3/2}$),  " + f"$V_\\mathrm{{cone}}$ = {V_CONE:.1f} mL")
    ax_r.set_ylabel("resíduo (Hz)")
    ax2.set_ylabel(r"$f_0$ (Hz)")
    ax2.set_xlim(0)
    ax2.set_ylim(0)
    for a_ in (ax2, ax_r):
        a_.grid(alpha=0.3)
    ax2.legend(loc="upper left", fontsize=8)
    estilo_tl.virgulas(fig2)
    fig2.tight_layout()
    nome = "f0_vs_invsqrtV" + ("_Aref" if args.ref else "") + ("_sem" + "_".join(f"{v:g}" for v in args.sem) if args.sem else "")
    fig2.savefig(os.path.join(AQUI, nome + ".png"), dpi=150)
    print(f"\nFiguras: IL_por_volume{'_Aref' if args.ref else ''}.png, {nome}.png")


if __name__ == "__main__":
    main()
