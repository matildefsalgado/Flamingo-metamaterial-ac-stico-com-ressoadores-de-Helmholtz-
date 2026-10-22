import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import welch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import estilo_tl  # noqa: E402  (estilo comum dos gráficos de TL e modelo, igual ao da separação)

# Parte A: ressonância local vs Bragg. 5 seringas com o volume que põe f0 ≈ fB (fusão),
# comparadas com as 5 seringas separadas (14,6 mL) da pasta numero_elementos.
# TL(f) = 20*log10(A_ref / A_medida), cada medição com a referência da SUA sessão.
# As sessões têm taxas de amostragem diferentes (44,1 e 48 kHz): os espectros são
# interpolados para uma grelha comum de 10 Hz.

# ---------- Ficheiros ----------
REF = 'fusao_ref.wav'                                  # êmbolos no fundo, lã na ponta (esta pasta)
MEDIDAS = ['fusao1.wav', 'fusao2.wav', 'fusao3.wav']   # 5 seringas a V, lã na ponta
# Separadas: 5 seringas a V_SEP (pasta numero_elementos, aula 3)
REF_SEP = 'A_ref_Aula3.wav'
MEDIDAS_SEP = ['dependencia_N5-1.wav', 'dependencia_N5-2.wav', 'dependencia_N5-3.wav']
# Teto: o gravado a 7/10 na pasta teto_sistema (outra sessão -> só indicativo)
TETO_REF, TETO_SILENCIO, TETO_FUGA = 'A_ref.wav', 'A_silencio.wav', 'A_fuga.wav'

# ---------- Parâmetros (SI) ----------
T = 22            # temperatura da sala (°C), a do registo (igual em todos os scripts)
U_T = 0.5         # incerteza da temperatura (°C)
a = 0.25          # espaçamento entre seringas (m)
U_A = 0.001       # incerteza do espaçamento (m), fita
U_V = 0.5e-6      # incerteza de leitura do volume (m³): metade da menor divisão (1 mL)
U_R = 0.1e-3      # incerteza dos raios do cone (m): craveira + parede
U_LEF = 2.5e-3    # incerteza de L_ef (m), da calibração (correcao_gargalo --sem 10)
RES = 10.0        # resolução em frequência (Hz): a incerteza de leitura é metade, 5 Hz
V = 4.1e-6        # volume de cada seringa na fusão (m³)
V_SEP = 14.6e-6   # volume das seringas separadas (m³)
NIVEIS = (6.0, 10.0)            # níveis para medir a largura das bandas (dB)
F = np.arange(200, 1201, 10.0)  # grelha comum (Hz)

sys.stdout.reconfigure(encoding='utf-8')
PASTA = os.path.dirname(os.path.abspath(__file__))
PASTA_SEP = os.path.join(PASTA, '..', 'numero_elementos')
PASTA_TETO = os.path.join(PASTA, '..', 'teto_sistema')


def espectro(nome, pasta=PASTA):
    # Sem normalizar (o TL é uma razão). Welch: pedaços de 0,1 s -> resolução 10 Hz.
    sr, data = wavfile.read(os.path.join(pasta, nome))
    data = data.astype(float)
    if data.ndim > 1:
        data = np.mean(data, axis=1)
    f, P = welch(data, sr, nperseg=sr // 10)
    return np.interp(F, f, np.sqrt(P))


def tl(ref, medidas, pasta=PASTA):
    """TL médio e desvio padrão das repetições."""
    A_ref = espectro(ref, pasta)
    reps = np.array([20 * np.log10(A_ref / espectro(m, pasta)) for m in medidas])
    return reps.mean(axis=0), reps.std(axis=0, ddof=1), reps


def pico(TL, f_min, f_max):
    m = np.flatnonzero((F >= f_min) & (F <= f_max))
    return m[np.argmax(TL[m])]


def largura(TL, i, nivel):
    """(f_esq, f_dir) da zona contínua à volta de i com TL >= nivel, interpolada; None se não chega."""
    if TL[i] < nivel:
        return None
    e, d = i, i
    while e > 0 and TL[e - 1] >= nivel:
        e -= 1
    while d < len(TL) - 1 and TL[d + 1] >= nivel:
        d += 1
    f_e = F[e] if e == 0 else np.interp(nivel, [TL[e - 1], TL[e]], [F[e - 1], F[e]])
    f_d = F[d] if d == len(TL) - 1 else np.interp(nivel, [TL[d + 1], TL[d]], [F[d + 1], F[d]])
    return f_e, f_d


# ---------- Previsões ----------
c = 331.3 * np.sqrt(1 + T / 273.15)
# Gargalo = bico cónico da seringa. Valores acertados aos dados de correcao_gargalo
# (medidas a olho ajustadas o mínimo para o modelo bater certo).
L = 0.035              # comprimento do cone, no eixo (m)
D_BASE = 0.0132        # diâmetro exterior na base (m)
PAREDE = 0.0012        # grossura da parede na base (m)
D_PONTA = 0.0041       # diâmetro interior na ponta (m)
R_B, R_P = D_BASE / 2 - PAREDE, D_PONTA / 2
S = np.pi * R_B * R_P  # área do cilindro equivalente ao cone: ∫dx/S(x) = L/(π R_B R_P)
dL = 0.73 * (R_B + R_P)   # correção de extremidade, cada ponta com o seu raio (m)
V_CONE = np.pi * L / 3 * (R_B**2 + R_B * R_P + R_P**2)   # ar no cone (m³): soma-se ao êmbolo


def f0_de(v):
    return c / (2 * np.pi) * np.sqrt(S / ((v + V_CONE) * (L + dL)))


def u_f0_de(v):
    """Incerteza do f0 previsto: f0 ∝ c·√(S/(V_total·L_ef)), tudo em quadratura."""
    u_S_rel = np.hypot(U_R / R_B, U_R / R_P)
    return f0_de(v) * np.sqrt((0.5 * U_V / (v + V_CONE))**2 + (0.5 * u_S_rel)**2
                              + (0.5 * U_LEF / (L + dL))**2 + u_rel_c**2)


u_rel_c = 0.5 * U_T / (T + 273.15)
f0_prev, f0_sep_prev = f0_de(V), f0_de(V_SEP)
u_f0_prev = u_f0_de(V)
fB_prev = c / (2 * a)
u_fB_prev = fB_prev * np.hypot(U_A / a, u_rel_c)
print(f"c = {c:.1f} m/s | fB previsto = {fB_prev:.0f} ± {u_fB_prev:.0f} Hz | "
      f"f0 previsto: {f0_prev:.0f} ± {u_f0_prev:.0f} Hz ({V * 1e6:.1f} mL), "
      f"{f0_sep_prev:.0f} Hz ({V_SEP * 1e6:.1f} mL)")

# ---------- TL ----------
TL_fus, sd_fus, reps_fus = tl(REF, MEDIDAS)
# As gravações N5 da aula 3 foram apagadas: sem elas, a parte das separadas é saltada
ha_sep = all(os.path.exists(os.path.join(PASTA_SEP, m)) for m in [REF_SEP] + MEDIDAS_SEP)
if ha_sep:
    TL_sep, sd_sep, _ = tl(REF_SEP, MEDIDAS_SEP, PASTA_SEP)
else:
    print("(faltam as gravações N5 da aula 3 em numero_elementos: a parte das separadas é saltada)")
teto = 20 * np.log10(espectro(TETO_REF, PASTA_TETO)
                     / np.maximum(espectro(TETO_SILENCIO, PASTA_TETO), espectro(TETO_FUGA, PASTA_TETO)))

# ---------- Picos ----------
# Ressonância local da fusão: perto do f0 previsto. Bragg: no TL das separadas, acima de fB
# (a banda de Bragg sai deslocada para cima de c/2a; janela larga para a apanhar).
i_f0 = pico(TL_fus, f0_prev - 150, f0_prev + 50)
f0_med = F[i_f0]
# Incerteza do f0 medido: dispersão do pico entre repetições (desvio padrão da média)
# e metade da resolução, em quadratura
f0_reps = np.array([F[pico(r, f0_prev - 150, f0_prev + 50)] for r in reps_fus])
u_f0 = np.hypot(f0_reps.std(ddof=1) / np.sqrt(len(f0_reps)), RES / 2)
u_TL = sd_fus[i_f0] / np.sqrt(len(reps_fus))
print(f"\nFusão ({V * 1e6:.1f} mL): ressonância local a {f0_med:.0f} ± {u_f0:.0f} Hz "
      f"(picos das repetições: {', '.join(f'{x:.0f}' for x in f0_reps)} Hz), "
      f"TL = {TL_fus[i_f0]:.1f} ± {u_TL:.1f} dB")
for nivel in NIVEIS:
    lim = largura(TL_fus, i_f0, nivel)
    if lim:
        print(f"   zona contínua com TL >= {nivel:.0f} dB: {lim[0]:.0f}–{lim[1]:.0f} Hz "
              f"(largura {lim[1] - lim[0]:.0f} Hz)")
if ha_sep:
    i_B = pico(TL_sep, fB_prev, fB_prev + 250)
    u_B = 20.0                      # a banda de Bragg é larga
    fB_med = F[i_B]
    print(f"Separadas ({V_SEP * 1e6:.1f} mL): banda de Bragg com máximo a {fB_med:.0f} Hz, "
          f"TL = {TL_sep[i_B]:.1f} dB  (previsto c/2a = {fB_prev:.0f} Hz, "
          f"desvio {100 * (fB_med / fB_prev - 1):+.0f}%)")
    delta = (fB_med - f0_med) / fB_med
    u_delta = np.hypot(u_f0 / fB_med, f0_med * u_B / fB_med**2)
    print(f"\nDessintonia contra fB medido: δ = {100 * delta:.0f} ± {100 * u_delta:.0f} %")

# Dessintonia contra o fB previsto: δ = 1 - f0/fB, propagada em quadratura
delta_p = 1 - f0_med / fB_prev
u_delta_p = (f0_med / fB_prev) * np.hypot(u_f0 / f0_med, u_fB_prev / fB_prev)
print(f"Dessintonia contra fB previsto: δ = (fB - f0)/fB = {100 * delta_p:.1f} ± {100 * u_delta_p:.1f} %"
      f"   (z contra 0 = {delta_p / u_delta_p:.1f})")
z = abs(f0_med - f0_prev) / np.hypot(u_f0, u_f0_prev)
print(f"f0 medido vs previsto ({V * 1e6:.1f} mL): {f0_med:.0f} ± {u_f0:.0f} vs "
      f"{f0_prev:.0f} ± {u_f0_prev:.0f} Hz ({100 * (f0_med / f0_prev - 1):+.0f} %, z = {z:.1f})")

# ---------- Pico de Bragg e modelo ----------
# Pico de Bragg: o maior acima da ressonância local (pelo menos 40 Hz acima) e até 1,2 fB
# (o mesmo critério de separacao.py)
jan_brg = (f0_med + 40, 1.2 * fB_prev)
i_brg = pico(TL_fus, *jan_brg)
brg_reps = np.array([F[pico(r, *jan_brg)] for r in reps_fus])
u_brg = np.hypot(brg_reps.std(ddof=1) / np.sqrt(len(brg_reps)), RES / 2)
print(f"Pico de Bragg: {F[i_brg]:.0f} ± {u_brg:.0f} Hz, TL = {TL_fus[i_brg]:.1f} dB")
# Modelo de matrizes de transferência com o f0 medido e Q = 13 (do N1 da aula 3, como na separação)
Q = 13.2
D_TUBO, RHO = 0.025, 1.2
v_eq = S * c**2 / ((2 * np.pi * f0_med)**2 * (L + dL)) - V_CONE
F_m = np.arange(F[0], F[-1] + 0.1, 2.0)
M5 = estilo_tl.tl_modelo(5, v_eq, Q, F_m, c, S, L + dL, V_CONE, a, np.pi * D_TUBO**2 / 4, RHO)

# ---------- Gráfico (estilo comum: estilo_tl.py, igual ao da separação) ----------
fig, ax = plt.subplots(figsize=(11, 6))
estilo_tl.grafico_tl(
    ax, F, TL_fus, sd_fus,
    titulo=f'Fusão: 5 seringas a {V * 1e6:.1f} mL (13/10), média de {len(reps_fus)}, faixa = desvio padrão',
    rotulo=f'5 seringas a {V * 1e6:.1f} mL', f0_prev=f0_prev, u_f0_prev=u_f0_prev,
    fB=fB_prev, u_fB=u_fB_prev, f_loc=f0_med, u_loc=u_f0, TL_loc=TL_fus[i_f0],
    f_brg=F[i_brg], u_brg=u_brg, TL_brg=TL_fus[i_brg],
    F_m=F_m, modelo=M5, rotulo_modelo=f'Modelo, 5 seringas (Q = {Q:.0f})')
estilo_tl.virgulas(fig)
fig.tight_layout()
# Sem as separadas não se escreve por cima da figura completa que já existe
NOME_FIG = 'tl_parteA.png' if ha_sep else 'tl_parteA_so_fusao.png'
fig.savefig(os.path.join(PASTA, NOME_FIG), dpi=150)
print(f"\nFigura: {NOME_FIG}")
plt.show()
