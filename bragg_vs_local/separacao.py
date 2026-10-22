import os
import re
import sys
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import welch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import estilo_tl  # noqa: E402  (estilo comum dos gráficos de TL e modelo)

# Separação da ressonância local e do efeito de Bragg.
# 5 seringas a V = 5 mL, espaçadas de a = 0,25 m: f0 fica abaixo de fB = c/2a, por isso os dois
# picos deviam aparecer separados no TL(f). Na fusão (4,1 mL) estavam a ~110 Hz um do outro;
# a 5 mL a ressonância local desce e o pico de Bragg deve ficar no mesmo sítio.
# TL(f) = 20*log10(A_ref / A_medida), média das repetições.
#
# Gravações (nesta pasta, todas da MESMA sessão):
#   ref_<n>.wav         êmbolos todos no fundo, lã na ponta              (referência)
#   sep5mL_<n>.wav      5 seringas a 5 mL, lã na ponta                   (medida)
#   ctrl5mL_<n>.wav     só a seringa 1 a 5 mL, as outras no fundo        (controlo, opcional)
# Uma seringa sozinha não tem Bragg: um pico que aparece com 5 seringas e não com 1 é Bragg.
# Para comparar, usa também a fusão (tl_parteA.py, nesta pasta, 4,1 mL), se lá estiver.
# Modelo: matrizes de transferência (tubo sem perdas, seringas como ramos laterais de Helmholtz,
# terminação sem reflexão), com o f0 medido e um Q tirado do controlo (ou Q = 13 sem controlo).

# ---------- Ficheiros ----------
PADRAO_REF = r'ref'                    # ref_1.wav, ref16_10.wav, ...
PADRAO_MED = r'(?:sep5mL|sep_bandas)'  # sep5mL_1.wav, sep_bandas1.wav, ...
PADRAO_CTRL = r'ctrl5mL'               # ctrl5mL_1.wav, ...
REF_FUS = 'fusao_ref.wav'
MEDIDAS_FUS = ['fusao1.wav', 'fusao2.wav', 'fusao3.wav']

# ---------- Parâmetros (SI) ----------
T = 22            # temperatura da sala (°C), a do registo (igual em todos os scripts)
U_T = 0.5         # incerteza da temperatura (°C)
a = 0.25          # espaçamento entre seringas (m)
U_A = 0.001       # incerteza do espaçamento (m), fita
U_V = 0.5e-6      # incerteza de leitura do volume (m³): metade da menor divisão (1 mL)
U_R = 0.1e-3      # incerteza dos raios do cone (m): craveira + parede
U_LEF = 2.5e-3    # incerteza de L_ef (m), da calibração (correcao_gargalo --sem 10)
RES = 10.0        # resolução em frequência (Hz): a incerteza de leitura é metade, 5 Hz
V = 5e-6          # volume de cada seringa (m³)
V_FUS = 4.1e-6    # volume da fusão, para comparar (m³)
Q_SEM_CTRL = 13.2 # Q usado se não houver controlo (tirado do N1 da aula 3)
D_TUBO = 0.025    # diâmetro interior do tubo de PVC (m)
RHO = 1.2         # densidade do ar (kg/m³)
F = np.arange(100, 2001, 10.0)   # grelha comum (Hz): varrimento de 100 a 2000 Hz

sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
PASTA = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(description='Separação ressonância local / Bragg')
ap.add_argument('--pasta', default=PASTA, help='pasta com as gravações (por omissão, esta)')
args = ap.parse_args()
PASTA_DADOS = args.pasta
PASTA_FUS = PASTA                                   # a fusão (tl_parteA.py) está nesta mesma pasta


def ficheiros(prefixo):
    """Gravações <prefixo><números>.wav, com _ ou - (ex.: ref_1, ref16_10, sep_bandas2; aceita ' (1)')."""
    p = re.compile(rf'^{prefixo}[_-]?[\d_-]*(?: \(\d+\))?\.wav$', re.IGNORECASE)
    return sorted(f for f in os.listdir(PASTA_DADOS) if p.match(f))


def espectro(nome, pasta):
    # Sem normalizar (o TL é uma razão). Welch: pedaços de 0,1 s -> resolução 10 Hz.
    sr, data = wavfile.read(os.path.join(pasta, nome))
    data = data.astype(float)
    if data.ndim > 1:
        data = np.mean(data, axis=1)
    f, P = welch(data, sr, nperseg=sr // 10)
    return np.interp(F, f, np.sqrt(P))


def tl(refs, medidas, pasta):
    """TL médio, desvio padrão e repetições. Com várias referências usa a média delas."""
    A_ref = np.mean([espectro(r, pasta) for r in refs], axis=0)
    reps = np.array([20 * np.log10(A_ref / espectro(m, pasta)) for m in medidas])
    sd = reps.std(axis=0, ddof=1) if len(reps) > 1 else np.zeros(len(F))
    return reps.mean(axis=0), sd, reps


def pico(TL, f_min, f_max):
    m = np.flatnonzero((F >= f_min) & (F <= f_max))
    return m[np.argmax(TL[m])]


def f_pico(reps, TLmedio, f_min, f_max):
    """Frequência do pico, com incerteza: dispersão entre repetições (s/√n) e metade da resolução."""
    i = pico(TLmedio, f_min, f_max)
    fr = np.array([F[pico(r, f_min, f_max)] for r in reps])
    u_rep = fr.std(ddof=1) / np.sqrt(len(fr)) if len(fr) > 1 else 0.0
    return i, np.hypot(u_rep, RES / 2), fr


# ---------- Previsões ----------
c = 331.3 * np.sqrt(1 + T / 273.15)
# Gargalo = bico cónico da seringa (valores acertados aos dados de correcao_gargalo)
L = 0.035              # comprimento do cone, no eixo (m)
D_BASE = 0.0132        # diâmetro exterior na base (m)
PAREDE = 0.0012        # grossura da parede na base (m)
D_PONTA = 0.0041       # diâmetro interior na ponta (m)
R_B, R_P = D_BASE / 2 - PAREDE, D_PONTA / 2
S = np.pi * R_B * R_P  # área do cilindro equivalente ao cone
dL = 0.73 * (R_B + R_P)   # correção de extremidade (m)
L_EF = L + dL
V_CONE = np.pi * L / 3 * (R_B**2 + R_B * R_P + R_P**2)   # ar no cone (m³)
S_TUBO = np.pi * D_TUBO**2 / 4
u_rel_c = 0.5 * U_T / (T + 273.15)


def f0_de(v):
    return c / (2 * np.pi) * np.sqrt(S / ((v + V_CONE) * L_EF))


def u_f0_de(v):
    """Incerteza do f0 previsto, tudo em quadratura."""
    u_S_rel = np.hypot(U_R / R_B, U_R / R_P)
    return f0_de(v) * np.sqrt((0.5 * U_V / (v + V_CONE))**2 + (0.5 * u_S_rel)**2
                              + (0.5 * U_LEF / L_EF)**2 + u_rel_c**2)


def tl_modelo(N, v, Q, f):
    """TL de N ressoadores espaçados de a (modelo de matrizes de transferência, em estilo_tl.py)."""
    return estilo_tl.tl_modelo(N, v, Q, f, c, S, L_EF, V_CONE, a, S_TUBO, RHO)


def v_equivalente(f0):
    """Volume que, no modelo, dá este f0."""
    return S * c**2 / ((2 * np.pi * f0)**2 * L_EF) - V_CONE


f0_prev, u_f0_prev = f0_de(V), u_f0_de(V)
fB_prev = c / (2 * a)
u_fB_prev = fB_prev * np.hypot(U_A / a, u_rel_c)
print(f"c = {c:.1f} m/s")
print(f"Previsto ({V * 1e6:.1f} mL): f0 = {f0_prev:.0f} ± {u_f0_prev:.0f} Hz   "
      f"fB = c/2a = {fB_prev:.0f} ± {u_fB_prev:.0f} Hz   (f0/fB = {f0_prev / fB_prev:.2f})")
print(f"   (na fusão o f0 medido ficou 8 % abaixo do previsto: a 5 mL espera-se ~{0.924 * f0_prev:.0f} Hz)")

# ---------- TL medido ----------
refs, meds, ctrls = ficheiros(PADRAO_REF), ficheiros(PADRAO_MED), ficheiros(PADRAO_CTRL)
if not refs or not meds:
    raise SystemExit(f"Faltam gravações em {PASTA_DADOS}: preciso de {PADRAO_REF}_1.wav, ... "
                     f"e {PADRAO_MED}_1.wav, ... (ver README.md).")
print(f"\nReferência: {', '.join(refs)}\nMedidas:    {', '.join(meds)}"
      + (f"\nControlo:   {', '.join(ctrls)}" if ctrls else "\nControlo:   (não há: sem ele não se prova que o 2.º pico é Bragg)"))
TL5, sd5, reps5 = tl(refs, meds, PASTA_DADOS)
n5 = len(reps5)

# Ressonância local: o maior pico entre 0,75 f0 previsto e um pouco abaixo de fB
jan_loc = (0.75 * f0_prev, fB_prev - 30)
i_loc, u_loc, loc_reps = f_pico(reps5, TL5, *jan_loc)
f_loc = F[i_loc]
z_loc = abs(f_loc - f0_prev) / np.hypot(u_loc, u_f0_prev)
print(f"\nRessonância local: {f_loc:.0f} ± {u_loc:.0f} Hz (repetições: {', '.join(f'{x:.0f}' for x in loc_reps)} Hz), "
      f"TL = {TL5[i_loc]:.1f} ± {sd5[i_loc] / np.sqrt(n5):.1f} dB   vs previsto {f0_prev:.0f} ± {u_f0_prev:.0f} Hz "
      f"({100 * (f_loc / f0_prev - 1):+.0f} %, z = {z_loc:.1f})")

# Pico de Bragg: o maior acima da ressonância local (pelo menos 40 Hz acima) e até 1,2 fB
jan_brg = (f_loc + 40, 1.2 * fB_prev)
i_brg, u_brg, brg_reps = f_pico(reps5, TL5, *jan_brg)
f_brg = F[i_brg]
i_min = i_loc + np.argmin(TL5[i_loc:i_brg + 1])
z_fB = abs(f_brg - fB_prev) / np.hypot(u_brg, u_fB_prev)
print(f"Mínimo entre os dois: {F[i_min]:.0f} Hz, TL = {TL5[i_min]:.1f} dB")
print(f"Pico de Bragg: {f_brg:.0f} ± {u_brg:.0f} Hz (repetições: {', '.join(f'{x:.0f}' for x in brg_reps)} Hz), "
      f"TL = {TL5[i_brg]:.1f} ± {sd5[i_brg] / np.sqrt(n5):.1f} dB   vs fB = {fB_prev:.0f} ± {u_fB_prev:.0f} Hz "
      f"({100 * (f_brg / fB_prev - 1):+.0f} %, z = {z_fB:.1f})")
print(f"Separação entre os picos: {f_brg - f_loc:.0f} ± {np.hypot(u_loc, u_brg):.0f} Hz; "
      f"o pico de Bragg está {TL5[i_brg] - TL5[i_min]:.1f} dB acima do mínimo")

# Controlo: 1 seringa ao mesmo volume, mesma sessão. No pico de Bragg deve dar ~0 dB.
Q = Q_SEM_CTRL
if ctrls:
    TL1, sd1, reps1 = tl(refs, ctrls, PASTA_DADOS)
    n1 = len(reps1)
    dif = TL5[i_brg] - TL1[i_brg]
    u_dif = np.hypot(sd5[i_brg] / np.sqrt(n5), sd1[i_brg] / np.sqrt(n1))
    u_dif = max(u_dif, 0.1)
    print(f"\nControlo (1 seringa) à frequência do pico de Bragg: TL = {TL1[i_brg]:.1f} dB")
    print(f"   5 seringas − 1 seringa = {dif:.1f} ± {u_dif:.1f} dB  (z contra 0 = {dif / u_dif:.1f}"
          f"{': o pico é da rede, não do tubo' if dif / u_dif > 3 else ': não se distingue do controlo'})")
    i_1 = pico(TL1, *jan_loc)
    print(f"   ressonância local do controlo: {F[i_1]:.0f} Hz, TL = {TL1[i_1]:.1f} dB")
    # Q escolhido para o pico do modelo com 1 seringa ter a altura do pico do controlo (bisseção)
    F_fino = np.linspace(F[i_1] - 60, F[i_1] + 60, 400)
    q_lo, q_hi = 0.5, 200.0
    for _ in range(60):
        Q = np.sqrt(q_lo * q_hi)
        if tl_modelo(1, v_equivalente(F[i_1]), Q, F_fino).max() > TL1[i_1]:
            q_hi = Q
        else:
            q_lo = Q

# ---------- Modelo ----------
F_m = np.arange(F[0], F[-1] + 0.1, 2.0)
M5 = tl_modelo(5, v_equivalente(f_loc), Q, F_m)
M1 = tl_modelo(1, v_equivalente(f_loc), Q, F_m)
maximos = np.flatnonzero((M5[1:-1] > M5[:-2]) & (M5[1:-1] > M5[2:])) + 1
maximos = maximos[(F_m[maximos] >= f_loc + 40) & (F_m[maximos] <= 1.2 * fB_prev)]
print(f"\nModelo (f0 = {f_loc:.0f} Hz medido, Q = {Q:.1f}"
      f"{', do controlo' if ctrls else ', do N1 da aula 3'}):")
if len(maximos):
    i_bm = maximos[np.argmax(M5[maximos])]
    z_m = abs(f_brg - F_m[i_bm]) / np.hypot(u_brg, u_fB_prev)
    print(f"   pico de Bragg a {F_m[i_bm]:.0f} Hz com {M5[i_bm]:.1f} dB  (medido {f_brg:.0f} ± {u_brg:.0f} Hz, z = {z_m:.1f})")
else:
    print("   o modelo não tem pico de Bragg separado com este f0 e este Q")

# ---------- Comparação com a fusão (4,1 mL) ----------
ha_fus = all(os.path.exists(os.path.join(PASTA_FUS, x)) for x in [REF_FUS] + MEDIDAS_FUS)
if ha_fus:
    TLf, sdf, repsf = tl([REF_FUS], MEDIDAS_FUS, PASTA_FUS)
    jl = (0.75 * f0_de(V_FUS), fB_prev - 30)
    il_f, ul_f, _ = f_pico(repsf, TLf, *jl)
    ib_f, ub_f, _ = f_pico(repsf, TLf, F[il_f] + 40, 1.2 * fB_prev)
    print(f"\nComparação com a fusão ({V_FUS * 1e6:.1f} mL, outra sessão):")
    print(f"   ressonância local: {F[il_f]:.0f} ± {ul_f:.0f} -> {f_loc:.0f} ± {u_loc:.0f} Hz "
          f"(desloca {f_loc - F[il_f]:+.0f} ± {np.hypot(ul_f, u_loc):.0f} Hz; previsto "
          f"{f0_de(V) - f0_de(V_FUS):+.0f} Hz)")
    print(f"   pico de Bragg:     {F[ib_f]:.0f} ± {ub_f:.0f} -> {f_brg:.0f} ± {u_brg:.0f} Hz "
          f"(desloca {f_brg - F[ib_f]:+.0f} ± {np.hypot(ub_f, u_brg):.0f} Hz; previsto 0)")

# ---------- Gráficos (estilo comum: estilo_tl.py, igual ao da fusão) ----------
comum = dict(fB=fB_prev, u_fB=u_fB_prev)
fig, ax = plt.subplots(figsize=(11, 6))
estilo_tl.grafico_tl(
    ax, F, TL5, sd5,
    titulo=f'Separação: 5 seringas a {V * 1e6:.0f} mL (16/10), média de {n5}, faixa = desvio padrão',
    rotulo=f'5 seringas a {V * 1e6:.0f} mL', f0_prev=f0_prev, u_f0_prev=u_f0_prev,
    f_loc=f_loc, u_loc=u_loc, TL_loc=TL5[i_loc], f_brg=f_brg, u_brg=u_brg, TL_brg=TL5[i_brg],
    F_m=F_m, modelo=M5, rotulo_modelo=f'Modelo, 5 seringas (Q = {Q:.0f})',
    controlo=(TL1, sd1, f'Controlo: 1 seringa a {V * 1e6:.0f} mL') if ctrls else None, **comum)
estilo_tl.virgulas(fig)
fig.tight_layout()
NOME_FIG = os.path.join(PASTA_DADOS, 'separacao.png')
fig.savefig(NOME_FIG, dpi=150)
print(f"\nFigura: {NOME_FIG}")

plt.show()
