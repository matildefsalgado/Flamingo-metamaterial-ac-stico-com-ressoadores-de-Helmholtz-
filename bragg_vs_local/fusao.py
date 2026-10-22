import os
import sys
import glob
import re
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import welch

# Fusão das bandas: 5 seringas com o volume que põe f0 = fB.
# Mede a largura da banda a 10 dB e compara com os picos separados (A_ref.wav e A_40mL.wav nesta pasta, 40 mL).

# ---------- Parâmetros (SI) ----------
T = 22            # temperatura da sala (°C), a do registo (igual em todos os scripts)
a = 0.25          # espaçamento entre seringas (m), medido no tubo
V_SEPARADO = 40e-6  # volume dos picos separados (m³)
NIVEL = 10.0      # nível para medir a largura da banda (dB)
F_MIN, F_MAX = 100, 2000

sys.stdout.reconfigure(encoding='utf-8')
PASTA = os.path.dirname(os.path.abspath(__file__))
PASTA_TETO = os.path.join(PASTA, '..', 'teto_sistema')
PASTA_SEP = PASTA
eps = 1e-12


def espectro(nome, pasta=PASTA):
    # Sem normalizar (o TL é uma razão). Welch: pedaços de 0,1 s -> resolução 10 Hz.
    sr, data = wavfile.read(os.path.join(pasta, nome))
    data = data.astype(float)
    if data.ndim > 1:
        data = np.mean(data, axis=1)
    f, P = welch(data, sr, nperseg=sr // 10)
    return sr, f, np.sqrt(P)


def largura(freq, TL, i_pico, nivel=NIVEL):
    """Largura da zona contínua à volta de i_pico onde TL >= nivel. Devolve (f_esq, f_dir) ou None."""
    if TL[i_pico] < nivel:
        return None
    i = i_pico
    while i > 0 and TL[i - 1] >= nivel:
        i -= 1
    j = i_pico
    while j < len(TL) - 1 and TL[j + 1] >= nivel:
        j += 1
    # Interpolação linear onde a curva cruza o nível
    f_esq = freq[i] if i == 0 else np.interp(nivel, [TL[i - 1], TL[i]], [freq[i - 1], freq[i]])
    f_dir = freq[j] if j == len(TL) - 1 else np.interp(nivel, [TL[j + 1], TL[j]], [freq[j + 1], freq[j]])
    return f_esq, f_dir


def pico_perto(freq, TL, f_alvo, janela=100):
    m = np.flatnonzero((freq > f_alvo - janela) & (freq < f_alvo + janela))
    return m[np.argmax(TL[m])]


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
L_ef = L + dL
fB = c / (2 * a)
V_sint = S * a**2 / (np.pi**2 * L_ef) - V_CONE   # volume a pôr na seringa


def f0_de(V_mL):
    return c / (2 * np.pi) * np.sqrt(S / ((V_mL * 1e-6 + V_CONE) * L_ef))


f0_sep = f0_de(V_SEPARADO * 1e6)
print(f"c = {c:.2f} m/s | fB = {fB:.0f} Hz | V sintonizado = {V_sint * 1e6:.1f} mL")

# ---------- Teto ----------
sr_t, freq, A_tref = espectro('A_ref.wav', PASTA_TETO)
_, _, A_tsil = espectro('A_silencio.wav', PASTA_TETO)
_, _, A_tfug = espectro('A_fuga.wav', PASTA_TETO)
idx = (freq >= F_MIN) & (freq <= F_MAX)
freq = freq[idx]
teto = 20 * np.log10(np.maximum(A_tref[idx], eps)
                     / np.maximum(np.maximum(A_tsil[idx], A_tfug[idx]), eps))


def calcular_TL(nome_ref, nome_med, pasta):
    sr1, _, A_ref = espectro(nome_ref, pasta)
    sr2, _, A_med = espectro(nome_med, pasta)
    assert sr1 == sr2 == sr_t, "Taxa de amostragem diferente da do teto (usar a mesma no Audacity)"
    return 20 * np.log10(np.maximum(A_ref[idx], eps) / np.maximum(A_med[idx], eps))


def relatorio(nome, TL, i):
    lim = largura(freq, TL, i)
    perto = TL[i] >= teto[i] - 3
    txt = f"{nome}: pico a {freq[i]:.0f} Hz, TL = {TL[i]:.1f} dB, teto = {teto[i]:.1f} dB"
    if perto:
        txt += "  (pico cortado pelo teto)"
    if lim is None:
        print(txt + f"\n   não chega a {NIVEL:.0f} dB -> sem largura")
        return None, lim
    larg = lim[1] - lim[0]
    print(txt + f"\n   largura a {NIVEL:.0f} dB = {larg:.0f} Hz  ({lim[0]:.0f} a {lim[1]:.0f} Hz)")
    if teto[(freq >= lim[0]) & (freq <= lim[1])].min() < NIVEL + 3:
        print(f"   atenção: dentro da banda o teto desce abaixo de {NIVEL + 3:.0f} dB")
    return larg, lim


# ---------- 1. Fusão (esta pasta) ----------
# Referência: o primeiro A_ref*.wav desta pasta.
# Medições: fusao_<V>mL_<n>.wav  (V com '-' ou ',' no lugar da vírgula: fusao_11-1mL_1.wav = 11,1 mL)
# (se não houver aqui, usa a da mesma aula que está em numero_elementos; caminho completo)
refs = (sorted(glob.glob(os.path.join(PASTA, 'A_ref_*.wav')))
        or sorted(glob.glob(os.path.join(PASTA, '..', 'numero_elementos', 'A_ref*.wav'))))
if not refs:
    raise SystemExit("Falta a referência A_ref*.wav (nem aqui nem em numero_elementos).")
nome_ref = os.path.abspath(refs[0])
padrao = re.compile(r'^fusao_(\d+(?:[-.,]\d+)?)mL_(\d+)\.wav$', re.IGNORECASE)
por_volume = defaultdict(list)       # V (mL) -> lista de TL (uma por repetição)
print(f"\nReferência: {nome_ref}")
for p in sorted(glob.glob(os.path.join(PASTA, 'fusao_*.wav'))):
    n = os.path.basename(p)
    m = padrao.match(n)
    if not m:
        print(f"(ignorado: {n}, o nome não é fusao_<V>mL_<n>.wav)")
        continue
    V_mL = float(re.sub('[-,]', '.', m.group(1)))
    TL = calcular_TL(nome_ref, n, PASTA)
    i = pico_perto(freq, TL, f0_de(V_mL))
    print(f"{n}: pico a {freq[i]:.0f} Hz, TL = {TL[i]:.1f} dB")
    por_volume[V_mL].append(TL)
if not por_volume:
    raise SystemExit("Não encontrei ficheiros fusao_<V>mL_<n>.wav nesta pasta.")

resultados = {}                      # V -> (TL médio, largura, limites)
for V_mL in sorted(por_volume):
    TL_fus = np.mean(por_volume[V_mL], axis=0)
    fp = f0_de(V_mL)
    print(f"\n=== {V_mL:g} mL (média de {len(por_volume[V_mL])}) ===")
    i_f0 = pico_perto(freq, TL_fus, fp)
    larg_fus, lim_fus = relatorio(f"Ressonância local (previsto {fp:.0f} Hz)", TL_fus, i_f0)
    relatorio(f"Bragg (previsto {fB:.0f} Hz)", TL_fus, pico_perto(freq, TL_fus, fB))
    print(f"Distância entre o pico medido e fB: {fB - freq[i_f0]:.0f} Hz")
    resultados[V_mL] = (TL_fus, larg_fus, lim_fus)
# Para a comparação usa-se o volume mais pequeno (o mais perto da fusão)
V_fus = min(resultados)
TL_fus, larg_fus, lim_fus = resultados[V_fus]

# ---------- 2. Picos separados (A_ref.wav e A_40mL.wav nesta pasta) ----------
tem_sep = all(os.path.exists(os.path.join(PASTA_SEP, n)) for n in ['A_ref.wav', 'A_40mL.wav'])
larg_sep = {}
if tem_sep:
    TL_sep = calcular_TL('A_ref.wav', 'A_40mL.wav', PASTA_SEP)
    print(f"\n=== Picos separados ({V_SEPARADO * 1e6:.0f} mL) ===")
    for nome, f_alvo in [('f0', f0_sep), ('fB', fB)]:
        larg_sep[nome] = relatorio(nome, TL_sep, pico_perto(freq, TL_sep, f_alvo))
else:
    print("\n(sem A_ref.wav e A_40mL.wav nesta pasta: não dá para comparar com os picos separados)")

# ---------- 3. Comparação ----------
if larg_fus and larg_sep:
    print("\n=== Comparação ===")
    for nome, (larg, _) in larg_sep.items():
        if larg:
            print(f"Fusão / pico {nome} separado = {larg_fus:.0f} / {larg:.0f} = {larg_fus / larg:.2f}×")

# ---------- 4. Gráfico ----------
plt.figure(figsize=(11, 6))
plt.plot(freq, teto, color='gray', linestyle='--', linewidth=1, label='Teto')
if tem_sep:
    plt.plot(freq, TL_sep, color='tab:orange', linewidth=1.2,
             label=f'Separados ({V_SEPARADO * 1e6:.0f} mL)')
cores = plt.cm.viridis(np.linspace(0, 0.85, len(resultados)))
for cor, (V_mL, (TL, larg, lim)) in zip(cores, sorted(resultados.items())):
    plt.plot(freq, TL, color=cor, linewidth=1.8,
             label=f'{V_mL:g} mL (média de {len(por_volume[V_mL])}, f0 previsto {f0_de(V_mL):.0f} Hz)')
    if lim:
        plt.axvspan(lim[0], lim[1], color=cor, alpha=0.12)
plt.axhline(NIVEL, color='black', linewidth=0.8, linestyle=':', label=f'{NIVEL:.0f} dB')
plt.axvline(fB, color='green', linestyle=':', label=f'fB = {fB:.0f} Hz')
plt.title('5 seringas: ressonância local e Bragg', fontsize=14)
plt.xlabel('Frequência (Hz)')
plt.ylabel('TL (dB)')
plt.grid(True, alpha=0.3)
plt.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(PASTA, 'fusao.png'), dpi=150)
plt.show()
