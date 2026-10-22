import os
import re
import sys
import glob
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import welch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import estilo_tl  # noqa: E402  (vírgula decimal nas figuras)

# Dependência com o número de seringas ativas N (todas com o mesmo volume sintonizado).
# Para cada N: TL_N(f) = 20*log10(A_ref / A_N), média das repetições.
# Depois: TL máximo vs N e largura da banda a 10 dB vs N, com barras de erro.
#
# Ficheiros nesta pasta:
#   A_ref.wav (ou A_ref_<qualquer coisa>.wav)  todos os êmbolos no fundo, lã na ponta
#   N1_1.wav, N1_2.wav, N1_3.wav, N2_1.wav, ...  (N = nº de seringas ativas, _n = repetição)
# O teto vem da pasta teto_sistema.

# ---------- Parâmetros (SI) ----------
T = 22            # temperatura da sala (°C), a do registo (igual em todos os scripts)
U_T = 0.5         # incerteza da temperatura (°C), igual em todos os scripts
a = 0.25          # espaçamento entre seringas (m)
U_A = 0.001       # incerteza do posicionamento das seringas (m)
U_D = 0.1e-3      # incerteza do diâmetro interior da ponta do gargalo (m), craveira
V = 14.6e-6       # volume posto em cada seringa ativa (m³)
U_V = 0.5e-6      # incerteza de leitura do volume (m³): metade da menor divisão (1 mL)
U_R = 0.1e-3      # incerteza dos raios do cone (m): craveira + parede
U_LEF = 2.5e-3    # incerteza de L_ef (m), da calibração (correcao_gargalo --sem 10)
# Todas as seringas têm o mesmo gargalo cónico (ver abaixo).
# {N: volume (m³) ou 'fusao' (volume que põe f0 = fB)} para os N com volume diferente de V
GARGALO = {}
NIVEL = 10.0      # nível para medir a largura da banda (dB)
F_MIN, F_MAX = 100, 2000

sys.stdout.reconfigure(encoding='utf-8')
PASTA = os.path.dirname(os.path.abspath(__file__))
PASTA_TETO = os.path.join(PASTA, '..', 'teto_sistema')
eps = 1e-12


def espectro(nome, pasta=PASTA):
    # Sem normalizar (o TL é uma razão). Pedaços de 0,1 s -> resolução 10 Hz.
    # Os pedaços com o sinal saturado (cortado no topo) são deitados fora.
    sr, data = wavfile.read(os.path.join(pasta, nome))
    data = data.astype(float)
    if data.ndim > 1:
        data = np.mean(data, axis=1)
    n = sr // 10
    topo = np.max(np.abs(data))
    saturado = np.abs(data) >= 0.999 * topo
    if saturado.sum() < 10:          # pico isolado, não é saturação
        saturado[:] = False
    P_soma, usados = 0, 0
    for k in range(0, len(data) - n + 1, n // 2):
        if saturado[k:k + n].any():
            continue
        f, P = welch(data[k:k + n], sr, nperseg=n)
        P_soma, usados = P_soma + P, usados + 1
    total = (len(data) - n) // (n // 2) + 1
    if usados < total:
        print(f"   {nome}: {total - usados} de {total} pedaços saturados ignorados")
    return sr, f, np.sqrt(P_soma / usados)


def largura(freq, TL, i_pico, nivel=NIVEL):
    """Zona contínua à volta de i_pico onde TL >= nivel.
    Devolve (f_esq, f_dir, inclinação esq, inclinação dir) ou None. Inclinações em dB/Hz."""
    if TL[i_pico] < nivel:
        return None
    i = i_pico
    while i > 0 and TL[i - 1] >= nivel:
        i -= 1
    j = i_pico
    while j < len(TL) - 1 and TL[j + 1] >= nivel:
        j += 1
    df = freq[1] - freq[0]
    # Interpolação linear onde a curva cruza o nível; inclinação = declive desse segmento
    if i == 0:
        f_esq, inc_esq = freq[i], np.nan
    else:
        f_esq = np.interp(nivel, [TL[i - 1], TL[i]], [freq[i - 1], freq[i]])
        inc_esq = (TL[i] - TL[i - 1]) / df
    if j == len(TL) - 1:
        f_dir, inc_dir = freq[j], np.nan
    else:
        f_dir = np.interp(nivel, [TL[j + 1], TL[j]], [freq[j + 1], freq[j]])
        inc_dir = (TL[j] - TL[j + 1]) / df
    return f_esq, f_dir, inc_esq, inc_dir


# ---------- Previsões e incertezas das frequências ----------
c = 331.3 * np.sqrt(1 + T / 273.15)
fB = c / (2 * a)
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
# c ∝ √(T em K)  ->  u_c/c = 0,5·u_T/(T + 273,15)
u_rel_T = 0.5 * U_T / (T + 273.15)
u_fB = fB * np.sqrt((U_A / a)**2 + u_rel_T**2)            # fB ∝ c/a


def gargalo(N):
    """(d, u_d, V, S, L_ef) das seringas usadas com N ativas (d = diâmetro da ponta)."""
    v = GARGALO.get(N, V)
    L_ef = L + dL
    if v == 'fusao':
        # f0 = fB  ->  V = S·a² / (π²·L_ef) - V_cone  (como em bragg_vs_local/fusao.py)
        v = S * a**2 / (np.pi**2 * L_ef) - V_CONE
    return D_PONTA, U_D, v, S, L_ef


def prever_f0(N):
    """f0 previsto e incerteza; (nan, nan) se o volume não é conhecido."""
    d, u_d, v, S, L_ef = gargalo(N)
    if v is None:
        return np.nan, np.nan
    f = c / (2 * np.pi) * np.sqrt(S / ((v + V_CONE) * L_ef))
    # f0 ∝ c·√(S/(V_total·L_ef)): em quadratura, metade da incerteza relativa de S, V_total e L_ef
    u_S_rel = np.hypot(U_R / R_B, U_R / R_P)
    return f, f * np.sqrt((0.5 * U_V / (v + V_CONE))**2 + (0.5 * u_S_rel)**2
                          + (0.5 * U_LEF / L_ef)**2 + u_rel_T**2)


f0, u_f0 = prever_f0(1)
print(f"c = {c:.2f} m/s | f0 (gargalo cónico, V = {V * 1e6:.1f} mL) = {f0:.0f} ± {u_f0:.0f} Hz | "
      f"fB = {fB:.0f} ± {u_fB:.0f} Hz")
for N, v in sorted(GARGALO.items()):
    f, uf = prever_f0(N)
    v_N = gargalo(N)[2]
    txt = (f"V = {v_N * 1e6:.1f} mL{' (ajuste f0 = fB)' if v == 'fusao' else ''}, "
           f"f0 = {f:.0f} ± {uf:.0f} Hz" if np.isfinite(f) else "volume desconhecido -> sem f0 previsto")
    print(f"   N = {N}: {txt}")

# Zona onde se procura o pico: só à volta do f0 (0,75 a 1,35 × f0), para todos os N medirem
# a ressonância local e não a banda de Bragg (com N = 4, 5 aparece um pico perto de fB).
F_LO, F_HI = 0.75 * f0, 1.35 * f0


def zona_de(N):
    f = f0s[N] if np.isfinite(f0s[N]) else f0
    return np.flatnonzero((freq >= 0.75 * f) & (freq <= 1.35 * f))

# ---------- Teto ----------
sr_t, freq, A_tref = espectro('A_ref.wav', PASTA_TETO)
_, _, A_tsil = espectro('A_silencio.wav', PASTA_TETO)
_, _, A_tfug = espectro('A_fuga.wav', PASTA_TETO)
idx = (freq >= F_MIN) & (freq <= F_MAX)
freq = freq[idx]
df = freq[1] - freq[0]
teto = 20 * np.log10(np.maximum(A_tref[idx], eps)
                     / np.maximum(np.maximum(A_tsil[idx], A_tfug[idx]), eps))

# ---------- Referências e gravações ----------
# Duas sessões, cada uma com a SUA referência (todos os êmbolos no fundo, lã na ponta):
#   aula 3 (12/10):  A_ref*.wav        + dependencia_N<n>-<rep>.wav
#   14/10:           referencia*.wav   + dependencia_N<n>-<rep>_nova.wav  (N = 4, 5 repetidos)
SESSOES = {'aula 3': 'A_ref*.wav', '14/10': 'referencia*.wav'}
A_refs = {}
for sessao, padrao_ref in SESSOES.items():
    refs = sorted(glob.glob(os.path.join(PASTA, padrao_ref)))
    if refs:
        sr_r, _, A_r = espectro(os.path.basename(refs[0]))
        assert sr_r == sr_t, "Taxa de amostragem diferente da do teto (usar a mesma no Audacity)"
        A_refs[sessao] = np.maximum(A_r[idx], eps)
        print(f"Referência ({sessao}): {os.path.basename(refs[0])}")
if not A_refs:
    raise SystemExit("Falta a referência (A_ref*.wav ou referencia*.wav) nesta pasta.")

# Nomes aceites: N1_1.wav, N1-1.wav, dependencia_N1-1.wav, dependencia_N4-2_nova.wav, ... (N4-1_nova (2).wav)
padrao = re.compile(r'^(?:dependencia_)?N(\d+)[_-](\d+)(_nova)?(?: \(\d+\))?\.wav$', re.IGNORECASE)
TL_por_N = defaultdict(list)
sessao_de = {}
for p in sorted(glob.glob(os.path.join(PASTA, '*.wav'))):
    nome = os.path.basename(p)
    if nome.lower().startswith(('a_ref', 'referencia')):
        continue
    m = padrao.match(nome)
    if not m:
        print(f"(ignorado: {nome}, o nome não é N<n>_<rep>.wav)")
        continue
    sessao = '14/10' if m.group(3) else 'aula 3'
    if sessao not in A_refs:
        print(f"(ignorado: {nome}, falta a referência da sessão {sessao})")
        continue
    N = int(m.group(1))
    if sessao_de.setdefault(N, sessao) != sessao:
        print(f"(ignorado: {nome}, N = {N} já tem gravações da sessão {sessao_de[N]})")
        continue
    sr, _, A = espectro(nome)
    assert sr == sr_t, f"{nome}: taxa de amostragem diferente da do teto"
    TL_por_N[N].append(20 * np.log10(A_refs[sessao] / np.maximum(A[idx], eps)))
A_ref = A_refs.get('aula 3', next(iter(A_refs.values())))

# N = 5: as gravações da pasta bragg_vs_local (fusao_<V>mL_<n>.wav) com o mesmo volume são exatamente "5 seringas ativas"
# (mesma aula, mesma referência). Só entram se não houver N5 nesta pasta.
if 5 not in TL_por_N:
    v_txt = f"{V * 1e6:g}".replace('.', '-')
    for p in sorted(glob.glob(os.path.join(PASTA, '..', 'bragg_vs_local', f'fusao_{v_txt}mL_*.wav'))):
        sr, _, A = espectro(os.path.basename(p), os.path.dirname(p))
        TL_por_N[5].append(20 * np.log10(A_ref / np.maximum(A[idx], eps)))
        print(f"N = 5 tirado de bragg_vs_local: {os.path.basename(p)}")
if not TL_por_N:
    raise SystemExit("Não encontrei ficheiros N<n>_<rep>.wav nesta pasta.")

# ---------- Cada N: TL máximo e largura ----------
f0s = {N: prever_f0(N)[0] for N in TL_por_N}
Ns, TLmax, u_TLmax, f_pico, perto_teto = [], [], [], [], []
larg, u_larg, inc = [], [], []
curvas, curvas_sd = {}, {}
print(f"\nBanda procurada entre {F_LO:.0f} e {F_HI:.0f} Hz (f0 de N = 1)")
for N in sorted(set(f0s) & set(GARGALO)):
    z = zona_de(N)
    print(f"   N = {N}: entre {freq[z[0]]:.0f} e {freq[z[-1]]:.0f} Hz")
print()
print("  N  rep   TL máx (dB)    f pico   largura a 10 dB (Hz)   inclinação esq / dir (dB/Hz)")
for N in sorted(TL_por_N):
    reps = np.array(TL_por_N[N])
    TLm = reps.mean(axis=0)
    curvas[N] = TLm
    curvas_sd[N] = reps.std(axis=0, ddof=1) if len(reps) > 1 else np.zeros_like(TLm)
    zona = zona_de(N)
    i = zona[np.argmax(TLm[zona])]

    # Incerteza do TL máximo: desvio padrão da média das repetições.
    # Com uma só repetição: a dispersão do TL fora da banda (ondulação da montagem).
    if len(reps) > 1:
        u_tl = np.std(reps[:, i], ddof=1) / np.sqrt(len(reps))
    else:
        fora = (freq < F_LO) | (freq > F_HI)
        u_tl = np.std(TLm[fora])

    # Largura de cada repetição -> dispersão; e da média -> valor
    lim = largura(freq, TLm, i)
    larg_reps = [lr[1] - lr[0] for lr in (largura(freq, t, i) for t in reps) if lr]
    if lim:
        w = lim[1] - lim[0]
        u_rep = np.std(larg_reps, ddof=1) / np.sqrt(len(larg_reps)) if len(larg_reps) > 1 else 0
        # Resolução: cada margem é lida numa grelha de 10 Hz -> 10/√12 por margem
        u_res = np.sqrt(2) * df / np.sqrt(12)
        # Volume (0,5 mL), diâmetro do gargalo e posição (1 mm) espalham o f0 e o fB de
        # cada seringa; como os erros das N seringas são independentes, a média cai com √N.
        # Sem f0 previsto (volume desconhecido) usa-se a frequência do pico.
        # A temperatura muda a escala de frequências toda (largura ∝ c).
        d_N, u_d_N, v_N = gargalo(N)[:3]
        f0_N = f0s[N] if np.isfinite(f0s[N]) else freq[i]
        u_rel_f0 = np.hypot(0.5 * U_V / v_N if v_N else 0, 0.5 * u_d_N / d_N)
        u_w = np.sqrt(u_rep**2 + u_res**2
                      + ((u_rel_f0 * f0_N)**2 + (U_A / a * fB)**2) / N
                      + (w * u_rel_T)**2)
    else:
        w, u_w = np.nan, np.nan
        lim = (np.nan, np.nan, np.nan, np.nan)

    Ns.append(N)
    TLmax.append(TLm[i])
    u_TLmax.append(u_tl)
    f_pico.append(freq[i])
    perto_teto.append(TLm[i] >= teto[i] - 3)
    larg.append(w)
    u_larg.append(u_w)
    inc.append((lim[2], lim[3]))
    aviso = "  <-- perto do teto" if perto_teto[-1] else ""
    print(f"  {N}  {len(reps)}    {TLm[i]:5.1f} ± {u_tl:4.1f}   {freq[i]:5.0f} Hz   "
          f"{w:6.0f} ± {u_w:4.0f}  ({lim[0]:.0f}–{lim[1]:.0f})     "
          f"{lim[2]:+.2f} / {lim[3]:+.2f}{aviso}")
    if not np.isfinite(f0s[N]):
        # Volume na seringa que põe o f0 no pico medido: V = S c² / ((2π f)² L_ef) - V_cone
        d_N, u_d_N, _, S_N, Lef_N = gargalo(N)
        V_ef = S_N * c**2 / ((2 * np.pi * freq[i])**2 * Lef_N) - V_CONE
        u_Vef = V_ef * np.sqrt((u_d_N / d_N)**2 + (2 * df / np.sqrt(12) / freq[i])**2
                               + (2 * u_rel_T)**2)
        print(f"       (o pico a {freq[i]:.0f} Hz corresponde a "
              f"V ≈ {V_ef * 1e6:.1f} ± {u_Vef * 1e6:.1f} mL se for o f0)")

Ns, TLmax, u_TLmax = np.array(Ns), np.array(TLmax), np.array(u_TLmax)
larg, u_larg = np.array(larg), np.array(u_larg)
perto_teto = np.array(perto_teto)

# ---------- O crescimento do TL máximo é linear? ----------
# Ajuste TL = k·N + b (reta com ordenada na origem) só com os pontos longe do teto.
# Não se força b = 0: com N = 0 o TL é 0, mas a 1.ª seringa pode não seguir a reta.
bons = ~perto_teto
if bons.sum() < 2:
    # O teto de manhã pode não servir para esta sessão (ganho/montagem diferentes):
    # nesse caso ajusta-se com todos os pontos e avisa-se.
    print("\nAVISO: quase todos os pontos estão 'perto do teto' de manhã -> ajuste com todos os N")
    bons = np.ones_like(perto_teto)
k = b = np.nan
if bons.sum() >= 2:
    x, y = Ns[bons].astype(float), TLmax[bons]
    w_fit = 1 / np.maximum(u_TLmax[bons], 0.1)**2
    # Mínimos quadrados pesados: matriz normal 2x2
    Sw, Sx, Sy = w_fit.sum(), (w_fit * x).sum(), (w_fit * y).sum()
    Sxx, Sxy = (w_fit * x * x).sum(), (w_fit * x * y).sum()
    D = Sw * Sxx - Sx**2
    k = (Sw * Sxy - Sx * Sy) / D
    b = (Sxx * Sy - Sx * Sxy) / D
    u_k, u_b = np.sqrt(Sw / D), np.sqrt(Sxx / D)
    res = y - (k * x + b)
    if len(x) > 2:
        # Se os pontos se espalham mais do que as barras de erro, a incerteza cresce com isso
        chi2_red = np.sum(w_fit * res**2) / (len(x) - 2)
        u_k, u_b = u_k * max(1, np.sqrt(chi2_red)), u_b * max(1, np.sqrt(chi2_red))
    else:
        chi2_red = np.nan
    print(f"\nAjuste TL máx = k·N + b (N usados = {Ns[bons].tolist()}):")
    print(f"   k = {k:.2f} ± {u_k:.2f} dB por seringa   b = {b:.1f} ± {u_b:.1f} dB   χ²/ν = {chi2_red * (len(x) - 2):.1f}/{len(x) - 2} = {chi2_red:.2f}")
    if perto_teto.any() and not bons.all():
        print(f"   pontos perto do teto (fora do ajuste): N = {Ns[perto_teto].tolist()}")
else:
    print("\nPoucos pontos longe do teto para testar se o crescimento é linear.")

# ---------- Gráficos ----------
fig = plt.figure(figsize=(14, 11))
gs = fig.add_gridspec(3, 2, height_ratios=[3, 2.2, 0.9])
ax_tl = fig.add_subplot(gs[0, :])   # curvas TL(f) na linha de cima, a toda a largura
ax_max = fig.add_subplot(gs[1, 0])
ax_res = fig.add_subplot(gs[2, 0], sharex=ax_max)   # resíduos do ajuste
ax_larg = fig.add_subplot(gs[1:, 1])

cores = plt.cm.viridis(np.linspace(0, 0.9, len(curvas)))
for cor, (N, TLm) in zip(cores, curvas.items()):
    ax_tl.fill_between(freq, TLm - curvas_sd[N], TLm + curvas_sd[N], color=cor, alpha=0.12, linewidth=0)
    ax_tl.plot(freq, TLm, color=cor, linewidth=1.4, label=f'N = {N}')
# O teto de 7/10 já não se desenha: é de outra sessão (ver teto_sistema).
ax_tl.axhline(NIVEL, color='black', linewidth=0.8, linestyle=':', label=f'TL = {NIVEL:.0f} dB (largura da banda)')
ax_tl.axvline(f0, color='red', linestyle=':', label=f'$f_0$ previsto = {f0:.0f} Hz')
for f in sorted({round(f0s[N]) for N in GARGALO if N in f0s and np.isfinite(f0s[N])}):
    ax_tl.axvline(f, color='orange', linestyle=':', label=f'$f_0$ = {f} Hz (outro volume)')
ax_tl.axvline(fB, color='green', linestyle=':', label=f'$f_B$ = {fB:.0f} Hz')
ax_tl.set_xlim(F_MIN, min(F_MAX, 2 * max(f0, fB)))
ax_tl.set_xlabel('Frequência (Hz)')
ax_tl.set_ylabel('TL (dB)')
ax_tl.set_title('TL(f) para cada número de seringas N (média de 3, faixa = desvio padrão)')
ax_tl.grid(True, alpha=0.3)
ax_tl.legend(fontsize=8)

# Um símbolo por sessão (cada uma com a sua referência). O "perto do teto" não se desenha:
# o teto é de 7/10, outra sessão, e marcava quase todos os pontos.
sess = np.array([sessao_de.get(N, 'aula 3') for N in Ns])
for nome_s, fmt, cor in (('aula 3', 'o', 'black'), ('14/10', 's', 'tab:blue')):
    m_s = sess == nome_s
    if m_s.any():
        ax_max.errorbar(Ns[m_s], TLmax[m_s], yerr=u_TLmax[m_s], fmt=fmt, color=cor, capsize=4,
                        label=f'Sessão {nome_s}')
        ax_larg.errorbar(Ns[m_s & np.isfinite(larg)], larg[m_s & np.isfinite(larg)],
                         yerr=u_larg[m_s & np.isfinite(larg)], fmt=fmt, color=cor, capsize=4,
                         label=f'Sessão {nome_s}')
if np.isfinite(k):
    NN = np.linspace(Ns.min(), Ns.max(), 50)
    ax_max.plot(NN, k * NN + b, color='red',
                label=f'ajuste TL = k·N + b\nk = {k:.1f} ± {u_k:.1f} dB,  b = {b:.1f} ± {u_b:.1f} dB\n'
                      f'χ²/ν = {chi2_red * (bons.sum() - 2):.1f}/{bons.sum() - 2}')
    for nome_s, fmt, cor in (('aula 3', 'o', 'black'), ('14/10', 's', 'tab:blue')):
        m_s = (sess == nome_s) & bons
        ax_res.errorbar(Ns[m_s], TLmax[m_s] - (k * Ns[m_s] + b), yerr=u_TLmax[m_s], fmt=fmt,
                        color=cor, capsize=4)
    ax_res.axhline(0, color='red', linewidth=1)
ax_res.set_xlabel('N (seringas ativas)')
ax_res.set_ylabel('resíduo (dB)')
ax_res.grid(True, alpha=0.3)
plt.setp(ax_max.get_xticklabels(), visible=False)
ax_max.set_xticks(Ns)
ax_max.set_ylabel('TL máximo (dB)')
ax_max.set_title('Atenuação máxima vs N')
ax_max.grid(True, alpha=0.3)
ax_max.legend(fontsize=8, loc='upper left')

ax_larg.legend()
ax_larg.set_xlabel('N (seringas ativas)')
ax_larg.set_xticks(Ns)
ax_larg.set_ylabel(f'Largura a {NIVEL:.0f} dB (Hz)')
ax_larg.set_title(f'Largura da banda vs N')
ax_larg.grid(True, alpha=0.3)

estilo_tl.virgulas(fig)
plt.tight_layout()
plt.savefig(os.path.join(PASTA, 'tl_vs_N.png'), dpi=150)
plt.show()
