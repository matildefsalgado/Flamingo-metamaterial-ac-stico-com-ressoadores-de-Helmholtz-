import os
import re
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import butter, sosfiltfilt, hilbert

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import estilo_tl  # noqa: E402  (cores comuns)

# Decaimento da ressonância de uma seringa: o Q pelo tempo, em vez de pela largura do pico.
# Tom contínuo a ~f0, cortado de repente: a seringa continua a vibrar à sua frequência natural
# e a amplitude cai como exp(−t/τ), com τ = Q/(π f0). Em dB a queda é uma reta:
#     declive = −20/(τ·ln 10) dB/s   ->   τ = −8,686/declive,   Q = π·f0·τ
# O próprio sistema (altifalante, tubo, lã) também demora a calar-se (ida e volta no tubo ~10 ms),
# por isso grava-se o mesmo com os êmbolos todos no fundo (controlo). O decaimento medido nunca é
# mais rápido do que o do sistema: se τ_seringa não for claramente maior do que τ_controlo,
# o Q que sai é só um limite superior.
#
# Gravações (nesta pasta, mesma sessão, mic na posição 2, lã na ponta):
#   dec_seringa_<n>.wav    seringa 1 a 14,6 mL, as outras no fundo; tom ~2 s ligado / ~1 s desligado, 5×
#   dec_controlo_<n>.wav   igual, com todos os êmbolos no fundo
# Cada gravação pode ter vários cortes: o código encontra-os sozinho.

# ---------- Ficheiros ----------
# Também aceita os nomes usados a 22/10: 'decaimento 1.wav' (seringa, tom a 430 Hz) e 'controlo (1).wav'.
# 'decaimento.wav' (tom a 600 Hz, frequência errada) e 'A_referencia.wav' (varrimento) não entram.
PADRAO_SER = r'(?:dec_seringa|decaimento 1)'
PADRAO_CTRL = r'(?:dec_controlo|controlo)'

# ---------- Parâmetros ----------
F0 = 430.0          # f0 da seringa 1 a 14,6 mL (Hz): medido na aula 3 (430) e a 14/10 (435)
U_F0 = 5.0          # incerteza de f0 (Hz): metade da resolução de 10 Hz
Q_LARG = 13.2       # Q pela largura do pico, N1 da aula 3 (12/10)
Q_DIST, U_Q_DIST = 15.1, 1.1   # Q pela largura do pico, seringa 1 sozinha a 14/10 (distancia_altifalante)
BANDA = 150.0       # meia largura do filtro passa-banda à volta de f0 (Hz): tira ruído e o estalo do corte
SUAVE = 1e-3        # suavização da envolvente (s): muito menor que τ (~10 ms)
INICIO_DB = 3.0     # o ajuste começa 3 dB abaixo do nível com o tom ligado (fora da transição)
FIM_DB = 25.0       # e acaba 25 dB abaixo, ou 10 dB acima do ruído, o que vier primeiro
MARGEM_RUIDO = 10.0
MIN_LIGADO = 0.5    # um corte só conta se o tom esteve ligado pelo menos 0,5 s antes (s)
MIN_DESLIGADO = 0.3 # e ficou desligado pelo menos 0,3 s depois (s)

sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
PASTA = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(description='Decaimento da ressonância (Q pelo tempo)')
ap.add_argument('--pasta', default=PASTA, help='pasta com as gravações (por omissão, esta)')
ap.add_argument('--f0', type=float, default=F0, help=f'frequência do tom / f0 (Hz), por omissão {F0:.0f}')
ap.add_argument('--banda', type=float, default=BANDA,
                help=f'meia largura do filtro à volta de f0 (Hz), por omissão {BANDA:.0f}')
ap.add_argument('--teste', action='store_true',
                help='corre com sinais sintéticos (τ conhecido) em vez das gravações')
args = ap.parse_args()
PASTA_DADOS = args.pasta
F0 = args.f0
BANDA = args.banda


def ficheiros(prefixo):
    """Gravações <prefixo><números>.wav, com _ ou - (aceita o ' (1)' que o Windows acrescenta)."""
    p = re.compile(rf'^{prefixo}[_-]?[\d_-]*(?: \(\d+\))?\.wav$', re.IGNORECASE)
    return sorted(f for f in os.listdir(PASTA_DADOS) if p.match(f))


def ler(nome):
    sr, x = wavfile.read(os.path.join(PASTA_DADOS, nome))
    x = x.astype(float)
    if x.ndim > 1:
        x = x.mean(axis=1)
    return x - x.mean(), sr


def media_movel(v, n):
    """Média móvel centrada de n amostras (com somas acumuladas, rápida em gravações longas)."""
    if n <= 1:
        return v
    c = np.cumsum(np.concatenate([[0.0], v]))
    m = (c[n:] - c[:-n]) / n
    return np.concatenate([np.full(n // 2, m[0]), m, np.full(n - 1 - n // 2, m[-1])])


def envolvente_db(x, sr):
    """Envolvente em dB: passa-banda à volta de f0 (fase zero), Hilbert, média de 1 ms."""
    sos = butter(2, [F0 - BANDA, F0 + BANDA], btype='bandpass', fs=sr, output='sos')
    env = np.abs(hilbert(sosfiltfilt(sos, x)))
    env = media_movel(env, int(SUAVE * sr))
    # Chão a −120 dB do máximo: o Audacity grava zeros exatos nos silêncios e log(0) estraga as medianas
    return 20 * np.log10(np.maximum(env, 1e-6 * env.max()))


def achar_cortes(L, sr):
    """Índices onde o tom passa de ligado a desligado (envolvente suavizada a 50 ms, limiar 20 dB abaixo do tom)."""
    Lg = 20 * np.log10(media_movel(10**(L / 20), int(0.05 * sr)) + 1e-12)
    alto, baixo = np.percentile(Lg, 90), np.percentile(Lg, 10)
    if alto - baixo < 15:
        return []
    ligado = Lg > alto - 20
    bordas = np.flatnonzero(ligado[:-1] & ~ligado[1:]) + 1
    n_on, n_off = int(MIN_LIGADO * sr), int(MIN_DESLIGADO * sr)
    return [b for b in bordas
            if b - n_on >= 0 and b + n_off <= len(L)
            and ligado[b - n_on:b].mean() > 0.9 and (~ligado[b:b + n_off]).mean() > 0.9]


def ajustar(L, sr, b):
    """Reta à queda em dB depois de um corte. Devolve τ (s) e o que é preciso para o gráfico."""
    ms = sr // 1000
    L_on = np.median(L[b - 250 * ms:b - 50 * ms])
    L_ruido = np.median(L[b + 150 * ms:b + 300 * ms])
    # O corte verdadeiro: a última amostra ainda a menos de 1 dB do nível ligado
    zona = np.arange(b - 60 * ms, b + 100 * ms)
    acima = zona[L[zona] > L_on - 1]
    if len(acima) == 0:
        return None
    k0 = acima[-1]
    L_fim = max(L_on - FIM_DB, L_ruido + MARGEM_RUIDO)
    if L_on - L_fim < INICIO_DB + 5:
        return None
    depois = np.arange(k0, b + 300 * ms)
    i0 = depois[np.argmax(L[depois] < L_on - INICIO_DB)]
    i1 = i0 + np.argmax(L[i0:b + 300 * ms] < L_fim)
    if i1 - i0 < 10:
        return None
    t = np.arange(i0, i1) / sr
    p = np.polyfit(t, L[i0:i1], 1)
    return dict(tau=-20 / (np.log(10) * p[0]), declive=p[0], L_on=L_on, L_ruido=L_ruido,
                faixa=L_on - L_fim, k0=k0, i0=i0, i1=i1, p=p)


def analisar(gravacoes):
    """gravacoes = [(nome, x, sr)]. Devolve a lista de ajustes de todos os cortes."""
    ajustes = []
    for nome, x, sr in gravacoes:
        L = envolvente_db(x, sr)
        cortes = achar_cortes(L, sr)
        bons = [a for a in (ajustar(L, sr, b) for b in cortes) if a is not None]
        print(f"   {nome}: {len(cortes)} cortes, {len(bons)} ajustados"
              + (f"  (τ = {', '.join(f'{a['tau'] * 1e3:.1f}' for a in bons)} ms)" if bons else ""))
        for a in bons:
            a.update(nome=nome, L=L, sr=sr)
            print(f"      queda ajustada de −{INICIO_DB:.0f} a −{a['faixa']:.0f} dB, "
                  f"ruído {a['L_on'] - a['L_ruido']:.0f} dB abaixo do tom")
        ajustes += bons
    return ajustes


def resumo(ajustes, nome):
    taus = np.array([a['tau'] for a in ajustes])
    m = taus.mean()
    u = taus.std(ddof=1) / np.sqrt(len(taus)) if len(taus) > 1 else 0.0
    print(f"{nome}: τ = {m * 1e3:.2f} ± {u * 1e3:.2f} ms  ({len(taus)} cortes, "
          f"T60 = {6.908 * m * 1e3:.0f} ms)")
    return m, u


def sintetico(tau, f, sr=48000, n=5, ruido_db=-50, seed=0):
    """Tom de f0 cortado n vezes (2 s ligado / 1 s desligado); depois do corte decai com τ."""
    rng = np.random.default_rng(seed)
    bloco = np.arange(int(3 * sr)) / sr
    ligado = bloco < 2
    env = np.where(ligado, 1.0, np.exp(-np.clip(bloco - 2, 0, None) / tau))
    x = np.tile(env * np.sin(2 * np.pi * f * bloco), n)
    return x + 10**(ruido_db / 20) * rng.standard_normal(len(x)), sr


# ---------- Limite do método: o filtro e a suavização também "tocam" ----------
x_f, sr_f = sintetico(1e-6, F0, n=2, ruido_db=-90)
L_f = envolvente_db(x_f, sr_f)
lim = [ajustar(L_f, sr_f, b) for b in achar_cortes(L_f, sr_f)]
tau_filtro = np.mean([a['tau'] for a in lim if a is not None])
print(f"f0 = {F0:.0f} ± {U_F0:.0f} Hz | filtro {F0 - BANDA:.0f}–{F0 + BANDA:.0f} Hz | "
      f"τ mínimo que o método consegue ver (corte instantâneo): {tau_filtro * 1e3:.2f} ms")

# ---------- Gravações ----------
if args.teste:
    TAU_S, TAU_C = 10e-3, 3e-3
    print(f"\nTESTE com sinais sintéticos: seringa τ = {TAU_S * 1e3:.0f} ms, controlo τ = {TAU_C * 1e3:.0f} ms")
    grav_s = [('sintético seringa', *sintetico(TAU_S, F0 - 8, seed=1))]
    grav_c = [('sintético controlo', *sintetico(TAU_C, F0 - 8, seed=2))]
else:
    ns, nc = ficheiros(PADRAO_SER), ficheiros(PADRAO_CTRL)
    if not ns:
        raise SystemExit(f"Faltam gravações em {PASTA_DADOS}: preciso de {PADRAO_SER}_1.wav, ... "
                         f"(e {PADRAO_CTRL}_1.wav, ... para o controlo). Ver LEIA-ME.")
    grav_s = [(n, *ler(n)) for n in ns]
    grav_c = [(n, *ler(n)) for n in nc]

print("\nSeringa:")
aj_s = analisar(grav_s)
print("Controlo:" if grav_c else "Controlo: (não há: sem ele não se sabe que parte do decaimento é da seringa)")
aj_c = analisar(grav_c) if grav_c else []
if not aj_s:
    raise SystemExit("Nenhum corte ajustado nas gravações da seringa (o tom foi mesmo cortado de repente?).")

print()
tau_s, u_s = resumo(aj_s, "Seringa ")
Q = np.pi * F0 * tau_s
u_Q = Q * np.hypot(u_s / tau_s, U_F0 / F0)
if aj_c:
    tau_c, u_c = resumo(aj_c, "Controlo")
    razao = tau_s / tau_c
    u_dif = max(np.hypot(u_s, u_c), 1e-5)
    print(f"τ_seringa / τ_controlo = {razao:.2f}   (diferença {(tau_s - tau_c) * 1e3:.2f} ± {u_dif * 1e3:.2f} ms, "
          f"z = {(tau_s - tau_c) / u_dif:.1f})")
    if razao < 2:
        print("   AVISO: o decaimento da seringa não é claramente mais lento que o do sistema:\n"
              "   o Q abaixo é só um LIMITE SUPERIOR.")
for a in aj_s + aj_c:
    if a['tau'] < 3 * tau_filtro:
        print(f"   AVISO: {a['nome']}: τ = {a['tau'] * 1e3:.1f} ms está perto do limite do método "
              f"({tau_filtro * 1e3:.1f} ms)")
        break

print(f"\nQ = π·f0·τ = {Q:.1f} ± {u_Q:.1f}" + ("   (limite superior)" if aj_c and razao < 2 else ""))
print(f"   vs largura do pico, aula 3 (12/10): Q = {Q_LARG:.1f}  ({100 * (Q / Q_LARG - 1):+.0f} %)")
print(f"   vs largura do pico, 14/10:          Q = {Q_DIST:.1f} ± {U_Q_DIST:.1f}  "
      f"(z = {abs(Q - Q_DIST) / np.hypot(u_Q, U_Q_DIST):.1f})")
if args.teste:
    print(f"\nTESTE: recuperado seringa {tau_s * 1e3:.2f} ms (verdadeiro {TAU_S * 1e3:.0f}), "
          f"controlo {tau_c * 1e3:.2f} ms (verdadeiro {TAU_C * 1e3:.0f})")

# ---------- Gráficos ----------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw=dict(width_ratios=[3, 2]))
for ajs, cor, rot in [(aj_s, estilo_tl.AZUL, 'seringa 1 a 14,6 mL'), (aj_c, estilo_tl.LARANJA, 'controlo (êmbolos no fundo)')]:
    for j, a in enumerate(ajs):
        ms = a['sr'] / 1000
        idx = np.arange(int(a['k0'] - 20 * ms), int(a['k0'] + 150 * ms))
        t = (idx - a['k0']) / ms
        ax1.plot(t, a['L'][idx] - a['L_on'], color=cor, alpha=0.5, linewidth=0.8,
                 label=rot if j == 0 else None)
        ti = np.arange(a['i0'], a['i1'])
        ax1.plot((ti - a['k0']) / ms, np.polyval(a['p'], ti / a['sr']) - a['L_on'], color=cor,
                 linestyle='--', linewidth=1.5)
ax1.axhline(-INICIO_DB, color='#999', linewidth=0.8, linestyle=':')
ax1.set_xlim(-20, 150)
ax1.set_ylim(-50, 5)
ax1.set_xlabel('tempo depois do corte (ms)')
ax1.set_ylabel('nível relativo ao tom ligado (dB)')
txt = f'τ seringa = {tau_s * 1e3:.1f} ± {u_s * 1e3:.1f} ms,  Q = {Q:.1f} ± {u_Q:.1f}'
if aj_c and razao < 2:
    txt += ' (limite superior)'
if aj_c:
    txt += f'\nτ controlo = {tau_c * 1e3:.1f} ± {u_c * 1e3:.1f} ms'
ax1.set_title(f'Decaimento depois de cada corte (f0 = {F0:.0f} Hz); tracejado = ajuste')
ax1.text(0.97, 0.95, txt, transform=ax1.transAxes, ha='right', va='top', fontsize=10,
         bbox=dict(facecolor='white', edgecolor='#ccc'))
ax1.legend(loc='lower left')
ax1.grid(alpha=0.3)

a0 = aj_s[0]
L0, sr0 = a0['L'], a0['sr']
passo = max(1, sr0 // 2000)
ax2.plot(np.arange(0, len(L0), passo) / sr0, L0[::passo] - a0['L_on'], color=estilo_tl.AZUL, linewidth=0.6)
for a in aj_s:
    if a['nome'] == a0['nome']:
        ax2.axvline(a['k0'] / sr0, color=estilo_tl.LARANJA, linewidth=0.8)
ax2.set_ylim(-60, 10)
ax2.set_xlabel('tempo (s)')
ax2.set_ylabel('nível relativo ao tom ligado (dB)')
ax2.set_title('seringa: gravação inteira e cortes encontrados (laranja)')
ax2.grid(alpha=0.3)
estilo_tl.virgulas(fig)
fig.tight_layout()
NOME_FIG = os.path.join(PASTA if args.teste else PASTA_DADOS, 'decaimento_teste.png' if args.teste else 'decaimento.png')
fig.savefig(NOME_FIG, dpi=150)
print(f"\nFigura: {NOME_FIG}")
plt.show()
