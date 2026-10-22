"""
Ressonância local de UMA seringa em função da distância ao altifalante.

Volume fixo, só uma seringa ativa de cada vez (as outras com o êmbolo no fundo).
A única coisa que muda entre gravações é a distância d da seringa ativa ao altifalante.

Para cada gravação:
    TL(f) = 20 log10( A_ref(f) / A(f) )      (A_ref = êmbolos todos no fundo, MESMA sessão)
e daí:
    f0     = frequência do máximo de TL (perto do f0 previsto)
    TL_max = valor desse máximo
    Q      = f0 / (largura do pico a TL_max - 3 dB)

O que se espera:
  - f0 NÃO deve depender de d: é uma propriedade da seringa (volume + gargalo).
    Testa-se com a média ponderada e o chi2 (f0 constante?) e com o declive df0/dd (compatível com 0?).
  - TL_max PODE depender de d: se a ponta do tubo refletir, há ondas estacionárias e a seringa
    funciona melhor num máximo de pressão do que num nodo. Esse padrão repete-se a cada
    lambda/2 = c / (2 f0) (o script escreve este valor para comparar).

Ficheiros (nesta pasta, todos da MESMA sessão, sem mexer em volume/ganho/microfone):
    A_ref.wav          todos os êmbolos no fundo, lã na ponta (gravar no INÍCIO)
    A_ref_fim.wav      igual, gravado no FIM (controlo: tem de coincidir com A_ref)  [opcional]
    A_silencio.wav     5 s com o altifalante desligado (teto da própria sessão)       [opcional]
    d<cm>_<rep>.wav    uma seringa a V, à distância d (em cm) do altifalante; rep = 1, 2, 3
                       ex.: d30_1.wav, d55.6_2.wav (também aceita vírgula: d55,6_2.wav)
    ou distancia_<k>-<rep>.wav   seringa k ativa (k = 1 a 5); d = D1 + (k - 1)·A_ESP
    (a referência também pode chamar-se referencia*.wav)

Uso:  python distancia_altifalante.py
      python distancia_altifalante.py --pasta OUTRA_PASTA     (para testar com outros ficheiros)
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import welch
from scipy.ndimage import uniform_filter1d

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import estilo_tl  # noqa: E402  (vírgula decimal nas figuras)

# ---------- Parâmetros (SI) ----------
T = 22             # temperatura da sala (°C), a do registo (igual em todos os scripts)
U_T = 0.5          # incerteza da temperatura (°C)
V = 14.6e-6        # volume FIXO da seringa ativa (m³), o mesmo da série de N
U_V = 0.5e-6       # incerteza de leitura do volume (m³): metade da menor divisão (1 mL)
U_R = 0.1e-3       # incerteza dos raios do cone (m): craveira + parede
U_LEF = 2.5e-3     # incerteza de L_ef (m), da calibração (correcao_gargalo --sem 10)
U_D = 0.002        # incerteza da distância medida com a fita (m)
A_ESP = 0.25       # espaçamento entre seringas (m), para os ficheiros distancia_<k>-<rep>.wav
D1 = 0.20          # distância da seringa 1 ao altifalante (m); None = eixo relativo à seringa 1
F_MIN, F_MAX = 150.0, 1500.0
NPERSEG = 2**15                        # resolução ~1,5 Hz
LARGURAS_HZ = [20, 30, 40, 50, 60, 70]  # suavização em frequência (o tubo dá uma ondulação de ~30 Hz)
JANELA = 0.35                          # procura o pico em f0_previsto ± 35 %

# Gargalo = bico cónico da seringa. Valores acertados aos dados de correcao_gargalo.
L = 0.035              # comprimento do cone, no eixo (m)
D_BASE = 0.0132        # diâmetro exterior na base (m)
PAREDE = 0.0012        # grossura da parede na base (m)
D_PONTA = 0.0041       # diâmetro interior na ponta (m)
R_B, R_P = D_BASE / 2 - PAREDE, D_PONTA / 2
S = np.pi * R_B * R_P  # área do cilindro equivalente ao cone
dL = 0.73 * (R_B + R_P)   # correção de extremidade (m)
V_CONE = np.pi * L / 3 * (R_B**2 + R_B * R_P + R_P**2)   # ar no cone (m³): soma-se ao êmbolo

sys.stdout.reconfigure(encoding='utf-8')


def espectro(caminho):
    """Espectro de amplitude (Welch). Sem normalizar: o TL é uma razão."""
    sr, x = wavfile.read(caminho)
    x = x.astype(float)
    if x.ndim > 1:
        x = np.mean(x, axis=1)
    f, P = welch(x, sr, nperseg=NPERSEG)
    return f, np.sqrt(P)


def na_grelha(f_dest, caminho):
    """Espectro interpolado para a grelha f_dest (as taxas de amostragem podem ser diferentes)."""
    f, A = espectro(caminho)
    return np.interp(f_dest, f, A)


def suavizar(f, y, largura):
    return uniform_filter1d(y, max(1, int(round(largura / (f[1] - f[0])))))


def largura_3dB(f, TL, i):
    """Largura do pico à volta de i, a TL[i] - 3 dB, interpolada. None se não desce 3 dB."""
    nivel = TL[i] - 3
    e, d = i, i
    while e > 0 and TL[e - 1] >= nivel:
        e -= 1
    while d < len(TL) - 1 and TL[d + 1] >= nivel:
        d += 1
    if e == 0 or d == len(TL) - 1:
        return None
    f_e = np.interp(nivel, [TL[e - 1], TL[e]], [f[e - 1], f[e]])
    f_d = np.interp(nivel, [TL[d + 1], TL[d]], [f[d + 1], f[d]])
    return f_d - f_e


def media_ponderada(y, u):
    w = 1 / u**2
    m = np.sum(w * y) / np.sum(w)
    chi2 = np.sum(w * (y - m)**2)
    return m, 1 / np.sqrt(np.sum(w)), chi2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pasta', default=os.path.dirname(os.path.abspath(__file__)))
    pasta = ap.parse_args().pasta

    # ---------- Previsões ----------
    c = 331.3 * np.sqrt(1 + T / 273.15)
    u_rel_c = 0.5 * U_T / (T + 273.15)
    f0_prev = c / (2 * np.pi) * np.sqrt(S / ((V + V_CONE) * (L + dL)))
    # f0 ∝ c·√(S/(V_total·L_ef)): em quadratura, metade da incerteza relativa de S, V_total e L_ef
    u_S_rel = np.hypot(U_R / R_B, U_R / R_P)
    u_f0_prev = f0_prev * np.sqrt((0.5 * U_V / (V + V_CONE))**2 + (0.5 * u_S_rel)**2
                                  + (0.5 * U_LEF / (L + dL))**2 + u_rel_c**2)
    meia_onda = c / (2 * f0_prev)
    print(f"c = {c:.1f} m/s | V = {V * 1e6:.1f} mL | f0 previsto = {f0_prev:.0f} ± {u_f0_prev:.0f} Hz | "
          f"λ/2 a f0 = {meia_onda * 100:.1f} cm (período esperado do padrão de ondas estacionárias)")

    # ---------- Gravações ----------
    padrao_cm = re.compile(r'^d([\d.,]+)_(\d+)\.wav$', re.IGNORECASE)
    padrao_k = re.compile(r'^distancia_(\d+)-(\d+)\.wav$', re.IGNORECASE)
    medidas = {}
    for p in sorted(glob.glob(os.path.join(pasta, '*.wav'))):
        nome = os.path.basename(p)
        m = padrao_cm.match(nome)
        if m:
            d = float(m.group(1).replace(',', '.')) / 100       # cm -> m
        elif padrao_k.match(nome):
            # seringa k: distância relativa à seringa 1 (ou absoluta, se D1 for conhecido)
            d = (D1 or 0) + (int(padrao_k.match(nome).group(1)) - 1) * A_ESP
        else:
            continue
        medidas.setdefault(round(d, 4), []).append(p)
    refs = [os.path.join(pasta, 'A_ref.wav')] + sorted(glob.glob(os.path.join(pasta, 'referencia*.wav')))
    refs = [r for r in refs if os.path.exists(r)]
    if not refs:
        raise SystemExit("Falta A_ref.wav ou referencia*.wav (êmbolos no fundo, mesma sessão).")
    caminho_ref = refs[0]
    print(f"Referência: {os.path.basename(caminho_ref)}")
    if not medidas:
        raise SystemExit("Não encontrei gravações d<cm>_<rep>.wav nem distancia_<k>-<rep>.wav.")

    f, A_ref = espectro(caminho_ref)
    banda = (f >= F_MIN) & (f <= F_MAX)
    zona = banda & (f >= f0_prev * (1 - JANELA)) & (f <= f0_prev * (1 + JANELA))
    # "fora": longe do pico, onde a seringa quase não atua -> TL deve dar ~0 dB
    fora = banda & ((f < f0_prev * 0.6) | (f > f0_prev * 1.6))

    # ---------- Controlos da sessão ----------
    caminho_fim = os.path.join(pasta, 'A_ref_fim.wav')
    if os.path.exists(caminho_fim):
        dif = suavizar(f, 20 * np.log10(na_grelha(f, caminho_fim) / A_ref), 30)[banda]
        ok = abs(dif.mean()) < 1 and dif.std() < 1.5
        print(f"Controlo A_ref_fim vs A_ref: desvio médio {dif.mean():+.1f} dB, ondulação {dif.std():.1f} dB"
              + ("  -> OK" if ok else "  -> ATENÇÃO: a montagem mudou durante a sessão"))
    else:
        print("(sem A_ref_fim.wav: não dá para confirmar que nada mudou durante a sessão)")
    teto = None
    caminho_sil = os.path.join(pasta, 'A_silencio.wav')
    if os.path.exists(caminho_sil):
        teto = suavizar(f, 20 * np.log10(A_ref / na_grelha(f, caminho_sil)), 30)

    # ---------- Cada distância ----------
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    cores = plt.cm.viridis(np.linspace(0, 0.9, len(medidas)))
    linhas = []
    print("\n d (cm)  f0 por repetição (Hz)        f0 médio      TL_max (dB)    Q          TL fora (dB)")
    for cor, d in zip(cores, sorted(medidas)):
        f0s, sist, picos, Qs, foras = [], [], [], [], []
        TL_media = 0
        TL_reps = []
        for p in medidas[d]:
            TL_bruto = 20 * np.log10(A_ref / na_grelha(f, p))
            # f0: máximo com várias suavizações (a ondulação do tubo puxa o pico); média e dispersão
            est = [f[np.argmax(np.where(zona, suavizar(f, TL_bruto, W), -np.inf))] for W in LARGURAS_HZ]
            f0s.append(np.mean(est))
            sist.append(np.std(est))
            TL = suavizar(f, TL_bruto, 30)
            i = np.argmax(np.where(zona, TL, -np.inf))
            picos.append(TL[i])
            w3 = largura_3dB(f, TL, i)
            Qs.append(f[i] / w3 if w3 else np.nan)
            foras.append(TL[fora].mean())
            TL_media = TL_media + TL / len(medidas[d])
            TL_reps.append(TL)
        n = len(f0s)
        f0s, picos, Qs = np.array(f0s), np.array(picos), np.array(Qs)
        sd = lambda v: np.nanstd(v, ddof=1) / np.sqrt(np.sum(np.isfinite(v))) if np.sum(np.isfinite(v)) > 1 else np.nan
        u_f0 = np.sqrt(np.nan_to_num(sd(f0s))**2 + np.mean(np.square(sist)) + 1.0**2)   # +1 Hz resolução
        u_pico = sd(picos) if n > 1 else np.nan
        # Sem pico: o máximo cai na borda da janela de procura ou mal sai do zero -> f0 não é medível
        f_zona = f[zona]
        sem_pico = picos.mean() < 3 or any(abs(x - f_zona[0]) < 2 or abs(x - f_zona[-1]) < 2 for x in f0s)
        linhas.append((d, f0s.mean(), u_f0, picos.mean(), u_pico, np.nanmean(Qs), sd(Qs), np.mean(foras), sem_pico))
        aviso = "  <-- sem pico de ressonância (fica fora dos ajustes)" if sem_pico else ""
        if abs(np.mean(foras)) > 2:
            aviso += "  <-- TL fora da banda longe de 0: referência não corresponde"
        if teto is not None:
            i = np.argmax(np.where(zona, TL_media, -np.inf))
            if TL_media[i] > teto[i] - 3:
                aviso += "  <-- perto do teto"
        print(f" {d * 100:6.1f}  {', '.join(f'{x:6.1f}' for x in f0s):27s}  {f0s.mean():6.1f} ± {u_f0:4.1f}  "
              f"{picos.mean():5.1f} ± {np.nan_to_num(u_pico):3.1f}   {np.nanmean(Qs):4.1f} ± {np.nan_to_num(sd(Qs)):3.1f}"
              f"   {np.mean(foras):+5.1f}{aviso}")
        sd_TL = np.std(TL_reps, axis=0, ddof=1) if len(TL_reps) > 1 else np.zeros_like(TL_media)
        ax[0, 0].fill_between(f[banda], (TL_media - sd_TL)[banda], (TL_media + sd_TL)[banda], color=cor,
                              alpha=0.12, linewidth=0)
        ax[0, 0].plot(f[banda], TL_media[banda], color=cor, label=f'd = {d * 100:g} cm')

    D = np.array([l[0] for l in linhas])
    F0, U_F0 = np.array([l[1] for l in linhas]), np.array([l[2] for l in linhas])
    TLM, U_TLM = np.array([l[3] for l in linhas]), np.array([l[4] for l in linhas])
    QQ, U_Q = np.array([l[5] for l in linhas]), np.array([l[6] for l in linhas])
    OK = ~np.array([l[8] for l in linhas])

    # ---------- f0 depende de d? ----------
    D_ok, F0_ok, U_ok = D[OK], F0[OK], U_F0[OK]
    if len(D_ok) >= 2:
        m, u_m, chi2 = media_ponderada(F0_ok, U_ok)
        dof = len(D_ok) - 1
        print(f"\nf0 médio (ponderado) = {m:.1f} ± {u_m:.1f} Hz   (previsto {f0_prev:.0f} ± {u_f0_prev:.0f} Hz)")
        print(f"   f0 constante? chi2/dof = {chi2:.1f}/{dof}" + ("  -> sim" if chi2 / dof < 2 else
              "  -> os f0 espalham-se mais do que as barras de erro"))
        if len(D_ok) >= 3:
            # reta f0 = a + k·d, pesada; k compatível com 0 -> f0 não depende da distância
            w = 1 / U_ok**2
            Sw, Sx, Sy = w.sum(), (w * D_ok).sum(), (w * F0_ok).sum()
            Sxx, Sxy = (w * D_ok * D_ok).sum(), (w * D_ok * F0_ok).sum()
            det = Sw * Sxx - Sx**2
            k = (Sw * Sxy - Sx * Sy) / det
            a0 = (Sxx * Sy - Sx * Sxy) / det
            u_k = np.sqrt(Sw / det)
            print(f"   declive df0/dd = {k / 100:+.2f} ± {u_k / 100:.2f} Hz/cm  "
                  f"(|k|/u = {abs(k) / u_k:.1f}; < 2 -> compatível com 0)")
        z = abs(m - f0_prev) / np.hypot(u_m, u_f0_prev)
        print(f"   medido vs previsto: z = {z:.1f}" + ("  (acordo)" if z < 2 else "  (desacordo)"))
    print(f"\nTL_max: de {TLM.min():.1f} a {TLM.max():.1f} dB. Se oscilar com d, comparar o período com "
          f"λ/2 = {meia_onda * 100:.1f} cm (ondas estacionárias no tubo).")

    # ---------- Gráficos ----------
    # O teto de 7/10 já não se desenha: é de outra sessão (ver teto_sistema).
    ax[0, 0].axvline(f0_prev, color='red', linestyle=':', label=f'$f_0$ previsto = {f0_prev:.0f} Hz')
    ax[0, 0].axhline(0, color='#bbb', linewidth=0.8)
    ax[0, 0].set_xlim(F_MIN, min(F_MAX, 2.5 * f0_prev))
    ax[0, 0].set_xlabel('Frequência (Hz)')
    ax[0, 0].set_ylabel('TL (dB)')
    ax[0, 0].set_title(f'TL(f) de uma seringa a {V * 1e6:g} mL, por distância (média de 3, faixa = desvio padrão)')
    ax[0, 0].grid(alpha=0.3)
    ax[0, 0].legend(fontsize=8)

    dd = D * 100
    eixo_d = 'Distância ao altifalante (cm)' if D1 is not None else 'Posição relativa à seringa 1 (cm)'
    ax[0, 1].errorbar(dd[OK], F0[OK], xerr=U_D * 100, yerr=U_F0[OK], fmt='o', color='black', capsize=4,
                      label='medido (média de 3 ± σ)')
    if (~OK).any():
        ax[0, 1].errorbar(dd[~OK], F0[~OK], xerr=U_D * 100, yerr=U_F0[~OK], fmt='o', mfc='white', color='gray',
                          capsize=4, label='sem pico (fora do ajuste)')
    ax[0, 1].axhspan(f0_prev - u_f0_prev, f0_prev + u_f0_prev, color='red', alpha=0.15, label=f'previsto {f0_prev:.0f} ± {u_f0_prev:.0f} Hz')
    if len(D_ok) >= 2:
        ax[0, 1].axhline(m, color='black', linestyle='--', linewidth=1,
                         label=f'média ponderada = {m:.0f} ± {u_m:.0f} Hz (χ²/ν = {chi2:.2f}/{dof})')
        if len(D_ok) >= 3:
            xx = np.array([D_ok.min(), D_ok.max()])
            ax[0, 1].plot(xx * 100, a0 + k * xx, color='tab:blue', linewidth=1.2,
                          label=f'reta: d$f_0$/dd = {k / 100:.1f} ± {u_k / 100:.1f} Hz/cm')
    ax[0, 1].set_xlabel(eixo_d)
    ax[0, 1].set_ylabel('$f_0$ (Hz)')
    ax[0, 1].set_title('Frequência de ressonância vs distância')
    ax[0, 1].grid(alpha=0.3)
    ax[0, 1].legend(fontsize=8)

    # Pontos sem linhas a ligá-los; as posições sem pico ficam a cinzento (o Q delas não existe)
    ax[1, 0].errorbar(dd[OK], TLM[OK], xerr=U_D * 100, yerr=np.nan_to_num(U_TLM[OK]), fmt='o',
                      color='tab:blue', capsize=4, label='medido (média de 3 ± σ)')
    if (~OK).any():
        ax[1, 0].errorbar(dd[~OK], TLM[~OK], xerr=U_D * 100, yerr=np.nan_to_num(U_TLM[~OK]), fmt='o',
                          mfc='white', color='gray', capsize=4, label='sem pico (máximo da janela)')
    ax[1, 0].axhline(0, color='#bbb', linewidth=0.8)
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].set_xlabel(eixo_d)
    ax[1, 0].set_ylabel('TL máximo (dB)')
    ax[1, 0].set_title(f'Atenuação máxima vs distância (λ/2 = {meia_onda * 100:.0f} cm)')
    ax[1, 0].grid(alpha=0.3)

    ax[1, 1].errorbar(dd[OK], QQ[OK], xerr=U_D * 100, yerr=np.nan_to_num(U_Q[OK]), fmt='o', color='tab:green',
                      capsize=4, label='medido (média de 3 ± σ); posições sem pico omitidas')
    ax[1, 1].set_xlim(ax[1, 0].get_xlim())
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].set_xlabel(eixo_d)
    ax[1, 1].set_ylabel('Q = $f_0$ / largura a −3 dB')
    ax[1, 1].set_title('Fator de qualidade do pico vs distância')
    ax[1, 1].grid(alpha=0.3)

    estilo_tl.virgulas(fig)
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, 'f0_vs_distancia.png'), dpi=150)
    print("\nFigura: f0_vs_distancia.png")


if __name__ == '__main__':
    main()
