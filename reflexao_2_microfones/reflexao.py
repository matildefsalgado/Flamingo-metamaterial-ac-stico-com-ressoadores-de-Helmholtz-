"""
Reflexão de UMA seringa e onda estacionária no tubo, com UM só microfone (mono).

O microfone ocupa duas posições, em gravações separadas:
    posição 1 = montante (antes da 1.ª seringa): ouve a onda incidente + a refletida
    posição 2 = jusante (depois da última seringa): ouve a onda transmitida
Duas terminações do tubo: lã (anecoica) e tampão de LEGO (refletora).

Ficheiros (nesta pasta, 48 kHz, mono; tudo sem mexer no gerador nem no volume do PC):
    la_pos2_ref.wav        lã, mic na posição 2, êmbolos todos no fundo        (varrimento)
    la_pos2_seringa.wav    lã, mic na posição 2, seringa 1 a V                  (varrimento)
    la_pos1_ref.wav        lã, mic na posição 1, êmbolos no fundo               (varrimento)
    la_pos1_seringa.wav    lã, mic na posição 1, seringa 1 a V                  (varrimento)
    lego_pos2_ref.wav, lego_pos2_seringa.wav, lego_pos1_ref.wav, lego_pos1_seringa.wav
                           o mesmo, com o tampão de LEGO no lugar da lã
    tom_la.wav, tom_lego.wav   seno contínuo em f0, mic na posição 1, ~5 s     [opcional]
Repetições: la_pos2_seringa_1.wav, la_pos2_seringa_2.wav, ... (faz-se a média; as que faltam
são ignoradas). Se faltar um par, essa parte da análise é saltada.

O que calcula:
  1. Lã, posição 2:  TL(f) = 20 log10(A_ref / A_seringa);  f0 = máximo de TL;
                     |T| = A_seringa / A_ref em f0;  |R| = sqrt(1 - |T|²)  (sem perdas).
  2. Modelo:  Z_c = rho c / S_tubo,
              Z_b = R_a + j w rho L_ef / (pi R_B R_P) + rho c² / (j w (V + V_cone)),
              R_teórico = -Z_c / (2 Z_b + Z_c),   T_teórico = 2 Z_b / (2 Z_b + Z_c).
     Sem perdas (R_a = 0) o modelo dá |R| = 1 em f0. O atrito R_a ajusta-se às ABAS do pico de TL
     (sem usar o topo) e depois prevê-se a altura do pico -> comparação independente com o medido.
  3. LEGO, posição 2: o mesmo TL, para ver o efeito da terminação.
  4. Posição 1: TL = 20 log(A_ref / A_seringa) (lã e LEGO) -> negativo onde a reflexão da seringa se soma.
  5. Onda estacionária (varrimentos de referência na posição 1):
     TL do LEGO em relação à lã, 20 log(A_lã / A_LEGO), vs f  e vs  d/lambda = d f / c.  Espaçamento entre máximos Δf (autocorrelação)
     comparado com c / (2 d).  |R_tampão| a partir da ondulação (pico a vale).
  6. Tons: TL do LEGO em relação à lã à frequência do seno (um ponto da curva anterior, como verificação).
Fase/desfasamento: NÃO se calcula (um só microfone, gerador sem sincronização com o PC).

Uso:  python reflexao.py
      python reflexao.py --pasta OUTRA_PASTA     (para testar com outros ficheiros)
"""
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import welch
from scipy.ndimage import uniform_filter1d
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import estilo_tl  # noqa: E402  (vírgula decimal nas figuras)

# ---------- Parâmetros (SI) ----------
T = 22             # temperatura da sala (°C), a do registo (igual em todos os scripts)
U_T = 0.5          # incerteza da temperatura (°C)
V = 4.0e-6         # volume da seringa 1 (m³), medido: 4,0 mL
U_V = 0.5e-6       # incerteza de leitura do volume (m³): metade da menor divisão (1 mL)
U_R = 0.1e-3       # incerteza dos raios do cone (m): craveira + parede
U_LEF = 2.5e-3     # incerteza de L_ef (m), da calibração (correcao_gargalo --sem 10)
D_TUBO = 0.025     # diâmetro INTERIOR do tubo (m), medido: 2,5 cm
D_MIC_PAREDE = 1.555  # (medido com a fita: 155,5 cm) distância do furo da posição 1 à face interior do tampão de LEGO (m); None = estimar
F_MIN, F_MAX = 150.0, 2000.0
NPERSEG = 2**15                         # resolução ~1,5 Hz
LARGURAS_HZ = [20, 30, 40, 50, 60, 70]  # suavização em frequência (o tubo dá ondulação)
JANELA = 0.35                           # procura o pico em f0_previsto ± 35 %

# Gargalo = bico cónico da seringa. Valores acertados aos dados de correcao_gargalo.
L = 0.035              # comprimento do cone, no eixo (m)
D_BASE = 0.0132        # diâmetro exterior na base (m)
PAREDE = 0.0012        # grossura da parede na base (m)
D_PONTA = 0.0041       # diâmetro interior na ponta (m)
R_B, R_P = D_BASE / 2 - PAREDE, D_PONTA / 2
S = np.pi * R_B * R_P  # área do cilindro equivalente ao cone
dL = 0.73 * (R_B + R_P)   # correção de extremidade (m)
L_EF = L + dL
V_CONE = np.pi * L / 3 * (R_B**2 + R_B * R_P + R_P**2)   # ar no cone (m³): soma-se ao êmbolo

sys.stdout.reconfigure(encoding='utf-8')
c = 331.3 * np.sqrt(1 + T / 273.15)
RHO = 101325 / (287.05 * (T + 273.15))     # densidade do ar (kg/m³)
S_TUBO = np.pi * (D_TUBO / 2)**2
Z_C = RHO * c / S_TUBO
f0_prev = c / (2 * np.pi) * np.sqrt(S / ((V + V_CONE) * L_EF))
# f0 ∝ c·√(S/(V_total·L_ef)): em quadratura, metade da incerteza relativa de S, V_total e L_ef
u_f0_prev = f0_prev * np.sqrt((0.5 * U_V / (V + V_CONE))**2 + (0.5 * np.hypot(U_R / R_B, U_R / R_P))**2
                              + (0.5 * U_LEF / L_EF)**2 + (0.5 * U_T / (T + 273.15))**2)


def espectro(caminho):
    sr, x = wavfile.read(caminho)
    x = x.astype(float)
    if x.ndim > 1:
        x = x.mean(axis=1)
    f, P = welch(x, sr, nperseg=NPERSEG)
    return f, np.sqrt(P)


# Outras gravações que servem para a mesma coisa. ref16_10.wav (16/10, 44,1 kHz) é a referência da
# separação de bandas: êmbolos todos no fundo, lã, mic na posição 2, mesmo dia que posicao2_la_4ml.
ALIASES = {'la_pos2_ref': ['ref16_10'], 'la_pos2_seringa': ['pos2la4ml_FINAL']}


def nome_alternativo(nome):
    """la_pos2_ref -> posicao2_la_0ml ; lego_pos1_seringa -> posicao1_lego_4ml (nomes usados na aula)."""
    term, pos, tipo = nome.split('_')
    return f"posicao{pos[-1]}_{term}_{'0ml' if tipo == 'ref' else '4*ml'}"


def carregar(pasta, nome, f_grelha=None):
    """Espectros de nome.wav e nome_<n>.wav (ou do nome alternativo). Devolve (f, amplitudes) ou (None, [])."""
    ps = set()
    for n in (nome, nome_alternativo(nome), *ALIASES.get(nome, [])):
        ps |= set(glob.glob(os.path.join(pasta, n + '.wav'))) | set(glob.glob(os.path.join(pasta, n + '_*.wav')))
    ps = sorted(ps)
    if not ps:
        return None, []
    f = f_grelha
    lista = []
    for p in ps:
        fp, A = espectro(p)
        if f is None:
            f = fp
        lista.append(np.interp(f, fp, A))
    return f, lista


def suavizar(f, y, largura):
    return uniform_filter1d(y, max(1, int(round(largura / (f[1] - f[0])))))


def T_modelo(f, R_a):
    w = 2 * np.pi * f
    Z_b = R_a + 1j * w * RHO * L_EF / S + RHO * c**2 / (1j * w * (V + V_CONE))
    return 2 * Z_b / (2 * Z_b + Z_C)


def R_modelo(f, R_a):
    w = 2 * np.pi * f
    Z_b = R_a + 1j * w * RHO * L_EF / S + RHO * c**2 / (1j * w * (V + V_CONE))
    return -Z_C / (2 * Z_b + Z_C)


def analisar_TL(f, refs, sers, zona):
    """TL médio, f0 (várias suavizações), TL no pico e incertezas."""
    A_ref = np.mean(refs, axis=0)
    reps = np.array([20 * np.log10(A_ref / A) for A in sers])
    TL_bruto = reps.mean(axis=0)
    est = [f[np.argmax(np.where(zona, suavizar(f, TL_bruto, W), -np.inf))] for W in LARGURAS_HZ]
    f0 = np.mean(est)
    TL = suavizar(f, TL_bruto, 30)
    i = np.argmax(np.where(zona, TL, -np.inf))
    if len(reps) > 1:
        picos = [suavizar(f, r, 30)[i] for r in reps]
        u_TL = np.std(picos, ddof=1) / np.sqrt(len(picos))
        f0s = [f[np.argmax(np.where(zona, suavizar(f, r, 30), -np.inf))] for r in reps]
        u_rep = np.std(f0s, ddof=1) / np.sqrt(len(f0s))
    else:
        u_TL, u_rep = np.nan, 0.0
    u_f0 = np.sqrt(np.std(est)**2 + u_rep**2 + 1.0)
    return TL, f0, u_f0, i, u_TL, len(reps)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pasta', default=os.path.dirname(os.path.abspath(__file__)))
    pasta = ap.parse_args().pasta

    print(f"c = {c:.1f} m/s | rho = {RHO:.3f} kg/m³ | V = {V * 1e6:.1f} mL | "
          f"f0 previsto = {f0_prev:.0f} ± {u_f0_prev:.0f} Hz | Z_c = {Z_C:.3g} Pa·s/m³ (D_tubo = {D_TUBO * 1e3:.0f} mm)")

    # grelha comum: a do primeiro ficheiro que existir
    f = None
    for nome in ('la_pos2_ref', 'la_pos1_ref', 'lego_pos1_ref', 'lego_pos2_ref'):
        f, _ = carregar(pasta, nome)
        if f is not None:
            break
    if f is None:
        raise SystemExit("Não encontrei gravações (la_pos2_ref.wav, la_pos2_seringa.wav, ...).")
    banda = (f >= F_MIN) & (f <= F_MAX)
    zona = banda & (f >= f0_prev * (1 - JANELA)) & (f <= f0_prev * (1 + JANELA))
    G = {k: carregar(pasta, k, f)[1] for k in
         ('la_pos2_ref', 'la_pos2_seringa', 'la_pos1_ref', 'la_pos1_seringa',
          'lego_pos2_ref', 'lego_pos2_seringa', 'lego_pos1_ref', 'lego_pos1_seringa')}
    for k, v in G.items():
        print(f"   {k:18s} {len(v)} gravação(ões)" if v else f"   {k:18s} (falta)")

    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    resumo = []

    # ---------- 1. Lã, posição 2: TL, f0, |T|, |R| ----------
    res_la = None
    if G['la_pos2_ref'] and G['la_pos2_seringa']:
        TL, f0, u_f0, i, u_TL, n = analisar_TL(f, G['la_pos2_ref'], G['la_pos2_seringa'], zona)
        Tm = 10**(-TL[i] / 20)
        Rm = np.sqrt(max(0.0, 1 - Tm**2))
        u_Tm = Tm * np.log(10) / 20 * u_TL if np.isfinite(u_TL) else np.nan
        u_Rm = Tm * u_Tm / Rm if Rm > 0 and np.isfinite(u_Tm) else np.nan
        res_la = (TL, f0, i)
        print(f"\n1. Lã, posição 2 ({n} rep.): f0 = {f0:.1f} ± {u_f0:.1f} Hz  (previsto {f0_prev:.0f} ± "
              f"{u_f0_prev:.0f} Hz, {100 * (f0 / f0_prev - 1):+.0f} %)")
        print(f"   TL(f0) = {TL[i]:.1f}" + (f" ± {u_TL:.1f}" if np.isfinite(u_TL) else "") + " dB"
              f"  ->  |T| = {Tm:.3f}" + (f" ± {u_Tm:.3f}" if np.isfinite(u_Tm) else "")
              + f",  |R| = sqrt(1 - |T|²) = {Rm:.3f}" + (f" ± {u_Rm:.3f}" if np.isfinite(u_Rm) else ""))
        resumo += [f"f0 (lã) = {f0:.0f} ± {u_f0:.0f} Hz  (prev. {f0_prev:.0f})",
                   f"TL(f0) = {TL[i]:.1f} dB,  |T| = {Tm:.3f},  |R| = {Rm:.3f}"]

        # ---------- 2. Modelo de impedâncias ----------
        # O modelo é deslocado para o f0 medido (escala de frequência) e R_a ajusta-se às abas
        # do pico (3 dB <= TL <= TL_max - 3 dB); a altura do pico fica como previsão.
        escala = f0_prev / f0
        abas = zona & (TL >= 3) & (TL <= TL[i] - 3)
        if abas.sum() >= 4:
            def custo(logR):
                TLm = -20 * np.log10(np.abs(T_modelo(f[abas] * escala, 10**logR)))
                return np.sum((TLm - TL[abas])**2)
            logR = minimize_scalar(custo, bounds=(3, 9), method='bounded').x
            R_a = 10**logR
            R_prev_pico = np.abs(R_modelo(f0_prev, R_a))
            TL_prev_pico = -20 * np.log10(np.abs(T_modelo(f0_prev, R_a)))
            print(f"\n2. Modelo: sem perdas |R_teórico(f0)| = 1 (reflexão total).")
            print(f"   Atrito ajustado às abas do pico: R_a = {R_a:.3g} Pa·s/m³ (R_a/Z_c = {R_a / Z_C:.2f})")
            print(f"   Previsão no pico: |R_teórico| = Z_c/(2R_a + Z_c) = {R_prev_pico:.3f}, "
                  f"TL = {TL_prev_pico:.1f} dB   vs medido |R| = {Rm:.3f}, TL = {TL[i]:.1f} dB")
            # Balanço de energia: com atrito, parte da energia fica no gargalo (|R|² + |T|² < 1),
            # e sqrt(1 - |T|²) passa a sobrestimar |R|.
            T_pico = np.abs(T_modelo(f0_prev, R_a))
            absorvido = 1 - R_prev_pico**2 - T_pico**2
            print(f"   Balanço no pico (modelo): |R|² = {R_prev_pico**2:.2f}, |T|² = {T_pico**2:.2f}, "
                  f"absorvido = {absorvido:.2f}")
            if absorvido > 0.1:
                print(f"   -> {100 * absorvido:.0f} % da energia é absorvida: sqrt(1 - |T|²) = {Rm:.3f} "
                      f"sobrestima |R|; o valor do modelo ({R_prev_pico:.3f}) é mais realista.")
            resumo += [f"modelo: R_a = {R_a:.2g} Pa·s/m³ (R_a/Z_c = {R_a / Z_C:.2f})",
                       f"|R| modelo = {R_prev_pico:.3f}  vs  sqrt(1-|T|²) = {Rm:.3f}",
                       f"absorvido no gargalo (modelo) = {100 * absorvido:.0f} %"]
            perto = zona
            Rmed = np.sqrt(np.clip(1 - 10**(-TL / 10), 0, 1))
            ax[0, 1].plot(f[perto], Rmed[perto], color='black', label='medido: sqrt(1 - |T|²)')
            ax[0, 1].plot(f[perto], np.abs(R_modelo(f[perto] * escala, 0.0)), color='tab:red', linestyle='--',
                          label='teórico sem perdas')
            ax[0, 1].plot(f[perto], np.abs(R_modelo(f[perto] * escala, R_a)), color='tab:red',
                          label=f'teórico, R_a ajustado ({R_a / Z_C:.2f} Z_c)')
        else:
            print("\n2. Modelo: o pico tem poucos pontos nas abas para ajustar o atrito.")
        ax[0, 1].axvline(f0, color='gray', linestyle=':')
        ax[0, 1].set_xlabel('Frequência (Hz)')
        ax[0, 1].set_ylabel('|R|')
        ax[0, 1].set_ylim(0, 1.05)
        ax[0, 1].set_title('Coeficiente de reflexão da seringa (lã, posição 2)\n'
                           'modelo deslocado para o f0 medido')
        ax[0, 1].grid(alpha=0.3)
        ax[0, 1].legend(fontsize=8)
        ax[0, 0].plot(f[banda], TL[banda], color='tab:blue', label='lã')
        ax[0, 0].axvline(f0, color='tab:blue', linestyle=':')
    else:
        print("\n1-2. Faltam la_pos2_ref / la_pos2_seringa: sem |T| nem |R|.")

    # ---------- 3. LEGO, posição 2 ----------
    if G['lego_pos2_ref'] and G['lego_pos2_seringa']:
        TLg, f0g, u_f0g, ig, u_TLg, ng = analisar_TL(f, G['lego_pos2_ref'], G['lego_pos2_seringa'], zona)
        print(f"\n3. LEGO, posição 2 ({ng} rep.): f0 = {f0g:.1f} ± {u_f0g:.1f} Hz, TL(f0) = {TLg[ig]:.1f} dB"
              + (f"   (lã: {res_la[1]:.1f} Hz, {res_la[0][res_la[2]]:.1f} dB)" if res_la else ""))
        resumo.append(f"LEGO, pos. 2: f0 = {f0g:.0f} Hz, TL = {TLg[ig]:.1f} dB")
        ax[0, 0].plot(f[banda], TLg[banda], color='tab:orange', label='LEGO')
        ax[0, 0].axvline(f0g, color='tab:orange', linestyle=':')
    ax[0, 0].axvline(f0_prev, color='red', linestyle='--', linewidth=1, label=f'f0 previsto = {f0_prev:.0f} Hz')
    ax[0, 0].axhline(0, color='#bbb', linewidth=0.8)
    ax[0, 0].set_xlim(F_MIN, min(F_MAX, 2.5 * f0_prev))
    ax[0, 0].set_xlabel('Frequência (Hz)')
    ax[0, 0].set_ylabel('TL = 20 log(A_ref / A_seringa)  (dB)')
    ax[0, 0].set_title(f'Perda de transmissão de 1 seringa a {V * 1e6:g} mL (posição 2)')
    ax[0, 0].grid(alpha=0.3)
    ax[0, 0].legend(fontsize=8)

    # ---------- 4. Posição 1: reflexão vista a montante ----------
    for term, cor in (('la', 'tab:blue'), ('lego', 'tab:orange')):
        if G[f'{term}_pos1_ref'] and G[f'{term}_pos1_seringa']:
            # TL na posição 1 (antes da seringa): negativo onde o som refletido pela seringa se soma
            g = suavizar(f, 20 * np.log10(np.mean(G[f'{term}_pos1_ref'], axis=0)
                                          / np.mean(G[f'{term}_pos1_seringa'], axis=0)), 30)
            ax[0, 2].plot(f[banda], g[banda], color=cor, label='lã' if term == 'la' else 'LEGO')
            j = np.argmax(np.where(zona, np.abs(g), -np.inf))
            print(f"\n4. Posição 1, {'lã' if term == 'la' else 'LEGO'}: maior desvio do TL perto de "
                  f"f0 = {g[j]:+.1f} dB a {f[j]:.0f} Hz  (a montante soma-se a onda refletida pela seringa)")
    ax[0, 2].axvline(f0_prev, color='red', linestyle='--', linewidth=1)
    ax[0, 2].axhline(0, color='#bbb', linewidth=0.8)
    ax[0, 2].set_xlim(F_MIN, min(F_MAX, 2.5 * f0_prev))
    ax[0, 2].set_xlabel('Frequência (Hz)')
    ax[0, 2].set_ylabel('TL = 20 log(A_ref / A_seringa)  (dB)')
    ax[0, 2].set_title('Perda de transmissão na posição 1 (antes da seringa)')
    ax[0, 2].grid(alpha=0.3)
    if ax[0, 2].lines:
        ax[0, 2].legend(fontsize=8)

    # ---------- 5. Onda estacionária: A_LEGO / A_lã nos varrimentos de referência, posição 1 ----------
    if G['la_pos1_ref'] and G['lego_pos1_ref']:
        # TL do LEGO com a lã como referência: 20 log(A_lã / A_LEGO)
        razao = 20 * np.log10(np.mean(G['la_pos1_ref'], axis=0) / np.mean(G['lego_pos1_ref'], axis=0))
        y = suavizar(f, razao, 10)[banda]
        fb = f[banda]
        # Δf: primeiro máximo da autocorrelação (sem a tendência lenta)
        yd = y - suavizar(fb, y, 400)
        ac = np.correlate(yd, yd, mode='full')[len(yd) - 1:]
        ac /= ac[0]
        df = fb[1] - fb[0]
        k0 = int(30 / df)
        # primeiro máximo local claro (> 0,3); os seguintes são múltiplos de Δf
        picos_ac = [j for j in range(k0, len(ac) - 1) if ac[j] >= ac[j - 1] and ac[j] >= ac[j + 1]]
        claros = [j for j in picos_ac if ac[j] > 0.3]
        k = claros[0] if claros else max(picos_ac, key=lambda j: ac[j])
        delta_f = k * df
        d_est = c / (2 * delta_f)
        # ondulação pico a vale -> |R_tampão|: máximo - mínimo em janelas de largura Δf, mediana
        n_j = max(2, int(round(delta_f / df)))
        ond = np.median([y[j:j + n_j].max() - y[j:j + n_j].min() for j in range(0, len(y) - n_j, n_j // 2)])
        R_tampao = (10**(ond / 20) - 1) / (10**(ond / 20) + 1)
        print(f"\n5. Onda estacionária (TL do LEGO em relação à lã, posição 1, êmbolos no fundo):")
        print(f"   espaçamento entre máximos Δf = {delta_f:.0f} Hz  ->  d = c/(2Δf) = {d_est * 100:.1f} cm")
        if D_MIC_PAREDE:
            print(f"   d medido com a fita = {D_MIC_PAREDE * 100:.1f} cm  ->  Δf teórico = c/(2d) = "
                  f"{c / (2 * D_MIC_PAREDE):.0f} Hz  ({100 * (delta_f / (c / (2 * D_MIC_PAREDE)) - 1):+.0f} %)")
        else:
            print("   (D_MIC_PAREDE = None: medir com a fita e pôr no topo do código para comparar)")
        print(f"   ondulação pico a vale = {ond:.1f} dB  ->  |R_tampão| = {R_tampao:.2f}")
        resumo += [f"onda estacionária: Δf = {delta_f:.0f} Hz -> d = {d_est * 100:.0f} cm"
                   + (f" (fita: {D_MIC_PAREDE * 100:.0f} cm)" if D_MIC_PAREDE else ""),
                   f"|R_tampão| (ondulação {ond:.1f} dB) = {R_tampao:.2f}"]
        d_uso = D_MIC_PAREDE or d_est
        ax[1, 0].plot(fb, y, color='tab:purple')
        ax[1, 0].set_xlabel('Frequência (Hz)')
        ax[1, 0].set_title(f'Onda estacionária: TL do LEGO em relação à lã (posição 1)\nΔf = {delta_f:.0f} Hz')
        ax[1, 1].plot(d_uso * fb / c, y, color='tab:purple')
        for n in np.arange(0, d_uso * fb[-1] / c + 0.5, 0.5):
            ax[1, 1].axvline(n, color='gray', linestyle=':', linewidth=0.8)
        ax[1, 1].set_xlabel(f'd/λ = d·f/c   (d = {d_uso * 100:.1f} cm, '
                            + ('fita' if D_MIC_PAREDE else 'estimado de Δf') + ')')
        ax[1, 1].set_title('A mesma curva contra d/λ (linhas em 0, ½, 1, 3/2, ...)')
        for a_ in (ax[1, 0], ax[1, 1]):
            a_.set_ylabel('TL = 20 log(A_lã / A_LEGO)  (dB)')
            a_.axhline(0, color='#bbb', linewidth=0.8)
            a_.grid(alpha=0.3)

        # ---------- 6. Tons ----------
        # tom_la.wav / tom_la_<n>.wav, ou os nomes usados na aula (655_1.wav, 655_1_lego.wav, ...)
        p_la = sorted(set(glob.glob(os.path.join(pasta, 'tom_la*.wav')))
                      | set(glob.glob(os.path.join(pasta, '[0-9][0-9][0-9]_[0-9].wav'))))
        p_lego = sorted(set(glob.glob(os.path.join(pasta, 'tom_lego*.wav')))
                        | set(glob.glob(os.path.join(pasta, '[0-9][0-9][0-9]_[0-9]_lego.wav'))))
        if p_la and p_lego:
            ft, _ = espectro(p_la[0])
            Als = [np.interp(ft, *espectro(p)) for p in p_la]
            Ags = [np.interp(ft, *espectro(p)) for p in p_lego]
            j = np.argmax(np.where(ft > F_MIN, np.mean(Als, axis=0), 0))
            dB_la = np.array([20 * np.log10(A[j]) for A in Als])
            dB_lego = np.array([20 * np.log10(A[j]) for A in Ags])
            r_tom = dB_la.mean() - dB_lego.mean()   # TL do LEGO em relação à lã
            # incerteza: desvio-padrão da média de cada conjunto, em quadratura
            u = lambda v: np.std(v, ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
            u_tom = np.hypot(u(dB_la), u(dB_lego))
            r_var = np.interp(ft[j], fb, y)
            print(f"\n6. Tom a {ft[j]:.0f} Hz ({len(p_la)} lã, {len(p_lego)} LEGO): "
                  f"TL LEGO/lã = {r_tom:+.1f} ± {u_tom:.1f} dB  (varrimento a essa frequência: {r_var:+.1f} dB)")
            print(f"   lã:   " + ', '.join(f"{v:.1f}" for v in dB_la) + " dB")
            print(f"   LEGO: " + ', '.join(f"{v:.1f}" for v in dB_lego) + " dB")
            resumo.append(f"tom {ft[j]:.0f} Hz: {r_tom:+.1f} ± {u_tom:.1f} dB (varrimento {r_var:+.1f} dB)")
            ax[1, 0].errorbar(ft[j], r_tom, yerr=u_tom, fmt='o', color='black', capsize=3,
                              label='tom (seno contínuo)')
            ax[1, 1].errorbar(d_uso * ft[j] / c, r_tom, yerr=u_tom, fmt='o', color='black', capsize=3,
                              label='tom (seno contínuo)')
            ax[1, 0].legend(fontsize=8)
            ax[1, 1].legend(fontsize=8)
    else:
        print("\n5. Faltam la_pos1_ref / lego_pos1_ref: sem onda estacionária.")
        ax[1, 0].axis('off')
        ax[1, 1].axis('off')

    ax[1, 2].axis('off')
    ax[1, 2].text(0, 1, 'Resumo\n\n' + '\n'.join(resumo) +
                  '\n\nFase não medida: 1 microfone,\ngerador sem sincronização com o PC.',
                  va='top', family='monospace', fontsize=9)
    # o gráfico da posição 1 sai da figura (só confirma a reflexão; os números ficam no terminal)
    ax[0, 2].set_visible(False)
    estilo_tl.virgulas(fig)
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, 'reflexao.png'), dpi=150)

    # cada gráfico também num ficheiro à parte, para se poder referir a eles separadamente
    separados = {'reflexao_a_TL_posicao2.png': ax[0, 0], 'reflexao_b_coef_reflexao.png': ax[0, 1],
                 'reflexao_c_onda_estacionaria_f.png': ax[1, 0], 'reflexao_d_onda_estacionaria_d_lambda.png': ax[1, 1]}
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for nome, a_ in separados.items():
        if a_.axison:
            caixa = a_.get_tightbbox(rend).transformed(fig.dpi_scale_trans.inverted()).padded(0.1)
            fig.savefig(os.path.join(pasta, nome), dpi=150, bbox_inches=caixa)
    print("\nFiguras: reflexao.png e, em separado, " + ', '.join(separados))


if __name__ == '__main__':
    main()
