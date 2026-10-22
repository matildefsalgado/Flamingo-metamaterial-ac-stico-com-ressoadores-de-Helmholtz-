import os
import numpy as np
from scipy.io import wavfile
from scipy.signal import welch

# Estilo comum dos gráficos de TL(f) da fusão (bragg_vs_local/tl_parteA.py) e da separação
# (bragg_vs_local/separacao.py): mesmas escalas, cores e anotações, para se poderem comparar.
# Também tem o modelo de matrizes de transferência e o teto, usados pelos dois scripts.

XLIM = (300, 1100)     # Hz
YLIM = (-15, 35)       # dB
AZUL, LARANJA, VERDE, CINZA_FB, CINZA_TETO = '#2a78d6', '#eb6834', '#2e9e5b', '#555', '#999'
PASTA_TETO = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'teto_sistema')


def virgulas(fig):
    """Vírgula decimal em toda a figura (eixos, títulos, legendas, anotações), como no texto do relatório.
    Chamar mesmo antes de guardar."""
    import re
    from matplotlib.ticker import FuncFormatter
    from matplotlib.text import Text
    pt = FuncFormatter(lambda v, _: f'{v:g}'.replace('.', ',').replace('-', '−'))
    for ax in fig.axes:
        for eixo in (ax.xaxis, ax.yaxis):
            if eixo.get_scale() == 'linear':
                eixo.set_major_formatter(pt)
    for t in fig.findobj(Text):
        s = t.get_text()
        novo = re.sub(r'(?<=\d)\.(?=\d)', ',', s)
        if novo != s:
            t.set_text(novo)


def tl_modelo(N, v, Q, f, c, S, L_ef, V_cone, a, S_tubo, rho=1.2):
    """TL de N ressoadores de Helmholtz espaçados de a num tubo sem perdas, terminação sem reflexão.
    A resistência de cada ressoador é a que dá o fator de qualidade Q."""
    w = 2 * np.pi * f
    k = w / c
    Zc = rho * c / S_tubo
    M = rho * L_ef / S                    # massa acústica do gargalo
    C = (v + V_cone) / (rho * c**2)       # compliância da cavidade
    R = np.sqrt(M / C) / Q
    Zr = R + 1j * w * M + 1 / (1j * w * C)
    out = np.empty(len(f))
    for i in range(len(f)):
        ramo = np.array([[1, 0], [1 / Zr[i], 1]])
        tubo = np.array([[np.cos(k[i] * a), 1j * Zc * np.sin(k[i] * a)],
                         [1j * np.sin(k[i] * a) / Zc, np.cos(k[i] * a)]])
        Tm = ramo
        for _ in range(N - 1):
            Tm = Tm @ tubo @ ramo
        out[i] = 20 * np.log10(abs(Tm[0, 0] + Tm[0, 1] / Zc + Zc * Tm[1, 0] + Tm[1, 1]) / 2)
    return out


def teto(F):
    """Teto de 7/10 (teto_sistema): 20·log(A_ref / max(A_silêncio, A_fuga)). Outra sessão: só indicativo."""
    def esp(nome):
        sr, x = wavfile.read(os.path.join(PASTA_TETO, nome))
        x = x.astype(float)
        if x.ndim > 1:
            x = x.mean(axis=1)
        f, P = welch(x, sr, nperseg=sr // 10)
        return np.interp(F, f, np.sqrt(P))
    return 20 * np.log10(esp('A_ref.wav') / np.maximum(esp('A_silencio.wav'), esp('A_fuga.wav')))


def grafico_tl(ax, F, TL, sd, *, titulo, rotulo, f0_prev, u_f0_prev, fB, u_fB,
               f_loc, u_loc, TL_loc, f_brg=None, u_brg=None, TL_brg=None,
               F_m=None, modelo=None, rotulo_modelo=None, controlo=None, mostrar_teto=False):
    """Desenha um TL(f) com o estilo comum. controlo = (TL, sd, rótulo) de 1 seringa, opcional."""
    ax.axvspan(f0_prev - u_f0_prev, f0_prev + u_f0_prev, color=VERDE, alpha=0.12, linewidth=0)
    ax.axvline(f0_prev, color=VERDE, linestyle='--', linewidth=1)
    ax.text(f0_prev - 6, YLIM[1] - 1, f'f0 previsto\n{f0_prev:.0f} ± {u_f0_prev:.0f} Hz',
            fontsize=9, color=VERDE, va='top', ha='right')
    ax.axvspan(fB - u_fB, fB + u_fB, color=CINZA_FB, alpha=0.15, linewidth=0)
    ax.axvline(fB, color=CINZA_FB, linestyle='--', linewidth=1)
    ax.text(fB + 6, YLIM[1] - 1, f'fB = c/2a\n{fB:.0f} ± {u_fB:.0f} Hz', fontsize=9, color='#333', va='top')
    if mostrar_teto:
        ax.plot(F, teto(F), color=CINZA_TETO, linestyle='--', linewidth=1, label='Teto (7/10, indicativo)')
    if modelo is not None:
        ax.plot(F_m, modelo, color=AZUL, linewidth=1.2, linestyle=':', label=rotulo_modelo)
    if controlo is not None:
        TLc, sdc, rot_c = controlo
        ax.fill_between(F, TLc - sdc, TLc + sdc, color=LARANJA, alpha=0.15, linewidth=0)
        ax.plot(F, TLc, color=LARANJA, linewidth=1.5, label=rot_c)
    ax.fill_between(F, TL - sd, TL + sd, color=AZUL, alpha=0.15, linewidth=0)
    ax.plot(F, TL, color=AZUL, linewidth=2, label=rotulo)
    ax.plot(f_loc, TL_loc, 'o', color=AZUL, markersize=6, markeredgecolor='black')
    ax.annotate(f'ressonância local\n{f_loc:.0f} ± {u_loc:.0f} Hz', (f_loc, TL_loc),
                xytext=(f_loc - 190, TL_loc + 3), fontsize=9, color='#333',
                arrowprops=dict(arrowstyle='-', color='#888'))
    if f_brg is not None:
        ax.plot(f_brg, TL_brg, 'o', color=AZUL, markersize=6, markeredgecolor='black')
        ax.annotate(f'Bragg\n{f_brg:.0f} ± {u_brg:.0f} Hz', (f_brg, TL_brg),
                    xytext=(f_brg + 90, TL_brg + 9), fontsize=9, color='#333',
                    arrowprops=dict(arrowstyle='-', color='#888'))
    ax.axhline(0, color='#bbb', linewidth=0.8)
    ax.set_xlabel('Frequência (Hz)')
    ax.set_ylabel('TL = 20 log10(A_ref / A)  (dB)')
    ax.set_title(titulo)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.grid(True, alpha=0.25)
    for lado in ('top', 'right'):
        ax.spines[lado].set_visible(False)
    ax.legend(loc='lower right', frameon=False, fontsize=9)
