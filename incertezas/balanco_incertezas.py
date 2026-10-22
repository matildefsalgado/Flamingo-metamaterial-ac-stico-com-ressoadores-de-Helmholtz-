import os
import re
import sys
import numpy as np

# Balanço de incertezas para o Apêndice C do relatório.
# Para cada grandeza principal: a contribuição de cada fonte (em Hz ou mm e em %) e a total,
# somadas em quadratura (fontes independentes). Escreve uma tabela no terminal e o ficheiro
# tabela_incertezas.tex (booktabs) para o Overleaf.
#
# As previsões são calculadas aqui com a mesma geometria dos outros scripts.
# Os valores MEDIDOS vêm dos outros scripts (ver a fonte ao lado de cada um);
# se voltarem a correr esses scripts e os números mudarem, atualizem-nos aqui.

sys.stdout.reconfigure(encoding='utf-8')
PASTA = os.path.dirname(os.path.abspath(__file__))

# ---------- Valores de entrada e incertezas (iguais em todos os scripts) ----------
T, U_T = 22, 0.5               # temperatura (°C), a do registo
a, U_A = 0.25, 0.001           # espaçamento entre seringas (m), fita
U_V = 0.5e-6                   # leitura do volume (m³): metade da menor divisão (1 mL)
L = 0.035                      # comprimento do cone (m)
D_BASE, PAREDE, D_PONTA = 0.0132, 0.0012, 0.0041
U_R = 0.1e-3                   # raios do cone (m): craveira + parede
U_LEF = 2.5e-3                 # L_ef (m), da calibração (correcao_gargalo --sem 10)
RES = 10.0                     # resolução em frequência das análises de TL (Hz)

# ---------- Valores medidos (copiados dos outros scripts) ----------
F0_REPS_FUSAO = [630, 620, 620]   # bragg_vs_local/tl_parteA.py: picos das 3 repetições (Hz)
F0_FUSAO = 620.0                  # idem, pico do TL médio (Hz)
M_CAL, U_M_CAL = 1.6025, 0.0239   # correcao_gargalo.py --sem 10: declive (Hz·m^1,5), já × √(χ²/ν)
CHI2_CAL = (5.0, 3)               # idem, χ² e ν

# ---------- Geometria ----------
c = 331.3 * np.sqrt(1 + T / 273.15)
R_B, R_P = D_BASE / 2 - PAREDE, D_PONTA / 2
S = np.pi * R_B * R_P
dL = 0.73 * (R_B + R_P)
L_ef = L + dL
V_CONE = np.pi * L / 3 * (R_B**2 + R_B * R_P + R_P**2)
r_c = 0.5 * U_T / (T + 273.15)                 # δc/c
r_S = np.hypot(U_R / R_B, U_R / R_P)           # δS/S


def f0_de(v):
    return c / (2 * np.pi) * np.sqrt(S / ((v + V_CONE) * L_ef))


def bloco(nome, valor, unidade, fontes):
    """fontes: lista de (descrição, incerteza relativa). Devolve (linhas, total relativo)."""
    total = np.sqrt(sum(r**2 for _, r in fontes))
    linhas = [(desc, valor * r, 100 * r) for desc, r in fontes]
    return {'nome': nome, 'valor': valor, 'unidade': unidade, 'linhas': linhas,
            'u': valor * total, 'rel': 100 * total}


blocos = []

# f0 previsto: f0 ∝ c·√(S/(V_total·L_ef))  ->  expoente ½ em S, V_total e L_ef
for v in (4.1e-6, 14.6e-6):
    vt = v + V_CONE
    blocos.append(bloco(f'f0 previsto, V = {v * 1e6:.1f} mL', f0_de(v), 'Hz', [
        (f'temperatura (±{U_T} °C)', r_c),
        (f'volume (±{U_V * 1e6:.1f} mL em {vt * 1e6:.1f} mL)', 0.5 * U_V / vt),
        (f'área S (raios ±{U_R * 1e3:.1f} mm)', 0.5 * r_S),
        (f'L_ef (±{U_LEF * 1e3:.1f} mm, calibração)', 0.5 * U_LEF / L_ef),
    ]))

# fB = c/2a
fB = c / (2 * a)
blocos.append(bloco('fB previsto = c/2a', fB, 'Hz', [
    (f'temperatura (±{U_T} °C)', r_c),
    (f'espaçamento a (±{U_A * 1e3:.0f} mm)', U_A / a),
]))

# f0 medido na fusão: repetições + resolução
reps = np.array(F0_REPS_FUSAO, dtype=float)
u_rep = reps.std(ddof=1) / np.sqrt(len(reps))
b_f0 = bloco('f0 medido, fusão (4,1 mL)', F0_FUSAO, 'Hz', [
    (f'repetições (s/√{len(reps)})', u_rep / F0_FUSAO),
    (f'resolução (metade de {RES:.0f} Hz)', RES / 2 / F0_FUSAO),
])
blocos.append(b_f0)

# Dessintonia δ = 1 - f0/fB: δ_abs = (f0/fB)·√((δf0/f0)² + (δfB/fB)²)
delta = 1 - F0_FUSAO / fB
q = F0_FUSAO / fB
b_fB = blocos[2]
fontes_d = [(f'f0 medido (±{b_f0["u"]:.0f} Hz)', b_f0['u'] / F0_FUSAO),
            (f'fB previsto (±{b_fB["u"]:.0f} Hz)', b_fB['u'] / fB)]
total_d = q * np.sqrt(sum(r**2 for _, r in fontes_d))
blocos.append({'nome': 'dessintonia δ = 1 − f0/fB', 'valor': 100 * delta, 'unidade': '%',
               'linhas': [(d, 100 * q * r, None) for d, r in fontes_d],
               'u': 100 * total_d, 'rel': None})

# L_ef calibrado: L_ef = c² S / (4π² m²)
L_cal = c**2 * S / (4 * np.pi**2 * M_CAL**2)
blocos.append(bloco('L_ef calibrado', L_cal * 1e3, 'mm', [
    (f'declive m (ajuste, χ²/ν = {CHI2_CAL[0]:.1f}/{CHI2_CAL[1]})', 2 * U_M_CAL / M_CAL),
    (f'área S (raios ±{U_R * 1e3:.1f} mm)', r_S),
    (f'temperatura (±{U_T} °C)', 2 * r_c),
]))


# ---------- Terminal ----------
def fmt(x, casas=1):
    return f'{x:.{casas}f}'


for b in blocos:
    print(f"\n{b['nome']} = {fmt(b['valor'])} {b['unidade']}")
    for desc, u, r in b['linhas']:
        rel = f'  ({r:.2f} %)' if r is not None else ''
        print(f"   {desc:45s} ± {fmt(u, 2):>7s} {b['unidade']}{rel}")
    rel = f"  ({b['rel']:.1f} %)" if b['rel'] is not None else ''
    dom = max(b['linhas'], key=lambda l: abs(l[1]))[0]
    print(f"   {'TOTAL (quadratura)':45s} ± {fmt(b['u'], 2):>7s} {b['unidade']}{rel}   domina: {dom}")
print(f"\nNota: o L do cone ({L * 1e3:.0f} mm) entra sem incerteza (valor fixo do modelo).")


# ---------- LaTeX ----------
def tex(s):
    s = re.sub('√' + r'(\d+)', lambda m: r'$\sqrt{' + m.group(1) + '}$', s)
    return (s.replace('%', r'\%').replace('±', r'$\pm$').replace('δ', r'$\delta$')
            .replace('χ²/ν', r'$\chi^2/\nu$').replace('°C', r'$^\circ$C')
            .replace('L_ef', r'$L_\mathrm{ef}$').replace('f0', r'$f_0$').replace('fB', r'$f_B$')
            .replace('−', '$-$'))


linhas = [r'\begin{table}[h]', r'\centering', r'\small',
          r'\caption{Balanço de incertezas: contribuição de cada fonte e total em quadratura.}',
          r'\label{tab:incertezas}', r'\begin{tabular}{lrr}', r'\toprule',
          r'Fonte & Contribuição & Relativa \\', r'\midrule']
for b in blocos:
    linhas.append(r'\multicolumn{3}{l}{\textbf{' + tex(b['nome']) + f" = {fmt(b['valor'])} {tex(b['unidade'])}" + r'}} \\')
    for desc, u, r in b['linhas']:
        rel = f'{r:.2f} ' + r'\%' if r is not None else '--'
        linhas.append(f"\\quad {tex(desc)} & {fmt(u, 2)} {tex(b['unidade'])} & {rel} \\\\")
    rel = f"{b['rel']:.1f} " + r'\%' if b['rel'] is not None else '--'
    linhas.append(f"\\quad \\textit{{total}} & \\textit{{{fmt(b['u'], 2)} {tex(b['unidade'])}}} & \\textit{{{rel}}} \\\\")
    linhas.append(r'\addlinespace')
linhas += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
with open(os.path.join(PASTA, 'tabela_incertezas.tex'), 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(linhas) + '\n')
print("Tabela LaTeX: tabela_incertezas.tex (precisa de \\usepackage{booktabs})")


# ---------- PNG (só para ver a tabela; o .tex continua a ser o do relatório) ----------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

filas = []   # (texto da fonte, contribuição, relativa, tipo)
for b in blocos:
    filas.append((f"{b['nome']} = {fmt(b['valor'])} {b['unidade']}", '', '', 'titulo'))
    for desc, u, r in b['linhas']:
        filas.append(('    ' + desc, f"{fmt(u, 2)} {b['unidade']}", f'{r:.2f} %' if r is not None else '–', 'fonte'))
    filas.append(('    total (quadratura)', f"{fmt(b['u'], 2)} {b['unidade']}",
                  f"{b['rel']:.1f} %" if b['rel'] is not None else '–', 'total'))

alt = 0.32
fig, ax = plt.subplots(figsize=(8.5, alt * (len(filas) + 2)))
ax.axis('off')
ax.set_xlim(0, 1)
ax.set_ylim(len(filas) + 1.5, -0.5)
X_FONTE, X_CONTR, X_REL = 0.01, 0.78, 0.99
ax.text(X_FONTE, 0, 'Fonte', weight='bold', va='center')
ax.text(X_CONTR, 0, 'Contribuição', weight='bold', va='center', ha='right')
ax.text(X_REL, 0, 'Relativa', weight='bold', va='center', ha='right')
ax.axhline(-0.45, color='black', linewidth=1.2)
ax.axhline(0.5, color='black', linewidth=0.7)
for i, (txt, contr, rel, tipo) in enumerate(filas, start=1):
    y = i + 0.2
    if tipo == 'titulo':
        if i > 1:
            ax.axhline(y - 0.55, color='#ccc', linewidth=0.6)
        ax.text(X_FONTE, y, txt, weight='bold', va='center')
        continue
    estilo = 'italic' if tipo == 'total' else 'normal'
    ax.text(X_FONTE, y, txt, style=estilo, va='center')
    ax.text(X_CONTR, y, contr, style=estilo, va='center', ha='right')
    ax.text(X_REL, y, rel, style=estilo, va='center', ha='right')
ax.axhline(len(filas) + 0.75, color='black', linewidth=1.2)
ax.set_title('Balanço de incertezas: contribuição de cada fonte e total em quadratura', fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(PASTA, 'tabela_incertezas.png'), dpi=150)
print("Imagem: tabela_incertezas.png")


# ---------- Tabela C.1: incertezas de entrada (de onde vêm) ----------
# (grandeza, instrumento, menor divisão, incerteza, tipo, onde entra)
ENTRADAS = [
    ('V (êmbolo)', 'seringa', '1 mL', f'±{U_V * 1e6:.1f} mL', 'B', 'f0 previsto, calibração'),
    ('R_B, R_P, parede', 'craveira', '0.05 mm', f'±{U_R * 1e3:.1f} mm', 'B', 'S, ΔL'),
    ('a, distâncias', 'fita', '1 mm', f'±{U_A * 1e3:.0f} mm', 'B', 'fB, distância'),
    ('T', 'termómetro', '–', f'±{U_T} °C (estimativa)', 'B', 'c'),
    ('f (pico)', 'análise espectral', f'{RES:.0f} Hz', f'±{RES / 2:.0f} Hz', 'B', 'f0 medido'),
    ('TL, f0', '3 repetições', '–', 's/√3', 'A', 'valores medidos'),
]
CABEC = ('Grandeza', 'Instrumento', 'Menor divisão', 'Incerteza', 'Tipo', 'Onde entra')

def tex_ent(s):
    return tex(s.replace('R_B', 'R_B').replace('ΔL', r'$\Delta L$')).replace('R_B', r'$R_B$').replace('R_P', r'$R_P$')

linhas = [r'\begin{table}[h]', r'\centering', r'\small',
          r'\caption{Incertezas de entrada: instrumento, menor divisão e incerteza usada. '
          r'As de leitura são metade da menor divisão.}',
          r'\label{tab:entradas}', r'\begin{tabular}{llllcl}', r'\toprule',
          ' & '.join(CABEC) + r' \\', r'\midrule']
for fila in ENTRADAS:
    linhas.append(' & '.join(tex_ent(x) for x in fila) + r' \\')
linhas += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
with open(os.path.join(PASTA, 'tabela_entradas.tex'), 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(linhas) + '\n')

fig, ax = plt.subplots(figsize=(10, 0.42 * (len(ENTRADAS) + 3)))
ax.axis('off')
ax.set_xlim(0, 1)
ax.set_ylim(len(ENTRADAS) + 1, -0.6)
XS = [0.01, 0.19, 0.37, 0.51, 0.73, 0.79]
for x, h in zip(XS, CABEC):
    ax.text(x, 0, h, weight='bold', va='center')
ax.axhline(-0.5, color='black', linewidth=1.2)
ax.axhline(0.5, color='black', linewidth=0.7)
for i, fila in enumerate(ENTRADAS, start=1):
    for x, txt in zip(XS, fila):
        ax.text(x, i, txt, va='center')
ax.axhline(len(ENTRADAS) + 0.5, color='black', linewidth=1.2)
ax.set_title('Incertezas de entrada (leitura = metade da menor divisão)', fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(PASTA, 'tabela_entradas.png'), dpi=150)
print("Tabela de entradas: tabela_entradas.tex e tabela_entradas.png")
