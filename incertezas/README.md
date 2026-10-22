# Balanço de incertezas

**O que calcula:** a contribuição de cada fonte de incerteza (T, V, raios do cone, L_ef, a, resolução) para f0 previsto, fB, f0 medido na fusão, dessintonia δ e L_ef, somadas em quadratura.

**Script:** `balanco_incertezas.py`. Não lê gravações. Gera `tabela_entradas.tex/.png` (tabela C.1) e `tabela_incertezas.tex/.png` (tabela C.2). Os `.tex` são os que vão para o Overleaf e precisam de `\usepackage{booktabs}`.

**Dados:** os valores medidos estão copiados no topo do script, com indicação do script de onde vêm. Se esses scripts mudarem, é preciso atualizá-los aqui.

**Incertezas de entrada:** T = 22 ± 0,5 °C, V ± 0,5 mL, R ± 0,1 mm, L_ef ± 2,5 mm, a ± 1 mm e f ± 5 Hz.
