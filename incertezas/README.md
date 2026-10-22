# Balanço de incertezas

**O que calcula:** a contribuição de cada fonte de incerteza (T, V, raios do cone, L_ef, a, resolução) para f0 previsto, fB, f0 medido na fusão, dessintonia δ e L_ef, somadas em quadratura.

**Script:** `balanco_incertezas.py`. Não lê gravações. Gera `tabela_entradas.tex/.png` (tabela C.1) e `tabela_incertezas.tex/.png` (tabela C.2). Os `.tex` são os que vão para o Overleaf e precisam de `\usepackage{booktabs}`.

**Dados:** os valores medidos estão copiados no topo do script, com indicação do script de onde vêm. Se esses scripts mudarem, é preciso atualizá-los aqui.

**Incertezas de entrada:** T = 22 ± 0,5 °C, V ± 0,5 mL, R ± 0,1 mm, L_ef ± 2,5 mm, a ± 1 mm e f ± 5 Hz.

## Passo a passo do código

```text
BALANÇO DE INCERTEZAS (Apêndice C do relatório)
=================================================

Uso:  python balanco_incertezas.py

Para cada grandeza principal escreve a contribuição de cada fonte (em Hz ou mm e em %)
e a total em quadratura, e diz qual domina:
  - f0 previsto a 4,1 mL (fusão) e a 14,6 mL (série de N, distância)
  - fB previsto = c/2a
  - f0 medido na fusão
  - dessintonia δ = 1 − f0/fB
  - L_ef calibrado

Também escreve tabela_incertezas.tex, para colar no Overleaf (precisa de \usepackage{booktabs}).

Incertezas de entrada (iguais em todos os scripts):
  T     22 ± 0,5 °C
  V     ± 0,5 mL       metade da menor divisão da seringa (1 mL)
  R     ± 0,1 mm       raios do cone: craveira + parede
  L_ef  ± 2,5 mm       da calibração (correcao_gargalo.py --sem 10)
  a     ± 1 mm         fita
  f     ± 5 Hz         metade da resolução de 10 Hz das análises de TL

Os valores MEDIDOS (picos da fusão, declive da calibração) estão copiados no topo do script,
com a indicação do script de onde vêm. Se esses scripts voltarem a correr e os números
mudarem, atualizem-nos lá.

Também gera (tabelas do Apêndice C):
  C.1  tabela_entradas.tex   / tabela_entradas.png    incertezas de entrada (instrumento, divisão, tipo)
  C.2  tabela_incertezas.tex / tabela_incertezas.png  propagação, fonte a fonte
Os .png são só para ver; no Overleaf usam-se os .tex.
```
