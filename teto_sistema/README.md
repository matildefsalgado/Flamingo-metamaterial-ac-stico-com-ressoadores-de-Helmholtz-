# Teto do sistema

**O que mede:** a maior atenuação que a montagem consegue medir em cada frequência, limitada pelo ruído da sala e pelo som que chega ao microfone por fora do tubo. Compara a terminação de lã com a de LEGO.

**Script:** `teto.py`, que lê as quatro gravações e gera `teto.png` com as curvas do teto, `20·log10(A_ref / max(A_silencio, A_fuga))`, para a lã e para o LEGO.

**Gravações (7/10):** `A_ref.wav` (tubo aberto, lã), `A_silencio.wav` (gerador desligado), `A_fuga.wav` (tubo tapado à saída do altifalante, lã) e `A_fuga_lego.wav` (o mesmo, com LEGO).

**Parâmetros:** Welch com segmentos de 0,1 s (resolução de 10 Hz), entre 100 e 2000 Hz. Este teto é de outra sessão: nos outros scripts serve só de indicação.

Detalhes em `READ-ME.txt`.
