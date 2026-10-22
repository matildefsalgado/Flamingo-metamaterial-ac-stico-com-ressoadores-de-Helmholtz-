# Reflexão de uma seringa e onda estacionária

**O que mede:** o coeficiente de reflexão |R| de uma seringa a 4,0 mL e a onda estacionária no tubo, comparando a terminação de lã com o tampão de LEGO. Usa um só microfone em duas posições: a montante (1) e a jusante (2).

**Script:** `reflexao.py`. Lê as gravações e gera `reflexao.png`, com os quatro painéis também em separado (`reflexao_a…d_*.png`). No terminal escreve |R| medido e teórico, o espaçamento Δf da ondulação e o |R| do tampão.

**Gravações:**
- 14/10: `posicao1_la_0ml.wav`, `posicao1_la_4ml.wav`, `posicao1_lego_0ml.wav`, `posicao1_lego_4ml.wav`, `posicao2_lego_0ml.wav`, `posicao2_lego_4ml.wav`
- 16/10: `ref16_10.wav` e `pos2la4ml_FINAL.wav` (lã, posição 2), tons `655_<n>.wav` e `655_<n>_lego.wav`

**Parâmetros:** D_MIC_PAREDE = 155,5 cm, tubo de 25 mm de diâmetro interior e resolução de cerca de 1,5 Hz. Há uma só gravação por varrimento, e a posição do microfone nos tons está por confirmar. A análise da onda estacionária está em `anexo_onda_estacionaria.txt`.

Detalhes em `LEIA-ME.txt`.
