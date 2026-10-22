# Ressonância em função da distância ao altifalante

**O que mede:** f0, TL máximo e Q de uma só seringa ativa (14,6 mL), mudando qual das 5 está ativa, ou seja, a distância ao altifalante.

**Script:** `distancia_altifalante.py`. Lê as gravações e gera `f0_vs_distancia.png`, com TL(f) por distância, f0 vs d, TL máx vs d e Q vs d. No terminal escreve o χ² de "f0 constante" e o declive de f0 vs d.

**Gravações (14/10):** `referencia_14-10.wav` (êmbolos no fundo) e `distancia_<k>-<rep>.wav` (seringa k ativa, 3 repetições).

**Parâmetros:** V = 14,6 mL, espaçamento a = 0,25 m, resolução de cerca de 1,5 Hz e pico procurado em f0 previsto ± 35 %.

## Passo a passo do código

```text
RESSONÂNCIA LOCAL EM FUNÇÃO DA DISTÂNCIA AO ALTIFALANTE
=======================================================

A IDEIA
-------
Uma só seringa ativa de cada vez, sempre com o MESMO volume. As outras ficam com o
êmbolo no fundo. Muda-se só QUAL seringa está ativa, ou seja, a distância d ao altifalante.
Pergunta: f0, a atenuação máxima e o Q mudam com d?
  - f0 não deve mudar: depende só da seringa (volume e gargalo).
  - TL máximo pode mudar: se a ponta do tubo refletir, há ondas estacionárias e a seringa
    bloqueia mais num máximo de pressão do que num nodo. O padrão repete-se a cada
    λ/2 = c/(2 f0) (≈ 43 cm a 14,6 mL). As seringas estão a 25 cm umas das outras, ou seja
    ~0,6 × λ/2: o TL máximo deve variar de seringa para seringa sem seguir uma tendência simples.


PROTOCOLO (tudo seguido, SEM mexer em volume, ganho ou microfone)
-----------------------------------------------------------------
1. A_ref.wav           todos os êmbolos no fundo, lã na ponta
2. A_silencio.wav      5 s com o altifalante desligado (teto desta sessão)
3. Seringa 1 a V, as outras no fundo: 3 gravações  -> d<cm>_1.wav, d<cm>_2.wav, d<cm>_3.wav
4. Voltar a pôr a seringa 1 no fundo, seringa 2 a V: 3 gravações
5. ... até à seringa 5
6. A_ref_fim.wav       igual à 1 (controlo: tem de coincidir com A_ref)

Nome dos ficheiros: d + distância em cm + _ + repetição. Exemplos:
    d32_1.wav     d57,6_2.wav     d83.2_3.wav
A distância mede-se com a fita, da membrana do altifalante ao eixo do gargalo da seringa.
Gravar a 48 kHz, como o resto.

Volume: V = 14,6 mL no topo do código (f0 ≈ 399 Hz), o mesmo da série de N, para poder
comparar com o N = 1. Se usarem outro, MUDAR V no código antes de correr.

Validar logo no laboratório: a coluna "TL fora" da tabela tem de dar ~0 dB e o controlo
A_ref_fim tem de dizer OK. Se não, a montagem mexeu e convém repetir ali mesmo.


O QUE O CÓDIGO FAZ (distancia_altifalante.py), PASSO A PASSO
------------------------------------------------------------
1. Calcula o f0 previsto (gargalo cónico, V + ar do cone) e λ/2 a essa frequência.
   Incerteza do previsto, em quadratura: (δf0/f0)² = (δc/c)² + ¼(δS/S)² + ¼(δV/V_total)² + ¼(δL_ef/L_ef)²,
   com δT = 0,5 °C, δV = 0,5 mL (metade de 1 mL), δR = 0,1 mm e δL_ef = 2,5 mm (calibração).
2. Calcula o espectro de cada gravação (Welch, resolução ~1,5 Hz).
3. Compara A_ref_fim com A_ref (desvio médio e ondulação) -> diz se a montagem mudou.
4. Para cada gravação:  TL = 20 × log10( A_ref / A_medida )   [dB]
5. f0 = máximo de TL perto do previsto (± 35 %), procurado com várias suavizações
   (20 a 70 Hz, por causa da ondulação do tubo): f0 = média, incerteza = dispersão.
6. TL_max = altura desse pico; Q = f0 / largura do pico a TL_max − 3 dB.
7. "TL fora" = TL médio longe do pico (< 0,6 f0 e > 1,6 f0): deve dar ~0. Avisa se passar 2 dB
   (referência não corresponde às medições). Avisa também se o pico estiver perto do teto.
8. Média das 3 repetições para cada d.
9. f0 depende de d? Média ponderada + chi2 (f0 constante?) e reta f0 = a + k·d
   (k compatível com 0?). Compara a média com o f0 previsto (z < 2 = acordo).
10. Gráficos em f0_vs_distancia.png:
       - TL(f) para cada distância;
       - f0 vs d, com a faixa do previsto;
       - TL máximo vs d (comparar uma eventual oscilação com λ/2);
       - Q vs d.

Uso:  python distancia_altifalante.py
```
