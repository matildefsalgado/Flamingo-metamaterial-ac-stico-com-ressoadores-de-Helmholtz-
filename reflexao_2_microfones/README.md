# Reflexão de uma seringa e onda estacionária

**O que mede:** o coeficiente de reflexão |R| de uma seringa a 4,0 mL e a onda estacionária no tubo, comparando a terminação de lã com o tampão de LEGO. Usa um só microfone em duas posições: a montante (1) e a jusante (2).

**Script:** `reflexao.py`. Lê as gravações e gera `reflexao.png`, com os quatro painéis também em separado (`reflexao_a…d_*.png`). No terminal escreve |R| medido e teórico, o espaçamento Δf da ondulação e o |R| do tampão.

**Gravações:**
- 14/10: `posicao1_la_0ml.wav`, `posicao1_la_4ml.wav`, `posicao1_lego_0ml.wav`, `posicao1_lego_4ml.wav`, `posicao2_lego_0ml.wav`, `posicao2_lego_4ml.wav`
- 16/10: `ref16_10.wav` e `pos2la4ml_FINAL.wav` (lã, posição 2), tons `655_<n>.wav` e `655_<n>_lego.wav`

**Parâmetros:** D_MIC_PAREDE = 155,5 cm, tubo de 25 mm de diâmetro interior e resolução de cerca de 1,5 Hz. Há uma só gravação por varrimento, e a posição do microfone nos tons está por confirmar. A análise da onda estacionária está em `anexo_onda_estacionaria.txt`.

## Passo a passo do código

```text
REFLEXÃO DE UMA SERINGA E ONDA ESTACIONÁRIA (1 MICROFONE)
=========================================================

MONTAGEM
--------
Altifalante na ponta fechada do tubo, ligado ao gerador de funções (varrimento 100–2000 Hz, ~10 s).
Um só microfone ligado ao PC, Audacity a 48 kHz, MONO, melhoramentos de áudio do Windows desligados.
O microfone muda de furo entre gravações:
  posição 1 = montante (antes da 1.ª seringa): ouve a onda incidente + a refletida
  posição 2 = jusante (depois da última seringa): ouve a onda transmitida
Seringa 1 a V = 4,0 mL (próximo dos 4,1 mL da fusão). As outras sempre com o êmbolo no fundo.
NÃO mexer no gerador, no volume do PC nem no ganho durante TODAS as gravações.


AS GRAVAÇÕES (nesta pasta)
--------------------------
Com a LÃ na ponta (solta do lado das seringas, cada vez mais compacta até à ponta):
  la_pos2_ref.wav        mic na posição 2, êmbolo da seringa 1 no fundo
  la_pos2_seringa.wav    mic na posição 2, seringa 1 a 4,0 mL
  la_pos1_ref.wav        mic na posição 1, êmbolo no fundo
  la_pos1_seringa.wav    mic na posição 1, seringa 1 a 4,0 mL
Com o TAMPÃO DE LEGO (selado com plasticina) no lugar da lã:
  lego_pos2_ref.wav, lego_pos2_seringa.wav, lego_pos1_ref.wav, lego_pos1_seringa.wav
Tons (seno contínuo no gerador, à frequência f0 que o script dá a partir dos varrimentos),
mic na posição 1, ~5 s cada:                                                    [opcional]
  tom_la.wav             com a lã
  tom_lego.wav           com o LEGO
Repetições (recomendado, 3 de cada): la_pos2_seringa_1.wav, la_pos2_seringa_2.wav, ...
O código faz a média das que encontrar. Se faltar um par, essa parte é saltada.

MEDIR COM A FITA (e pôr no topo do código):
  D_MIC_PAREDE   distância do furo da posição 1 à face interior do tampão de LEGO: 155,5 cm (medido 16/10)
                 -> Δf previsto c/2d = 111 Hz; medido 100 Hz (-10 %), que bate melhor com o tubo
                 inteiro (c/2L = 98 Hz para L = 1,75 m): a onda vai e volta entre o altifalante e o tampão.
  D_TUBO         diâmetro INTERIOR do tubo: 25 mm (medido)


O QUE O CÓDIGO FAZ (reflexao.py), PASSO A PASSO
-----------------------------------------------
1. Lã, posição 2:  TL = 20 × log10( A_ref / A_seringa ).  O máximo perto do previsto é o f0.
   Incerteza do f0 previsto, em quadratura: (δf0/f0)² = (δc/c)² + ¼(δS/S)² + ¼(δV/V_total)² + ¼(δL_ef/L_ef)²,
   com δT = 0,5 °C, δV = 0,5 mL (metade de 1 mL), δR = 0,1 mm e δL_ef = 2,5 mm (calibração).
   Em f0:  |T| = A_seringa / A_ref,   |R| = √(1 − |T|²)   (supõe que não há perdas).
2. Modelo das impedâncias:
       Z_c = ρ c / S_tubo
       Z_b = R_a + jωρ L_ef / (π R_B R_P) + ρ c² / ( jω (V + V_cone) )
       R_teórico = −Z_c / (2 Z_b + Z_c),     T_teórico = 2 Z_b / (2 Z_b + Z_c)
   Sem atrito (R_a = 0) o modelo dá |R| = 1 em f0. O atrito R_a é ajustado às ABAS do pico de TL
   (sem usar o topo) e depois o modelo prevê o |R| e o TL no pico. Assim a comparação com o
   medido é independente.
   Balanço de energia: com atrito, |R|² + |T|² < 1 (parte fica no gargalo), e √(1 − |T|²)
   SOBRESTIMA o |R|. O código diz quanta energia é absorvida.
3. LEGO, posição 2: o mesmo TL, para comparar com a lã (efeito da terminação).
4. Posição 1: TL = 20·log(A_ref / A_seringa) (lã e LEGO): fica negativo onde o som refletido pela seringa se soma.
5. Onda estacionária, com os varrimentos de REFERÊNCIA na posição 1:
       TL do LEGO em relação à lã, 20·log(A_lã / A_LEGO), em função de f   e de  d/λ = d·f/c
   Dividir pela lã tira a resposta do altifalante e do microfone; fica só a reflexão no tampão.
   Espaçamento entre máximos Δf (autocorrelação) -> d = c/(2Δf), comparado com a fita.
   Ondulação pico a vale Δ (dB) -> |R_tampão| = (10^(Δ/20) − 1) / (10^(Δ/20) + 1).
6. Tons: TL do LEGO em relação à lã à frequência do seno; deve cair em cima da curva do passo 5.
7. Figura reflexao.png: TL lã/LEGO na posição 2, |R| medido vs teórico, onda estacionária vs f e vs d/λ,
   resumo. O gráfico da posição 1 já não se desenha (só confirma a reflexão; os números ficam no terminal).
   Os 4 gráficos também saem em separado: reflexao_a_TL_posicao2.png, reflexao_b_coef_reflexao.png,
   reflexao_c_onda_estacionaria_f.png, reflexao_d_onda_estacionaria_d_lambda.png.

A fase/desfasamento NÃO é calculada: com um só microfone as gravações não são simultâneas
e o gerador não está sincronizado com o PC.

Uso:  python reflexao.py
```
