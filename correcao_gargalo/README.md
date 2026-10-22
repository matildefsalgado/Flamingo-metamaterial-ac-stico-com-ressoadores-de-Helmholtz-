# Correção do gargalo

**O que mede:** f0 de uma seringa para vários volumes (10 a 50 mL). Do ajuste de f0 = m/√(V + V_cone) tira-se o comprimento efetivo do gargalo, L_ef, e a correção ΔL = L_ef − L.

**Script:** `correcao_gargalo.py --sem 10`. Lê as gravações e gera `IL_por_volume.png` e `f0_vs_invsqrtV_sem10.png`, e escreve no terminal L_ef, ΔL e χ²/ν.
- `--sem 10` tira os 10 mL do ajuste. É a versão do relatório.
- `--ref` usa a `A_ref_7outopenlab.wav` como referência em vez da mediana. Só serve para estimar o sistemático.

**Gravações (7/10):** `fo_<V>mL_<rep>.wav` (V = 10, 20, 30, 40, 50 mL, 3 repetições cada) e `A_ref_7outopenlab.wav` (só com `--ref`).

**Parâmetros:** a referência é a **mediana dos 15 espectros**, porque a A_ref desse dia tinha flutuações de potência. Resolução de cerca de 1,5 Hz e δV = ±0,5 mL. Resultado: L_ef = 40,7 ± 2,5 mm e ΔL = 5,7 mm. Consoante a referência escolhida, ΔL varia entre 2,4 e 5,7 mm.

## Passo a passo do código

```text
CORREÇÃO DO GARGALO
===================

O SCRIPT QUE CONTA PARA O RELATÓRIO
-----------------------------------
    python correcao_gargalo.py --sem 10
(o f0_vs_V.py é a versão antiga, com a A_ref_7outopenlab como referência; fica só para comparação)

AS GRAVAÇÕES (nesta pasta)
--------------------------
fo_<V>mL_<rep>.wav        Só a seringa 1 com volume V (10, 20, 30, 40, 50 mL), as outras no fundo.
                          3 repetições por volume.
A_ref_7outopenlab.wav     Êmbolos no fundo, lã na ponta. NÃO serve de referência (gravada noutra
                          altura, ~10 dB acima e com ondulação de ±4 dB). Só com --ref.


O QUE O CÓDIGO FAZ (correcao_gargalo.py), PASSO A PASSO
-------------------------------------------------------
1. Espectro de cada gravação inteira (Welch, resolução ~1,5 Hz).
2. Referência = MEDIANA dos 15 espectros: cada volume só altera 3 dos 15, por isso a mediana
   fica com a resposta do altifalante, dos modos do tubo e do microfone do mesmo dia.
3. IL = 10 × log10( P_ref / P )  [dB] para cada gravação.
4. f0 = máximo de IL entre 150 e 1500 Hz, com várias suavizações (20 a 70 Hz) por causa da
   ondulação de ~30 Hz do tubo: f0 = média, u_sist = desvio padrão entre suavizações.
5. Por volume: média das 3 repetições. Incerteza de cada ponto, em quadratura:
       desvio padrão da média (s/√3)  +  u_sist  +  volume (½ × 0,5 mL / V × f0)  +  1 Hz
   O ±0,5 mL é metade da menor divisão da seringa (1 mL).
6. Ajuste ponderado f0 = m / √(V + V_cone)  (reta pela origem, pesos 1/σ²).
   Escreve χ²/ν. Se χ²/ν > 1, a incerteza de m é multiplicada por √(χ²/ν).
7. L_ef = c² S / (4π² m²), com incerteza de m, de S (raios ±0,1 mm) e de c (T ±0,5 °C)
   em quadratura;  ΔL = L_ef − L, comparado com 0,73 (R_B + R_P).
8. Figuras: IL_por_volume.png e f0_vs_invsqrtV[_sem10].png.

Resultado (--sem 10): L_ef = 40,7 ± 2,5 mm, ΔL = 5,7 mm (teoria 5,0–6,3 mm), χ²/ν = 5,0/3.
Os 10 mL ficam fora do ajuste (f0 muito abaixo da reta). Sistemático da referência, discutido
à parte: construir a mediana sem o próprio volume dá ΔL = 4,7 mm, sem o volume e os vizinhos
(±10 mL) dá 2,4 mm (o ponto frágil são os 30 mL); a A_ref_7outopenlab dá 3,3 mm.
```
