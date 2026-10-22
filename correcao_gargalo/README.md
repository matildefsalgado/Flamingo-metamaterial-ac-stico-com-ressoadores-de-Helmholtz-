# Correção do gargalo

**O que mede:** f0 de uma seringa para vários volumes (10 a 50 mL). Do ajuste de f0 = m/√(V + V_cone) tira-se o comprimento efetivo do gargalo, L_ef, e a correção ΔL = L_ef − L.

**Script:** `correcao_gargalo.py --sem 10`. Lê as gravações e gera `IL_por_volume.png` e `f0_vs_invsqrtV_sem10.png`, e escreve no terminal L_ef, ΔL e χ²/ν.
- `--sem 10` tira os 10 mL do ajuste. É a versão do relatório.
- `--ref` usa a `A_ref_7outopenlab.wav` como referência em vez da mediana. Só serve para estimar o sistemático.

**Gravações (7/10):** `fo_<V>mL_<rep>.wav` (V = 10, 20, 30, 40, 50 mL, 3 repetições cada) e `A_ref_7outopenlab.wav` (só com `--ref`).

**Parâmetros:** a referência é a **mediana dos 15 espectros**, porque a A_ref desse dia tinha flutuações de potência. Resolução de cerca de 1,5 Hz e δV = ±0,5 mL. Resultado: L_ef = 40,7 ± 2,5 mm e ΔL = 5,7 mm. Consoante a referência escolhida, ΔL varia entre 2,4 e 5,7 mm.

Detalhes em `READ-ME.txt`.
