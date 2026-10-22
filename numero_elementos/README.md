# Dependência com o número de seringas

**O que mede:** TL máximo e largura da banda a 10 dB em função do número de seringas ativas, de N = 1 a 5, todas a 14,6 mL.

**Script:** `tl_vs_N.py`. Lê as gravações e gera `tl_vs_N.png`, com TL(f) para cada N, TL máx vs N com o ajuste linear e a largura vs N. No terminal escreve a tabela, a inclinação k (dB por seringa) e o χ²/ν.

**Gravações (duas sessões, cada uma com a sua referência):**
- aula 3 (12/10): `A_ref_Aula3.wav` + `dependencia_N1..N3-<rep>.wav`
- 14/10: `referencia_14-10 (1).wav` + `dependencia_N4/N5-<rep>_nova*.wav`

O teto vem de `../teto_sistema/`.

**Parâmetros:** grelha de 10 Hz e pico procurado entre 0,75 e 1,35 × f0. As duas sessões não têm medição de ponte, e o ajuste dá χ²/ν = 3,6.

Detalhes em `READ-ME.txt`.
