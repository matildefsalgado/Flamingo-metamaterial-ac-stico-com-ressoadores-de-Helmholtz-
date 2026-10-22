# Ressonância local vs. Bragg

**O que mede:** a relação entre a ressonância local das seringas (f0) e a banda de Bragg da rede (fB = c/2a). Quando f0 ≈ fB as bandas fundem-se; com mais volume separam-se.

**Scripts:**
- `tl_parteA.py`: fusão, 5 seringas a 4,1 mL. Gera `tl_parteA_so_fusao.png` e escreve f0, dessintonia δ e larguras a 6 e 10 dB.
- `separacao.py`: separação, 5 seringas a 5 mL. Gera `separacao.png` e escreve o pico local, o pico de Bragg e a comparação com a fusão.
- `fusao.py`: fusão a 11,1 mL, tentativa anterior à calibração do gargalo. Gera `fusao.png`.

Os dois primeiros partilham o estilo e o modelo de matrizes de transferência de `../estilo_tl.py`.

**Gravações (todas nesta pasta):**
- fusão 4,1 mL (13/10): `fusao_ref.wav`, `fusao1.wav`, `fusao2.wav`, `fusao3.wav`
- separação 5 mL (16/10): `ref16_10.wav`, `sep_bandas1.wav`, `sep_bandas2.wav`, `sep_bandas3.wav`
- fusão 11,1 mL (12/10): `fusao_11-1mL_1/2/3.wav`, com a referência `../numero_elementos/A_ref_Aula3.wav`

**Parâmetros:** grelha de 10 Hz (±5 Hz nas frequências), a = 0,25 m, Q = 13 no modelo quando não há controlo e δV = ±0,5 mL. Falta o controlo com uma só seringa a 5 mL, que seria `ctrl5mL_<n>.wav`.

Detalhes em `LEIA-ME.txt`.
