# Decaimento da ressonância

**O que mede:** o tempo de decaimento τ da seringa 1 (14,6 mL) depois de se cortar um tom a 430 Hz, e daí Q = π·f0·τ. Compara com um controlo, com todos os êmbolos no fundo.

**Script:** `decaimento.py`. Lê as gravações e gera `decaimento.png`, e escreve τ da seringa e do controlo, a razão entre os dois e Q.
- `--teste` corre o script com sinais sintéticos de τ conhecido.
- `--f0`, `--banda` e `--pasta` mudam a frequência, o filtro e a pasta das gravações.

**Gravações (22/10):** `decaimento 1.wav` (seringa) e `controlo (1).wav` (controlo).

**Parâmetros:** filtro de f0 ± 150 Hz, envolvente de Hilbert suavizada a 1 ms e ajuste de −3 a −25 dB. Resultado: τ da seringa ≈ τ do controlo ≈ 19 ms, por isso o decaimento é do tubo e só se pode dizer que Q < 26.

Detalhes em `READ-ME.txt`.
