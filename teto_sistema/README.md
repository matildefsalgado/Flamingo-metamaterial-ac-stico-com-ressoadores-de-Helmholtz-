# Teto do sistema

**O que mede:** a maior atenuação que a montagem consegue medir em cada frequência, limitada pelo ruído da sala e pelo som que chega ao microfone por fora do tubo. Compara a terminação de lã com a de LEGO.

**Script:** `teto.py`, que lê as quatro gravações e gera `teto.png` com as curvas do teto, `20·log10(A_ref / max(A_silencio, A_fuga))`, para a lã e para o LEGO.

**Gravações (7/10):** `A_ref.wav` (tubo aberto, lã), `A_silencio.wav` (gerador desligado), `A_fuga.wav` (tubo tapado à saída do altifalante, lã) e `A_fuga_lego.wav` (o mesmo, com LEGO).

**Parâmetros:** Welch com segmentos de 0,1 s (resolução de 10 Hz), entre 100 e 2000 Hz. Este teto é de outra sessão: nos outros scripts serve só de indicação.

## Passo a passo do código

```text
TETO DO SISTEMA
===============

AS GRAVAÇÕES
------------
A_ref.wav        Tubo aberto, lã na ponta, gerador ligado.
                 É o som "normal" que chega ao microfone.
A_silencio.wav   Gerador desligado.
                 É o ruído da sala e da eletrónica.
A_fuga.wav       Tubo tapado à saída do buzzer, lã na ponta.
                 É o som que chega ao microfone por fora do tubo.
A_fuga_lego.wav  Igual à anterior, mas com parede de LEGO em vez de lã.
                 Serve para comparar com a fuga com lã.


O QUE O CÓDIGO FAZ (teto.py), PASSO A PASSO
----------------------------------------
1. Lê as gravações.
2. Para cada frequência (de 10 em 10 Hz, de 100 a 2000 Hz) calcula a amplitude
   do som nessa frequência (média de muitos pedaços de 0,1 s).
3. Em cada frequência escolhe o maior entre o silêncio e a fuga.
4. Calcula o teto:
       Teto = 20 × log10( A_ref / maior(A_silencio, A_fuga) )     [dB]
5. Faz o passo 3 e 4 duas vezes: com A_fuga.wav (lã) e com A_fuga_lego.wav (LEGO).
6. Escreve uma tabela com os dois tetos em algumas frequências.
7. Desenha as duas curvas no mesmo gráfico e guarda-o em teto.png.
```
