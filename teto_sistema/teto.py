import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import welch

# Pasta onde está este script (os .wav estão aqui também)
PASTA = os.path.dirname(os.path.abspath(__file__))


def ler_audio(nome):
    samplerate, data = wavfile.read(os.path.join(PASTA, nome))
    data = data.astype(float)
    # Se for estéreo, converte para mono
    if data.ndim > 1:
        data = np.mean(data, axis=1)
    return samplerate, data


def espectro(data, samplerate):
    # NÃO normalizar: o teto é uma razão entre gravações,
    # por isso as amplitudes têm de ficar na mesma escala.
    # Welch: média de muitos pedaços de 0,1 s -> curva lisa, resolução 10 Hz,
    # e as frequências coincidem mesmo que as gravações tenham durações diferentes.
    xf, P = welch(data, samplerate, nperseg=samplerate // 10)
    return xf, np.sqrt(P)


def calcular_teto(nome_fuga, nome_ref='A_ref.wav', nome_silencio='A_silencio.wav'):
    """Teto TL_max(f) = 20*log10(A_ref / max(A_silencio, A_fuga)), entre 100 e 2000 Hz."""
    sr_ref, ref = ler_audio(nome_ref)
    sr_sil, sil = ler_audio(nome_silencio)
    sr_fug, fug = ler_audio(nome_fuga)
    assert sr_ref == sr_sil == sr_fug, "As gravações têm taxas de amostragem diferentes"

    frequencias, A_ref = espectro(ref, sr_ref)
    _, A_silencio = espectro(sil, sr_ref)
    _, A_fuga = espectro(fug, sr_ref)

    # Intervalo do varrimento (100 Hz a 2000 Hz)
    indices = (frequencias >= 100) & (frequencias <= 2000)

    eps = 1e-12
    A_ref_c = np.maximum(A_ref[indices], eps)
    ruido_maximo = np.maximum(np.maximum(A_silencio[indices], A_fuga[indices]), eps)
    return frequencias[indices], 20 * np.log10(A_ref_c / ruido_maximo)


# Teto 1: terminação de lã | Teto 2: terminação com parede de LEGO
freq, teto_la = calcular_teto('A_fuga.wav')
_, teto_lego = calcular_teto('A_fuga_lego.wav')

print(" f (Hz)   teto lã   teto LEGO   diferença")
for f0 in [300, 500, 600, 670, 700, 800, 1000, 1200, 1340, 1500, 1800, 2000]:
    i = np.argmin(np.abs(freq - f0))
    print(f"{f0:6d} {teto_la[i]:9.1f} {teto_lego[i]:11.1f} {teto_lego[i] - teto_la[i]:11.1f}")

# Gráfico
plt.figure(figsize=(10, 6))
plt.plot(freq, teto_la, label='Teto: terminação de lã', color='red', linewidth=1)
plt.plot(freq, teto_lego, label='Teto: parede de LEGO', color='blue', linewidth=1)
plt.title('Teto de medição', fontsize=14)
plt.xlabel('Frequência (Hz)')
plt.ylabel('TL máximo mensurável (dB)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(PASTA, 'teto.png'), dpi=150)
plt.show()
