"""Analiz 9: Füceyre demirleme alanında radar tespitlerinin mekânsal yoğunluğu (inen yörünge).
Girdi: veri/sar/tespit_konum_fuceyre.csv (sar_tespit_konum.py), veri/sar/fuceyre_sabit_maske.npy (kara ve sabit hedefler).
Çıktı: pdf/gorsel/yogunluk_{once,kapanma2}.png (aynı renk ölçeği), sonuclar/a9_yogunluk.json
Birim: 500 m'lik karede sahne başına ortalama tespit (gemi) sayısı. Yalnız Sentinel-1C ve 1D görüntüleri (1A farklı eşik üretiyor).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from PIL import Image
from pyproj import Transformer

import analiz_ortak as A

KUTU = (56.37, 25.03, 56.62, 25.37)
tr = Transformer.from_crs(4326, 32640, always_xy=True)
x0, y0 = tr.transform(KUTU[0], KUTU[1])
x1, y1 = tr.transform(KUTU[2], KUTU[3])
HUC = 500.0
nx, ny = int(np.ceil((x1 - x0) / HUC)), int(np.ceil((y1 - y0) / HUC))
D = pd.read_csv(A.VERI / "sar" / "tespit_konum_fuceyre.csv", parse_dates=["tarih"])
D = D[D.sahne.str[:3].isin(["S1C", "S1D"])]  # uydu kuralı: analiz7_sar.py başındaki not
X, Y = tr.transform(D.boylam.values, D.enlem.values)
D["i"] = ((y1 - np.asarray(Y)) // HUC).astype(int)
D["j"] = ((np.asarray(X) - x0) // HUC).astype(int)
D = D[(D.i >= 0) & (D.i < ny) & (D.j >= 0) & (D.j < nx)]
DON = {"once": ("2025-07-01", "2026-02-27"), "kapanma2": ("2026-07-13", "2026-09-22")}
yog, ozet = {}, {}
for k, (a, b) in DON.items():
    d = D[(D.tarih >= a) & (D.tarih <= b)]
    n = d.sahne.nunique()
    G = np.zeros((ny, nx))
    np.add.at(G, (d.i.values, d.j.values), 1)
    yog[k] = G / n
    ozet[k] = {"sahne": int(n), "sahne_basi_tespit": round(len(d) / n, 1)}
# kara maskesi: 20 m'lik sabit maskeyi 500 m'ye indir (karenin yarısından fazlası sabitse kara say)
M = np.load(A.VERI / "sar" / "fuceyre_sabit_maske.npy").astype(float)
f = int(HUC / 20)
h, w = (M.shape[0] // f) * f, (M.shape[1] // f) * f
Mk = M[:h, :w].reshape(h // f, f, w // f, f).mean(axis=(1, 3)) > 0.5
Mk = np.pad(Mk, ((0, max(0, ny - Mk.shape[0])), (0, max(0, nx - Mk.shape[1]))), constant_values=False)[:ny, :nx]
# Natural Earth karası (20 m) ile radar maskesinin birleşimi; kıyıya 1,5 km'den yakın hücreler (liman yapıları) hesap dışı
import sys
sys.argv = ["x"]
exec(open(A.HERE / "sar_gemi_sayimi.py", encoding="utf-8").read().split('if __name__ == "__main__":')[0].split("BURASI = ")[0])
from sar_gemi_sayimi import kara_maskesi  # noqa: E402
import rasterio.transform
T20 = rasterio.transform.from_bounds(x0, y0, x1, y1, M.shape[1], M.shape[0])
NEk = kara_maskesi(KUTU, 32640, T20, M.shape).astype(float)
NEk = NEk[:h, :w].reshape(h // f, f, w // f, f).mean(axis=(1, 3)) > 0.5
NEk = np.pad(NEk, ((0, max(0, ny - NEk.shape[0])), (0, max(0, nx - NEk.shape[1]))), constant_values=False)[:ny, :nx]
Mk = Mk | NEk
from scipy import ndimage as ndi
kiyi = ndi.distance_transform_edt(~Mk) * HUC / 1000 < 1.5
for k in yog:
    yog[k] = ndi.gaussian_filter(np.where(kiyi | Mk, 0, yog[k]), 1.0)
ust = max(np.percentile(yog["kapanma2"][yog["kapanma2"] > 0.01], 98), 1e-6)
RAMPA = np.array([[238, 242, 245], [246, 223, 211], [239, 191, 166], [224, 145, 111], [207, 106, 68], [180, 70, 31]]) / 255.0
def renk(v):
    t = np.clip(v / ust, 0, 1) * (len(RAMPA) - 1)
    i = np.floor(t).astype(int)
    i2 = np.minimum(i + 1, len(RAMPA) - 1)
    w_ = (t - i)[..., None]
    return RAMPA[i] * (1 - w_) + RAMPA[i2] * w_
for k in DON:
    rgb = renk(yog[k])
    rgb[kiyi & ~Mk] = np.array([244, 242, 237]) / 255.0
    rgb[Mk] = np.array([220, 212, 192]) / 255.0
    im = Image.fromarray((rgb * 255).astype(np.uint8)).resize((nx * 20, ny * 20), Image.BILINEAR)
    im.save(A.HERE / "pdf" / "gorsel" / f"yogunluk_{k}.png")
# açık denizdeki bölgelere göre artış: kıyıya uzaklık bantları (km)
from scipy import ndimage as ndi
uzak = ndi.distance_transform_edt(~Mk) * HUC / 1000
bant = {}
for a, b in ((1.5, 5), (5, 10), (10, 30)):
    m = (uzak >= a) & (uzak < b) & ~Mk
    bant[f"{a}-{b} km"] = {k: round(float(yog[k][m].sum()), 1) for k in DON}
ozet["kiyiya_uzaklik"] = bant
ozet["olcek_ust"] = round(float(ust), 3)
(A.SONUC / "a9_yogunluk.json").write_text(json.dumps(ozet, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps(ozet, ensure_ascii=False))
