"""Hürmüz geçiş koridorunun iki radar görüntüsü: savaş öncesi ve ikinci kapanma.

Seçim kuralı (elle seçim yok): iki görüntü de aynı uydudan (Sentinel-1C; 1A farklı eşik üretiyor, bkz. analiz7_sar.py).
Her dönemde gemi sayısı o dönemin 1C medyanına en yakın sahne alınır; eşitlikte en yeni tarih. Böylece görüntü dönemin
tipik gününü gösterir, en çarpıcı gününü değil.
Girdi: veri/sar/gemi_sayimi_hurmuz_40m.csv ve hurmuz_40m_sabit_maske.npy (sar_gemi_sayimi.py, SAR_RES=40, inen yörünge).
Çıktı: pdf/gorsel/hurmuz_radar_once.jpg, pdf/gorsel/hurmuz_radar_kapanma2.jpg, sonuclar/a10_hurmuz_radar.json
Kullanım: SAR_RES=40 SAR_YORUNGE=descending SAR_EK=_40m python sar_hurmuz_gorsel.py
"""
from __future__ import annotations

import json
import os

os.environ.setdefault("SAR_RES", "40")
os.environ.setdefault("SAR_YORUNGE", "descending")
os.environ.setdefault("SAR_EK", "_40m")

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from pyproj import Transformer
from rasterio.features import rasterize
from shapely import clip_by_rect
from shapely.geometry import shape
from shapely.ops import transform as shp_transform

import sar_gemi_sayimi as G

BURASI = G.BURASI
BBOX = G.KONUMLAR["hurmuz"]
DON = {"once": ("2025-07-01", "2026-02-27"), "kapanma2": ("2026-07-13", "2026-09-22")}


def kara(epsg, transform, shape_):
    """Görüntü için kara (tampon yok)."""
    tr = Transformer.from_crs(4326, epsg, always_xy=True)
    geoms = []
    for f in json.loads(G.NE.read_text(encoding="utf-8"))["features"]:
        g = shape(f["geometry"])
        x0, y0, x1, y1 = g.bounds
        if x1 < BBOX[0] - 1 or x0 > BBOX[2] + 1 or y1 < BBOX[1] - 1 or y0 > BBOX[3] + 1:
            continue
        g = clip_by_rect(g, BBOX[0] - 0.2, BBOX[1] - 0.2, BBOX[2] + 0.2, BBOX[3] + 0.2)
        if not g.is_empty:
            geoms.append(shp_transform(tr.transform, g))
    return rasterize([(g, 1) for g in geoms], out_shape=shape_, transform=transform, fill=0).astype(bool)


def ciz(d, sabit, karam, r, dosya):
    # deniz: lacivert tonlu gri; -27 dB koyu, -12 dB ve üstü açık. Kara düz renk, radar dokusu yüzde 12.
    x = np.clip((np.nan_to_num(d, nan=-30) + 27) / 15, 0, 1) ** 0.8
    lac = np.array([11, 26, 46]) / 255
    rgb = lac + (np.array([0.93, 0.95, 0.97]) - lac) * x[..., None]
    rgb[karam] = np.array([0x5b, 0x63, 0x6f]) / 255 * 0.88 + rgb[karam] * 0.12
    im = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8))
    dr = ImageDraw.Draw(im)
    for (y, xx) in r["merkez"]:  # 85 mm genişlikte yaklaşık 1,2 mm yarıçaplı halka
        dr.ellipse([xx - 24, y - 24, xx + 24, y + 24], outline=(214, 178, 88), width=5)
    im.convert("RGB").save(BURASI / "pdf" / "gorsel" / dosya, quality=88, optimize=True)


def main():
    S = pd.read_csv(G.CIKTI / "gemi_sayimi_hurmuz_40m.csv", parse_dates=["tarih"])
    S = S[S.sahne.str.startswith("S1C")]
    sabit = np.load(G.CIKTI / "hurmuz_40m_sabit_maske.npy")
    epsg0 = json.loads((G.CIKTI / "hurmuz_40m_sabit_maske.json").read_text())["epsg"]
    sahneler = {s["id"]: s for s in G.sahneler(BBOX, "2025-07-01", "2026-09-22")}
    token = G.sas()
    out = {}
    for ad, (a, b) in DON.items():
        g = S[(S.tarih >= a) & (S.tarih <= b)].copy()
        med = float(g.gemi.median())
        g["uzak"] = (g.gemi - med).abs()
        sec = g.sort_values(["uzak", "tarih"], ascending=[True, False]).iloc[0]
        arr, epsg, T = G.oku(sahneler[sec.sahne], BBOX, token)
        assert epsg == epsg0 and arr.shape == sabit.shape
        d = G.db(arr)
        r = G.say(d, sabit)
        karam = kara(epsg, T, d.shape)
        ciz(d, sabit, karam, r, f"hurmuz_radar_{ad}.jpg")
        out[ad] = {"tarih": sec.tarih.strftime("%Y-%m-%d"), "sahne": sec.sahne, "gemi": int(r["n"]), "gemi_csv": int(sec.gemi),
                   "donem_medyan_1c": med, "donem_sahne": int(len(g)), "boyut_px": [d.shape[1], d.shape[0]]}
        print(ad, out[ad], flush=True)
    (BURASI / "sonuclar" / "a10_hurmuz_radar.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
