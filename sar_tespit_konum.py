"""Füceyre demirleme alanında radar tespitlerinin konumlarını kaydeder (yoğunluk haritası için).
sar_gemi_sayimi.py'deki okuma, maske ve eşik mantığını aynen kullanır; yalnız inen yörünge sahneleri.
Çıktı: veri/sar/tespit_konum_fuceyre.csv (tarih, sahne, boylam, enlem, alan_piksel)."""
import csv
import sys
import numpy as np
from pyproj import Transformer
sys.argv = ["x"]
exec(open("sar_gemi_sayimi.py", encoding="utf-8").read().split('if __name__ == "__main__":')[0])

ad = "fuceyre"
bbox = KONUMLAR[ad]
sabit = np.load(CIKTI / f"{ad}_sabit_maske.npy")
epsg0 = json.loads((CIKTI / f"{ad}_sabit_maske.json").read_text())["epsg"]
token = sas()
S = [s for s in sahneler(bbox, "2025-07-01", "2026-09-22") if s["properties"].get("sat:orbit_state") == "descending"]
f = CIKTI / "tespit_konum_fuceyre.csv"
bitmis = set()
if f.exists():
    with open(f) as fh:
        bitmis = {r["sahne"] for r in csv.DictReader(fh)}
yeni = not f.exists()
with open(f, "a", newline="") as fh:
    wr = csv.writer(fh)
    if yeni:
        wr.writerow(["tarih", "sahne", "boylam", "enlem", "alan_piksel"])
    for i, s in enumerate(S):
        if s["id"] in bitmis:
            continue
        if i % 15 == 0:
            token = sas()
        try:
            a, epsg, T = oku(s, bbox, token)
        except Exception as e:
            token = sas()
            try:
                a, epsg, T = oku(s, bbox, token)
            except Exception as e2:
                print("okunamadı", s["id"], repr(e2)[:80], flush=True)
                continue
        if epsg != epsg0 or a.shape != sabit.shape:
            continue
        d = db(a)
        kapsama = float((np.isfinite(d) & ~sabit).sum() / max((~sabit).sum(), 1))
        if kapsama < 0.9:
            continue
        r = say(d, sabit)
        if r is None:
            continue
        inv = Transformer.from_crs(epsg, 4326, always_xy=True)
        for (yy, xx), alan in zip(r["merkez"], r["alan"]):
            X, Y = T * (xx + 0.5, yy + 0.5)
            lon, lat = inv.transform(X, Y)
            wr.writerow([s["properties"]["datetime"][:10], s["id"], round(lon, 5), round(lat, 5), int(alan)])
        fh.flush()
        print(s["properties"]["datetime"][:10], r["n"], flush=True)
print("BİTTİ", flush=True)
