"""Sentinel-1 radar (RTC, Microsoft Planetary Computer) görüntülerinden gemi sayımı.

Neden: AIS'e dayanan geçiş verisi (PortWatch) gemilerin kendi yayınına bağlı; radar gemiyi yayın olmasa da görür. Radar bulut ve
geceden etkilenmez; gemi, koyu deniz yüzeyinde parlak bir nokta olarak görünür.

Yöntem:
  1. Her konum için bir kutu (boylam/enlem) tanımlanır; sahneler Planetary Computer STAC'tan bulunur.
  2. VH polarizasyonu 20 m çözünürlükte okunur (VH'de deniz yüzeyi yansıması VV'ye göre daha zayıf,
     gemi/deniz karşıtlığı daha yüksek), dB'ye çevrilir.
  3. Sabit parlak hedefler (kara, liman yapıları, platformlar, şamandıralar) maskelenir: savaş öncesi
     (1 Temmuz 2025-27 Şubat 2026) sahnelerin piksel bazında medyanı eşiği aşan pikseller + Natural Earth karası
     (300 m tampon).
  4. Eşik: max(-18 dB, deniz medyanı + 9 dB). Eşiği aşan 8-komşulu kümelerden alanı >= 3 piksel (20 m'de
     >= 1.200 m2) olanlar gemi sayılır; alanı >= 12 piksel olanlar "büyük gemi" (kabaca 150 m ve üstü).
  5. Kutunun deniz alanının %90'ından azını kapsayan sahneler atlanır.
Sınırlar: gemi tipi ayırt edilmez (tanker, konteyner, savaş gemisi aynı sayılır); yan yana demirli gemiler tek
küme görünebilir; rüzgârlı günlerde yanlış alarm artabilir. Sayım bir düzey göstergesi değil, eğilim göstergesidir.
Çıktı: veri/sar/gemi_sayimi_<konum>.csv (sahne başına) ve veri/sar/kontrol/*.png (tespit kontrol görüntüleri).
"""
from __future__ import annotations

import json
import sys
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.windows import from_bounds
from scipy import ndimage as ndi
from shapely import clip_by_rect
from shapely.geometry import shape
from shapely.ops import transform as shp_transform

BURASI = Path(__file__).resolve().parent
CIKTI = BURASI / "veri" / "sar"
(CIKTI / "kontrol").mkdir(parents=True, exist_ok=True)
NE = BURASI / "veri" / "harita" / "ne_10m_admin_0_countries.geojson"
STAC = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
RES = float(__import__("os").environ.get("SAR_RES", "20"))
YORUNGE = __import__("os").environ.get("SAR_YORUNGE", "hepsi")  # hepsi | descending | ascending
EK = __import__("os").environ.get("SAR_EK", "")  # çıktı dosya adı eki, ör. _40m

KONUMLAR = {
    # ad: (batı, güney, doğu, kuzey)
    "fuceyre": (56.37, 25.03, 56.62, 25.37),       # Füceyre demirleme alanı (Hürmüz dışı, Umman Körfezi)
    "hurmuz": (56.05, 26.30, 56.75, 26.75),        # Hürmüz geçiş koridoru (Musandam-Larak arası)
    "rastanura": (49.95, 26.45, 50.45, 26.95),     # Ras Tanura ve Cuayme terminalleri (Körfez içi)
    "yanbu": (37.85, 23.75, 38.35, 24.20),         # Yanbu limanı ve demirleme alanı (Kızıldeniz)
    "babulmendep": (43.15, 12.40, 43.70, 12.90),   # Bâbülmendep boğazı (Perim çevresi)
    "basra": (48.70, 29.40, 49.20, 29.85),         # Basra açık deniz terminalleri (ABOT/Mina al-Bakr, Khor al-Amaya) ve demirleme alanı
}


def sas() -> str:
    j = json.load(urllib.request.urlopen("https://planetarycomputer.microsoft.com/api/sas/v1/token/sentinel-1-rtc", timeout=60))
    return j["token"]


def sahneler(bbox, bas, bit):
    out, body = [], {"collections": ["sentinel-1-rtc"], "bbox": list(bbox), "datetime": f"{bas}T00:00:00Z/{bit}T23:59:59Z", "limit": 250}
    url = STAC
    while url:
        req = urllib.request.Request(url, data=json.dumps(body).encode(), headers={"Content-Type": "application/json"})
        j = json.load(urllib.request.urlopen(req, timeout=180))
        out += j["features"]
        nxt = [l for l in j.get("links", []) if l.get("rel") == "next"]
        if nxt and nxt[0].get("body"):
            body = nxt[0]["body"]
        else:
            url = None
    return sorted(out, key=lambda f: f["properties"]["datetime"])


def kara_maskesi(bbox, epsg, transform, shape_):
    tr = Transformer.from_crs(4326, epsg, always_xy=True)
    j = json.loads(NE.read_text(encoding="utf-8"))
    geoms = []
    for f in j["features"]:
        g = shape(f["geometry"])
        x0, y0, x1, y1 = g.bounds
        if x1 < bbox[0] - 1 or x0 > bbox[2] + 1 or y1 < bbox[1] - 1 or y0 > bbox[3] + 1:
            continue
        g = clip_by_rect(g, bbox[0] - 0.2, bbox[1] - 0.2, bbox[2] + 0.2, bbox[3] + 0.2)
        if not g.is_empty:
            geoms.append(shp_transform(tr.transform, g).buffer(300))
    if not geoms:
        return np.zeros(shape_, bool)
    return rasterize([(g, 1) for g in geoms], out_shape=shape_, transform=transform, fill=0).astype(bool)


def oku(item, bbox, token):
    epsg = int(str(item["properties"].get("proj:epsg") or item["properties"]["proj:code"]).split(":")[-1])
    tr = Transformer.from_crs(4326, epsg, always_xy=True)
    x0, y0 = tr.transform(bbox[0], bbox[1])
    x1, y1 = tr.transform(bbox[2], bbox[3])
    href = item["assets"]["vh"]["href"] + "?" + token
    with rasterio.open(href) as src:
        w = from_bounds(x0, y0, x1, y1, src.transform)
        nx, ny = int(round((x1 - x0) / RES)), int(round((y1 - y0) / RES))
        a = src.read(1, window=w, out_shape=(ny, nx), boundless=True, fill_value=0, resampling=Resampling.average)
    transform = rasterio.transform.from_bounds(x0, y0, x1, y1, nx, ny)
    return a.astype("float32"), epsg, transform


def db(a):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(a > 0, 10 * np.log10(a), np.nan)


def say(d, sabit, esik_taban=-18.0, fark=9.0):
    deniz = np.isfinite(d) & ~sabit
    if deniz.sum() == 0:
        return None
    med = float(np.nanmedian(d[deniz]))
    esik = max(esik_taban, med + fark)
    b = (d > esik) & deniz
    lab, n = ndi.label(b, structure=np.ones((3, 3)))
    if n == 0:
        return {"n": 0, "n_buyuk": 0, "deniz_medyan_db": med, "esik_db": esik, "merkez": []}
    alan = ndi.sum(b, lab, range(1, n + 1))
    ok = alan >= (3 if RES <= 20 else 2)
    buyuk = alan >= (12 if RES <= 20 else 3)
    merkez = ndi.center_of_mass(b, lab, np.arange(1, n + 1)[ok])
    return {"n": int(ok.sum()), "n_buyuk": int(buyuk.sum()), "deniz_medyan_db": med, "esik_db": esik, "merkez": merkez,
            "alan": alan[ok].tolist()}


MASKE_ESIK_DB, MASKE_GENISLETME, MASKE_ORNEK, SAVAS_ONCESI_SON = -17.0, 2, 24, "2026-02-28"


def maske_kur(ad, S, bbox, token):
    """Sabit hedef maskesi. Savaş öncesi sahneler (S tarih sıralı) arasından yaklaşık 24 örnek alınır (her adim'inci sahne); ilk okunan
    sahnenin EPSG'si ve ızgarası esas alınır, ızgarası uymayan sahne atlanır. Piksel bazında medyanı -17 dB'yi aşan ya da hiç veri olmayan
    pikseller iki adım genişletilir ve Natural Earth karası (300 m tampon) eklenir. Dönüş: maske, EPSG ve künye (kullanılan ve atlanan
    sahneler, ızgara, parametreler); künye maskenin sonradan aynı sahnelerle yeniden kurulabilmesi için yazılır."""
    onces = [s for s in S if s["properties"]["datetime"][:10] < SAVAS_ONCESI_SON]
    adim = max(1, len(onces) // MASKE_ORNEK)
    yigin, epsg0, T0, kullanilan, atlanan = [], None, None, [], []
    for s in onces[::adim]:
        try:
            a, epsg, T = oku(s, bbox, token)
        except Exception as e:  # ağ hatası vb.
            print("  okunamadı", s["id"], repr(e)[:80], flush=True)
            atlanan.append({"sahne": s["id"], "neden": "okunamadı: " + repr(e)[:80]})
            continue
        if epsg0 is None:
            epsg0, T0 = epsg, T
        if epsg != epsg0 or a.shape != (yigin[0].shape if yigin else a.shape):
            atlanan.append({"sahne": s["id"], "goreli_yorunge": s["properties"].get("sat:relative_orbit"),
                            "neden": f"ızgara uymuyor (EPSG {epsg}, {a.shape[0]}x{a.shape[1]})"})
            continue
        yigin.append(np.where(a > 0, db(a), np.nan))
        kullanilan.append({"sahne": s["id"], "yorunge": s["properties"].get("sat:orbit_state"),
                           "goreli_yorunge": s["properties"].get("sat:relative_orbit")})
    med = np.nanmedian(np.stack(yigin), axis=0)
    sabit = (med > MASKE_ESIK_DB) | ~np.isfinite(med)
    sabit = ndi.binary_dilation(sabit, iterations=MASKE_GENISLETME) | kara_maskesi(bbox, epsg0, T0, med.shape)
    kunye = {"konum": ad, "kutu_bgdk": list(bbox), "cozunurluk_m": RES, "yorunge_secimi": YORUNGE, "dosya_eki": EK,
             "epsg": epsg0, "sekil_satir_sutun": list(med.shape), "donusum_affine": [round(float(v), 6) for v in tuple(T0)[:6]],
             "aday_savas_oncesi_sahne": len(onces), "ornek_adimi": adim, "kullanilan_sahne": kullanilan, "atlanan_sahne": atlanan,
             "parametre": {"medyan_esigi_db": MASKE_ESIK_DB, "genisletme_adimi": MASKE_GENISLETME, "hedef_ornek": MASKE_ORNEK,
                           "savas_oncesi": f"< {SAVAS_ONCESI_SON}", "kara": "Natural Earth 10m admin 0, 300 m tampon", "kanal": "VH"},
             "maskelenen_pay": round(float(sabit.mean()), 5)}
    return sabit, epsg0, kunye


def maske_kunyesi(ad, bas="2025-07-01", bit="2026-09-22"):
    """Diskteki maskeyi aynı kuralla yeniden kurar, bit düzeyinde karşılaştırır ve künyeyi yazar (sayımlara dokunmaz). Ortam değişkenleri
    (SAR_RES, SAR_YORUNGE, SAR_EK) maskenin kurulduğu çalıştırmadakiyle aynı olmalı; run_radar.sh'teki komutlara bakın."""
    import hashlib
    bbox = KONUMLAR[ad]
    S = sahneler(bbox, bas, bit)
    if YORUNGE != "hepsi":
        S = [s for s in S if s["properties"].get("sat:orbit_state") == YORUNGE]
    sabit, epsg0, kunye = maske_kur(ad, S, bbox, sas())
    mf, mj = CIKTI / f"{ad}{EK}_sabit_maske.npy", CIKTI / f"{ad}{EK}_sabit_maske.json"
    disk = np.load(mf)
    ayni = disk.shape == sabit.shape and bool(np.array_equal(disk, sabit)) and json.loads(mj.read_text())["epsg"] == epsg0
    kunye["maske_sha256"] = hashlib.sha256(mf.read_bytes()).hexdigest()
    kunye["yeniden_kurulum"] = {"tarih": datetime.now().strftime("%Y-%m-%d"), "diskteki_maskeyle_ayni": ayni,
                                "farkli_piksel": None if disk.shape != sabit.shape else int((disk != sabit).sum()),
                                "kod_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "sahne_araligi": [bas, bit]}
    (CIKTI / f"{ad}{EK}_sabit_maske_kunye.json").write_text(json.dumps(kunye, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"{ad}{EK}: {len(kunye['kullanilan_sahne'])} sahne, EPSG {epsg0}, diskteki maskeyle aynı: {ayni}")


def konum_isle(ad, bas="2025-07-01", bit="2026-09-22", kontrol_sayisi=6, en_fazla=None):
    import csv
    bbox = KONUMLAR[ad]
    token = sas()
    S = sahneler(bbox, bas, bit)
    if YORUNGE != "hepsi":
        S = [s for s in S if s["properties"].get("sat:orbit_state") == YORUNGE]
    if en_fazla:
        S = S[:en_fazla]
    print(f"{ad}: {len(S)} sahne", flush=True)
    # 1) savaş öncesi medyanı ile sabit hedef maskesi (daha önce hesaplandıysa diskten okunur)
    mf, mj = CIKTI / f"{ad}{EK}_sabit_maske.npy", CIKTI / f"{ad}{EK}_sabit_maske.json"
    if mf.exists() and mj.exists():
        sabit = np.load(mf)
        epsg0 = json.loads(mj.read_text())["epsg"]
        print(f"  sabit maske diskten okundu (EPSG {epsg0})", flush=True)
    else:
        sabit, epsg0, kunye = maske_kur(ad, S, bbox, token)
        np.save(mf, sabit)
        mj.write_text(json.dumps({"epsg": epsg0}))
        kunye["maske_sha256"] = __import__("hashlib").sha256(mf.read_bytes()).hexdigest()
        (CIKTI / f"{ad}{EK}_sabit_maske_kunye.json").write_text(json.dumps(kunye, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"  sabit maske: denizin %{100 * (1 - sabit.mean()):.0f}'ı kullanılabilir", flush=True)
    # 2) sahne başına sayım
    f = CIKTI / f"gemi_sayimi_{ad}{EK}.csv"
    yeni = not f.exists()
    bitmis = set()
    if not yeni:  # kaldığı yerden devam: işlenmiş sahneleri atla
        import csv as _csv
        with open(f) as fh0:
            bitmis = {r["sahne"] for r in _csv.DictReader(fh0)}
        print(f"  {len(bitmis)} sahne daha önce işlenmiş, atlanacak", flush=True)
    with open(f, "a", newline="") as fh:
        wr = csv.writer(fh)
        if yeni:
            wr.writerow(["konum", "tarih", "saat", "sahne", "yorunge", "kapsama", "gemi", "buyuk_gemi", "deniz_medyan_db", "esik_db"])
        kontrol = 0
        for i, s in enumerate(S):
            if s["id"] in bitmis:
                continue
            if i % 15 == 0:
                token = sas()  # SAS belirteci yaklaşık bir saatte doluyor
            try:
                a, epsg, T = oku(s, bbox, token)
            except Exception as e:
                if "403" in repr(e):  # belirteç dolmuş: yenile ve bir kez daha dene
                    token = sas()
                    try:
                        a, epsg, T = oku(s, bbox, token)
                    except Exception as e2:
                        print("  okunamadı", s["id"], repr(e2)[:80], flush=True)
                        continue
                else:
                    print("  okunamadı", s["id"], repr(e)[:80], flush=True)
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
            p = s["properties"]
            wr.writerow([ad, p["datetime"][:10], p["datetime"][11:16], s["id"], p.get("sat:orbit_state"), round(kapsama, 3),
                         r["n"], r["n_buyuk"], round(r["deniz_medyan_db"], 2), round(r["esik_db"], 2)])
            fh.flush()
            if kontrol < kontrol_sayisi and (i % max(1, len(S) // kontrol_sayisi) == 0):
                kontrol_goruntusu(ad, p["datetime"][:10], d, sabit, r)
                kontrol += 1
            if i % 20 == 0:
                print(f"  {p['datetime'][:10]} gemi={r['n']} büyük={r['n_buyuk']} medyan={r['deniz_medyan_db']:.1f}dB", flush=True)


def goreli_yorunge(ad, bas="2025-07-01", bit="2026-09-22"):
    """Sahne -> göreli yörünge (STAC sat:relative_orbit). Aynı yöndeki iki göreli yörünge kutuya farklı geliş açısıyla bakar; aynı seride
    birleşince modelde göreli yörünge kuklası gerekir. Çıktı: veri/sar/goreli_yorunge_<konum>.csv"""
    import csv
    with open(CIKTI / f"goreli_yorunge_{ad}.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["sahne", "yorunge", "goreli_yorunge"])
        for s in sahneler(KONUMLAR[ad], bas, bit):
            wr.writerow([s["id"], s["properties"].get("sat:orbit_state"), s["properties"].get("sat:relative_orbit")])


def kontrol_goruntusu(ad, tarih, d, sabit, r):
    from PIL import Image, ImageDraw
    x = np.clip((np.nan_to_num(d, nan=-30) + 30) / 30, 0, 1)
    rgb = np.stack([x, x, x], -1)
    rgb[sabit] = rgb[sabit] * 0.35 + np.array([0.35, 0.25, 0.1]) * 0.65
    im = Image.fromarray((rgb * 255).astype(np.uint8))
    dr = ImageDraw.Draw(im)
    for (y, xx) in r["merkez"]:
        dr.ellipse([xx - 6, y - 6, xx + 6, y + 6], outline=(255, 200, 0), width=2)
    im.save(CIKTI / "kontrol" / f"{ad}{EK}_{tarih}.png")


if __name__ == "__main__":
    if sys.argv[1] == "--goreli":  # yalnız göreli yörünge eşlemesi: python sar_gemi_sayimi.py --goreli basra
        goreli_yorunge(sys.argv[2])
        sys.exit()
    if sys.argv[1] == "--maske-kunye":  # maskeyi yeniden kur, diskteki ile karşılaştır, künyeyi yaz: SAR_... python sar_gemi_sayimi.py --maske-kunye basra
        maske_kunyesi(sys.argv[2])
        sys.exit()
    ad = sys.argv[1]
    bas = sys.argv[2] if len(sys.argv) > 2 else "2025-07-01"
    bit = sys.argv[3] if len(sys.argv) > 3 else "2026-09-22"
    en_fazla = int(sys.argv[4]) if len(sys.argv) > 4 else None
    konum_isle(ad, bas, bit, en_fazla=en_fazla)
