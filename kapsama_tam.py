"""Analiz 16 (tam): alev ölçüsünün gözlenebilirliği, bütün geceler ve Tablo 3'ün on bölgesi. Yöntem kapsama_pilot.py ile aynı: NOAA-21
VJ214IMG (375 m yangın maskesi, QA bit 22, yangın pikselleri) ve CLDMSK_L2_VIIRS_NOAA21 (bulut maskesi, konum, bakış açısı).

Disk ve kesinti güvenliği (kullanıcı isteği, 26 Eylül 2026: SSD'de yer az, konum değişebilir):
  - Aynı anda en fazla PENCERE granül çifti iner; her çift işlenir işlenmez silinir. Boş alan MIN_BOS_GB altına düşerse süreç durur.
  - Her çiftin sonucu hemen kendi dosyasına yazılır (veri/firms/kapsama/parca/<anahtar>.csv.gz ve <anahtar>_ates.csv.gz); var olan
    parçalar atlanır. Süreç her an kesilebilir, yeniden başlatınca kaldığı yerden sürer.
  - --sure DAKIKA: o kadar süre sonra yeni çift başlatmaz, eldekileri bitirip temiz çıkar ("DURAKLAMA"). Hepsi bitince "TAMAMLANDI".
  - Arka arkaya ARDISIK_HATA indirme hatası (ağ yok) olursa temiz çıkar ("AĞ SORUNU").
  - Her çıkışta veri/firms/kapsama/ARA_NOTLAR.md dosyasına ilerleme, hız, tahmini kalan süre ve ara sonuçlar eklenir; durum.json güncellenir.
Granül listesi bir kez CMR'den alınır (veri/firms/kapsama/granuller.csv). Earthdata token'ı .env'de (EARTHDATA_TOKEN; değer yazdırılmaz).
Kullanım: python kapsama_tam.py --sure 60        (bir saatlik parti)
          python kapsama_tam.py --ozet            (indirmeden, eldeki parçalarla ara özet)
          python kapsama_tam.py --yeniden 10      (bitmiş 10 çifti rastgele seçip yeniden indirir, işler, kayıtlı parçalarla bayt düzeyinde
                                                   karşılaştırır; yarısı 27 Eylül 2026 öncesi yazılmış parçalardan, yarısı sonrasından)
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import time
import urllib.parse
import urllib.request
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import analiz15_alev as AL
import analiz_ortak as A
import kapsama_pilot as KP

KLASOR = A.VERI / "firms" / "kapsama"
PARCA = KLASOR / "parca"
GECICI = Path(os.environ.get("KAPSAMA_GECICI", KLASOR / "_gecici"))
BOLGELER = AL.TABLO
PENCERE = 5          # aynı anda inen en fazla çift (her biri ~53 MB + ~1 MB; en çok ~270 MB geçici dosya)
ISCI = 8
MIN_BOS_GB = 5.0
ARDISIK_HATA = 8
ZAMAN = ("2024-09-01", "2026-09-24")


def hucreler() -> pd.DataFrame:
    H = pd.read_csv(A.VERI / "firms" / "kalici_hucreler.csv")
    H = H[H["bolge"].isin(BOLGELER)].copy()
    H["anahtar"] = H["hx"].astype(np.int64) * 100000 + H["hy"].astype(np.int64)
    return H


def kutu(H: pd.DataFrame) -> tuple:
    return (round(float(H["lon"].min()) - 0.1, 2), round(float(H["lat"].min()) - 0.1, 2),
            round(float(H["lon"].max()) + 0.1, 2), round(float(H["lat"].max()) + 0.1, 2))


def granul_listesi(K: tuple) -> pd.DataFrame:
    """Bütün dönemin gece granülleri, iki ürün eşlenmiş. Bir kez alınır ve saklanır."""
    f = KLASOR / "granuller.csv"
    if f.exists():
        return pd.read_csv(f)
    def cek(kol):
        q = urllib.parse.urlencode({"collection_concept_id": kol, "bounding_box": ",".join(map(str, K)), "day_night_flag": "night",
                                    "temporal": f"{ZAMAN[0]}T00:00:00Z,{ZAMAN[1]}T23:59:59Z", "page_size": 2000})
        out, sa = {}, None
        while True:
            req = urllib.request.Request("https://cmr.earthdata.nasa.gov/search/granules.json?" + q,
                                         headers={"CMR-Search-After": sa} if sa else {})
            with urllib.request.urlopen(req, timeout=120) as r:
                j = json.loads(r.read()); sa = r.headers.get("CMR-Search-After"); hits = int(r.headers["CMR-Hits"])
            for e in j["feed"]["entry"]:
                for ln in e.get("links", []):
                    h = ln.get("href", "")
                    m = __import__("re").search(r"\.A(\d{7})\.(\d{4})\.", h)
                    if h.startswith("https://") and h.endswith(".nc") and m:
                        out[f"A{m.group(1)}.{m.group(2)}"] = (h, float(e.get("granule_size", 0) or 0))
            if not j["feed"]["entry"] or len(out) >= hits or not sa:
                return out
    a, b = cek(KP.KOL["ates"]), cek(KP.KOL["bulut"])
    ortak = sorted(set(a) & set(b))
    G = pd.DataFrame({"anahtar": ortak, "url_ates": [a[k][0] for k in ortak], "url_bulut": [b[k][0] for k in ortak],
                      "mb": [a[k][1] + b[k][1] for k in ortak]})
    G["tarih"] = [str(datetime.strptime(k[1:8], "%Y%j").date()) for k in G["anahtar"]]
    KLASOR.mkdir(parents=True, exist_ok=True)
    G.to_csv(f, index=False)
    pd.DataFrame({"anahtar": sorted(set(a) ^ set(b))}).to_csv(KLASOR / "esi_olmayan_granuller.csv", index=False)
    return G


def bos_gb() -> float:
    return shutil.disk_usage(str(GECICI)).free / 1e9


def isle_ve_yaz(k: str, fa, fb, anahtarlar: np.ndarray, hedef: Path = PARCA) -> tuple[int, float]:
    pa, pb = fa.result(), fb.result()
    try:
        hc, fp = KP.isle(pa, pb, anahtarlar)
        mb = (pa.stat().st_size + pb.stat().st_size) / 1e6
    finally:
        pa.unlink(missing_ok=True)
        pb.unlink(missing_ok=True)
    tarih = str(datetime.strptime(k[1:8], "%Y%j").date())
    (hc.assign(tarih=tarih, gecis=k) if len(hc) else pd.DataFrame({"anahtar": []})).to_csv(hedef / f"{k}.csv.gz", index=False)
    if len(fp):
        fp.assign(tarih=tarih, gecis=k).to_csv(hedef / f"{k}_ates.csv.gz", index=False)
    return len(hc), mb


def ara_ozet(H: pd.DataFrame) -> dict:
    """Eldeki parçalardan dönem x bölge geçerli gözlem ve bulut payları (hücre-gece düzeyi, pilotla aynı tanım)."""
    dosyalar = [f for f in PARCA.glob("A*.csv.gz") if not f.name.endswith("_ates.csv.gz")]
    if not dosyalar:
        return {}
    R = pd.concat([pd.read_csv(f) for f in dosyalar], ignore_index=True)
    R = R.dropna(subset=["n_m"])
    if R.empty:
        return {}
    n_i = 4 * R["n_m"]
    gec = R[[f"fm_{c}" for c in KP.GECERLI]].sum(axis=1) / n_i
    bul = R["fm_4"] / n_i
    R["durum"] = np.select([gec >= 0.5, bul >= 0.5], ["geçerli", "bulutlu"], "belirsiz")
    hg = R.groupby(["tarih", "anahtar"])["durum"].agg(lambda x: "geçerli" if (x == "geçerli").any() else
                                                      ("bulutlu" if (x == "bulutlu").all() else "belirsiz")).reset_index()
    hg = hg.merge(H[["anahtar", "bolge"]], on="anahtar")
    hg["donem"] = [AL.donem(pd.Timestamp(t)) or "geçiş" for t in hg["tarih"]]
    out = {}
    for (b, d), x in hg.groupby(["bolge", "donem"]):
        out.setdefault(b, {})[d] = {"hucre_gece": int(len(x)), "gecerli": round(float((x["durum"] == "geçerli").mean()), 3),
                                    "bulutlu": round(float((x["durum"] == "bulutlu").mean()), 3)}
    return out


def not_yaz(durum: dict, ozet: dict | None):
    f = KLASOR / "ARA_NOTLAR.md"
    yeni = not f.exists()
    with open(f, "a", encoding="utf-8") as g:
        if yeni:
            g.write("# Bulut ve kapsama tam çalıştırması: ara notlar\n\nHer parti sonunda eklenir. Devam: "
                    "`cd MetalTakip_Rapor && KAPSAMA_GECICI=<geçici klasör> python kapsama_tam.py --sure 60`. Parçalar "
                    "`veri/firms/kapsama/parca/` altında; var olanlar atlanır.\n\n")
        g.write(f"## {durum['zaman']} · {durum['sonuc']}\n\n")
        g.write(f"- Çift: {durum['biten']} / {durum['toplam']} (%{100 * durum['biten'] / durum['toplam']:.1f}); bu partide {durum['parti_cift']} çift, "
                f"{durum['parti_gb']:.1f} GB, {durum['parti_dk']:.0f} dakika ({durum['hiz_cift_dk']:.1f} çift/dk).\n")
        g.write(f"- Tahmini kalan süre: {durum['kalan_saat']:.1f} saat. Boş disk: {durum['bos_gb']:.0f} GB. Hatalı çift (bu parti): {durum['parti_hata']}.\n")
        if ozet:
            g.write("- Ara sonuç (hücre-gece geçerli gözlem payı / bulutlu payı):\n")
            for b in ("Irak (güney)", "Irak (kuzey)", "Katar", "Bahreyn"):
                if b in ozet:
                    g.write(f"  - {b}: " + "; ".join(f"{d} {v['gecerli']:.2f} / {v['bulutlu']:.2f} ({v['hucre_gece']})"
                                                       for d, v in sorted(ozet[b].items())) + "\n")
        g.write("\n")


def yeniden(n: int, anahtarlar: np.ndarray, G: pd.DataFrame, tohum: int = 2026) -> list[dict]:
    """Yeniden üretim denetimi: bitmiş n çift (yarısı 27 Eylül 2026 00:00'dan önce, yarısı sonra yazılmış) yeniden indirilir ve işlenir;
    gzip içeriği kayıtlı parçayla bayt düzeyinde karşılaştırılır. Sonuç ARA_NOTLAR.md'ye eklenir."""
    esik = datetime(2026, 9, 27).timestamp()
    bitmis = sorted(f for f in PARCA.glob("A*.csv.gz") if not f.name.endswith("_ates.csv.gz"))
    rng = np.random.default_rng(tohum)
    sec = []
    for grup in ([f for f in bitmis if f.stat().st_mtime < esik], [f for f in bitmis if f.stat().st_mtime >= esik]):
        m = min(len(grup), n // 2 if not sec else n - len(sec))
        sec += [f.name.split(".csv")[0] for f in rng.choice(grup, size=m, replace=False)] if m else []
    hedef = GECICI / "yeniden"
    hedef.mkdir(parents=True, exist_ok=True)
    sonuc = []
    with ThreadPoolExecutor(4) as havuz:
        for k in sec:
            r = G[G["anahtar"] == k].iloc[0]
            isle_ve_yaz(k, havuz.submit(KP.indir, r["url_ates"], hedef), havuz.submit(KP.indir, r["url_bulut"], hedef), anahtarlar, hedef)
            satir = {"anahtar": k}
            for ek in ("", "_ates"):
                a, b = PARCA / f"{k}{ek}.csv.gz", hedef / f"{k}{ek}.csv.gz"
                satir["hucre" if not ek else "ates"] = ("ikisi de yok" if not a.exists() and not b.exists() else
                                                        "aynı" if a.exists() and b.exists() and gzip.open(a).read() == gzip.open(b).read() else "FARKLI")
                b.unlink(missing_ok=True)
            print(satir, flush=True)
            sonuc.append(satir)
    with open(KLASOR / "ARA_NOTLAR.md", "a", encoding="utf-8") as g:
        g.write(f"## {datetime.now().strftime('%Y-%m-%d %H:%M')} · YENİDEN ÜRETİM DENETİMİ\n\n")
        g.write(f"- {len(sonuc)} çift yeniden indirildi ve işlendi; hücre dosyası aynı: {sum(x['hucre'] == 'aynı' for x in sonuc)}, "
                f"ateş dosyası aynı ya da ikisinde de yok: {sum(x['ates'] != 'FARKLI' for x in sonuc)}.\n")
        g.write("- Çiftler: " + ", ".join(f"{x['anahtar']} ({x['hucre']}/{x['ates']})" for x in sonuc) + "\n\n")
    return sonuc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sure", type=float, default=60.0, help="dakika")
    ap.add_argument("--ozet", action="store_true")
    ap.add_argument("--yeniden", type=int, default=0, help="yeniden üretim denetimi: çift sayısı")
    arg = ap.parse_args()
    H = hucreler()
    K = kutu(H)
    KP.KUTU = K  # KP.isle bu kutuyu kullanır
    PARCA.mkdir(parents=True, exist_ok=True)
    GECICI.mkdir(parents=True, exist_ok=True)
    G = granul_listesi(K)
    bitmis = {f.name.split(".csv")[0] for f in PARCA.glob("A*.csv.gz") if not f.name.endswith("_ates.csv.gz")}
    kalan = G[~G["anahtar"].isin(bitmis)].reset_index(drop=True)
    print(f"kutu {K} | çift {len(G)} | biten {len(bitmis)} | kalan {len(kalan)}", flush=True)
    if arg.ozet:
        print(json.dumps(ara_ozet(H), ensure_ascii=False, indent=1))
        return
    if arg.yeniden:
        yeniden(arg.yeniden, H["anahtar"].to_numpy(), G)
        return
    for f in GECICI.glob("*.nc"):  # önceki kesintiden kalan yarım dosyalar (--ozet ve --yeniden çalışan partiye dokunmasın diye burada)
        f.unlink()
    anahtarlar = H["anahtar"].to_numpy()
    t0, parti_cift, parti_mb, hata, ardisik = time.time(), 0, 0.0, [], 0
    sonuc = "TAMAMLANDI"
    kuyruk: deque = deque()
    it = iter(kalan.itertuples(index=False))
    with ThreadPoolExecutor(ISCI) as havuz:
        def doldur():
            while len(kuyruk) < PENCERE and (time.time() - t0) / 60 < arg.sure:
                r = next(it, None)
                if r is None:
                    return
                kuyruk.append((r.anahtar, havuz.submit(KP.indir, r.url_ates, GECICI), havuz.submit(KP.indir, r.url_bulut, GECICI)))
        doldur()
        while kuyruk:
            if bos_gb() < MIN_BOS_GB:
                sonuc = "DİSK DOLU"
                break
            k, fa, fb = kuyruk.popleft()
            try:
                _, mb = isle_ve_yaz(k, fa, fb, anahtarlar)
                parti_cift += 1
                parti_mb += mb
                ardisik = 0
            except Exception as e:  # noqa: BLE001
                hata.append({"anahtar": k, "hata": repr(e)[:160], "zaman": datetime.now().isoformat(timespec="seconds")})
                ardisik += 1
                for fut in (fa, fb):  # yarım kalan indirmeyi temizle
                    try:
                        fut.result().unlink(missing_ok=True)
                    except Exception:  # noqa: BLE001
                        pass
                if ardisik >= ARDISIK_HATA:
                    sonuc = "AĞ SORUNU"
                    break
            if parti_cift and parti_cift % 25 == 0:
                print(f"  {parti_cift} çift | {(time.time() - t0) / 60:.0f} dk | {parti_mb / 1e3:.1f} GB", flush=True)
            doldur()
        if sonuc == "TAMAMLANDI" and next(it, None) is not None:
            sonuc = "DURAKLAMA"
        for _, fa, fb in kuyruk:  # erken çıkışta kuyruktakileri bekle ve sil
            for fut in (fa, fb):
                try:
                    fut.result().unlink(missing_ok=True)
                except Exception:  # noqa: BLE001
                    pass
    for f in GECICI.glob("*.nc"):
        f.unlink()
    if hata:
        pd.DataFrame(hata).to_csv(KLASOR / "hatalar.csv", mode="a", index=False, header=not (KLASOR / "hatalar.csv").exists())
    biten = len({f.name.split(".csv")[0] for f in PARCA.glob("A*.csv.gz") if not f.name.endswith("_ates.csv.gz")})
    dk = (time.time() - t0) / 60
    hiz = parti_cift / dk if dk > 0 else 0.0
    durum = {"zaman": datetime.now().strftime("%Y-%m-%d %H:%M"), "sonuc": sonuc, "biten": biten, "toplam": int(len(G)),
             "parti_cift": parti_cift, "parti_gb": parti_mb / 1e3, "parti_dk": dk, "hiz_cift_dk": hiz,
             "kalan_saat": (len(G) - biten) / hiz / 60 if hiz else float("nan"), "bos_gb": bos_gb(), "parti_hata": len(hata)}
    (KLASOR / "durum.json").write_text(json.dumps(durum, ensure_ascii=False, indent=1), encoding="utf-8")
    not_yaz(durum, ara_ozet(H))
    print(f"SONUÇ: {sonuc} | {biten}/{len(G)} | parti {parti_cift} çift, {parti_mb / 1e3:.1f} GB, {dk:.0f} dk | kalan ~{durum['kalan_saat']:.1f} saat",
          flush=True)


if __name__ == "__main__":
    main()
