"""Analiz 16 (pilot): alev ölçüsünün gözlenebilirliği. NOAA-21 VIIRS Level-2 dosyalarıyla kalıcı alev hücrelerinin her geçişteki gözlem
durumu (açık, bulutlu, işlenemedi), bakış açısı ve bowtie işareti. Üçüncü çapraz denetimin bulut/kapsama protokolünün ilk adımı.

Kaynaklar (Earthdata Login token'ı depo kökündeki .env dosyasında, EARTHDATA_TOKEN; değer yazdırılmaz):
  VJ214IMG v002 (LP DAAC): 375 m yangın maskesi. Sınıflar: 0 işlenemedi, 1 bowtie silmesi, 2 güneş parlaması, 3 su, 4 bulut, 5 açık kara,
    6 sınıflanamadı, 7-9 yangın (düşük, orta, yüksek güven). "algorithm QA" bit 22: kalan bowtie tekrarı (NASA kılavuzu, Tablo 2).
    Yangın pikselleri: FP_latitude, FP_longitude, FP_power (FRP), FP_ViewZenAng (bakış açısı), FP_line, FP_sample.
  CLDMSK_L2_VIIRS_NOAA21 v1 (LAADS): 750 m bulut maskesi (Integer_Cloud_Mask: 0 bulutlu, 1 muhtemelen bulutlu, 2 muhtemelen açık,
    3 kesin açık, -1 sonuç yok) ve konum (enlem, boylam, sensor_zenith x 0,01 derece). I-bant pikseli (2i+a, 2j+b), M-bant pikseli (i, j)
    içindedir (6432 x 6400 = 2 x 3216 x 3200).
Hücre: analiz15_alev.py'nin 0,01 derecelik kalıcı alev hücreleri (veri/firms/kalici_hucreler.csv). Bir hücrenin bir geçişteki durumu, konumu
  hücreye düşen M-bant piksellerinin içindeki I-bant piksellerinden: pikselin en az yarısı 3, 5 ya da 7-9 ise "geçerli", en az yarısı 4 ise
  "bulutlu", değilse "belirsiz".
Pilot geceler (tohum 2026): Katar'ın ilk kapanmada FIRMS'te hiç tespit olmayan geceleri; Bahreyn'in aynı türden gecelerinden 10'u; savaş
  öncesi, ilk ve ikinci kapanmadan 12'şer rastgele gece.
Dosyalar geçici klasöre (KAPSAMA_GECICI ortam değişkeni ya da sistem geçici klasörü) indirilir, işlenir, silinir. Çıktı: veri/firms/kapsama_pilot/hucre_gecis.csv.gz, sonuclar/a16_kapsama_pilot.json
Kullanım: python kapsama_pilot.py   (indirme ve işleme; yaklaşık 10 GB, 15 dakika)
          python kapsama_pilot.py --yalniz-ozet   (kayıtlı pilot dosyalarından özetleri yeniden hesaplar)
"""
from __future__ import annotations

import json
import os
import re
import tempfile
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

import analiz15_alev as AL
import analiz_ortak as A

KOK = Path(__file__).resolve().parent.parent
BOLGELER = ["Irak (güney)", "Irak (kuzey)", "Katar", "Kuveyt", "Bahreyn", "Batı İran (çoğu Huzistan)"]
KUTU = (42.5, 24.0, 52.0, 37.5)  # B, G, D, K; main() seçilen bölgelerin hücrelerinden 0,1 derece paylı olarak yeniden kurar
KOL = {"ates": "C2831626262-LPCLOUD", "bulut": "C3206162112-LAADS"}
CIKTI = A.VERI / "firms" / "kapsama_pilot"
GECERLI, BULUT, EKSIK = (3, 5, 7, 8, 9), (4,), (0, 1)
TOHUM = 2026


def token() -> str:
    for s in (KOK / ".env").read_text(encoding="utf-8").splitlines():
        if s.startswith("EARTHDATA_TOKEN="):
            return s.split("=", 1)[1].strip()
    raise SystemExit("EARTHDATA_TOKEN .env dosyasında yok")


TOK = token()


class _Yonlendir(urllib.request.HTTPRedirectHandler):
    """Earthdata alan adlarına yönlendirmede token'ı yeniden ekler; imzalı S3 adresine eklemez."""
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        yeni = super().redirect_request(req, fp, code, msg, headers, newurl)
        if yeni is not None and ("earthdata.nasa.gov" in newurl or "earthdatacloud.nasa.gov" in newurl):
            yeni.add_unredirected_header("Authorization", f"Bearer {TOK}")
        return yeni


ACICI = urllib.request.build_opener(_Yonlendir())


def indir(url: str, klasor: Path, deneme: int = 3) -> Path:
    hedef = klasor / url.rsplit("/", 1)[1]
    for k in range(deneme):
        try:
            req = urllib.request.Request(url, headers={"Authorization": f"Bearer {TOK}"})
            with ACICI.open(req, timeout=120) as r, open(hedef, "wb") as f:
                while blok := r.read(1 << 20):
                    f.write(blok)
            return hedef
        except Exception:  # noqa: BLE001
            if k == deneme - 1:
                raise
            time.sleep(5 * (k + 1))
    return hedef


def granuller(kol: str, gun: str) -> dict[str, str]:
    """Bir UTC günündeki gece granülleri: {"AYYYYDDD.HHMM": indirme adresi}."""
    q = urllib.parse.urlencode({"collection_concept_id": kol, "bounding_box": ",".join(map(str, KUTU)), "day_night_flag": "night",
                                "temporal": f"{gun}T00:00:00Z,{gun}T23:59:59Z", "page_size": 50})
    with urllib.request.urlopen("https://cmr.earthdata.nasa.gov/search/granules.json?" + q, timeout=120) as r:
        j = json.loads(r.read())
    out = {}
    yd = pd.Timestamp(gun).strftime("%Y%j")  # granül başlangıcı sorgu gününde olmalı (gece yarısını aşan granül iki güne düşmesin)
    for e in j["feed"]["entry"]:
        for ln in e.get("links", []):
            h = ln.get("href", "")
            m = re.search(r"\.A(\d{7})\.(\d{4})\.", h)
            if h.startswith("https://") and h.endswith(".nc") and m and m.group(1) == yd:
                out[f"A{m.group(1)}.{m.group(2)}"] = h
    return out


def geceler_sec(S: pd.DataFrame) -> pd.DataFrame:
    """Pilot geceleri ve seçilme nedenleri."""
    rng = np.random.default_rng(TOHUM)
    donem = pd.Series([AL.donem(t) for t in S.index], index=S.index)
    gozlenen = S.notna().any(axis=1)
    sec = {}
    k1 = gozlenen & (donem == "1. kapanma")
    for t in S.index[k1 & (S["Katar"] == 0)]:
        sec[t] = "Katar sıfır, ilk kapanma"
    bah = [t for t in S.index[k1 & (S["Bahreyn"] == 0)] if t not in sec]
    for t in rng.choice(bah, size=min(10, len(bah)), replace=False):
        sec[pd.Timestamp(t)] = "Bahreyn sıfır, ilk kapanma"
    for d in ("Savaş öncesi", "1. kapanma", "2. kapanma"):
        aday = [t for t in S.index[gozlenen & (donem == d)] if t not in sec]
        for t in rng.choice(aday, size=12, replace=False):
            sec[pd.Timestamp(t)] = f"rastgele, {d}"
    return pd.DataFrame({"tarih": list(sec), "neden": list(sec.values())}).sort_values("tarih").reset_index(drop=True)


def isle(ates: Path, bulut: Path, anahtarlar: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Bir granül çifti: hücre başına gözlem durumu ve yangın pikselleri."""
    with h5py.File(bulut, "r") as b:
        lat, lon = b["geolocation_data/latitude"][:], b["geolocation_data/longitude"][:]
        sz = b["geolocation_data/sensor_zenith"][:].astype(float)
        icm = b["geophysical_data/Integer_Cloud_Mask"][:]
    m = (lat >= KUTU[1]) & (lat <= KUTU[3]) & (lon >= KUTU[0]) & (lon <= KUTU[2])
    ii, jj = np.nonzero(m)
    anahtar = np.floor(lon[ii, jj] / AL.HUCRE).astype(np.int64) * 100000 + np.floor(lat[ii, jj] / AL.HUCRE).astype(np.int64)
    s = np.isin(anahtar, anahtarlar)
    ii, jj, anahtar = ii[s], jj[s], anahtar[s]
    with h5py.File(ates, "r") as a:
        fm = a["fire mask"][:]
        fp = {k: a[k][:] for k in ("FP_latitude", "FP_longitude", "FP_power", "FP_ViewZenAng", "FP_line", "FP_sample", "FP_confidence")}
        satirlar = np.unique(fp["FP_line"]) if len(fp["FP_line"]) else np.array([], dtype=int)
        qa_satir = {int(r): a["algorithm QA"][int(r), :] for r in satirlar}
    hucre = pd.DataFrame()
    if len(ii):
        sayim = {}
        for c in range(10):
            sayim[f"fm_{c}"] = sum((fm[2 * ii + x, 2 * jj + y] == c).astype(np.int16) for x in (0, 1) for y in (0, 1))
        for c in (-1, 0, 1, 2, 3):
            sayim[f"cm_{c}".replace("-", "y")] = (icm[ii, jj] == c).astype(np.int16)
        zen = np.where(sz[ii, jj] < -30000, np.nan, sz[ii, jj] * 0.01)
        hucre = pd.DataFrame({"anahtar": anahtar, "n_m": 1, **sayim, "zenit": zen})
        hucre = hucre.groupby("anahtar").agg({**{k: "sum" for k in ["n_m", *sayim]}, "zenit": "median"}).reset_index()
    fpd = pd.DataFrame()
    if len(fp["FP_line"]):
        fpd = pd.DataFrame({"anahtar": np.floor(fp["FP_longitude"] / AL.HUCRE).astype(np.int64) * 100000
                            + np.floor(fp["FP_latitude"] / AL.HUCRE).astype(np.int64),
                            "frp": fp["FP_power"], "vza": fp["FP_ViewZenAng"], "guven": fp["FP_confidence"],
                            "bowtie": [(int(qa_satir[int(r)][int(c)]) >> 22) & 1 for r, c in zip(fp["FP_line"], fp["FP_sample"])]})
        fpd = fpd[fpd["anahtar"].isin(anahtarlar)]
    return hucre, fpd


def ozet(R: pd.DataFrame, P: pd.DataFrame, G: pd.DataFrame, H: pd.DataFrame, S: pd.DataFrame, D: pd.DataFrame) -> dict:
    n_i = 4 * R["n_m"]
    R = R.assign(gecerli=R[[f"fm_{c}" for c in GECERLI]].sum(axis=1) / n_i, bulut=R["fm_4"] / n_i,
                 eksik=R[[f"fm_{c}" for c in EKSIK]].sum(axis=1) / n_i, cm_acik=(R["cm_2"] + R["cm_3"]) / R["n_m"])
    R["durum"] = np.select([R["gecerli"] >= 0.5, R["bulut"] >= 0.5], ["geçerli", "bulutlu"], "belirsiz")
    # hücre-gece: en az bir geçişte geçerli > hepsi bulutlu > diğer; hiç satır yoksa kapsanmadı
    hg = R.groupby(["tarih", "anahtar"])["durum"].agg(lambda x: "geçerli" if (x == "geçerli").any() else
                                                      ("bulutlu" if (x == "bulutlu").all() else "belirsiz")).reset_index()
    tum = G[["tarih"]].merge(H[["anahtar", "bolge"]], how="cross")
    hg = tum.merge(hg, on=["tarih", "anahtar"], how="left").fillna({"durum": "kapsanmadı"})
    bg = hg.groupby(["tarih", "bolge"])["durum"].value_counts(normalize=True).unstack(fill_value=0.0).reset_index()
    for c in ("geçerli", "bulutlu", "belirsiz", "kapsanmadı"):
        if c not in bg:
            bg[c] = 0.0
    firms = S.stack(future_stack=True).rename("firms_frp").reset_index().rename(columns={"level_0": "tarih", "level_1": "bolge"})
    firms.columns = ["tarih", "bolge", "firms_frp"]
    bg = bg.merge(firms, on=["tarih", "bolge"], how="left").merge(G, on="tarih")
    bg["donem"] = [AL.donem(t) or "geçiş" for t in bg["tarih"]]
    out = {"geceler": int(len(G)), "granul_cifti": int(R.groupby(["tarih", "gecis"]).ngroups), "hucre": {b: int((H["bolge"] == b).sum()) for b in BOLGELER},
           "bolge_donem": {}, "sifir_geceler": {}, "bowtie": {}, "firms_l2": {}}
    for (b, d), x in bg.groupby(["bolge", "donem"]):
        out["bolge_donem"].setdefault(b, {})[d] = {"gece": int(len(x)), **{c: round(float(x[c].mean()), 3) for c in ("geçerli", "bulutlu", "belirsiz", "kapsanmadı")}}
    for b in ("Katar", "Bahreyn"):
        z = bg[(bg["bolge"] == b) & (bg["firms_frp"] == 0)]
        out["sifir_geceler"][b] = {"gece": int(len(z)), "gecerli_pay_medyan": round(float(z["geçerli"].median()), 3) if len(z) else None,
                                   "gecerli_en_az_yuzde80": int((z["geçerli"] >= 0.8).sum()), "bulutlu_en_az_yuzde50": int((z["bulutlu"] >= 0.5).sum()),
                                   "geceler": [{"tarih": str(t.date()), **{c: round(float(v), 2) for c, v in r.items()}}
                                               for t, r in z.set_index("tarih")[["geçerli", "bulutlu", "belirsiz", "kapsanmadı"]].iterrows()]}
    if len(P):
        P = P.merge(H[["anahtar", "bolge"]], on="anahtar")
        out["bowtie"] = {"piksel": int(len(P)), "bowtie_piksel_payi": round(float(P["bowtie"].mean()), 4),
                         "bowtie_frp_payi": round(float(P.loc[P["bowtie"] == 1, "frp"].sum() / P["frp"].sum()), 4),
                         "bolge": {b: round(float(x["bowtie"].mean()), 4) for b, x in P.groupby("bolge")}}
        # FIRMS (NRT) ile L2 standart ürün: aynı gecelerde bölge başına toplam FRP (bütün geçişler, bowtie dahil) karşılaştırması
        l2 = P.groupby(["tarih", "bolge"])["frp"].sum().rename("l2_frp_toplam").reset_index()
        out["firms_l2"]["not"] = "FIRMS tarafı ana ölçü (geçiş içi toplam, en küçük track); L2 tarafı bütün geçişlerin toplamı: düzey değil, ilişki karşılaştırılır"
        k = bg.merge(l2, on=["tarih", "bolge"], how="left").fillna({"l2_frp_toplam": 0.0})
        out["firms_l2"]["bolge_korelasyon"] = {b: round(float(x[["firms_frp", "l2_frp_toplam"]].corr().iloc[0, 1]), 3)
                                               for b, x in k.groupby("bolge") if x["firms_frp"].notna().sum() > 5}
    return out, bg


def ozet_ek(bg: pd.DataFrame, P: pd.DataFrame, H: pd.DataFrame, G: pd.DataFrame) -> dict:
    """Rastgele gecelerde dönem karşılaştırması, açık gecelerde FIRMS FRP, bowtie payının dönemlere göre dağılımı ve geçiş kuralının
    gerçek bakış açısıyla sınanması (yangın piksellerinin FP_ViewZenAng değeri)."""
    out = {}
    r = bg[bg["neden"].str.startswith("rastgele")]
    out["rastgele_geceler"] = {b: {d: {"gece": int(len(y)), "gecerli": round(float(y["geçerli"].mean()), 3), "bulutlu": round(float(y["bulutlu"].mean()), 3)}
                                   for d, y in x.groupby("donem")} for b, x in r.groupby("bolge")}
    ac = bg[bg["geçerli"] >= 0.8]
    out["acik_geceler"] = {b: {d: {"gece": int(len(y)), "firms_frp_ort": round(float(y["firms_frp"].mean()), 1),
                                   "firms_frp_medyan": round(float(y["firms_frp"].median()), 1)} for d, y in x.groupby("donem")}
                           for b, x in ac.groupby("bolge")}
    if len(P):
        P = P.merge(H[["anahtar", "bolge"]], on="anahtar")
        P["donem"] = [AL.donem(t) or "geçiş" for t in P["tarih"]]
        pay = lambda x: round(float(x.loc[x["bowtie"] == 1, "frp"].sum() / x["frp"].sum()), 3)
        out["bowtie_frp_payi_donem"] = {d: pay(x) for d, x in P.groupby("donem")}
        out["bowtie_frp_payi_bolge_donem"] = {b: {d: pay(y) for d, y in x.groupby("donem")} for b, x in P.groupby("bolge")}
        # geçiş kuralı: FIRMS'te track ve scan kurallarının seçtiği geçiş, L2'de bakış açısı en küçük geçişle aynı mı (10 dakika içinde)
        P["dk"] = P["gecis"].str[-4:].astype(int).map(lambda x: (x // 100) * 60 + x % 100)
        L = P.groupby(["tarih", "anahtar", "dk"])["vza"].median().reset_index()
        L = L[L.groupby(["tarih", "anahtar"])["dk"].transform(lambda x: x.max() - x.min()) > 20]
        sv = L.sort_values("vza").groupby(["tarih", "anahtar"]).first().reset_index().rename(columns={"dk": "dk_vza"})
        D = AL.yukle(["VIIRS_NOAA21_NRT"])
        D = D[D["tarih"].isin(G["tarih"])]
        D["anahtar"] = D["hx"].astype(np.int64) * 100000 + D["hy"].astype(np.int64)
        F = D[D["anahtar"].isin(H["anahtar"])].groupby(["tarih", "anahtar", "gecis"]).agg(track=("track", "mean"), scan=("scan", "mean"),
                                                                                         dk=("dk", "min")).reset_index()
        F = F[F.groupby(["tarih", "anahtar"])["gecis"].transform("size") > 1]
        st = F.sort_values(["track", "scan", "gecis"]).groupby(["tarih", "anahtar"])["dk"].first().rename("dk_track")
        ss = F.sort_values(["scan", "gecis"]).groupby(["tarih", "anahtar"])["dk"].first().rename("dk_scan")
        K = sv.merge(st.reset_index(), on=["tarih", "anahtar"]).merge(ss.reset_index(), on=["tarih", "anahtar"])
        out["gecis_kurali_bakis_acisi"] = {"hucre_gece": int(len(K)),
                                           "track_en_kucuk_acili_gecisi_secti": round(float(((K["dk_track"] - K["dk_vza"]).abs() <= 10).mean()), 4),
                                           "scan_en_kucuk_acili_gecisi_secti": round(float(((K["dk_scan"] - K["dk_vza"]).abs() <= 10).mean()), 4)}
    return out


def main():
    import sys
    D = AL.yukle(["VIIRS_NOAA21_NRT"])
    H = pd.read_csv(A.VERI / "firms" / "kalici_hucreler.csv")
    S = AL.gece_serisi(D, H)[BOLGELER]
    H = H[H["bolge"].isin(BOLGELER)].copy()
    if "--yalniz-ozet" in sys.argv:  # indirme yapmadan kayıtlı pilot dosyalarından özet
        H["anahtar"] = H["hx"].astype(np.int64) * 100000 + H["hy"].astype(np.int64)
        R = pd.read_csv(CIKTI / "hucre_gecis.csv.gz", parse_dates=["tarih"])
        P = pd.read_csv(CIKTI / "yangin_pikselleri.csv.gz", parse_dates=["tarih"])
        G = pd.read_csv(CIKTI / "geceler.csv", parse_dates=["tarih"])
        onceki = json.loads((A.SONUC / "a16_kapsama_pilot.json").read_text(encoding="utf-8"))
        out, bg = ozet(R, P, G, H, S, D)
        out["sorunlar"] = onceki.get("sorunlar", [])
        out.update(ozet_ek(bg, P, H, G))
        (A.SONUC / "a16_kapsama_pilot.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
        print(json.dumps({k: out[k] for k in ("rastgele_geceler", "acik_geceler", "bowtie_frp_payi_donem", "gecis_kurali_bakis_acisi")},
                         ensure_ascii=False, indent=1))
        return
    H["anahtar"] = H["hx"].astype(np.int64) * 100000 + H["hy"].astype(np.int64)
    anahtarlar = H["anahtar"].to_numpy()
    global KUTU  # açık denizdeki hücreler dahil (Katar'ın bazı hücreleri 52,5 D'de)
    KUTU = (round(H["lon"].min() - 0.1, 2), round(H["lat"].min() - 0.1, 2), round(H["lon"].max() + 0.1, 2), round(H["lat"].max() + 0.1, 2))
    print("kutu:", KUTU)
    G = geceler_sec(S)
    print(f"pilot: {len(G)} gece;", G["neden"].value_counts().to_dict())
    CIKTI.mkdir(parents=True, exist_ok=True)
    satir, piksel, kayit = [], [], []
    with tempfile.TemporaryDirectory(dir=os.environ.get("KAPSAMA_GECICI")) as gecici, ThreadPoolExecutor(4) as havuz:
        gk = Path(gecici)
        isler = []
        for t in G["tarih"]:
            gun = str(t.date())
            ga, gb = granuller(KOL["ates"], gun), granuller(KOL["bulut"], gun)
            for k in sorted(set(ga) & set(gb)):
                isler.append((t, k, havuz.submit(indir, ga[k], gk), havuz.submit(indir, gb[k], gk)))
            for k in sorted(set(ga) ^ set(gb)):
                kayit.append({"tarih": gun, "gecis": k, "sorun": "eşi yok"})
        print(f"{len(isler)} granül çifti kuyrukta")
        t0 = time.time()
        for n, (t, k, fa, fb) in enumerate(isler, 1):
            try:
                pa, pb = fa.result(), fb.result()
                hc, fp = isle(pa, pb, anahtarlar)
                pa.unlink(missing_ok=True)
                pb.unlink(missing_ok=True)
            except Exception as e:  # noqa: BLE001
                kayit.append({"tarih": str(t.date()), "gecis": k, "sorun": repr(e)[:120]})
                continue
            if len(hc):
                satir.append(hc.assign(tarih=t, gecis=k))
            if len(fp):
                piksel.append(fp.assign(tarih=t, gecis=k))
            if n % 20 == 0:
                print(f"  {n}/{len(isler)} çift | {time.time() - t0:.0f} sn")
    R = pd.concat(satir, ignore_index=True)
    P = pd.concat(piksel, ignore_index=True) if piksel else pd.DataFrame()
    R.to_csv(CIKTI / "hucre_gecis.csv.gz", index=False)
    if len(P):
        P.to_csv(CIKTI / "yangin_pikselleri.csv.gz", index=False)
    G.to_csv(CIKTI / "geceler.csv", index=False)
    out, bg = ozet(R, P, G, H, S, D)
    bg.to_csv(CIKTI / "bolge_gece.csv", index=False)
    out["sorunlar"] = kayit
    out.update(ozet_ek(bg, P, H, G))
    (A.SONUC / "a16_kapsama_pilot.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k not in ("sifir_geceler",)}, ensure_ascii=False, indent=1)[:4000])
    print("sıfır geceler:", {b: {k: v for k, v in x.items() if k != "geceler"} for b, x in out["sifir_geceler"].items()})


if __name__ == "__main__":
    main()
