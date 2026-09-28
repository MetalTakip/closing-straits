"""Analiz 15: Körfez'de gaz alevleri (NASA FIRMS, VIIRS 375 m gece tespitleri).

Soru: Hürmüz kapanınca ihracatı boğaza bağlı sahalarda tespit edilen alev etkinliği azaldı mı? İhracatı boğazı atlayan bölgeler
(Umman, ihracatı Ceyhan'a boru hattıyla giden kuzey Irak) karşılaştırma için.
Yöntem (sürüm 4, 26 Eylül 2026 üçüncü çapraz denetiminden sonra; sürüm 3'te geçiş en küçük ortalama "scan" ile seçiliyordu):
  1. Gece tespitleri (daynight = N). Tarih, kaydın UTC tarihidir; yerel güneş gecesine göre yeniden gruplama sağlamlıkta. Ana seri
     NOAA-21 NRT, 1 Eylül 2024'ten (firms_cek.py bu tarihten çeker); ikinci uydu Suomi NPP (Haziran 2026'ya kadar SP, sonra NRT; 9-12 Mart ve
     28 Nisan-2 Haziran 2026 arasında hiç gece tespiti yok).
  2. Geçiş: bir gecenin tespit zamanları 20 dakikadan uzun boşlukla ayrılıyorsa ayrı uydu geçişleridir. Aynı dakikadaki birden çok kayıt
     aynı geçişin farklı pikselleridir (NASA kılavuzu: 750 m FRP eşleşen 375 m pikseller arasında paylaştırılabilir).
  3. Kalıcı alev hücreleri: 0,01 derecelik ızgarada savaş öncesi gecelerin (1 Eylül 2024-27 Şubat 2026) en az %5'inde tespit görülen
     hücreler (NOAA-21'den seçilir, iki uyduda da aynı liste). Bu ölçüt hücrenin gaz alevi olduğunu tek başına kanıtlamaz; Suomi NPP SP
     kayıtlarının tür alanı yardımcı kontrol olarak raporlanır.
  4. Bölge: Natural Earth ülke sınırları (denizdeki platformlar en yakın kıyı ülkesine, 1 dereceye kadar); çalışma kutusu (43,5-60 D,
     17-37,5 K) içindeki kısım. Irak enlem 32,5'in güneyi ve 34'ün kuzeyi olarak ikiye; İran'da batı kesim (boylam < 51, enlem >= 29; kalıcı
     hücrelerin %56'sı ve savaş öncesi alev gücünün %72'si Huzistan eyaletinde, Natural Earth eyalet sınırlarıyla, 27 Eylül 2026) ve
     Asaluyeh-Güney Pars (51,5-53 D, 27-28,2 K) ayrı.
  5. Ölçü ("tespit edilen FRP", MW): her hücre-gece için önce geçiş içindeki piksellerin FRP'si toplanır; birden çok geçiş varsa
     tespit edilen piksellerin ortalama "track" boyu (tarama yönüne dik piksel boyu) en küçük olan geçiş alınır; sonra bölgedeki hücreler
     toplanır. "track" bakış açısıyla düzenli artar; "scan" (tarama yönündeki boy) VIIRS'in piksel birleştirme bölgeleri yüzünden
     testere dişi çizer, bu yüzden en küçük "scan" nadire en yakın geçişi garanti etmez (üçüncü denetim, Z01). Kural geçişleri FRP'ye göre
     sıralamaz, ama yalnız tespit içeren geçişler arasından seçer (tespit olmayan geçiş FIRMS dosyasında yok).
     Kutunun hiçbir yerinde tespit olmayan geceler eksik sayılır. Bölgede tespit olmayan gece sıfırdır: bu, geçerli bir gözlemde alevin
     söndüğü anlamına gelmez; FIRMS noktaları bulut, yerel kapsama ve gözlem kalitesi bilgisi taşımaz (FRP'nin tarama örtüşmesi QA biti de
     FIRMS dosyalarında yok).
  6. Ana model: PPML (Poisson sözde en çok olabilirlik), FRP düzeyi = dönem + takvim ayı kuklaları; HAC standart hata (Bartlett, 14
     gecikme, gözlenen gece sırasıyla, küçük örneklem düzeltmesi yok). Sağlamlık: geçiş birleştirme kuralı (en küçük ortalama scan = sürüm 3,
     en küçük ortalama piksel alanı, en yüksek geçiş, ortalama geçiş, piksel en yükseği, bütün kayıtların toplamı), gece başına geçiş sayısı
     kuklaları, yerel güneş gecesi, log(1 + FRP), HAC 7 ve
     28, Suomi NPP ve iki uydunun ortak günleri, altı eşik/hücre kurgusu, güney-kuzey farkı (bölge x dönem PPML, iki uyduda), Holm
     (tablodaki 10 bölge x 2 kapanma; açılış dönemi, Irak'ın orta kesimi, "Diğer", İran'ın geri kalanı ve güney-kuzey farkı aile dışında).
  7. Karşılaştırmalar: 2026 aylarının 2025'in aynı günlerine göre değişimi (aylık medyan; eksik ayda aynı gün penceresi); SOMO'nun
     çıkış noktası bazında aylık ihracatı (somo_cek.py -> veri/rafineri/irak_somo_cikis_aylik.csv; betimsel korelasyon, p değeri yok).
Sınırlar: bulut ve yerel gözlenebilirlik ayrılamıyor; alev üretimle doğrusal değildir; FRP güçtür, yakılan gaz hacmi değildir.
Çıktı: sonuclar/a15_alev.json, veri/firms/kalici_hucreler.csv
"""
from __future__ import annotations

import calendar
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from shapely import STRtree, points
from shapely.geometry import shape
from statsmodels.stats.multitest import multipletests

import analiz_ortak as A

F = A.VERI / "firms"
ONCE_SON = pd.Timestamp("2026-02-27")
DON = {"1. kapanma": ("2026-03-02", "2026-06-13"), "Açılış": ("2026-06-18", "2026-07-10"), "2. kapanma": ("2026-07-13", "2026-09-22")}
ULKE = {"Iraq": "Irak", "Kuwait": "Kuveyt", "Iran": "İran", "Saudi Arabia": "Suudi Arabistan", "Qatar": "Katar",
        "United Arab Emirates": "BAE", "Bahrain": "Bahreyn", "Oman": "Umman"}
TABLO = ["Irak (güney)", "Kuveyt", "Katar", "Bahreyn", "Batı İran (çoğu Huzistan)", "İran (Asaluyeh)", "Suudi Arabistan", "BAE", "Umman", "Irak (kuzey)"]
HUCRE = 0.01
KALICI_PAY = 0.05
GECIS_ARALIK = 20  # dakika
KURAL = "track"   # ana ölçünün geçiş kuralı
KURALLAR = ("track", "scan", "alan", "en_yuksek", "ortalama", "piksel_en_yuksek", "toplam")
K = ("1. kapanma", "2. kapanma")


def gecis_ata(D: pd.DataFrame) -> pd.DataFrame:
    """Her gecenin tespit zamanlarını GECIS_ARALIK dakikadan uzun boşluklarda böler: 0, 1, 2... geçiş kimliği."""
    parca = []
    for t, grp in D.groupby("tarih"):
        z = np.sort(grp["dk"].unique())
        parca.append(pd.DataFrame({"tarih": t, "dk": z, "gecis": np.concatenate([[0], np.cumsum(np.diff(z) > GECIS_ARALIK)])}))
    return D.drop(columns=[c for c in ("gecis",) if c in D]).merge(pd.concat(parca), on=["tarih", "dk"], how="left")


def yukle(kaynaklar: list[str]) -> pd.DataFrame:
    D = pd.concat([pd.read_csv(F / f"{k}.csv.gz") for k in kaynaklar], ignore_index=True)
    D = D[D["daynight"] == "N"].copy()
    D["tarih"] = pd.to_datetime(D["acq_date"])  # kaydın UTC tarihi
    D["dk"] = (D["acq_time"] // 100) * 60 + D["acq_time"] % 100
    D = gecis_ata(D)
    D["hx"] = np.floor(D["longitude"] / HUCRE).astype(int)
    D["hy"] = np.floor(D["latitude"] / HUCRE).astype(int)
    return D


def yerel_gece(D: pd.DataFrame) -> pd.DataFrame:
    """Kayıtları yerel güneş gecesine göre tarihler (UTC + boylam/15 saat, 12 saat geri: gecenin başladığı gün); geçişleri yeniden atar."""
    D = D.copy()
    zaman = pd.to_datetime(D["acq_date"]) + pd.to_timedelta(D["dk"], unit="m") + pd.to_timedelta(D["longitude"] / 15, unit="h")
    D["tarih"] = (zaman - pd.Timedelta(hours=12)).dt.normalize()
    return gecis_ata(D)


def kalici_hucreler(D: pd.DataFrame) -> pd.DataFrame:
    O = D[D["tarih"] <= ONCE_SON]
    gece = O["tarih"].nunique()
    H = O.groupby(["hx", "hy"])["tarih"].nunique().rename("gece").reset_index()
    H = H[H["gece"] >= KALICI_PAY * gece].copy()
    H["lon"], H["lat"] = (H["hx"] + 0.5) * HUCRE, (H["hy"] + 0.5) * HUCRE
    return H


def bolge_ata(H: pd.DataFrame) -> pd.DataFrame:
    gj = json.loads((A.VERI / "harita" / "ne_10m_admin_0_countries.geojson").read_text(encoding="utf-8"))
    geo, ad = [], []
    for f in gj["features"]:
        p = f["properties"]
        n = p.get("ADMIN") or p.get("NAME")
        if n in ULKE:
            geo.append(shape(f["geometry"]))
            ad.append(ULKE[n])
    agac = STRtree(geo)
    nok = points(H["lon"].to_numpy(), H["lat"].to_numpy())
    ulke = []
    for pt in nok:
        ic = agac.query(pt, predicate="within")
        if len(ic):
            ulke.append(ad[ic[0]])
            continue
        j = agac.nearest(pt)
        ulke.append(ad[j] if geo[j].distance(pt) <= 1.0 else "Diğer")
    H = H.copy()
    H["ulke"] = ulke
    b = H["ulke"].copy()
    irak = H["ulke"] == "Irak"
    b[irak & (H["lat"] < 32.5)] = "Irak (güney)"
    b[irak & (H["lat"] >= 34.0)] = "Irak (kuzey)"
    b[irak & (H["lat"] >= 32.5) & (H["lat"] < 34.0)] = "Irak (orta)"
    iran = H["ulke"] == "İran"
    b[iran & (H["lon"] < 51.0) & (H["lat"] >= 29.0)] = "Batı İran (çoğu Huzistan)"
    b[iran & H["lon"].between(51.5, 53.0) & H["lat"].between(27.0, 28.2)] = "İran (Asaluyeh)"
    b[iran & ~b.isin(["Batı İran (çoğu Huzistan)", "İran (Asaluyeh)"])] = "İran (diğer)"
    H["bolge"] = b
    return H


def hucre_gece(D: pd.DataFrame, H: pd.DataFrame, kural: str = KURAL) -> pd.DataFrame:
    """Hücre-gece başına tespit edilen FRP. Geçiş kuralları: "track" (ana; geçiş içi toplam, ortalama track boyu en küçük geçiş; eşitlikte
    küçük scan, sonra önceki geçiş), "scan" (sürüm 3: ortalama scan boyu en küçük geçiş), "alan" (ortalama piksel alanı scan x track en
    küçük geçiş), "en_yuksek" ve "ortalama" (geçiş içi toplamların gecelik en yükseği ve ortalaması), "piksel_en_yuksek" (sürüm 2: tek
    pikselin en yüksek değeri), "toplam" (sürüm 1: bütün kayıtların toplamı)."""
    Kk = D.merge(H[["hx", "hy", "bolge"]], on=["hx", "hy"], how="inner")
    g = ["tarih", "hx", "hy", "bolge"]
    if kural == "toplam":
        return Kk.groupby(g, as_index=False)["frp"].sum()
    if kural == "piksel_en_yuksek":
        return Kk.groupby(g, as_index=False)["frp"].max()
    Kk = Kk.assign(alan=Kk["scan"] * Kk["track"])
    cg = Kk.groupby(g + ["gecis"], as_index=False).agg(frp=("frp", "sum"), scan=("scan", "mean"), track=("track", "mean"), alan=("alan", "mean"))
    sira = {"track": ["track", "scan", "gecis"], "scan": ["scan", "gecis"], "alan": ["alan", "gecis"]}
    if kural in sira:
        return cg.sort_values(sira[kural]).groupby(g, as_index=False).first()[g + ["frp"]]
    return cg.groupby(g, as_index=False)["frp"].agg("max" if kural == "en_yuksek" else "mean")


def gece_serisi(D: pd.DataFrame, H: pd.DataFrame, kural: str = KURAL) -> pd.DataFrame:
    S = hucre_gece(D, H, kural).groupby(["tarih", "bolge"])["frp"].sum().unstack(fill_value=0.0)
    geceler = D.groupby("tarih").size()
    tum = pd.date_range(D["tarih"].min(), D["tarih"].max(), freq="D")
    S = S.reindex(tum, fill_value=0.0)
    S[~S.index.isin(geceler.index)] = np.nan  # hiç tespit yoksa veri yok (uydu ya da işleme boşluğu)
    return S


def donem(t: pd.Timestamp) -> str | None:
    if t <= ONCE_SON:
        return "Savaş öncesi"
    for d, (a, b) in DON.items():
        if pd.Timestamp(a) <= t <= pd.Timestamp(b):
            return d
    return None


def tasarim(y: pd.Series, ek: pd.DataFrame | None = None):
    y = y.dropna()
    d = pd.Series([donem(t) for t in y.index], index=y.index)
    y, d = y[d.notna()], d[d.notna()]
    X = pd.get_dummies(d).drop(columns="Savaş öncesi").astype(float)
    X = X.join(pd.get_dummies(y.index.month, prefix="ay", drop_first=True).set_index(y.index).astype(float))
    if ek is not None:
        X = X.join(ek.reindex(y.index).fillna(0.0))
    return y, d, sm.add_constant(X)


def etki(m, k: str) -> dict:
    b, s = m.params[k], m.bse[k]
    return {"etki_yuzde": round(100 * (np.exp(b) - 1), 3), "alt": round(100 * (np.exp(b - 1.96 * s) - 1), 3),
            "ust": round(100 * (np.exp(b + 1.96 * s) - 1), 3), "p": float(m.pvalues[k])}


def regresyon(y: pd.Series, yontem: str = "ppml", gecikme: int = 14, ek: pd.DataFrame | None = None) -> dict:
    """yontem "ppml" (ana): FRP düzeyi için Poisson sözde en çok olabilirlik; etki ortalama tespit edilen FRP'deki yüzde değişim.
    "log": log(1 + FRP) için OLS (dönüştürülmüş ölçek; MW biriminde). İkisinde de HAC (gecikme gözlenen gece sayısıdır).
    ek: ek kontrol değişkenleri (ör. gece başına geçiş sayısı kuklaları)."""
    y, d, X = tasarim(y, ek)
    if yontem == "ppml":
        m = sm.GLM(y, X, family=sm.families.Poisson()).fit(cov_type="HAC", cov_kwds={"maxlags": gecikme})
    else:
        m = sm.OLS(np.log1p(y), X).fit(cov_type="HAC", cov_kwds={"maxlags": gecikme})
    out = {k: etki(m, k) for k in DON if k in m.params}
    out["n"] = int(m.nobs)
    out["n_donem"] = {k: int((d == k).sum()) for k in ("Savaş öncesi", *DON)}
    return out


def fark_testi(S: pd.DataFrame, a: str = "Irak (güney)", b: str = "Irak (kuzey)") -> dict:
    """a bölgesinin b'ye göre göreli değişimi: 100 x [exp(beta_a - beta_b) - 1]. İki bölgenin gece serileri alt alta, bölge x dönem ve
    bölge x ay etkileşimli PPML; standart hata aynı gecenin iki gözlemini birlikte ele alan HAC (hac-groupsum, 14 gece)."""
    Y = S[[a, b]].dropna()
    d = pd.Series([donem(t) for t in Y.index], index=Y.index)
    Y, d = Y[d.notna()], d[d.notna()]
    parca = []
    for bolge, g in ((a, 1.0), (b, 0.0)):
        X = pd.get_dummies(d).drop(columns="Savaş öncesi").astype(float)
        X = X.join(pd.get_dummies(Y.index.month, prefix="ay", drop_first=True).set_index(Y.index).astype(float))
        for c in list(X.columns):
            X[f"{c}_a"] = X[c] * g
        X["bolge_a"] = g
        X["y"] = Y[bolge].to_numpy()
        X["zaman"] = np.arange(len(Y))
        parca.append(X)
    U = pd.concat(parca).sort_values(["zaman", "bolge_a"])
    y, zaman = U.pop("y"), U.pop("zaman").to_numpy()
    m = sm.GLM(y, sm.add_constant(U), family=sm.families.Poisson()).fit(cov_type="hac-groupsum",
                                                                      cov_kwds={"time": zaman, "maxlags": 14})
    return {k: etki(m, f"{k}_a") for k in K}


def yillik(S: pd.DataFrame, b: str) -> dict:
    """2026 aylarının 2025'in aynı günlerine göre değişimi (aylık medyan). Eksik ayda (ör. Eylül 2026) iki yılda aynı gün penceresi."""
    son = S[b].dropna().index.max()
    out = {}
    for ay in range(3, 10):
        bit = min(calendar.monthrange(2026, ay)[1], son.day if (son.year, son.month) == (2026, ay) else 31)
        a26 = S.loc[f"2026-{ay:02d}-01":f"2026-{ay:02d}-{bit:02d}", b].median()
        a25 = S.loc[f"2025-{ay:02d}-01":f"2025-{ay:02d}-{min(bit, calendar.monthrange(2025, ay)[1]):02d}", b].median()
        out[f"2026-{ay:02d}"] = round(100 * (a26 / a25 - 1), 3) if a25 and a25 > 0 else None
    return out


def harita_hucreleri(D: pd.DataFrame, H: pd.DataFrame, adim: float = 0.1) -> list[dict]:
    """Harita için 0,1 derecelik hücreler: gece başına ortalama tespit edilen FRP (ana ölçü); kapanma dönemi, 2025'in aynı günleriyle
    karşılaştırılır (mevsim ayrılır)."""
    Kk = hucre_gece(D, H).merge(H[["hx", "hy", "lon", "lat"]], on=["hx", "hy"])
    Kk["gx"], Kk["gy"] = np.floor(Kk["lon"] / adim).astype(int), np.floor(Kk["lat"] / adim).astype(int)
    geceler = pd.Series(sorted(D["tarih"].unique()))
    pencere = {"once": ("2024-09-01", "2026-02-27"), "k1_2026": ("2026-03-02", "2026-06-13"), "k1_2025": ("2025-03-02", "2025-06-13"),
               "k2_2026": ("2026-07-13", "2026-09-22"), "k2_2025": ("2025-07-13", "2025-09-22")}
    ort = {}
    for ad, (a, b) in pencere.items():
        n = int(((geceler >= a) & (geceler <= b)).sum())
        W = Kk[(Kk["tarih"] >= a) & (Kk["tarih"] <= b)]
        ort[ad] = W.groupby(["gx", "gy"])["frp"].sum() / max(n, 1)
    T = pd.DataFrame(ort).fillna(0.0)
    bolge = Kk.groupby(["gx", "gy"])["bolge"].agg(lambda x: x.mode().iloc[0])
    T = T.join(bolge)
    T = T[T["once"] >= 2.0]  # gece başına en az 2 MW: gürültüyü ayıklar
    out = []
    for (gx, gy), r in T.iterrows():
        d1 = 100 * (r["k1_2026"] / r["k1_2025"] - 1) if r["k1_2025"] > 0.5 else None
        d2 = 100 * (r["k2_2026"] / r["k2_2025"] - 1) if r["k2_2025"] > 0.5 else None
        out.append({"lon": round((gx + 0.5) * adim, 3), "lat": round((gy + 0.5) * adim, 3), "bolge": r["bolge"],
                    "once_mw": round(float(r["once"]), 1), "d1": None if d1 is None else round(d1, 1), "d2": None if d2 is None else round(d2, 1)})
    return out


def saglamlik(ana: pd.DataFrame) -> dict:
    """Kalıcı nokta eşiği (%2, %5, %10) ve hücre boyu (0,01 ve 0,02 derece): altı kurguda ilk ve ikinci kapanma etkileri (ana model)."""
    global HUCRE, KALICI_PAY
    eski = (HUCRE, KALICI_PAY)
    sonuc = []
    for hucre in (0.01, 0.02):
        for pay in (0.02, 0.05, 0.10):
            HUCRE, KALICI_PAY = hucre, pay
            D = ana.copy()
            D["hx"], D["hy"] = np.floor(D["longitude"] / hucre).astype(int), np.floor(D["latitude"] / hucre).astype(int)
            H = bolge_ata(kalici_hucreler(D))
            S = gece_serisi(D, H)
            for b in ("Irak (güney)", "Irak (kuzey)", "Katar"):
                r = regresyon(S[b])
                sonuc.append({"hucre": hucre, "esik": pay, "bolge": b, "k1": r["1. kapanma"]["etki_yuzde"], "k2": r["2. kapanma"]["etki_yuzde"]})
    HUCRE, KALICI_PAY = eski
    R = pd.DataFrame(sonuc)
    return {"kurgu": 6, "tum": sonuc,
            **{b: {"k1_aralik": [float(g["k1"].min()), float(g["k1"].max())], "k2_aralik": [float(g["k2"].min()), float(g["k2"].max())]}
               for b, g in R.groupby("bolge")}}


def eksik_bloklar(S: pd.DataFrame, bas: str = "2024-09-01", bit: str = "2026-09-24") -> list[list[str]]:
    g = S.dropna(how="all").index
    eksik = [t for t in pd.date_range(bas, bit) if t not in g]
    bloklar = []
    for t in eksik:
        if bloklar and (t - pd.Timestamp(bloklar[-1][1])).days == 1:
            bloklar[-1][1] = str(t.date())
        else:
            bloklar.append([str(t.date()), str(t.date())])
    return bloklar


def cokluk(D: pd.DataFrame, H: pd.DataFrame) -> dict:
    """Kalıcı hücrelerde hücre-gece başına kayıt çokluğu: aynı geçişte birden çok piksel mi, birden çok geçiş mi? Döneme ve bölgeye göre."""
    Kk = D.merge(H[["hx", "hy", "bolge"]], on=["hx", "hy"], how="inner")
    g = Kk.groupby(["tarih", "hx", "hy", "bolge"]).agg(kayit=("frp", "size"), dakika=("dk", "nunique"), gecis=("gecis", "nunique"),
                                                         acilim=("dk", lambda x: x.max() - x.min())).reset_index()
    g["donem"] = [donem(t) or "geçiş" for t in g["tarih"]]
    coklu = g["kayit"] > 1
    ozet = {"hucre_gece": int(len(g)), "coklu_kayit_payi": round(float(coklu.mean()), 4),
            "coklu_icinde_ayni_dakika_payi": round(float(((g["dakika"] == 1) & coklu).sum() / coklu.sum()), 4),
            "coklu_gecis_payi": round(float((g["gecis"] > 1).mean()), 4), "acilim_60dk_ustu_payi": round(float((g["acilim"] > 60).mean()), 4)}
    tablo = {}
    for (b, d), x in g[g["bolge"].isin(TABLO)].groupby(["bolge", "donem"]):
        tablo.setdefault(b, {})[d] = {"ayni_geciste_coklu_piksel": round(float(((x["kayit"] > 1) & (x["gecis"] == 1)).mean()), 3),
                                      "coklu_gecis": round(float((x["gecis"] > 1).mean()), 3), "hucre_gece": int(len(x))}
    gs = D.groupby("tarih")["gecis"].nunique()
    dd = pd.Series([donem(t) or "geçiş" for t in gs.index], index=gs.index)
    ozet["gece_basina_gecis_ort"] = {k: round(float(v), 3) for k, v in gs.groupby(dd).mean().items()}
    return {"ozet": ozet, "bolge_donem": tablo}


def geometri(D: pd.DataFrame, H: pd.DataFrame) -> dict:
    """Piksel geometrisi: track (tarama yönüne dik boy) bakış açısıyla düzenli artar, scan (tarama yönündeki boy) VIIRS'in birleştirme
    bölgelerinde geri düşer. Kanıt: track'e göre scan medyanı (en az 500 kayıtlı 0,01 km'lik dilimler). Ayrıca çok geçişli hücre-gecelerde
    track ve scan kurallarının farklı geçiş seçtiği pay."""
    t = D.groupby(D["track"].round(2))["scan"].agg(["median", "size"])
    t = t[t["size"] >= 500]
    Kk = D.merge(H[["hx", "hy", "bolge"]], on=["hx", "hy"], how="inner")
    g = ["tarih", "hx", "hy"]
    cg = Kk.groupby(g + ["gecis"], as_index=False).agg(scan=("scan", "mean"), track=("track", "mean"))
    cok = cg.groupby(g)["gecis"].transform("size") > 1
    cg = cg[cok]
    sec_t = cg.sort_values(["track", "scan", "gecis"]).groupby(g)["gecis"].first()
    sec_s = cg.sort_values(["scan", "gecis"]).groupby(g)["gecis"].first()
    return {"track_scan_medyan": {f"{k:.2f}": round(float(v), 3) for k, v in t["median"].items()},
            "track_aralik": [float(t.index.min()), float(t.index.max())],
            "cok_gecisli_hucre_gece": int(len(sec_t)), "track_scan_farkli_secim": int((sec_t != sec_s).sum()),
            "track_scan_farkli_pay": round(float((sec_t != sec_s).mean()), 4)}


def snpp_tur(H: pd.DataFrame) -> dict | None:
    """Suomi NPP standart işleme (SP) kayıtlarının tür alanı: 0 bitki yangını varsayımı, 2 başka sabit kara kaynağı, 3 denizde.
    Kalıcı hücrelerin gaz alevi olup olmadığına yardımcı kontrol (kanıt değil)."""
    f = F / "VIIRS_SNPP_SP.csv.gz"
    if not f.exists():
        return None
    S = pd.read_csv(f)
    S = S[S["daynight"] == "N"].copy()
    S["hx"], S["hy"] = np.floor(S["longitude"] / HUCRE).astype(int), np.floor(S["latitude"] / HUCRE).astype(int)
    Kk = S.merge(H[["hx", "hy", "bolge"]], on=["hx", "hy"], how="left")
    Kk["kalici"] = Kk["bolge"].notna()
    pay = lambda x: {str(int(k)): round(float(v), 3) for k, v in x["type"].value_counts(normalize=True).sort_index().items()}
    return {"kalici": pay(Kk[Kk["kalici"]]), "diger": pay(Kk[~Kk["kalici"]]),
            "bolge": {b: pay(x) for b, x in Kk[Kk["bolge"].isin(TABLO)].groupby("bolge")}}


def somo_karsilastir(aylik: pd.DataFrame) -> dict:
    """SOMO çıkış noktası bazında aylık ihracat (somo_cek.py) ile güney ve kuzey Irak'ın aylık alev medyanı. Korelasyonlar betimseldir:
    aylar değiştirilebilir değildir (zaman bağımlılığı, rejim değişimi), bu yüzden p değeri verilmez."""
    so = pd.read_csv(A.VERI / "rafineri" / "irak_somo_cikis_aylik.csv")
    P = so.pivot_table(index="ay", columns="cikis", values="miktar_varil", aggfunc="sum")
    Pk = so[so.cikis == "Ceyhan"].pivot_table(index="ay", columns="kalem", values="miktar_varil", aggfunc="sum")
    ay_ = {}
    for ay, r in P.iterrows():
        gun = calendar.monthrange(int(ay[:4]), int(ay[5:]))[1]
        t = pd.Timestamp(f"{ay}-01")
        g = aylik.loc[t, "Irak (güney)"] if t in aylik.index else np.nan
        n = aylik.loc[t, "Irak (kuzey)"] if t in aylik.index else np.nan
        mv = lambda v: None if pd.isna(v) else round(float(v) / 1e6, 3)
        ay_[ay] = {"basra_milyon_varil": mv(r.get("Basra")), "basra_bin_varil_gun": None if pd.isna(r.get("Basra")) else round(r["Basra"] / gun / 1e3, 3),
                   "ceyhan_milyon_varil": mv(r.get("Ceyhan")), "ceyhan_bin_varil_gun": None if pd.isna(r.get("Ceyhan")) else round(r["Ceyhan"] / gun / 1e3, 3),
                   "ceyhan_kirkuk_milyon_varil": mv(Pk.loc[ay].get("Kirkuk")) if ay in Pk.index else None,
                   "ceyhan_krg_milyon_varil": mv(Pk.loc[ay].get("(KRG) Kirkuk")) if ay in Pk.index else None,
                   "khor_zubair_qaiyarah_milyon_varil": mv(r.get("Khor al-Zubair")),
                   "irak_guney_frp": None if pd.isna(g) else round(float(g), 2), "irak_kuzey_frp": None if pd.isna(n) else round(float(n), 2)}
    ortak = [a for a, v in ay_.items() if v["basra_bin_varil_gun"] is not None and v["irak_guney_frp"] is not None]
    x = np.array([ay_[a]["basra_bin_varil_gun"] for a in ortak])
    y = np.array([ay_[a]["irak_guney_frp"] for a in ortak])
    alt = lambda kume: float(stats.pearsonr([ay_[a]["basra_bin_varil_gun"] for a in kume], [ay_[a]["irak_guney_frp"] for a in kume])[0])
    once = [a for a in ortak if a <= "2026-02"]
    kap = [a for a in ortak if a >= "2026-03"]
    kor = {"aylar": ortak, "n": len(ortak), "pearson": float(stats.pearsonr(x, y)[0]), "spearman": float(stats.spearmanr(x, y)[0]),
           "savas_oncesi_aylari": once, "savas_oncesi_pearson": alt(once), "kapanma_aylari": kap, "kapanma_pearson": alt(kap) if len(kap) >= 3 else None}
    return {"aylik": ay_, "korelasyon": kor, "kaynak": "SOMO aylık rapor PDF'leri (Haziran 2025-Mart 2026, Ağustos 2026) ve ihracat "
                                                      "grafiği (Nisan-Haziran 2026); somo_cek.py"}


def erken_mart(S: pd.DataFrame) -> dict:
    """Kapanmanın ilk haftaları: bölge başına gece ortalaması tespit edilen FRP'nin savaştan hemen önceki altı haftaya (15 Ocak-27 Şubat 2026)
    göre değişimi (%). Payne Institute'un VIIRS Nightfire değişim tespitiyle (Zhizhin ve Bazilian, 24 Mart 2026; olaylar 1-10 Mart)
    karşılaştırmak içindir; betimsel, model yok."""
    S = S.assign(**{"İran (toplam)": S[[c for c in S if "İran" in c]].sum(axis=1, min_count=1),
                    "Irak (toplam)": S[[c for c in S if c.startswith("Irak")]].sum(axis=1, min_count=1)})
    taban = S.loc["2026-01-15":"2026-02-27"].mean()
    pencere = {"1-10 Mart": ("2026-03-01", "2026-03-10"), "11-31 Mart": ("2026-03-11", "2026-03-31"), "Nisan": ("2026-04-01", "2026-04-30"),
               "Mayıs": ("2026-05-01", "2026-05-31"), "1-13 Haziran": ("2026-06-01", "2026-06-13")}
    return {"taban": ["2026-01-15", "2026-02-27"],
            "bolgeler": {b: {"taban_mw": round(float(taban[b]), 2),
                             **{k: round(100 * (float(S.loc[a:z, b].mean()) / float(taban[b]) - 1), 3) for k, (a, z) in pencere.items()}}
                         for b in [*TABLO, "İran (toplam)", "Irak (toplam)"]}}


def main():
    ana = yukle(["VIIRS_NOAA21_NRT"])
    H = bolge_ata(kalici_hucreler(ana))
    S = gece_serisi(ana, H)
    Sk = {k: gece_serisi(ana, H, k) for k in KURALLAR if k != KURAL}
    Sy = gece_serisi(yerel_gece(ana), H)
    S2 = None
    if all((F / f"{k}.csv.gz").exists() for k in ("VIIRS_SNPP_SP", "VIIRS_SNPP_NRT")):
        S2 = gece_serisi(yukle(["VIIRS_SNPP_SP", "VIIRS_SNPP_NRT"]), H)
    ortak = S.dropna(how="all").index.intersection(S2.dropna(how="all").index) if S2 is not None else None
    gs = ana.groupby("tarih")["gecis"].nunique().clip(upper=4)
    gecis_kukla = pd.get_dummies(gs, prefix="gecis", drop_first=True).astype(float)
    aylik = S.resample("MS").median()
    out = {"meta": {"kaynak": "NASA FIRMS, VIIRS 375 m aktif ateş, gece; NOAA-21 NRT (ana), Suomi NPP SP+NRT (ikinci uydu)",
                    "donem": {"Savaş öncesi": ["2024-09-01", str(ONCE_SON.date())], **{k: list(v) for k, v in DON.items()}},
                    "kalici_hucre_esigi": f"savaş öncesi gecelerde tespit payı en az %{int(KALICI_PAY * 100)}", "hucre_derece": HUCRE,
                    "olcu": "tespit edilen FRP: geçiş içindeki piksellerin toplamı; birden çok geçişte tespit edilen piksellerin ortalama track "
                            "boyu en küçük geçiş (bakış açısının vekili; FRP'ye göre sıralamaz, yalnız tespitli geçişler arasından seçer); "
                            "bölgede toplam (MW)",
                    "gecis_tanimi": f"bir gecenin tespit zamanlarında {GECIS_ARALIK} dakikadan uzun boşluk yeni geçiş",
                    "ana_model": "PPML, dönem + takvim ayı kuklaları, HAC (Bartlett) 14 gözlenen gece, küçük örneklem düzeltmesi yok",
                    "gece_sayisi": int(S.notna().any(axis=1).sum()), "son_tarih": str(S.index.max().date()),
                    "gece_sayisi_snpp": int(S2.notna().any(axis=1).sum()) if S2 is not None else None,
                    "eksik_bloklar_noaa21": eksik_bloklar(S), "eksik_bloklar_snpp": eksik_bloklar(S2) if S2 is not None else None,
                    "ortak_gun": int(len(ortak)) if ortak is not None else None,
                    "cokluk": cokluk(ana, H), "geometri": geometri(ana, H), "snpp_tur": snpp_tur(H)},
           "bolgeler": {}}
    for b in S.columns:
        y = S[b]
        out["bolgeler"][b] = {
            "hucre": int((H["bolge"] == b).sum()),
            "once_medyan_frp": round(float(y[y.index <= ONCE_SON].median()), 2),
            "once_ortalama_frp": round(float(y[y.index <= ONCE_SON].mean()), 2),
            "etki": regresyon(y),
            "etki_log": regresyon(y, "log"),
            **{f"etki_kural_{k}": regresyon(v[b]) if b in v else None for k, v in Sk.items()},
            "etki_gecis_kukla": regresyon(y, ek=gecis_kukla),
            "etki_yerel_gece": regresyon(Sy[b]) if b in Sy else None,
            "etki_hac7": regresyon(y, gecikme=7), "etki_hac28": regresyon(y, gecikme=28),
            "etki_snpp": regresyon(S2[b]) if S2 is not None and b in S2 else None,
            "etki_ortak_noaa21": regresyon(y.loc[ortak]) if ortak is not None else None,
            "etki_ortak_snpp": regresyon(S2.loc[ortak, b]) if S2 is not None and b in S2 else None,
            "yillik": yillik(S, b),
            "sifir_gece_payi": {k: round(float((y[[donem(t) == k for t in y.index]].dropna() == 0).mean()), 3) for k in ("Savaş öncesi", *DON)},
            "aylik": [[str(t.date()), None if pd.isna(v) else round(float(v), 2)] for t, v in aylik[b].items()]}
    # Holm: tablodaki 10 bölge x 2 kapanma (ana model); açılış dönemi ve öteki bölgeler aile dışında
    p = [(b, k, out["bolgeler"][b]["etki"][k]["p"]) for b in TABLO for k in K]
    duz = multipletests([x[2] for x in p], method="holm")[1]
    for (b, k, _), h in zip(p, duz):
        out["bolgeler"][b]["etki"][k]["p_holm"] = float(h)
    out["fark_guney_kuzey"] = fark_testi(S)
    out["fark_guney_kuzey_snpp"] = fark_testi(S2) if S2 is not None else None
    out["somo"] = somo_karsilastir(aylik)
    out["erken_mart"] = erken_mart(S)
    out["saglamlik"] = saglamlik(ana)
    out["harita"] = harita_hucreleri(ana, H)
    (A.SONUC / "a15_alev.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    H.to_csv(F / "kalici_hucreler.csv", index=False)
    pd.set_option("display.width", 240)
    print(pd.DataFrame({b: {"hücre": v["hucre"], "öncesi medyan": v["once_medyan_frp"],
                            **{k: f"{e['etki_yuzde']:+.0f}% [{e['alt']:+.0f}, {e['ust']:+.0f}]" for k, e in v["etki"].items() if k in DON},
                            **{f"K1 {k}": f"{v[f'etki_kural_{k}']['1. kapanma']['etki_yuzde']:+.0f}" for k in KURALLAR if k != KURAL},
                            "K1 yerel": f"{v['etki_yerel_gece']['1. kapanma']['etki_yuzde']:+.0f}"}
                        for b, v in out["bolgeler"].items()}).T.to_string())
    print("güney-kuzey farkı:", {k: round(v["etki_yuzde"], 1) for k, v in out["fark_guney_kuzey"].items()},
          "| Suomi NPP:", {k: round(v["etki_yuzde"], 1) for k, v in (out["fark_guney_kuzey_snpp"] or {}).items()})
    print("çokluk:", out["meta"]["cokluk"]["ozet"])
    print("SOMO korelasyon:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out["somo"]["korelasyon"].items()
                               if k not in ("aylar", "savas_oncesi_aylari", "kapanma_aylari")})


if __name__ == "__main__":
    main()
