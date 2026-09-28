"""Analiz 17: iki kapanmanın farkı. Denizde (radar, AIS) ve karada (alev) göstergeler, ihracat (SOMO) ve fiziki fiyat (Dated Brent - ICE)
iki kapanmada nasıl ayrıştı? Astra'nın değerlendirme notundaki birinci öneri (27 Eylül 2026).

1. Dönem farkı: aynı modelde ikinci ve ilk kapanma katsayılarının farkı, aynı kovaryansla (Wald). Fark = 100 x [exp(b2 - b1) - 1]:
   ikinci kapanmadaki beklenen düzeyin ilk kapanmadakine göre yüzde farkı. İki ayrı güven aralığına bakılarak karar verilmez.
   - Radar: analiz7_sar.py'nin ana modeli (Poisson, dönem + uydu kuklaları, HC1), beş konum ve ana yörünge; Hürmüz'de büyük yankılar da.
   - Alev: analiz15_alev.py'nin ana modeli (PPML, dönem + takvim ayı kuklaları, HAC 14 gözlenen gece), Tablo 3'ün on bölgesi. Güney ve
     kuzey Irak için bulut denetimli sürümler (analiz16_kapsama.py): Level-2 ölçüsü, kapsama >= %90 olan geceler; hücre sabit etkili PPML.
2. Eşit uzunlukta başlangıç pencereleri: iki kapanmanın ilk 30 ve ilk 72 günü (ikinci kapanma veri kesiminde 72 gün). Her kapanma
   "ilk n gün" ve "kalan" diye iki kuklaya ayrılır; ilk n günlerin farkı aynı modelle sınanır. Radar sahneleri seyrek (yörünge başına
   yaklaşık 12 günde bir), bu yüzden 30 günlük pencere yalnız alevlerde.
3. Mekanizma ipuçları (betimsel): Hürmüz koridorunda radarın büyük yankı payı; PortWatch AIS geçişleri (Hürmüz, günlük; toplam ve tanker
   sayısı, tankerlerin tahmini yükü); SOMO Basra ihracatı (aylık, günlük hız); Dated Brent - ICE ön vade farkı ($/varil).
4. Şekil verisi: radar sahneleri, güney Irak alevi (7 gecelik ortalama), Basra aylık günlük hızı, Dated-ICE farkı ve AIS tanker geçişi
   (7 günlük ortalama); sıklıklar korunur, aylık seri günlüğe yayılmaz.
5. Basra testi: Ağustos'ta Basra'dan yüklenen petrol boğazı nasıl geçti? (a) İki serinin karşılaştırması (betimsel): PortWatch'un
   Hürmüz'ü geçen bütün tankerler için AIS'ten tahmin ettiği yük (capacity_tanker: su çekiminden tahmin edilen doluluk oranı x DWT,
   metrik ton; iki yön) ton başına 7,33 varille (Energy Institute ortalaması) ham petrol eşdeğerine çevrilir ve SOMO'nun Basra yüklemesiyle
   (konşimento tarihine göre) karşılaştırılır. Seri DWT toplamı değil, tahmini yüktür: oran taşınabilecek payın üst sınırı ya da görünmeyen
   petrolün payı değildir (delta denetimi D01). (b) Basra açık deniz terminalleri ve demirleme alanında radar sayımı (sar_gemi_sayimi.py,
   kutu 48,70-49,20 D, 29,40-29,85 K; iki yörünge ayrı seri, ayrı sabit hedef maskesi; ana seri çıkan yörünge 174). Dönem etkileri radar
   ana modeliyle; eğilim: gemi ~ gün + uydu (+ göreli yörünge) kuklası, OLS (HC1 ve HC3, t dağılımı), Poisson sağlamlık; bütün ikinci
   kapanma, yalnız Ağustos ve inen seride göreli yörünge başına. Karşılaştırma ölçeği varsayımsaldır: Ağustos yüklemesinin her yükü kutuya
   eklenen ayrı ve tespit edilebilir bir gemide kalsaydı sayı 30 günde yükleme kadar (VLCC ya da Suezmax eşdeğeri) artardı; gemi sayısı
   yük stoku değildir (kutudaki boş kapasite, çıkışlar ve yan yana gemiler bu bağı bozar).
6. Basra terminalinde AIS görünürlüğü (28 Eylül 2026): PortWatch'un günlük liman verisi (Daily_Ports_Data; portwatch_liman_cek.py) Basra
   petrol terminalinde tankerlerin yüklediği tahmini yükü (export_tanker, ton; su çekimi değişiminden) veriyor. Aylık ortalama 7,33 varil/ton
   ile ham petrol eşdeğerine çevrilir ve SOMO'nun Basra yüklemesine bölünür: savaştan önce bu oran kamuya açık AIS verisinin yüklemenin ne
   kadarını gördüğünü, Ağustos 2026'daki değeri ise görünürlüğün ne kadar düştüğünü gösterir. Betimseldir; üst sınır değildir.
Lisanslı fiyat verisi ve PortWatch içeren Excel ya da PortWatch liman dosyası yoksa (herkese açık kod paketi) PortWatch, AIS/Basra, liman ve
Dated-ICE bölümleri atlanır; radar, alev, Basra radarı, SOMO ve küçük akış hesabı aynen üretilir.
Çıktı: sonuclar/a17_iki_kapanma.json
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

import analiz15_alev as AL
import analiz16_kapsama as KP
import analiz_ortak as A

K1, K2 = "1. kapanma", "2. kapanma"
ONCE_R = "Tem 2025-Şub 2026"
DON_R = {ONCE_R: ("2025-07-01", "2026-02-27"), K1: ("2026-03-02", "2026-06-13"), "Açılış": ("2026-06-18", "2026-07-10"),
         K2: ("2026-07-13", "2026-09-22")}
RADAR = {"hurmuz": ("_40m", "descending", "Hürmüz koridoru"), "fuceyre": ("", "descending", "Füceyre"), "yanbu": ("", "ascending", "Yanbu"),
         "babulmendep": ("_40m_asc", "ascending", "Bâbülmendep"), "rastanura": ("_40m", "descending", "Ras Tanura")}
# Basra açıkları: savaş öncesi verisi a7'de yok (yeniden üretim denetimine girmez); inen yörünge sağlamlık serisi
BASRA = {"basra": ("_40m", "ascending", "Basra terminalleri (çıkan yörünge)"),
         "basra_inen": ("_40m_desc", "descending", "Basra terminalleri (inen yörünge)")}
CEVRIM = 7.33  # varil/ton, ham petrol (Energy Institute, Statistical Review of World Energy, yaklaşık çevrim katsayıları)
# PortWatch'un capacity_tanker alanı: AIS'e dayalı tahmini yük (metrik ton; geminin özellikleri ve bildirdiği su çekiminden tahmin edilen
# doluluk oranı x DWT). DWT toplamı değildir (IMF PortWatch veri sözlüğü). Excel'de sütunun adı 28 Eylül 2026'ya kadar yanlış olarak
# "Tanker kapasitesi (DWT)" idi; delta denetiminden sonra düzeltildi.
PW_YUK = "Tanker tahmini yükü (ton)"
VLCC, SUEZMAX = 2e6, 1e6  # yük, varil
TUTARLI = ("S1C", "S1D")
PENCERE = (30, 72)
BASLA = {K1: pd.Timestamp("2026-03-02"), K2: pd.Timestamp("2026-07-13")}


def fark(b: pd.Series, V: pd.DataFrame, i: str, j: str) -> dict:
    """j'nin i'ye göre farkı: 100 x [exp(b_j - b_i) - 1], %95 aralık ve Wald p (aynı kovaryansla)."""
    d = b[j] - b[i]
    s = float(np.sqrt(V.loc[i, i] + V.loc[j, j] - 2 * V.loc[i, j]))
    return {"fark_yuzde": round(100 * (np.exp(d) - 1), 3), "alt": round(100 * (np.exp(d - 1.96 * s) - 1), 3),
            "ust": round(100 * (np.exp(d + 1.96 * s) - 1), 3), "p": float(2 * stats.norm.sf(abs(d / s)))}


def etki(b: pd.Series, V: pd.DataFrame, k: str) -> dict:
    s = float(np.sqrt(V.loc[k, k]))
    return {"etki_yuzde": round(100 * (np.exp(b[k]) - 1), 3), "alt": round(100 * (np.exp(b[k] - 1.96 * s) - 1), 3),
            "ust": round(100 * (np.exp(b[k] + 1.96 * s) - 1), 3), "p": float(2 * stats.norm.sf(abs(b[k] / s)))}


def pencere_etiketi(t: pd.Timestamp, n: int, taban: str, donem) -> str | None:
    """Dönem etiketi, kapanmalar "ilk n gün" ve "kalan" diye ayrılarak."""
    d = donem(t)
    if d in (K1, K2):
        return f"{d} ilk {n}" if (t - BASLA[d]).days < n else f"{d} kalan"
    return d


# ---------------------------------------------------------------- radar
def radar_veri(k: str) -> pd.DataFrame:
    ek, yor, _ = RADAR[k] if k in RADAR else BASRA[k]
    S = pd.read_csv(A.VERI / "sar" / f"gemi_sayimi_{k.split('_')[0]}{ek}.csv", parse_dates=["tarih"])
    S = S[S["yorunge"] == yor].copy()
    S["uydu"] = S["sahne"].str[:3]
    gy = A.VERI / "sar" / f"goreli_yorunge_{k.split('_')[0]}.csv"  # yalnız Basra'da var (sar_gemi_sayimi.py --goreli)
    if gy.exists():
        S = S.merge(pd.read_csv(gy)[["sahne", "goreli_yorunge"]], on="sahne", how="left")
    return S.sort_values("tarih")


def radar_donem(t: pd.Timestamp) -> str | None:
    for d, (a, b) in DON_R.items():
        if pd.Timestamp(a) <= t <= pd.Timestamp(b):
            return d
    return None


def radar_model(g: pd.DataFrame, etiket: pd.Series, sutun: str = "gemi"):
    """analiz7_sar.py ana modeli: Poisson, dönem + uydu kuklaları, HC0 x n/(n-k) (HC1)."""
    g = g.assign(donem=etiket).dropna(subset=["donem"])
    X = pd.get_dummies(g["donem"]).drop(columns=ONCE_R).astype(float)
    if g["uydu"].nunique() > 1:
        X = X.join(pd.get_dummies(g["uydu"], prefix="uydu", drop_first=True).astype(float))
    if "goreli_yorunge" in g and g["goreli_yorunge"].nunique() > 1:  # aynı yönde iki göreli yörünge (Basra inen: 35 ve 108)
        X = X.join(pd.get_dummies(g["goreli_yorunge"], prefix="yor", drop_first=True).astype(float))
    X = sm.add_constant(X)
    m = sm.GLM(g[sutun].astype(float), X, family=sm.families.Poisson()).fit(cov_type="HC0")
    V = m.cov_params() * m.nobs / (m.nobs - len(m.params))
    return m.params, V, g


def radar_bolumu() -> dict:
    out = {}
    for k, (_, _, ad) in RADAR.items():
        S = radar_veri(k)
        b, V, g = radar_model(S, S["tarih"].map(radar_donem))
        r = {"ad": ad, "n": int(len(g)), K1: etki(b, V, K1), K2: etki(b, V, K2), "fark_2_1": fark(b, V, K1, K2),
             "n_donem": {d: int((g["donem"] == d).sum()) for d in DON_R}}
        if k == "hurmuz":
            b, V, _ = radar_model(S, S["tarih"].map(radar_donem), "buyuk_gemi")
            r["buyuk"] = {K1: etki(b, V, K1), K2: etki(b, V, K2), "fark_2_1": fark(b, V, K1, K2)}
        n = 72
        et = S["tarih"].map(lambda t: pencere_etiketi(t, n, ONCE_R, radar_donem))
        b, V, g = radar_model(S, et)
        i, j = f"{K1} ilk {n}", f"{K2} ilk {n}"
        if i in b and j in b:
            r[f"pencere_{n}"] = {i: etki(b, V, i), j: etki(b, V, j), "fark_2_1": fark(b, V, i, j),
                                 "n": {i: int((g["donem"] == i).sum()), j: int((g["donem"] == j).sum())}}
        out[k] = r
    # büyük yankı payı ve sahne medyanları (yalnız 1C ve 1D; Hürmüz koridoru)
    S = radar_veri("hurmuz")
    T = S[S["uydu"].isin(TUTARLI)].assign(donem=lambda d: d["tarih"].map(radar_donem)).dropna(subset=["donem"])
    out["hurmuz_betimsel"] = {d: {"sahne": int(len(x)), "gemi_medyan": float(x["gemi"].median()), "buyuk_medyan": float(x["buyuk_gemi"].median()),
                                  "buyuk_payi": round(float(x["buyuk_gemi"].sum() / x["gemi"].sum()), 3)} for d, x in T.groupby("donem")}
    return out


def birikim(S: pd.DataFrame, a: str = DON_R[K2][0], b: str = DON_R[K2][1]) -> dict:
    """Doğrusal eğilim, gemi/30 gün (varsayılan: bütün ikinci kapanma): gemi ~ gün + uydu (+ göreli yörünge) kuklası. OLS; aralıklar HC1
    ve HC3 standart hatasıyla, t dağılımından (küçük örneklemde HC3 daha temkinli). Poisson (aynı kuklalar, HC1): oransal eğilim, %/30 gün;
    sağlamlık. Eğilim gemi sayısınındır; yük stoku değildir."""
    g = S[(S["tarih"] >= a) & (S["tarih"] <= b)].copy()
    X = pd.DataFrame({"gun": (g["tarih"] - pd.Timestamp(a)).dt.days.astype(float)}, index=g.index)
    if g["uydu"].nunique() > 1:
        X = X.join(pd.get_dummies(g["uydu"], prefix="uydu", drop_first=True).astype(float))
    if "goreli_yorunge" in g and g["goreli_yorunge"].nunique() > 1:
        X = X.join(pd.get_dummies(g["goreli_yorunge"], prefix="yor", drop_first=True).astype(float))
    X = sm.add_constant(X)
    out = {"n": int(len(g)), "parametre": int(X.shape[1]), "ilk": str(g["tarih"].min().date()), "son": str(g["tarih"].max().date())}
    for sutun in ("gemi", "buyuk_gemi"):
        r = {}
        for cov, ek in (("HC1", ""), ("HC3", "_hc3")):
            m = sm.OLS(g[sutun].astype(float), X).fit(cov_type=cov)
            e, s, tq = 30 * m.params["gun"], 30 * m.bse["gun"], stats.t.ppf(0.975, m.df_resid)
            if not ek:
                r["egim_30_gun"] = round(float(e), 2)
            r.update({f"alt{ek}": round(float(e - tq * s), 2), f"ust{ek}": round(float(e + tq * s), 2),
                      f"p{ek}": float(2 * stats.t.sf(abs(e / s), m.df_resid))})
        out[sutun] = r
    m = sm.GLM(g["gemi"].astype(float), X, family=sm.families.Poisson()).fit(cov_type="HC0")
    V = m.cov_params() * m.nobs / (m.nobs - len(m.params))
    bb, ss = float(m.params["gun"]), float(np.sqrt(V.loc["gun", "gun"]))
    out["poisson_yuzde_30_gun"] = {"etki": round(100 * (np.exp(30 * bb) - 1), 2), "alt": round(100 * (np.exp(30 * (bb - 1.96 * ss)) - 1), 2),
                                   "ust": round(100 * (np.exp(30 * (bb + 1.96 * ss)) - 1), 2), "p": float(2 * stats.norm.sf(abs(bb / ss)))}
    return out


def basra_bolumu(so: dict) -> dict:
    """Basra açıklarında radar: dönem etkileri (radar ana modeli), ikinci-ilk kapanma farkı, betimsel medyanlar (1C ve 1D) ve birikme testi."""
    ag = so["2026-08"]["bin_varil_gun"] * 1e3 * 31
    out = {"beklenen_30_gun": {"vlcc": round(ag / VLCC * 30 / 31, 1), "suezmax": round(ag / SUEZMAX * 30 / 31, 1),
                               "varsayim": "Varsayımsal ölçek: SOMO'nun Ağustos 2026 Basra yüklemesinin her yükü (VLCC 2 milyon, Suezmax 1 milyon "
                                           "varil) kutuya eklenen ayrı ve tespit edilebilir bir gemide kalsaydı 30 günde eklenecek gemi. Kutudaki "
                                           "boş kapasite, öteki gemilerin çıkışı ve ay içi zamanlama modellenmez; gözlenen stok ya da karşıolgusal "
                                           "tahmin değildir."}}
    for k, (ek, yor, ad) in BASRA.items():
        if not (A.VERI / "sar" / f"gemi_sayimi_basra{ek}.csv").exists():
            continue
        S = radar_veri(k)
        b, V, g = radar_model(S, S["tarih"].map(radar_donem))
        T = S[S["uydu"].isin(TUTARLI)].assign(donem=lambda d: d["tarih"].map(radar_donem)).dropna(subset=["donem"])
        out[k] = {"ad": ad, "yorunge": yor, "sahne": int(len(S)), "sahne_1c1d": int(S["uydu"].isin(TUTARLI).sum()), "n": int(len(g)),
                  K1: etki(b, V, K1), K2: etki(b, V, K2),
                  "fark_2_1": fark(b, V, K1, K2), "n_donem": {d: int((g["donem"] == d).sum()) for d in DON_R},
                  "uydular": sorted(g["uydu"].unique().tolist()),
                  "goreli_yorunge": {str(int(y)): int(n) for y, n in g["goreli_yorunge"].value_counts().sort_index().items()},
                  "medyan_1c1d": {d: {"sahne": int(len(x)), "gemi": float(x["gemi"].median()), "buyuk": float(x["buyuk_gemi"].median())}
                                  for d, x in T.groupby("donem")},
                  "birikim": birikim(S[S["uydu"].isin(TUTARLI)]),
                  "birikim_agustos": birikim(S[S["uydu"].isin(TUTARLI)], "2026-08-01", "2026-08-31")}
        if S["goreli_yorunge"].nunique() > 1:  # inen seri: göreli yörünge başına (35 ve 108), bütün ikinci kapanma
            T2 = S[S["uydu"].isin(TUTARLI)]
            out[k]["birikim_goreli"] = {str(int(y)): birikim(T2[T2["goreli_yorunge"] == y]) for y in sorted(T2["goreli_yorunge"].dropna().unique())}
        bb, VV, _ = radar_model(S, S["tarih"].map(radar_donem), "buyuk_gemi")
        out[k]["buyuk"] = {K1: etki(bb, VV, K1), K2: etki(bb, VV, K2), "fark_2_1": fark(bb, VV, K1, K2)}
    return out


def ais_basra(P: pd.DataFrame, so: dict) -> dict:
    """PortWatch'un Hürmüz'ü geçen tankerler için AIS'ten tahmin ettiği günlük yük (ham petrol eşdeğeri) ile SOMO Basra yüklemesi, aylık."""
    m = P.resample("MS").agg({"Tanker": "mean", PW_YUK: "mean"})
    m["gun"] = P["Tanker"].resample("MS").size()
    m["gecis"] = P["Tanker"].resample("MS").sum()
    out = {}
    for t, r in m.loc["2025-10-01":].iterrows():
        ay = t.strftime("%Y-%m")
        yuk = float(r[PW_YUK]) * CEVRIM / 1e3
        bas = so.get(ay, {}).get("bin_varil_gun")
        out[ay] = {"gun": int(r["gun"]), "tanker_gun": round(float(r["Tanker"]), 2), "tanker_gecis": int(r["gecis"]),
                   "tahmini_yuk_bin_varil_gun": round(yuk, 1),
                   "basra_bin_varil_gun": bas, "oran_yuk_basra": round(yuk / bas, 3) if bas else None}
    return {"cevrim_varil_ton": CEVRIM,
            "tanim": "tahmini_yuk: PortWatch capacity_tanker (AIS'e dayalı tahmini tanker yükü, metrik ton; doluluk tahmini x DWT; iki yön, "
                     "bütün tankerler ve yük türleri) x 7,33 varil/ton; Basra: SOMO yüklemesi, konşimento tarihine göre. Oran iki tahmini "
                     "karşılaştırır; taşınabilecek payın üst sınırı değildir.",
            "aylik": out}


def pw_liman(so: dict) -> dict | None:
    """PortWatch liman verisi: Basra petrol terminali (port2479) ve Umm Kasr (port1341); aylık tanker uğrağı, tahmini yükleme (bin varil/gün,
    ham petrol eşdeğeri) ve Basra terminalinin SOMO Basra yüklemesine oranı. Dosya yoksa None."""
    f = A.VERI / "portwatch" / "irak_liman_gunluk.csv"
    if not f.exists():
        return None
    L = pd.read_csv(f, parse_dates=["date"])
    out = {"kaynak": "IMF PortWatch Daily_Ports_Data (portwatch_liman_cek.py)", "cevrim_varil_ton": CEVRIM, "aylik": {}}
    for pid, ad in (("port2479", "basra_terminali"), ("port1341", "umm_kasr")):
        x = L[L["portid"] == pid].set_index("date").sort_index()
        m = x.resample("MS").agg({"portcalls_tanker": "sum", "export_tanker": "mean"})
        m["gun"] = x["export_tanker"].resample("MS").size()
        for t, r in m.loc["2025-10-01":].iterrows():
            ay = t.strftime("%Y-%m")
            yuk = float(r["export_tanker"]) * CEVRIM / 1e3
            bas = so.get(ay, {}).get("bin_varil_gun") if ad == "basra_terminali" else None
            out["aylik"].setdefault(ay, {})[ad] = {"gun": int(r["gun"]), "tanker_ugragi": int(r["portcalls_tanker"]),
                                                   "tahmini_yukleme_bin_varil_gun": round(yuk, 1)}
            if ad == "basra_terminali":
                out["aylik"][ay][ad]["oran_somo_basra"] = round(yuk / bas, 4) if bas else None
    once = [out["aylik"][a]["basra_terminali"]["oran_somo_basra"] for a in ("2025-10", "2025-11", "2025-12", "2026-01", "2026-02")]
    out["savas_oncesi_oran"] = {"min": min(once), "max": max(once), "ortalama": round(sum(once) / len(once), 4)}
    out["son_tarih"] = str(L["date"].max().date())
    return out


# ---------------------------------------------------------------- alev
def alev_model(y: pd.Series, etiket):
    """analiz15_alev.py ana modeli: PPML, dönem + takvim ayı kuklaları, HAC 14 (gözlenen gece sırası)."""
    y = y.dropna()
    d = pd.Series([etiket(t) for t in y.index], index=y.index)
    y, d = y[d.notna()], d[d.notna()]
    X = pd.get_dummies(d).drop(columns="Savaş öncesi").astype(float)
    X = X.join(pd.get_dummies(y.index.month, prefix="ay", drop_first=True).set_index(y.index).astype(float))
    X = sm.add_constant(X)
    m = sm.GLM(y, X, family=sm.families.Poisson()).fit(cov_type="HAC", cov_kwds={"maxlags": 14})
    return m.params, m.cov_params(), d


def alev_ozet(y: pd.Series) -> dict:
    b, V, d = alev_model(y, AL.donem)
    r = {K1: etki(b, V, K1), K2: etki(b, V, K2), "fark_2_1": fark(b, V, K1, K2)}
    for n in PENCERE:
        b, V, d = alev_model(y, lambda t: pencere_etiketi(t, n, "Savaş öncesi", AL.donem))
        i, j = f"{K1} ilk {n}", f"{K2} ilk {n}"
        r[f"pencere_{n}"] = {i: etki(b, V, i), j: etki(b, V, j), "fark_2_1": fark(b, V, i, j),
                             "gece": {i: int((d == i).sum()), j: int((d == j).sum())}}
    return r


def hucre_fark(T: pd.DataFrame, b: str) -> dict:
    """analiz16_kapsama.hucre_model ile aynı tasarım (hücre sabit etkili PPML, Driscoll-Kraay); ikinci ve ilk kapanma farkı."""
    x = T[(T["bolge"] == b) & T["gecerli"]].copy()
    x["donem"] = KP.donem_sutunu(x["tarih"])
    x = x[x["donem"] != "geçiş"]
    X = pd.get_dummies(x["donem"]).reindex(columns=[K1, "Açılış", K2], fill_value=False).astype(float)
    X = X.join(pd.get_dummies(x["tarih"].dt.month, prefix="ay").astype(float).drop(columns="ay_1", errors="ignore"))
    X = X.loc[:, X.sum() > 0]
    gece = pd.Index(sorted(x["tarih"].unique()))
    r = KP.ppml_fe(x["frp"].to_numpy(float), pd.factorize(x["anahtar"])[0], gece.get_indexer(x["tarih"]), X.to_numpy(), 14)
    bb, V = pd.Series(r["b"], index=X.columns), pd.DataFrame(r["V"], index=X.columns, columns=X.columns)
    return {K1: etki(bb, V, K1), K2: etki(bb, V, K2), "fark_2_1": fark(bb, V, K1, K2)}


def alev_bolumu(S: pd.DataFrame) -> dict:
    out = {"bolgeler": {b: alev_ozet(S[b]) for b in AL.TABLO}}
    B = pd.read_csv(KP.KLASOR / "bolge_gece.csv", parse_dates=["tarih"])
    T = pd.read_csv(KP.KLASOR / "hucre_gece.csv.gz", parse_dates=["tarih"], usecols=["tarih", "anahtar", "bolge", "gecerli", "frp"])
    for b in ("Irak (güney)", "Irak (kuzey)"):
        out["bolgeler"][b]["l2_k90"] = alev_ozet(KP.seri(B, b, "l2", 0.9))
        out["bolgeler"][b]["hucre_l2"] = hucre_fark(T, b)
    return out


# ---------------------------------------------------------------- betimsel göstergeler
def portwatch() -> tuple[dict, pd.DataFrame]:
    P = pd.read_excel(A.XLSX, sheet_name="PortWatch_Gunluk", parse_dates=["Tarih"])
    P = P[P["Geçit"] == "Hürmüz"].set_index("Tarih").sort_index()
    donem = {"Oca-Şub 2026": ("2026-01-01", "2026-02-27"), "Tem 2025-Şub 2026": ("2025-07-01", "2026-02-27"), K1: DON_R[K1],
             "Açılış": DON_R["Açılış"], K2: DON_R[K2]}
    out = {}
    for d, (a, b) in donem.items():
        x = P.loc[a:b]
        out[d] = {"gun": int(len(x)), "gecis": round(float(x["Toplam geçiş"].mean()), 2), "tanker": round(float(x["Tanker"].mean()), 2),
                  "tanker_tahmini_yuk_mt": round(float(x[PW_YUK].mean()) / 1e6, 4),
                  "gecis_basina_tahmini_yuk_bin": round(float(x[PW_YUK].sum() / max(x["Tanker"].sum(), 1)) / 1e3, 1),
                  "son_tarih": str(x.index.max().date()) if len(x) else None}
    return out, P


def somo() -> dict:
    so = pd.read_csv(A.VERI / "rafineri" / "irak_somo_cikis_aylik.csv")
    b = so[so["cikis"] == "Basra"].groupby("ay").agg(varil=("miktar_varil", "sum"), gun=("gun", "first"))
    b["bin_varil_gun"] = b["varil"] / b["gun"] / 1e3
    k = so[so["cikis"] == "Basra"].groupby("ay")["kaynak"].first()
    return {ay: {"bin_varil_gun": round(float(r["bin_varil_gun"]), 1), "kaynak": "grafik" if "grafi" in str(k[ay]).lower() else "PDF"}
            for ay, r in b.iterrows() if ay >= "2025-10"}


def dated_ice() -> tuple[dict, pd.Series]:
    F = A.fiyatlar(A.excel())
    s = (F["dated"] - F["ice_m1"]).dropna()
    s = s[s.index >= "2025-07-01"]
    # "Oca-Şub 2026": raporun fiyat bölümlerindeki savaş öncesi penceresi (s. 26'daki 1,7 dolar); uydu bölümü radar tabanını kullanır
    donem = {"Oca-Şub 2026": ("2026-01-01", "2026-02-27"), ONCE_R: DON_R[ONCE_R], K1: DON_R[K1], "Açılış": DON_R["Açılış"], K2: DON_R[K2]}
    out = {d: {"gun": int(len(s.loc[a:b])), "ortalama": round(float(s.loc[a:b].mean()), 2), "medyan": round(float(s.loc[a:b].median()), 2)}
           for d, (a, b) in donem.items()}
    out["son_tarih"] = str(s.index.max().date())
    return out, s


def main():
    H = pd.read_csv(A.VERI / "firms" / "kalici_hucreler.csv")
    ana = AL.yukle(["VIIRS_NOAA21_NRT"])
    S = AL.gece_serisi(ana, H)
    a15 = json.loads((A.SONUC / "a15_alev.json").read_text(encoding="utf-8"))
    a7 = json.loads((A.SONUC / "a7_sar_ozet.json").read_text(encoding="utf-8"))
    R, L = radar_bolumu(), alev_bolumu(S)
    # yeniden üretim denetimi: dönem etkileri a7 ve a15 ile aynı olmalı
    den = {"radar_en_buyuk_fark": max(abs(R[k][d]["etki_yuzde"] - a7[k]["yorunge"][RADAR[k][1]]["regresyon"][d]["etki_yuzde"])
                                      for k in RADAR for d in (K1, K2)),
           "alev_en_buyuk_fark": max(abs(L["bolgeler"][b][d]["etki_yuzde"] - a15["bolgeler"][b]["etki"][d]["etki_yuzde"])
                                     for b in AL.TABLO for d in (K1, K2))}
    excel = A.XLSX.exists()  # paket: lisanslı fiyat verisi ve PortWatch yok
    if excel:
        pw, P = portwatch()
        di, spread = dated_ice()
    # şekil verisi
    Hr = radar_veri("hurmuz")
    Hr = Hr[Hr["uydu"].isin(TUTARLI) & (Hr["tarih"] >= "2025-12-01")]
    g = S["Irak (güney)"]
    taban = float(g.loc[:AL.ONCE_SON].median())
    g7 = (100 * g / taban).rolling(7, min_periods=4).mean().loc["2025-12-01":]
    sekil = {"radar": [[str(t.date()), int(v)] for t, v in zip(Hr["tarih"], Hr["gemi"])],
             "alev_guney_7": [[str(t.date()), round(float(v), 2)] for t, v in g7.dropna().items()],
             "basra_aylik": {ay: v["bin_varil_gun"] for ay, v in somo().items()},
             "radar_once_medyan": a7["hurmuz"]["yorunge"]["descending"]["donem"][ONCE_R]["medyan"]}
    if excel:
        pw7 = P["Tanker"].rolling(7, min_periods=4).mean().loc["2025-12-01":]
        sekil["dated_ice"] = [[str(t.date()), round(float(v), 3)] for t, v in spread.loc["2025-12-01":].items()]
        sekil["ais_tanker_7"] = [[str(t.date()), round(float(v), 2)] for t, v in pw7.dropna().items()]
    # küçük bir tanker akışı koridordaki anlık gemi sayısına ne ekler? (Little yasası: anlık sayı = akış x geçiş süresi; yaklaşık hesap)
    so = somo()
    akis = so["2026-08"]["bin_varil_gun"] * 1e3
    km, knot = 70.0, 12.0  # koridor kutusunun doğu-batı boyu (0,7 derece, 26,5 K) ve tanker hızı (varsayım)
    saat = km / (knot * 1.852)
    kucuk = {"varsayim": "Ağustos 2026 Basra yüklemesi; VLCC 2 milyon, Suezmax 1 milyon varil; kutudan geçiş 70 km, 12 knot",
             "basra_varil_gun": akis, "gecis_saat": round(saat, 2),
             "anlik_gemi_vlcc": round(akis / 2e6 * saat / 24, 4), "anlik_gemi_suezmax": round(akis / 1e6 * saat / 24, 4),
             "gunluk_vlcc": round(akis / 2e6, 2)}
    # Basra testi
    ba = basra_bolumu(so)
    if excel:
        ais = ais_basra(P, so)
        kucuk["ais_agustos"] = ais["aylik"]["2026-08"]
    for k in BASRA:
        if k in ba:
            B = radar_veri(k)
            B = B[B["uydu"].isin(TUTARLI)]
            sekil[f"{k}_radar"] = [[str(t.date()), int(v)] for t, v in zip(B["tarih"], B["gemi"])]
            sekil[f"{k}_once_medyan"] = ba[k]["medyan_1c1d"][ONCE_R]["gemi"]
    out = {"meta": {"tanim": __doc__.strip(), "denetim": den}, "kucuk_akis": kucuk, "radar": R, "alev": L, "somo_basra": somo(), "basra": ba,
           "sekil": sekil}
    if excel:
        out.update({"portwatch_hurmuz": pw, "dated_ice": di, "ais_basra": ais})
    li = pw_liman(somo())
    if li:
        out["portwatch_basra_terminal"] = li
    (A.SONUC / "a17_iki_kapanma.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    # özet
    print("denetim:", den)
    for k, r in R.items():
        if k in RADAR:
            f = r["fark_2_1"]
            print(f"radar {r['ad']:16s} K1 {r[K1]['etki_yuzde']:+6.1f} K2 {r[K2]['etki_yuzde']:+6.1f} | K2/K1 {f['fark_yuzde']:+6.1f} "
                  f"[{f['alt']:+.0f}, {f['ust']:+.0f}] p={f['p']:.3f}" + (f" | ilk 72: {r['pencere_72']['fark_2_1']['fark_yuzde']:+.1f} "
                                                                           f"p={r['pencere_72']['fark_2_1']['p']:.3f}" if "pencere_72" in r else ""))
    for b, r in L["bolgeler"].items():
        f = r["fark_2_1"]
        ek = " | ".join(f"ilk {n}: {r[f'pencere_{n}']['fark_2_1']['fark_yuzde']:+.0f} p={r[f'pencere_{n}']['fark_2_1']['p']:.3f}" for n in PENCERE)
        print(f"alev {b:28s} K1 {r[K1]['etki_yuzde']:+6.1f} K2 {r[K2]['etki_yuzde']:+6.1f} | K2/K1 {f['fark_yuzde']:+7.1f} [{f['alt']:+.0f}, "
              f"{f['ust']:+.0f}] p={f['p']:.3f} | {ek}")
        for v in ("l2_k90", "hucre_l2"):
            if v in r:
                print(f"      {v:10s} K2/K1 {r[v]['fark_2_1']['fark_yuzde']:+7.1f} [{r[v]['fark_2_1']['alt']:+.0f}, {r[v]['fark_2_1']['ust']:+.0f}] "
                      f"p={r[v]['fark_2_1']['p']:.4f}")
    if excel:
        print("PortWatch Hürmüz:", json.dumps(pw, ensure_ascii=False))
    print("Hürmüz radar betimsel:", json.dumps(R["hurmuz_betimsel"], ensure_ascii=False))
    print("SOMO Basra:", json.dumps(out["somo_basra"], ensure_ascii=False))
    if excel:
        print("Dated-ICE:", json.dumps(di, ensure_ascii=False))
        print("AIS / Basra:", json.dumps(ais, ensure_ascii=False))
    for k in BASRA:
        if k in ba:
            r = ba[k]
            print(f"Basra {k}: n={r['n']} {r['n_donem']} K1 {r[K1]['etki_yuzde']:+.1f} [{r[K1]['alt']:+.0f}, {r[K1]['ust']:+.0f}] "
                  f"K2 {r[K2]['etki_yuzde']:+.1f} [{r[K2]['alt']:+.0f}, {r[K2]['ust']:+.0f}] | K2/K1 {r['fark_2_1']['fark_yuzde']:+.1f} "
                  f"p={r['fark_2_1']['p']:.3f}")
            print("   medyan 1C/1D:", r["medyan_1c1d"], "| birikim:", r["birikim"], "| Ağustos:", r["birikim_agustos"])
    print("beklenen (30 gün):", ba["beklenen_30_gun"])
    if li:
        print("Basra terminali (PortWatch / SOMO):", {a: v["basra_terminali"] for a, v in li["aylik"].items()}, "| savaş öncesi:", li["savas_oncesi_oran"])


if __name__ == "__main__":
    main()
