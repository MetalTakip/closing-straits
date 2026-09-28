"""Analiz 16: alev ölçüsünün gözlenebilirliği (tam çalıştırma: 1 Eylül 2024-24 Eylül 2026, Tablo 3'ün on bölgesi). İkinci çapraz denetimin
4. bölümündeki bulut ve kapsama protokolünü ve üçüncü denetimin Z01 notunu (evren tespit olmayan geçişleri de içermeli) uygular.

Girdi: kapsama_tam.py'nin parçaları (veri/firms/kapsama/parca/). Her NOAA-21 gece granül çifti (VJ214IMG + CLDMSK_L2_VIIRS_NOAA21) için
hücre başına yangın maskesi ve bulut maskesi sayımları, bakış açısı; hücredeki yangın pikselleri, FRP, QA bit 22 (kalan bowtie tekrarı).
FIRMS tarafı analiz15_alev.py'nin ana ölçüsü ve modeli.

Tanımlar (pilot ve denetim protokolüne göre, sonuçlara bakmadan yazıldı):
  Geçiş: aynı UTC gecesinin granülleri 20 dakikadan uzun boşlukla ayrılıyorsa ayrı uydu geçişleridir; aynı yörüngenin art arda granülleri tek
    geçiştir (granül sınırındaki hücrenin sayımları toplanır).
  Hücre-geçiş durumu (pilotla aynı): hücreye düşen M-bant piksellerindeki I-bant piksellerinin en az yarısı su, açık kara ya da yangın ise
    "geçerli", en az yarısı bulut ise "bulutlu", değilse "belirsiz" (işlenemedi, bowtie silmesi, güneş parlaması, sınıflanamadı).
  Hücre-gece durumu: "granül yok" (o gece işlenmiş granül çifti yok); "kapsanmadı" (granül var, hiçbir geçiş hücreyi görmedi); "bulutlu"
    (bütün geçişler bulutlu); "belirsiz" (geçerli geçiş yok, öteki durumlar); "geçerli, tespit yok"; "geçerli, tespit var". İlk dördü gözlenen
    sıfır yapılmaz. Geçerli tespit yokluğu da fiziksel sıfırın kesin kanıtı değildir.
  Geçiş seçimi (sonuçtan bağımsız): hücrenin geçerli olduğu geçişler arasında bakış açısı (bulut ürününün sensor_zenith değeri, hücredeki
    M-bant piksellerinin medyanı) en küçük olan; eşitlikte önceki geçiş. FRP'ye ve tespite bakmaz.
  L2 ölçüsü: seçilen geçişte hücredeki yangın piksellerinin FRP toplamı; QA bit 22 işaretli pikseller çıkarılır; tespit yoksa 0. Bölge-gece
    değeri, bölgenin geçerli hücrelerinin toplamı (MW).
  Kapsama: bölgenin sabit hücre listesinde o gece geçerli olan hücrelerin payı. Eşikler %80, %90, %95 (protokolün önerisi, evrensel değil).
Modeller (analiz15_alev.py'nin ana modeli: PPML, dönem + takvim ayı kuklaları, HAC 14 gözlenen gece, küçük örneklem düzeltmesi yok):
  firms_tum     FIRMS ana ölçüsü, bütün geceler (a15_alev.json'u birebir yeniden üretmeli; denetim).
  firms_kXX     FIRMS ana ölçüsü, bölgenin kapsamasının en az %XX olduğu geceler.
  l2_kXX        L2 ölçüsü, aynı geceler.
  l2bt_k90      L2 ölçüsü, bowtie işaretli pikseller dahil (%90).
  l2silmesiz_kXX, hucre_l2silmesiz: geçerlilik kuralında bowtie silmesi (yangın maskesi sınıf 1) paydadan çıkarılır; silme gözlem kaybı değil.
  firmsgecerli_k90: FIRMS'in hücre-gece değerlerinin yalnız geçerli hücrelerdeki toplamı (%90).
  l2siki_kXX, hucre_l2siki: REDDEDİLDİ, yalnız belge için. Bağımsız bulut ürünü de açık demeli; ama ürün alev piksellerini bulut sayıyor
    (bulut_urunu_frp: bulutlu payı FRP ile artıyor), yani bu seçim sonuca bağlı.
  hucre_l2, hucre_firms: hücre-gece PPML, hücre sabit etkili, yalnız geçerli hücre-geceler; takvim ayı kuklaları. Standart hata
    Driscoll-Kraay: skorlar gece içinde toplanır (aynı gecenin hücreleri bağımsız sayılmaz), Bartlett 14 gözlenen gece; statsmodels'ın
    hac-groupsum varsayılanı gibi küme düzeltmeli. FIRMS sürümünde bağımlı değişken FIRMS'in hücre-gece değeri (tespit yoksa 0).
  Güney-kuzey farkı (analiz15_alev.fark_testi): iki bölgenin de kapsamasının eşiği geçtiği gecelerde; FIRMS ve L2 ölçüleriyle.
  Holm: her model ailesinde tablodaki 10 bölge x 2 kapanma.
Tanılar: kompozisyon (bölge x dönem durum payları, bakış açısı, geçiş sayısı), ayni_gunler_2025_2026, sifir_nedenleri, tutarlilik (FIRMS-L2),
  gecis_kurali (Z01), bowtie, bulut_urunu_frp, gecerlilik_tespit, bulut_urunu_uyumu; meta.denetimler, meta.pilot_uyumu, meta.ppml_fe_denetim.
Çıktı: sonuclar/a16_kapsama.json; veri/firms/kapsama/hucre_gece_gecis.csv.gz (hücre x gece x geçiş), hucre_gece.csv.gz, bolge_gece.csv
Kullanım: python analiz16_kapsama.py
"""
from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

import analiz15_alev as AL
import analiz_ortak as A

KLASOR = A.VERI / "firms" / "kapsama"
PARCA = KLASOR / "parca"
JSON_CIKTI = A.SONUC / "a16_kapsama.json"
# parçaların birleşik hali (kod paketi bunları taşır; parça klasörü yoksa bunlar okunur)
BIRLESIK = {"hucre": KLASOR / "l2_hucre_granul.csv.gz", "ates": KLASOR / "l2_yangin_pikselleri.csv.gz", "granul": KLASOR / "l2_islenen_granuller.csv"}
GECERLI, BULUT = (3, 5, 7, 8, 9), (4,)  # kapsama_pilot.py ile aynı
ESIKLER = (0.80, 0.90, 0.95)
ZAMAN = ("2024-09-01", "2026-09-24")
DONEMLER = ("Savaş öncesi", "1. kapanma", "Açılış", "2. kapanma")
K = AL.K
TABLO = AL.TABLO
DURUMLAR = ("granül yok", "kapsanmadı", "bulutlu", "belirsiz", "geçerli, tespit yok", "geçerli, tespit var")
SAYIM = [f"fm_{c}" for c in range(10)] + ["cm_y1", "cm_0", "cm_1", "cm_2", "cm_3"]


def donem_sutunu(t: pd.Series) -> pd.Series:
    return pd.Series([AL.donem(x) or "geçiş" for x in t], index=t.index)


def parcalar() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Parçaları birleştirir: hücre x granül sayımları, yangın pikselleri ve işlenmiş granül listesi (boş parçalar dahil). Parça klasörü varsa
    birleşik dosyaları da yazar; yoksa (kod paketi) birleşik dosyaları okur. İki yol aynı satırları aynı sırayla verir."""
    hc = sorted(f for f in PARCA.glob("A*.csv.gz") if not f.name.endswith("_ates.csv.gz"))
    if hc:
        at = sorted(PARCA.glob("A*_ates.csv.gz"))
        with ThreadPoolExecutor(8) as h:
            R = pd.concat(list(h.map(pd.read_csv, hc)), ignore_index=True)
            P = pd.concat(list(h.map(pd.read_csv, at)), ignore_index=True)
        granuller = [f.name.split(".csv")[0] for f in hc]
        R.dropna(subset=["n_m"]).to_csv(BIRLESIK["hucre"], index=False)
        P.to_csv(BIRLESIK["ates"], index=False)
        pd.DataFrame({"anahtar": granuller}).to_csv(BIRLESIK["granul"], index=False)
    else:
        R, P = pd.read_csv(BIRLESIK["hucre"]), pd.read_csv(BIRLESIK["ates"])
        granuller = pd.read_csv(BIRLESIK["granul"])["anahtar"].tolist()
    R = R.dropna(subset=["n_m"])
    for c in ["n_m", *SAYIM]:
        R[c] = R[c].astype(np.int32)
    R["anahtar"] = R["anahtar"].astype(np.int64)
    P["anahtar"] = P["anahtar"].astype(np.int64)
    P["bowtie"] = P["bowtie"].astype(np.int8)
    return R, P, granuller


def gecis_ata(granuller: list[str]) -> pd.DataFrame:
    """Granül anahtarlarını (AYYYYDDD.HHMM) UTC gecesi içinde geçişlere ayırır: 20 dakikadan uzun boşluk yeni geçiş; kimlik ilk granül."""
    g = pd.DataFrame({"granul": sorted(set(granuller))})
    g["tarih"] = pd.to_datetime(g["granul"].str[1:8], format="%Y%j")
    g["dk"] = g["granul"].str[-4:-2].astype(int) * 60 + g["granul"].str[-2:].astype(int)
    g = g.sort_values(["tarih", "dk"]).reset_index(drop=True)
    yeni = (g["tarih"] != g["tarih"].shift()) | (g["dk"].diff() > AL.GECIS_ARALIK)
    g["gecis"] = g["granul"].where(yeni).ffill()
    return g


def hucre_gecis(R: pd.DataFrame, P: pd.DataFrame, G: pd.DataFrame) -> pd.DataFrame:
    """Hücre x gece x geçiş tablosu: sayımlar, geçerli ve bulut payı, durum, bakış açısı, yangın pikselleri (bowtie ayrı)."""
    R = R.rename(columns={"gecis": "granul"}).drop(columns="tarih").merge(G[["granul", "tarih", "gecis"]], on="granul")
    z = R["zenit"].notna()
    R["zn"] = np.where(z, R["zenit"] * R["n_m"], 0.0)
    R["nz"] = np.where(z, R["n_m"], 0)
    C = R.groupby(["tarih", "gecis", "anahtar"], as_index=False).agg(
        n_granul=("granul", "size"), n_m=("n_m", "sum"), **{c: (c, "sum") for c in SAYIM}, zn=("zn", "sum"), nz=("nz", "sum"))
    C["zenit"] = (C.pop("zn") / C.pop("nz").replace(0, np.nan)).round(2)
    n_i = 4 * C["n_m"]
    C["gecerli_pay"] = C[[f"fm_{c}" for c in GECERLI]].sum(axis=1) / n_i
    C["bulut_pay"] = C[[f"fm_{c}" for c in BULUT]].sum(axis=1) / n_i
    C["durum"] = np.select([C["gecerli_pay"] >= 0.5, C["bulut_pay"] >= 0.5], ["geçerli", "bulutlu"], "belirsiz")
    P = P.rename(columns={"gecis": "granul"}).drop(columns="tarih").merge(G[["granul", "tarih", "gecis"]], on="granul")
    P["frp_bt0"] = P["frp"].where(P["bowtie"] == 0, 0.0)
    F = P.groupby(["tarih", "gecis", "anahtar"], as_index=False).agg(
        frp=("frp_bt0", "sum"), frp_tum=("frp", "sum"), n_ates=("frp", "size"), n_bowtie=("bowtie", "sum"), vza=("vza", "median"))
    C = C.merge(F, on=["tarih", "gecis", "anahtar"], how="outer")
    # yangın pikseli hücrede, ama bu geçişte hücreye M-bant pikseli merkezi düşmedi: durumu bilinmiyor, geçerli sayılmaz
    C["durum"] = C["durum"].fillna("durum yok")
    # silmesiz durum (duyarlılık): bowtie silmesi (sınıf 1) gözlem kaybı değil, yinelenen tarama örtüşmesinin atılmasıdır; paydadan çıkarılır
    n_e = (4 * C["n_m"] - C["fm_1"]).where(lambda x: x > 0)
    C["durum_silmesiz"] = np.select([C[[f"fm_{c}" for c in GECERLI]].sum(axis=1) / n_e >= 0.5, C["fm_4"] / n_e >= 0.5], ["geçerli", "bulutlu"],
                                    "belirsiz")
    C.loc[C["n_m"] == 0, "durum_silmesiz"] = "durum yok"
    # sıkı durum (REDDEDİLDİ, yalnız belge için): yangın maskesi geçerli ve bağımsız bulut ürünü (Integer_Cloud_Mask 2-3) M-bant piksellerinin en
    # az yarısında açık. Bulut ürünü sıcak alev pikselini bulut sayıyor (bulutlu payı FRP ile birlikte artıyor: bulut_urunu_frp), bu yüzden
    # bu seçim sonuca bağlıdır ve ölçü için kullanılamaz.
    C["cm_acik"] = (C["cm_2"] + C["cm_3"]) / C["n_m"].where(C["n_m"] > 0)
    C["durum_siki"] = np.where((C["durum"] == "geçerli") & ~(C["cm_acik"] >= 0.5), "bulutlu", C["durum"])
    for c in ("frp", "frp_tum"):
        C[c] = C[c].fillna(0.0)
    for c in ("n_ates", "n_bowtie", "n_m", "n_granul", *SAYIM):
        C[c] = C[c].fillna(0).astype(np.int32)
    return C


def hucre_gece(C: pd.DataFrame, H: pd.DataFrame, granul_geceleri: pd.DatetimeIndex, durum: str = "durum") -> pd.DataFrame:
    """Hücre-gece: dört durum (ayrıntılı altı etiket), seçilen geçiş ve L2 ölçüsü; tam ızgara (bütün geceler x sabit hücre listesi).
    durum: hücre-geçiş durum sütunu ("durum" ana, "durum_siki" iki maskeli duyarlılık)."""
    C = C.copy()
    C["g"] = C[durum] == "geçerli"
    C["b"] = C[durum] == "bulutlu"
    C["k"] = C[durum] != "durum yok"
    C["t"] = C["frp"] > 0
    C["gt"] = C["g"] & C["t"]
    o = C.groupby(["tarih", "anahtar"], as_index=False).agg(n_gecis=("k", "sum"), n_gecerli=("g", "sum"), n_bulutlu=("b", "sum"),
                                                            n_tespitli=("t", "sum"), n_gecerli_tespitli=("gt", "sum"))
    V = C[C["g"]].sort_values(["zenit", "gecis"]).groupby(["tarih", "anahtar"], as_index=False).first()
    V = V[["tarih", "anahtar", "gecis", "zenit", "frp", "frp_tum", "n_ates", "n_bowtie", "gecerli_pay", "bulut_pay"]]
    izgara = pd.MultiIndex.from_product([pd.date_range(*ZAMAN), H["anahtar"].to_numpy()], names=["tarih", "anahtar"]).to_frame(index=False)
    T = izgara.merge(o, on=["tarih", "anahtar"], how="left").merge(V, on=["tarih", "anahtar"], how="left")
    for c in ("n_gecis", "n_gecerli", "n_bulutlu", "n_tespitli", "n_gecerli_tespitli"):
        T[c] = T[c].fillna(0).astype(np.int16)
    secili = T["gecis"].notna()
    T["durum"] = np.select(
        [~T["tarih"].isin(granul_geceleri), T["n_gecis"] == 0, secili & (T["frp"] > 0), secili, T["n_bulutlu"] == T["n_gecis"]],
        ["granül yok", "kapsanmadı", "geçerli, tespit var", "geçerli, tespit yok", "bulutlu"], "belirsiz")
    T["gecerli"] = secili
    return T.merge(H[["anahtar", "bolge"]], on="anahtar")


def firms_hucre(ana: pd.DataFrame, H: pd.DataFrame) -> pd.DataFrame:
    """FIRMS ana ölçüsünün hücre-gece değeri (analiz15_alev.hucre_gece, track kuralı) ve her FIRMS geçişinin zamanı."""
    FH = AL.hucre_gece(ana, H)
    FH["anahtar"] = FH["hx"].astype(np.int64) * 100000 + FH["hy"].astype(np.int64)
    return FH[["tarih", "anahtar", "frp"]].rename(columns={"frp": "frp_firms"})


def bolge_gece(T: pd.DataFrame, S: pd.DataFrame) -> pd.DataFrame:
    """Bölge-gece: kapsama, durum payları, L2 ölçüsü (geçerli hücrelerin toplamı), FIRMS ana ölçüsü ve FIRMS'in geçerli hücrelerdeki toplamı."""
    T = T.assign(l2=T["frp"].where(T["gecerli"]), l2bt=T["frp_tum"].where(T["gecerli"]), firms_g=T["frp_firms"].where(T["gecerli"]))
    g = T.groupby(["bolge", "tarih"])
    B = g.agg(n_hucre=("anahtar", "size"), kapsama=("gecerli", "mean"), zenit_medyan=("zenit", "median")).reset_index()
    for c in ("l2", "l2bt", "firms_g"):
        B[c] = g[c].sum(min_count=1).to_numpy()
    pay = pd.crosstab([T["bolge"], T["tarih"]], T["durum"], normalize="index").reindex(columns=list(DURUMLAR), fill_value=0.0)
    B = B.merge(pay.reset_index(), on=["bolge", "tarih"])
    f = S.stack(future_stack=True).rename("firms").reset_index()
    f.columns = ["tarih", "bolge", "firms"]
    B = B.merge(f, on=["tarih", "bolge"], how="left")
    B["donem"] = donem_sutunu(B["tarih"])
    return B


def seri(B: pd.DataFrame, b: str, sutun: str, esik: float | None = None, bolgeler: tuple[str, ...] | None = None) -> pd.Series:
    """Bölgenin gece serisi; esik verilirse kapsaması eşiğin altındaki geceler eksik (bolgeler verilirse hepsinin kapsaması eşiği geçmeli)."""
    x = B[B["bolge"] == b].set_index("tarih")
    y = x[sutun].astype(float).copy()
    if esik is not None:
        for c in bolgeler or (b,):
            kp = B[B["bolge"] == c].set_index("tarih")["kapsama"].reindex(y.index)
            y[~(kp >= esik)] = np.nan
    return y.reindex(pd.date_range(*ZAMAN))


def ppml_fe(y: np.ndarray, hucre: np.ndarray, zaman: np.ndarray, X: np.ndarray, gecikme: int = 14) -> dict:
    """Hücre sabit etkili PPML. Hücre etkisi kapalı biçimde yoğunlaştırılır (alfa_i = log toplam y_i - log toplam exp(x b)); Newton.
    Varyans: yoğunlaştırılmış skorlar (x - hücre içi mu ağırlıklı ortalama)(y - mu) gece içinde toplanır, Bartlett HAC (Driscoll-Kraay),
    küme düzeltmesi T/(T-1) x (n-1)/(n-k), k = katsayı + hücre sayısı (statsmodels hac-groupsum varsayılanı). Bütün y'si sıfır hücre düşer."""
    toplam = np.bincount(hucre, weights=y)
    dusen = int(((toplam == 0) & (np.bincount(hucre) > 0)).sum())
    tut = toplam[hucre] > 0
    y, X = y[tut], X[tut]
    hucre = np.unique(hucre[tut], return_inverse=True)[1]
    zaman = np.unique(zaman[tut], return_inverse=True)[1]
    Gn, Tn, n, k = hucre.max() + 1, zaman.max() + 1, len(y), X.shape[1]
    sy = np.bincount(hucre, weights=y, minlength=Gn)

    def durum(b):
        eta = X @ b
        m = eta.max()
        alfa = np.log(sy) - np.log(np.bincount(hucre, weights=np.exp(eta - m), minlength=Gn)) - m
        mu = np.exp(alfa[hucre] + eta)
        smu = np.bincount(hucre, weights=mu, minlength=Gn)
        xbar = np.column_stack([np.bincount(hucre, weights=mu * X[:, j], minlength=Gn) for j in range(k)]) / smu[:, None]
        Xt = X - xbar[hucre]
        return mu, Xt, (Xt * mu[:, None]).T @ Xt

    b = np.zeros(k)
    for _ in range(200):
        mu, Xt, Am = durum(b)
        adim = np.linalg.solve(Am, Xt.T @ (y - mu))
        b = b + adim
        if np.max(np.abs(adim)) < 1e-11:
            break
    mu, Xt, Am = durum(b)
    u = Xt * (y - mu)[:, None]
    s = np.column_stack([np.bincount(zaman, weights=u[:, j], minlength=Tn) for j in range(k)])
    Sm = s.T @ s
    for lag in range(1, gecikme + 1):
        Gm = s[lag:].T @ s[:-lag]
        Sm += (1 - lag / (gecikme + 1)) * (Gm + Gm.T)
    Ai = np.linalg.inv(Am)
    V = Ai @ Sm @ Ai * Tn / (Tn - 1) * (n - 1) / (n - k - Gn)
    return {"b": b, "se": np.sqrt(np.diag(V)), "V": V, "n": int(n), "hucre": int(Gn), "gece": int(Tn), "dusen_hucre": dusen}


def hucre_model(T: pd.DataFrame, b: str, sutun: str) -> dict | None:
    """Bölgenin geçerli hücre-gecelerinde hücre sabit etkili PPML; dönem ve takvim ayı kuklaları (Ocak ve savaş öncesi taban)."""
    x = T[(T["bolge"] == b) & T["gecerli"]].copy()
    x["donem"] = donem_sutunu(x["tarih"])
    x = x[x["donem"] != "geçiş"]
    if x[sutun].sum() <= 0 or x.loc[x["donem"] == "Savaş öncesi", "tarih"].nunique() < 30:
        return None
    sifir = [k for k in K if (x["donem"] == k).any() and x.loc[x["donem"] == k, sutun].sum() == 0]
    x = x[~x["donem"].isin(sifir)]
    X = pd.get_dummies(x["donem"]).reindex(columns=["1. kapanma", "Açılış", "2. kapanma"], fill_value=False).astype(float)
    X = X.join(pd.get_dummies(x["tarih"].dt.month, prefix="ay").astype(float).drop(columns="ay_1", errors="ignore"))
    X = X.loc[:, X.sum() > 0]
    gece = pd.Index(sorted(x["tarih"].unique()))
    r = ppml_fe(x[sutun].to_numpy(float), pd.factorize(x["anahtar"])[0], gece.get_indexer(x["tarih"]), X.to_numpy(), 14)
    out = {}
    for j, c in enumerate(X.columns):
        if c in AL.DON:
            bb, se = r["b"][j], r["se"][j]
            out[c] = {"etki_yuzde": round(100 * (np.exp(bb) - 1), 3), "alt": round(100 * (np.exp(bb - 1.96 * se) - 1), 3),
                      "ust": round(100 * (np.exp(bb + 1.96 * se) - 1), 3), "p": float(2 * stats.norm.sf(abs(bb / se)))}
    for k in sifir:
        out[k] = {"etki_yuzde": -100.0, "alt": None, "ust": None, "p": None, "not": "dönemde bütün hücre-geceler sıfır"}
    out.update({"n": r["n"], "hucre": r["hucre"], "gece": r["gece"], "dusen_hucre": r["dusen_hucre"],
                "n_donem": {d: int((x["donem"] == d).sum()) for d in DONEMLER},
                "gece_donem": {d: int(x.loc[x["donem"] == d, "tarih"].nunique()) for d in DONEMLER}})
    for k in K:
        if k in out and out["gece_donem"][k] < 5 and k not in sifir:
            out.pop(k)
    return out


def ppml_fe_denetim(T: pd.DataFrame, b: str = "Bahreyn") -> dict:
    """ppml_fe'yi statsmodels GLM (açık hücre kuklaları, hac-groupsum) ile karşılaştırır: katsayı ve standart hata aynı olmalı."""
    x = T[(T["bolge"] == b) & T["gecerli"]].copy()
    x["donem"] = donem_sutunu(x["tarih"])
    x = x[x["donem"] != "geçiş"]
    x = x[x.groupby("anahtar")["frp"].transform("sum") > 0].sort_values(["tarih", "anahtar"])
    X = pd.get_dummies(x["donem"]).reindex(columns=["1. kapanma", "Açılış", "2. kapanma"], fill_value=False).astype(float)
    X = X.join(pd.get_dummies(x["tarih"].dt.month, prefix="ay").astype(float).drop(columns="ay_1", errors="ignore"))
    X = X.loc[:, X.sum() > 0]
    gece = pd.Index(sorted(x["tarih"].unique()))
    zaman = gece.get_indexer(x["tarih"])
    r = ppml_fe(x["frp"].to_numpy(float), pd.factorize(x["anahtar"])[0], zaman, X.to_numpy(), 14)
    D = pd.get_dummies(x["anahtar"], prefix="h").astype(float)
    m = sm.GLM(x["frp"].to_numpy(float), np.column_stack([X.to_numpy(), D.to_numpy()]), family=sm.families.Poisson()).fit(
        cov_type="hac-groupsum", cov_kwds={"time": zaman, "maxlags": 14}, tol=1e-12, maxiter=200)
    k = X.shape[1]
    return {"bolge": b, "n": int(len(x)), "b_fark_en_buyuk": float(np.max(np.abs(r["b"] - m.params[:k]))),
            "se_oran_aralik": [float(np.min(r["se"] / m.bse[:k])), float(np.max(r["se"] / m.bse[:k]))]}


def reg(y: pd.Series) -> dict | None:
    """analiz15_alev.regresyon (ana model). Savaş öncesinde en az 30 gece ve kapanmalardan birinde en az 5 gece yoksa None; 5 geceden az
    gözlenen kapanma döneminin etkisi raporlanmaz (n_donem'de görünür). Bir dönemde bütün değerler sıfırsa etki -100 ve aralık yok."""
    y = y.dropna()
    d = pd.Series([AL.donem(t) for t in y.index], index=y.index)
    if (d == "Savaş öncesi").sum() < 30 or not any((d == k).sum() >= 5 for k in K):
        return None
    sifir = [k for k in K if (d == k).sum() >= 5 and y[d == k].sum() == 0]
    yy = y.where(~d.isin(sifir)).dropna() if sifir else y
    X = AL.tasarim(yy)[2].to_numpy()
    rk = np.linalg.matrix_rank(X)
    belirsiz = [] if rk == X.shape[1] else [k for k in K if k in AL.tasarim(yy)[2] and
                                             np.linalg.matrix_rank(np.delete(X, list(AL.tasarim(yy)[2].columns).index(k), axis=1)) == rk]
    r = AL.regresyon(yy)
    r["n_donem"] = {k: int((d == k).sum()) for k in DONEMLER}
    for k in K:
        if k in sifir:
            r[k] = {"etki_yuzde": -100.0, "alt": None, "ust": None, "p": None, "not": "dönemde bütün geceler sıfır"}
        elif r["n_donem"][k] < 5 or k in belirsiz:
            r.pop(k, None)
    if rk < X.shape[1]:
        r["uyari"] = f"tasarım matrisi tam ranklı değil ({rk}/{X.shape[1]}); belirlenemeyen dönem: {belirsiz or 'yok'}"
    return r


def fark_guvenli(Sf: pd.DataFrame) -> dict | None:
    """analiz15_alev.fark_testi; bir kapanma döneminde uygun gece yoksa None. Kullanılan gece sayısını ekler."""
    Y = Sf.dropna()
    d = [AL.donem(t) for t in Y.index]
    if d.count("Savaş öncesi") < 30 or not all(d.count(k) >= 5 for k in K):
        return None
    return {**AL.fark_testi(Sf), "gece": int(len(Y)), "n_donem": {d: int(sum(AL.donem(t) == d for t in Y.index)) for d in DONEMLER}}


def ozet_etki(r: dict | None) -> dict | None:
    if r is None:
        return None
    return {**{k: r[k] for k in K if k in r}, "n": r.get("n"), "n_donem": r.get("n_donem"), **({"gece_donem": r["gece_donem"]} if "gece_donem" in r else {}),
            **({"uyari": r["uyari"]} if "uyari" in r else {}),
            **({k: r[k] for k in ("hucre", "gece", "dusen_hucre")} if "hucre" in r else {})}


def holm(M: dict, aile: str):
    p = [(b, k, M[b][aile][k]["p"]) for b in TABLO if M[b].get(aile) for k in K if k in M[b][aile] and M[b][aile][k]["p"] is not None]
    if not p:
        return
    for (b, k, _), h in zip(p, multipletests([x[2] for x in p], method="holm")[1]):
        M[b][aile][k]["p_holm"] = float(h)


def kompozisyon(T: pd.DataFrame, B: pd.DataFrame) -> dict:
    """Bölge x dönem: hücre-gece durum payları, eşiği geçen gece payları, seçilen geçişin bakış açısı, geçiş sayısı dağılımı ve tespit olmayan
    geçiş seçimi (Z01: seçilen geçerli geçişte tespit yok, başka bir geçerli geçişte var)."""
    T = T[T["bolge"].isin(TABLO)].assign(donem=lambda d: donem_sutunu(d["tarih"]))
    B = B[B["bolge"].isin(TABLO)]
    out = {}
    for (b, d), x in T.groupby(["bolge", "donem"]):
        v = x[x["gecerli"]]
        kap = x[x["n_gecis"] > 0]
        bx = B[(B["bolge"] == b) & (B["donem"] == d)]
        cok = v[v["n_gecerli"] >= 2]
        out.setdefault(b, {})[d] = {
            "gece": int(x["tarih"].nunique()), "hucre_gece": int(len(x)),
            "durum": {s: round(float((x["durum"] == s).mean()), 4) for s in DURUMLAR},
            "kapsama_ort": round(float(bx["kapsama"].mean()), 4),
            **{f"gece_kapsama_{int(q * 100)}": round(float((bx["kapsama"] >= q).mean()), 4) for q in ESIKLER},
            "zenit_medyan": None if v.empty else round(float(v["zenit"].median()), 2),
            "zenit_c25_c75": None if v.empty else [round(float(v["zenit"].quantile(0.25)), 2), round(float(v["zenit"].quantile(0.75)), 2)],
            "gecis_sayisi_payi": {str(n): round(float(((kap["n_gecis"].clip(upper=3)) == n).mean()), 4) for n in (1, 2, 3)} if len(kap) else None,
            "gecerli_gecis_ort": round(float(v["n_gecerli"].mean()), 3) if len(v) else None,
            "coklu_gecerli_payi": round(float(len(cok) / len(v)), 4) if len(v) else None,
            "secilen_tespitsiz_baskasinda_tespit": round(float(((cok["frp"] == 0) & (cok["n_gecerli_tespitli"] > 0)).mean()), 4) if len(cok) else None,
        }
    return out


def ayni_gunler(T: pd.DataFrame, B: pd.DataFrame) -> dict:
    """Kapanma pencerelerinde 2026 ile 2025'in aynı günleri: bulutlu hücre-gece payı ve kapsaması en az %90 olan gece payı. Takvim ayı kuklaları
    2025'e göre olağan mevsim bulutunu ayırır; 2026'daki fazla bulutu ayıramaz."""
    out = {}
    for k in K:
        a, b = AL.DON[k]
        for b_ in TABLO:
            for y in (2025, 2026):
                bas, son = f"{y}{a[4:]}", f"{y}{b[4:]}"
                t = T[(T["bolge"] == b_) & T["tarih"].between(bas, son)]
                bb = B[(B["bolge"] == b_) & B["tarih"].between(bas, son)]
                out.setdefault(k, {}).setdefault(b_, {})[str(y)] = {"bulutlu": round(float((t["durum"] == "bulutlu").mean()), 4),
                                                                     "gece_kapsama_90": round(float((bb["kapsama"] >= 0.9).mean()), 4),
                                                                     "gece": int(bb["tarih"].nunique())}
    return out


def sifir_nedenleri(T: pd.DataFrame, B: pd.DataFrame) -> dict:
    """FIRMS sıfırlarının nedenleri. Bölge-gece: FIRMS'in bölgede tespit görmediği (ama kutuda verisi olan) gecelerde kapsama. Hücre-gece:
    FIRMS'te tespit olmayan kalıcı hücre-gecelerin L2 durumlarına dağılımı."""
    out = {"bolge_gece": {}, "hucre_gece": {}}
    for (b, d), x in B[B["bolge"].isin(TABLO) & (B["donem"] != "geçiş")].groupby(["bolge", "donem"]):
        z = x[x["firms"] == 0]
        out["bolge_gece"].setdefault(b, {})[d] = {
            "firms_sifir_gece": int(len(z)), "firms_gece": int(x["firms"].notna().sum()),
            "kapsama_medyan": None if z.empty else round(float(z["kapsama"].median()), 3),
            "kapsama_80_ustu": int((z["kapsama"] >= 0.8).sum()),
            "kapsama_50_alti": int((z["kapsama"] < 0.5).sum()),
            "bulutlu_min": None if z.empty else round(float(z["bulutlu"].min()), 3),
            "bulutlu_medyan": None if z.empty else round(float(z["bulutlu"].median()), 3),
            "l2_tespit_var": int((z["l2"] > 0).sum())}
    fg = set(B.loc[B["firms"].notna(), "tarih"])
    Tz = T[T["bolge"].isin(TABLO) & T["tarih"].isin(fg) & (T["frp_firms"] == 0)].assign(donem=lambda d: donem_sutunu(d["tarih"]))
    for (b, d), x in Tz[Tz["donem"] != "geçiş"].groupby(["bolge", "donem"]):
        out["hucre_gece"].setdefault(b, {})[d] = {"n": int(len(x)), **{s: round(float((x["durum"] == s).mean()), 4) for s in DURUMLAR}}
    return out


def tutarlilik(T: pd.DataFrame, B: pd.DataFrame, P: pd.DataFrame, ana: pd.DataFrame, H: pd.DataFrame) -> dict:
    """FIRMS (NRT) ile L2 (standart işleme) karşılaştırması: bölge-gece korelasyonu ve dönemlere göre düzey oranı (%90 kapsamalı geceler),
    kalıcı hücrelerdeki piksel sayıları; hücre-gece düzeyinde tespit uyumu."""
    out = {"bolge": {}}
    for b in TABLO:
        x = B[(B["bolge"] == b) & (B["kapsama"] >= 0.9) & B["firms"].notna() & (B["donem"] != "geçiş")]
        if len(x) < 10:
            continue
        out["bolge"][b] = {"gece": int(len(x)), "pearson_firms_l2": round(float(x[["firms", "l2"]].corr().iloc[0, 1]), 4),
                           "pearson_firmsgecerli_l2": round(float(x[["firms_g", "l2"]].corr().iloc[0, 1]), 4),
                           "oran_l2_firms": {d: (round(float(y["l2"].sum() / y["firms"].sum()), 4) if y["firms"].sum() > 0 else None)
                                             for d, y in x.groupby("donem")}}
    Hh = H[H["bolge"].isin(TABLO)]
    anahtar = set(Hh["anahtar"])
    Pk = P[P["anahtar"].isin(anahtar)]
    Fa = ana.assign(anahtar=ana["hx"].astype(np.int64) * 100000 + ana["hy"].astype(np.int64))
    Fa = Fa[Fa["anahtar"].isin(anahtar)]
    ortak = sorted(set(pd.to_datetime(Pk["tarih"])) & set(Fa["tarih"]))
    pl = Pk.assign(tarih=pd.to_datetime(Pk["tarih"])).groupby("tarih").size().reindex(ortak, fill_value=0)
    fl = Fa.groupby("tarih").size().reindex(ortak, fill_value=0)
    dd = pd.Series([AL.donem(t) or "geçiş" for t in ortak], index=ortak)
    out["piksel_sayisi"] = {"gece": len(ortak), "l2": int(pl.sum()), "firms": int(fl.sum()), "oran_l2_firms": round(float(pl.sum() / fl.sum()), 4),
                            "donem_oran": {d: round(float(pl[dd == d].sum() / fl[dd == d].sum()), 4) for d in DONEMLER if fl[dd == d].sum() > 0},
                            "gece_pearson": round(float(np.corrcoef(pl, fl)[0, 1]), 4)}
    v = T[T["gecerli"] & T["bolge"].isin(TABLO) & T["tarih"].isin(set(fl.index))]
    out["hucre_gece_tespit_uyumu"] = {"gecerli_hucre_gece": int(len(v)),
                                      "ikisinde_tespit": round(float(((v["frp"] > 0) & (v["frp_firms"] > 0)).mean()), 4),
                                      "yalniz_l2": round(float(((v["frp"] > 0) & (v["frp_firms"] == 0)).mean()), 4),
                                      "yalniz_firms": round(float(((v["frp"] == 0) & (v["frp_firms"] > 0)).mean()), 4),
                                      "ikisinde_yok": round(float(((v["frp"] == 0) & (v["frp_firms"] == 0)).mean()), 4)}
    return out


def gecis_kurali(C: pd.DataFrame, ana: pd.DataFrame, H: pd.DataFrame, G: pd.DataFrame) -> dict:
    """Z01 izlemesi: FIRMS'in track kuralının seçtiği geçiş, L2'de (bulut ürününün konum dosyasından) bakış açısı en küçük geçiş mi?
    (a) yalnız FIRMS'te tespitli geçişler arasında, (b) hücrenin geçerli olduğu bütün geçişler arasında (tespit olmayanlar dahil)."""
    Hh = H[H["bolge"].isin(TABLO)]
    Kk = ana.merge(Hh[["hx", "hy"]], on=["hx", "hy"])
    Kk["anahtar"] = Kk["hx"].astype(np.int64) * 100000 + Kk["hy"].astype(np.int64)
    cg = Kk.groupby(["tarih", "anahtar", "gecis"], as_index=False).agg(track=("track", "mean"), scan=("scan", "mean"), dk=("dk", "min"))
    cg = cg[cg.groupby(["tarih", "anahtar"])["gecis"].transform("size") > 1]
    # FIRMS geçişini L2 geçişine eşle: aynı UTC gecesinde başlangıcı FIRMS zamanından en çok 20 dakika önce olan en yakın L2 geçişi
    L = G.groupby("gecis", as_index=False).agg(tarih=("tarih", "first"), bas=("dk", "min"), son=("dk", "max"))
    e = cg.merge(L, on="tarih", how="inner")
    e = e[(e["dk"] >= e["bas"] - 2) & (e["dk"] <= e["son"] + 8)]
    e = e.sort_values("bas").groupby(["tarih", "anahtar", "gecis_x"], as_index=False).last().rename(columns={"gecis_x": "fgecis", "gecis_y": "gecis"})
    e = e.merge(C[["tarih", "gecis", "anahtar", "zenit", "durum"]], on=["tarih", "gecis", "anahtar"], how="left")
    e = e[e.groupby(["tarih", "anahtar"])["fgecis"].transform("size") > 1].dropna(subset=["zenit"])
    e = e[e.groupby(["tarih", "anahtar"])["fgecis"].transform("size") > 1]
    sec_t = e.sort_values(["track", "scan", "fgecis"]).groupby(["tarih", "anahtar"])["gecis"].first()
    sec_z = e.sort_values(["zenit", "gecis"]).groupby(["tarih", "anahtar"])["gecis"].first()
    Cg = C[C["durum"] == "geçerli"]
    sec_g = Cg.sort_values(["zenit", "gecis"]).groupby(["tarih", "anahtar"])["gecis"].first()
    ortak = sec_t.index.intersection(sec_g.index)
    return {"cok_gecisli_hucre_gece": int(len(sec_t)),
            "track_tespitli_gecisler_arasinda_en_kucuk_aciyi_secti": round(float((sec_t == sec_z.reindex(sec_t.index)).mean()), 4),
            "karsilastirilan_gecerli": int(len(ortak)),
            "track_gecerli_gecisler_arasinda_en_kucuk_aciyi_secti": round(float((sec_t.loc[ortak] == sec_g.loc[ortak]).mean()), 4)}


def bowtie(T: pd.DataFrame, P: pd.DataFrame, H: pd.DataFrame) -> dict:
    Pk = P.merge(H[H["bolge"].isin(TABLO)][["anahtar", "bolge"]], on="anahtar")
    Pk["donem"] = donem_sutunu(pd.to_datetime(Pk["tarih"]))
    pay = lambda x: round(float(x.loc[x["bowtie"] == 1, "frp"].sum() / x["frp"].sum()), 4) if x["frp"].sum() > 0 else None
    return {"piksel": int(len(Pk)), "bowtie_piksel_payi": round(float(Pk["bowtie"].mean()), 4), "bowtie_frp_payi": pay(Pk),
            "donem": {d: pay(x) for d, x in Pk.groupby("donem") if d != "geçiş"},
            "bolge_donem": {b: {d: pay(y) for d, y in x.groupby("donem") if d != "geçiş"} for b, x in Pk.groupby("bolge")}}


def bulut_urunu_frp(C: pd.DataFrame, Ht: pd.DataFrame) -> dict:
    """Bulut ürününün 'bulutlu' payı, yangın maskesinin geçerli dediği hücre-geçişlerde, tespit edilen FRP dilimine göre. Artıyorsa bulut ürünü
    alevleri bulut sayıyor demektir: onunla seçim yapmak sonuca bağlı olur."""
    x = C[C["anahtar"].isin(set(Ht["anahtar"])) & (C["durum"] == "geçerli") & (C["n_m"] > 0)]
    bul = (x["cm_0"] + x["cm_1"]) / x["n_m"] >= 0.5
    dilim = pd.cut(x["frp"], [-1, 0, 2, 5, 10, 20, 50, np.inf], labels=["0", "0-2", "2-5", "5-10", "10-20", "20-50", "50+"])
    return {str(k): {"bulutlu_pay": round(float(v.mean()), 4), "n": int(len(v))} for k, v in bul.groupby(dilim, observed=True)}


def gecerlilik_tespit(C: pd.DataFrame, Ht: pd.DataFrame) -> dict:
    """Geçerlilik kurallarının tespitle ilişkisi (hücre-geçiş): ana kural, silmesiz kural ve yangın piksellerini paydadan çıkaran kural
    (ikincisi alevli hücrenin aleyhine: alev pikseli açık bir pikseldir). Sınıf payları tespitsiz ve tespitli hücre-geçişlerde."""
    x = C[C["anahtar"].isin(set(Ht["anahtar"])) & (C["n_m"] > 0)]
    t = x["frp"] > 0
    ates = x["fm_7"] + x["fm_8"] + x["fm_9"]
    n_a = (4 * x["n_m"] - ates).where(lambda v: v > 0)
    yangin_disi = np.where((x["fm_3"] + x["fm_5"]) / n_a >= 0.5, "geçerli", "değil")
    n_i = 4 * x["n_m"]
    return {"gecerli_pay": {ad: {"tespitsiz": round(float((d[~t] == "geçerli").mean()), 4), "tespitli": round(float((d[t] == "geçerli").mean()), 4)}
                            for ad, d in (("ana", x["durum"]), ("silmesiz", x["durum_silmesiz"]), ("yangin_payda_disi", pd.Series(yangin_disi, index=x.index)))},
            "sinif_payi": {f"fm_{c}": {"tespitsiz": round(float(x.loc[~t, f"fm_{c}"].sum() / n_i[~t].sum()), 4),
                                       "tespitli": round(float(x.loc[t, f"fm_{c}"].sum() / n_i[t].sum()), 4)} for c in range(10)}}


def bulut_urunu_uyumu(C: pd.DataFrame) -> dict:
    """Yangın maskesinin bulut sınıfı ile bağımsız bulut ürününün (Integer_Cloud_Mask 0-1 bulutlu) uyumu, hücre-geçiş düzeyinde."""
    x = C[C["n_m"] > 0]
    cm_bulut = (x["cm_0"] + x["cm_1"]) / x["n_m"] >= 0.5
    fm_bulut = x["bulut_pay"] >= 0.5
    return {"hucre_gecis": int(len(x)), "ayni_karar": round(float((cm_bulut == fm_bulut).mean()), 4),
            "fm_bulutlu_cm_acik": round(float((fm_bulut & ~cm_bulut).mean()), 4), "fm_acik_cm_bulutlu": round(float((~fm_bulut & cm_bulut).mean()), 4),
            "pearson_bulut_payi": round(float(np.corrcoef(x["bulut_pay"], (x["cm_0"] + x["cm_1"]) / x["n_m"])[0, 1]), 4)}


def denetimler(R: pd.DataFrame, C: pd.DataFrame, granuller: list[str]) -> dict:
    """İç tutarlılık: her satırda I-bant sınıf sayımlarının toplamı 4 x M-bant pikseli, bulut ürünü sınıflarının toplamı M-bant pikseli olmalı;
    parça listesi granül listesiyle aynı olmalı."""
    Gl = pd.read_csv(KLASOR / "granuller.csv")
    fm = R[[f"fm_{c}" for c in range(10)]].sum(axis=1)
    cm = R[["cm_y1", "cm_0", "cm_1", "cm_2", "cm_3"]].sum(axis=1)
    hata = KLASOR / "hatalar.csv"
    return {"granul_listesi": int(len(Gl)), "islenen_granul_cifti": len(granuller),
            "listede_olup_islenmeyen": sorted(set(Gl["anahtar"]) - set(granuller)),
            "fm_toplam_4nm_disi_satir": int((fm != 4 * R["n_m"]).sum()), "cm_toplam_nm_disi_satir": int((cm != R["n_m"]).sum()),
            "hucre_granul_satiri": int(len(R)), "zenit_eksik_satir": int(R["zenit"].isna().sum()),
            "durumu_olmayan_yangin_hucre_gecisi": int((C["durum"] == "durum yok").sum()),
            "durumu_olmayan_frp_payi": round(float(C.loc[C["durum"] == "durum yok", "frp_tum"].sum() / C["frp_tum"].sum()), 5),
            "durumu_olmayan_frp_payi_donem": {d: round(float(x.loc[x["durum"] == "durum yok", "frp_tum"].sum() / x["frp_tum"].sum()), 5)
                                              for d, x in C.assign(donem=donem_sutunu(C["tarih"])).groupby("donem") if x["frp_tum"].sum() > 0},
            "hata_kaydi": int(len(pd.read_csv(hata))) if hata.exists() else 0,
            "esi_olmayan_granul": pd.read_csv(KLASOR / "esi_olmayan_granuller.csv")["anahtar"].tolist()}


def pilot_uyumu(R: pd.DataFrame, P: pd.DataFrame, H: pd.DataFrame) -> dict:
    """Tam çalıştırma, pilotun (kapsama_pilot.py, 59 gece, 6 bölge) granüllerinde pilotla aynı satırları üretmeli: aynı granül ve hücreler için
    bütün sayımlar, bakış açısı ve yangın pikselleri birebir karşılaştırılır."""
    d = A.VERI / "firms" / "kapsama_pilot"
    if not (d / "hucre_gecis.csv.gz").exists():
        return {"pilot": "yok"}
    Rp = pd.read_csv(d / "hucre_gecis.csv.gz")
    Pp = pd.read_csv(d / "yangin_pikselleri.csv.gz")
    hucre = set(H.loc[H["bolge"].isin(["Irak (güney)", "Irak (kuzey)", "Katar", "Kuveyt", "Bahreyn", "Batı İran (çoğu Huzistan)"]), "anahtar"])
    ortak = set(Rp["gecis"]) & set(R["gecis"])
    sut = ["gecis", "anahtar", "n_m", *SAYIM, "zenit"]
    a = Rp[Rp["gecis"].isin(ortak)][sut].sort_values(["gecis", "anahtar"]).reset_index(drop=True)
    b = R[R["gecis"].isin(ortak) & R["anahtar"].isin(hucre)][sut].sort_values(["gecis", "anahtar"]).reset_index(drop=True)
    ps = ["gecis", "anahtar", "frp", "vza", "guven", "bowtie"]
    pa = Pp[Pp["gecis"].isin(ortak)][ps].sort_values(ps).reset_index(drop=True)
    pb = P[P["gecis"].isin(ortak) & P["anahtar"].isin(hucre)][ps].sort_values(ps).reset_index(drop=True)
    ayni = len(a) == len(b) and bool((a[sut[:-1]].to_numpy() == b[sut[:-1]].to_numpy()).all()) and bool(np.allclose(a["zenit"], b["zenit"], equal_nan=True, rtol=0, atol=0))
    ayni_p = len(pa) == len(pb) and bool((pa[["gecis", "anahtar", "guven", "bowtie"]].to_numpy() == pb[["gecis", "anahtar", "guven", "bowtie"]].to_numpy()).all()) \
        and bool(np.array_equal(pa[["frp", "vza"]].to_numpy(), pb[["frp", "vza"]].to_numpy()))
    return {"pilot_granul": int(Rp["gecis"].nunique()), "ortak_granul": len(ortak), "hucre_satir": [int(len(a)), int(len(b))],
            "yangin_piksel": [int(len(pa)), int(len(pb))], "hucre_birebir": ayni, "yangin_birebir": ayni_p}


def main():
    t0 = time.time()
    H = pd.read_csv(A.VERI / "firms" / "kalici_hucreler.csv")
    H["anahtar"] = H["hx"].astype(np.int64) * 100000 + H["hy"].astype(np.int64)
    Ht = H[H["bolge"].isin(TABLO)].copy()
    ana = AL.yukle(["VIIRS_NOAA21_NRT"])
    S = AL.gece_serisi(ana, H)
    R, P, granuller = parcalar()
    print(f"parçalar: {len(granuller)} granül çifti, {len(R)} hücre-granül satırı, {len(P)} yangın pikseli | {time.time() - t0:.0f} sn", flush=True)
    G = gecis_ata(granuller)
    C = hucre_gecis(R, P, G)
    T = hucre_gece(C, Ht, pd.DatetimeIndex(G["tarih"].unique()))
    T = T.merge(firms_hucre(ana, H), on=["tarih", "anahtar"], how="left")
    T["frp_firms"] = T["frp_firms"].fillna(0.0)
    B = bolge_gece(T, S[[b for b in TABLO if b in S]])
    Tb = hucre_gece(C, Ht, pd.DatetimeIndex(G["tarih"].unique()), "durum_silmesiz")
    Tb = Tb.merge(firms_hucre(ana, H), on=["tarih", "anahtar"], how="left")
    Tb["frp_firms"] = Tb["frp_firms"].fillna(0.0)
    Bb = bolge_gece(Tb, S[[b for b in TABLO if b in S]])
    Ts = hucre_gece(C, Ht, pd.DatetimeIndex(G["tarih"].unique()), "durum_siki")
    Ts = Ts.merge(firms_hucre(ana, H), on=["tarih", "anahtar"], how="left")
    Ts["frp_firms"] = Ts["frp_firms"].fillna(0.0)
    Bs = bolge_gece(Ts, S[[b for b in TABLO if b in S]])
    print(f"tablolar: {len(C)} hücre-geçiş, {len(T)} hücre-gece | {time.time() - t0:.0f} sn", flush=True)

    a15 = json.loads((A.SONUC / "a15_alev.json").read_text(encoding="utf-8"))
    M = {}
    for b in TABLO:
        m = {"firms_tum": ozet_etki(reg(seri(B, b, "firms")))}
        for q in ESIKLER:
            e = int(q * 100)
            m[f"firms_k{e}"] = ozet_etki(reg(seri(B, b, "firms", q)))
            m[f"l2_k{e}"] = ozet_etki(reg(seri(B, b, "l2", q)))
        for q in ESIKLER:
            m[f"l2silmesiz_k{int(q * 100)}"] = ozet_etki(reg(seri(Bb, b, "l2", q)))
        for q in ESIKLER:
            m[f"l2siki_k{int(q * 100)}"] = ozet_etki(reg(seri(Bs, b, "l2", q)))
        m["firmsgecerli_k90"] = ozet_etki(reg(seri(B, b, "firms_g", 0.9)))
        m["l2bt_k90"] = ozet_etki(reg(seri(B, b, "l2bt", 0.9)))
        m["hucre_l2"] = hucre_model(T, b, "frp")
        m["hucre_firms"] = hucre_model(T, b, "frp_firms")
        m["hucre_l2silmesiz"] = hucre_model(Tb, b, "frp")
        m["hucre_l2siki"] = hucre_model(Ts, b, "frp")
        tum = m["firms_tum"]["n_donem"]
        for ad, r in m.items():
            if r and ad != "firms_tum" and "n_donem" in r:
                r["dislanan_gece_payi"] = {d: round(1 - r["n_donem"][d] / tum[d], 4) if tum[d] else None for d in DONEMLER} if "hucre" not in r else None
        M[b] = m
    for aile in M[TABLO[0]]:
        holm(M, aile)
    yeniden = {b: max(abs(M[b]["firms_tum"][k]["etki_yuzde"] - a15["bolgeler"][b]["etki"][k]["etki_yuzde"]) for k in K) for b in TABLO}
    fark = {}
    for q in (None, *ESIKLER):
        ad = "tum" if q is None else f"k{int(q * 100)}"
        Sf = pd.DataFrame({b: seri(B, b, "firms", q, ("Irak (güney)", "Irak (kuzey)")) for b in ("Irak (güney)", "Irak (kuzey)")})
        fark[f"firms_{ad}"] = fark_guvenli(Sf)
        if q is not None:
            Sl = pd.DataFrame({b: seri(B, b, "l2", q, ("Irak (güney)", "Irak (kuzey)")) for b in ("Irak (güney)", "Irak (kuzey)")})
            fark[f"l2_{ad}"] = fark_guvenli(Sl)
            Sb = pd.DataFrame({b: seri(Bb, b, "l2", q, ("Irak (güney)", "Irak (kuzey)")) for b in ("Irak (güney)", "Irak (kuzey)")})
            fark[f"l2silmesiz_{ad}"] = fark_guvenli(Sb)
            Ss = pd.DataFrame({b: seri(Bs, b, "l2", q, ("Irak (güney)", "Irak (kuzey)")) for b in ("Irak (güney)", "Irak (kuzey)")})
            fark[f"l2siki_{ad}"] = fark_guvenli(Ss)
    print(f"modeller | {time.time() - t0:.0f} sn", flush=True)

    out = {"meta": {"tanim": __doc__.split("Modeller")[0].strip(),
                    "donem": {"Savaş öncesi": ["2024-09-01", str(AL.ONCE_SON.date())], **{k: list(v) for k, v in AL.DON.items()}},
                    "hucre": {b: int((Ht["bolge"] == b).sum()) for b in TABLO}, "gece": int(len(pd.date_range(*ZAMAN))),
                    "granul_gecesi": int(G["tarih"].nunique()), "gecis": int(G["gecis"].nunique()),
                    "gece_basina_gecis": {str(k): int(v) for k, v in G.groupby("tarih")["gecis"].nunique().value_counts().sort_index().items()},
                    "granul_yok_geceler": [str(t.date()) for t in pd.date_range(*ZAMAN) if t not in set(G["tarih"])],
                    "firms_yok_geceler": [str(t.date()) for t in S.index[S.isna().all(axis=1)]],
                    "firms_tum_a15_en_buyuk_fark_puan": yeniden,
                    "ppml_fe_denetim": [ppml_fe_denetim(T, b) for b in ("Bahreyn", "Katar", "Irak (kuzey)", "Irak (güney)")],
                    "denetimler": denetimler(R, C, granuller),
                    "pilot_uyumu": pilot_uyumu(R, P, H)},
           "kompozisyon": kompozisyon(T, B), "kompozisyon_silmesiz": kompozisyon(Tb, Bb), "ayni_gunler_2025_2026": ayni_gunler(T, B),
           "sifir_nedenleri": sifir_nedenleri(T, B),
           "modeller": M, "fark_guney_kuzey": fark,
           "siki_uyari": "l2siki_* ve hucre_l2siki reddedildi: bağımsız bulut ürünü alev piksellerini bulut sayıyor (bulut_urunu_frp); "
                         "bu seçim sonuca bağlı, ölçü için kullanılamaz. Yalnız belge için tutuldu.",
           "bulut_urunu_frp": bulut_urunu_frp(C, Ht), "gecerlilik_tespit": gecerlilik_tespit(C, Ht),
           "tutarlilik": tutarlilik(T, B, P, ana, H), "gecis_kurali": gecis_kurali(C, ana, H, G), "bowtie": bowtie(T, P, H),
           "bulut_urunu_uyumu": bulut_urunu_uyumu(C)}
    JSON_CIKTI.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
    C.drop(columns=["cm_y1"]).to_csv(KLASOR / "hucre_gece_gecis.csv.gz", index=False, float_format="%.4g")
    T.drop(columns=["gecerli_pay", "bulut_pay"]).to_csv(KLASOR / "hucre_gece.csv.gz", index=False, float_format="%.4g")
    B.to_csv(KLASOR / "bolge_gece.csv", index=False, float_format="%.5g")
    print(f"bitti | {time.time() - t0:.0f} sn", flush=True)
    pd.set_option("display.width", 250)
    tab = {}
    for b in TABLO:
        tab[b] = {f"{ad} {k[:1]}": (f"{M[b][ad][k]['etki_yuzde']:+.0f}" if M[b].get(ad) and k in M[b][ad] else "-")
                  for ad in ("firms_tum", "firms_k90", "l2_k80", "l2_k90", "l2_k95", "l2silmesiz_k90", "hucre_l2", "hucre_firms") for k in K}
    print(pd.DataFrame(tab).T.to_string())
    print("güney-kuzey:", {k: ({d: round(v[d]["etki_yuzde"], 1) for d in K} if v else None) for k, v in fark.items()})
    print("a15 yeniden üretim, en büyük fark (puan):", yeniden)
    print("ppml_fe denetimi:", out["meta"]["ppml_fe_denetim"])
    print("pilot uyumu:", out["meta"]["pilot_uyumu"])


if __name__ == "__main__":
    main()
