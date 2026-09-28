"""Analiz 7: Sentinel-1 radar gemi sayımlarının özeti (aylık medyan, dönem medyanları ve dönem etkileri, yörünge yönüne göre).
Girdi: veri/sar/gemi_sayimi_<konum><ek>.csv (sar_gemi_sayimi.py). Çıktı: sonuclar/a7_sar_ozet.json

Uydu kuralı: Sentinel-1A görüntülerinde deniz yüzeyinin medyanı 1C ve 1D'dekinden yaklaşık 2,5 dB yüksek; görüntüye özgü eşik (deniz
medyanı + 9 dB) yükseldiği için 1A aynı dönemde daha az gemi sayıyor (Ras Tanura savaş öncesi medyanı: 1A 80, 1C 115). Bunun gürültü
tabanından mı, işleme farkından mı geldiğini ayıramıyoruz; mekanizma iddia edilmez, uydu farkı istatistiksel olarak denetlenir. 1A savaş
öncesinde ve ilk kapanmada var, ikinci kapanmada yok; 1D Mart 2026'da başlıyor. Dönem medyanları yalnız 1C ve 1D görüntüleriyle.
Ana model (sürüm 2, 25 Eylül 2026 çapraz denetiminden sonra): gemi sayısı = dönem + uydu kuklaları, Poisson (beklenen sayıdaki yüzde
değişim), bütün görüntüler; standart hata HC1: statsmodels GLM'de "HC1" çağrısı HC0 ile aynı sonucu verdiği için HC0 alınıp n/(n-k)
çarpanı açıkça uygulanır (ikinci denetim, 25 Eylül). log(1 + y) OLS modelinde statsmodels HC1 zaten bu çarpanı içerir. Sağlamlık ("saglamlik"): log(1 + gemi) OLS (sürüm 1'in ana modeli; dönüştürülmüş
ölçek), yalnız 1C ve 1D, yalnız 1C (bütün dönemlerde gözlenen tek uydu), büyük yankılar (20 m'de en az 12, 40 m'de en az 3 piksel;
alan eşiği, gemi boyunu kanıtlamaz), Newey-West (4 görüntü gecikmesi), aya göre kümelenmiş hata (15 küme; aralık ayrıca t(G-1) ile),
ikinci kapanma görüntülerini tek tek çıkarma. Holm düzeltmesi: beş ana konum x iki kapanma (Füceyre çıkan yörünge doğrulama sayılır).
Yörünge yönleri ayrı tutulur: geliş açısı farklı olduğu için eşik ve sayım düzeyleri farklıdır.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

import analiz_ortak as A

DON = {"Tem 2025-Şub 2026": ("2025-07-01", "2026-02-27"), "1. kapanma": ("2026-03-02", "2026-06-13"),
       "Açılış": ("2026-06-18", "2026-07-10"), "2. kapanma": ("2026-07-13", "2026-09-22")}
TUTARLI = ("S1C", "S1D")
# konum: (dosya eki, kullanılacak yörünge(ler), çözünürlük notu)
KAYNAK = {"fuceyre": ("", ("descending", "ascending"), "20 m, iki yörünge"),
          "hurmuz": ("_40m", ("descending",), "40 m, inen yörünge"),
          "rastanura": ("_40m", ("descending",), "40 m, inen yörünge"),
          "yanbu": ("", ("ascending",), "20 m, çıkan yörünge"),
          "babulmendep": ("_40m_asc", ("ascending",), "40 m, çıkan yörünge")}
ANA = [("hurmuz", "descending"), ("fuceyre", "descending"), ("yanbu", "ascending"), ("babulmendep", "ascending"), ("rastanura", "descending")]
TERIM = "C(donem, Treatment('Tem 2025-Şub 2026'))"


def donem_etiketi(t):
    for d, (a, b) in DON.items():
        if pd.Timestamp(a) <= t <= pd.Timestamp(b):
            return d
    return None


def donem_medyan(g):
    return {d: {"medyan": float(g[(g.tarih >= a) & (g.tarih <= b)].gemi.median()) if len(g[(g.tarih >= a) & (g.tarih <= b)]) else None,
                "n": int(len(g[(g.tarih >= a) & (g.tarih <= b)]))} for d, (a, b) in DON.items()}


def regresyon(g, sutun="gemi", yontem="poisson"):
    """yontem: "poisson" (ana, HC1), "log" (log(1 + sayı) OLS, HC1), "HAC" (Poisson, Newey-West 4), "ay" (Poisson, aya göre küme)."""
    g = g.dropna(subset=["donem"]).sort_values("tarih").copy()
    sag = f"{TERIM}" + (" + C(uydu)" if g.uydu.nunique() > 1 else "")
    if yontem == "log":
        g["y"] = np.log1p(g[sutun])
        m = smf.ols(f"y ~ {sag}", data=g).fit(cov_type="HC1")
    else:
        # GLM'de statsmodels "HC1" çağrısı HC0 ile aynı kovaryansı veriyor (küçük örneklem çarpanı uygulanmıyor); bu yüzden ana modelde
        # HC0 alınır ve HC1 çarpanı n/(n-k) açıkça uygulanır. Newey-West ve küme kovaryansı statsmodels'ın kendi tanımıyla.
        kov = {"poisson": ("HC0", {}), "HAC": ("HAC", {"maxlags": 4}),
               "ay": ("cluster", {"groups": pd.factorize(g.tarih.dt.to_period("M"))[0]})}[yontem]
        m = smf.glm(f"{sutun} ~ {sag}", data=g, family=sm.families.Poisson()).fit(cov_type=kov[0], cov_kwds=kov[1] or None)
    carpan = np.sqrt(m.nobs / (m.nobs - len(m.params))) if yontem in ("poisson", "log") else 1.0
    G = g.tarih.dt.to_period("M").nunique()
    out = {}
    for d in DON:
        k = f"{TERIM}[T.{d}]"
        if k in m.params:
            b, s = m.params[k], m.bse[k] * (carpan if yontem == "poisson" else 1.0)
            out[d] = {"etki_yuzde": round(100 * (np.exp(b) - 1), 3), "alt": round(100 * (np.exp(b - 1.96 * s) - 1), 3),
                      "ust": round(100 * (np.exp(b + 1.96 * s) - 1), 3), "p": float(2 * stats.norm.sf(abs(b / s)))}
            if yontem == "ay":  # az küme: t(G-1) kritik değeriyle aralık
                t = stats.t.ppf(0.975, G - 1)
                out[d]["alt_t"], out[d]["ust_t"] = round(100 * (np.exp(b - t * s) - 1), 3), round(100 * (np.exp(b + t * s) - 1), 3)
    out["n"] = int(m.nobs)
    out["n_donem"] = {d: int((g.donem == d).sum()) for d in DON}
    if yontem == "ay":
        out["kume"] = int(G)
    return out


def tek_cikar(g, donem="2. kapanma"):
    """İkinci kapanma görüntülerini tek tek çıkararak ana model: nokta tahmini aralığı ve aralığı sıfırı içeren deneme sayısı."""
    g = g.dropna(subset=["donem"])
    sonuc = []
    for i in g.index[g.donem == donem]:
        r = regresyon(g.drop(i))[donem]
        sonuc.append({"tarih": str(g.loc[i, "tarih"].date()), "gemi": int(g.loc[i, "gemi"]), **r})
    return {"denemeler": sonuc, "min": min(x["etki_yuzde"] for x in sonuc), "max": max(x["etki_yuzde"] for x in sonuc),
            "sifir_iceren": int(sum(x["alt"] <= 0 <= x["ust"] for x in sonuc)), "n": len(sonuc)}


out = {}
for k, (ek, yorler, notu) in KAYNAK.items():
    f = A.VERI / "sar" / f"gemi_sayimi_{k}{ek}.csv"
    if not f.exists():
        continue
    S = pd.read_csv(f, parse_dates=["tarih"])
    S = S[S.yorunge.isin(yorler)].copy()
    S["uydu"] = S.sahne.str[:3]
    S["donem"] = S.tarih.map(donem_etiketi)
    T = S[S.uydu.isin(TUTARLI)]
    tamam = bool(len(S) and S.tarih.max() >= pd.Timestamp("2026-09-01"))
    o = {"kaynak": notu, "sahne": int(len(S)), "sahne_tutarli": int(len(T)), "uydular": "Sentinel-1C ve 1D",
         "son_tarih": S.tarih.max().strftime("%Y-%m-%d") if len(S) else None, "tamam": tamam, "yorunge": {}}
    for yor in yorler:
        g, gt = S[S.yorunge == yor], T[T.yorunge == yor]
        if gt.empty:
            continue
        m = gt.set_index("tarih").gemi.resample("ME").median().dropna()
        g1c = g[g.uydu == "S1C"]
        o["yorunge"][yor] = {"aylik": [[d.strftime("%Y-%m-15"), float(v)] for d, v in m.items()],
                            "donem": donem_medyan(gt),
                            "uydu_donem": {u: donem_medyan(gu) for u, gu in g.groupby("uydu")},
                            "deniz_medyan_db": {u: round(float(gu.deniz_medyan_db.median()), 2) for u, gu in g.groupby("uydu")},
                            "regresyon": regresyon(g),
                            "saglamlik": {"log": regresyon(g, yontem="log"), "tutarli": regresyon(gt), "yalniz_1c": regresyon(g1c),
                                          "buyuk": regresyon(g, "buyuk_gemi"), "HAC": regresyon(g, yontem="HAC"), "ay": regresyon(g, yontem="ay"),
                                          "tek_cikar": tek_cikar(g)}}
    out[k] = o
# Holm: beş ana konum x iki kapanma (ana model)
aile = [(k, y, d) for k, y in ANA for d in ("1. kapanma", "2. kapanma")]
duz = multipletests([out[k]["yorunge"][y]["regresyon"][d]["p"] for k, y, d in aile], method="holm")[1]
for (k, y, d), h in zip(aile, duz):
    out[k]["yorunge"][y]["regresyon"][d]["p_holm"] = float(h)
(A.SONUC / "a7_sar_ozet.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
for k, v in out.items():
    for yor, y in v["yorunge"].items():
        d = y["donem"]
        print(f"{k:12s} {yor:10s} tutarlı n={v['sahne_tutarli']:3d}/{v['sahne']:3d}  " +
              "  ".join(f"{a}: {d[a]['medyan']} (n={d[a]['n']})" for a in d) +
              "  | Poisson: " + ", ".join(f"{a} {r['etki_yuzde']:+.0f}% [{r['alt']:+.0f}, {r['ust']:+.0f}]"
                                          for a, r in y["regresyon"].items() if a in ("1. kapanma", "2. kapanma")))
