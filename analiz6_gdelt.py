"""Analiz 6: GDELT haber ve olay endeksleri; "tehdit ucuz, gerçekleşme pahalı" testi; Türk medyası.

Girdi: veri/gdelt_bq/haber_endeksi_v2.csv (GKG AllNames eşleşmesi), olay_endeksi_v2.csv (Events, boğaz kutuları).
Test: Brent ICE ön vade günlük getirisi (vade devri günleri çıkarılmış) ve mutlak getirisi, aynı günün tehdit, eylem ve
haber ilgisi ölçülerine regresyonla bağlanır. Ölçüler log(1 + on binde pay), son 250 işlem gününe göre standartlaştırılır
(bir gün gecikmeli pencere). Hafta sonu haberleri bir sonraki işlem gününe taşınır (en yüksek değer). HAC (5 gecikme).
Çıktı: sonuclar/a6_gdelt_testi.json
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import statsmodels.api as sm

import analiz_ortak as A

S = A.excel()
r = A.getiriler(S)["Brent (ICE en yakın vadeli)"].rename("r")
O = pd.read_csv(A.VERI / "gdelt_bq" / "olay_endeksi_v2.csv", parse_dates=["gun"]).set_index("gun")
H = pd.read_csv(A.VERI / "gdelt_bq" / "haber_endeksi_v2.csv", parse_dates=["gun"]).set_index("gun")
X = pd.DataFrame({"tehdit": np.log1p(1e4 * O.h_tehdit / O.makale), "eylem": np.log1p(1e4 * O.h_eylem / O.makale),
                  "ilgi": np.log1p(1e4 * H.h_ad / H.belge)})
X = X.reindex(pd.date_range(X.index.min(), X.index.max()))
tg = r.dropna().index
X["islem"] = pd.Series(tg, index=tg).reindex(X.index).bfill()
Xg = X.groupby("islem")[["tehdit", "eylem", "ilgi"]].max()
Z = (Xg - Xg.rolling(250, min_periods=120).mean().shift(1)) / Xg.rolling(250, min_periods=120).std().shift(1)
D = pd.concat([r, Z], axis=1).dropna()
D["abs_r"] = D.r.abs()
out = {"regresyon": [], "sicrama": []}
for bag in ("r", "abs_r"):
    for per, (a, b) in {"2017-2025": ("2017", "2025"), "2026": ("2026", "2026")}.items():
        d = D.loc[a:b]
        m = sm.OLS(d[bag], sm.add_constant(d[["tehdit", "eylem", "ilgi"]])).fit(cov_type="HAC", cov_kwds={"maxlags": 5})
        out["regresyon"].append({"bagimli": bag, "donem": per, "n": int(len(d)), "r2": round(float(m.rsquared), 3),
                                 **{f"{k}_b": round(float(m.params[k]), 4) for k in ("tehdit", "eylem", "ilgi")},
                                 **{f"{k}_t": round(float(m.tvalues[k]), 3) for k in ("tehdit", "eylem", "ilgi")}})
for per, (a, b) in {"2017-2025": ("2017", "2025"), "2026": ("2026", "2026")}.items():
    d = D.loc[a:b]
    gruplar = {"yalnız tehdit sıçraması": d[(d.tehdit > 2) & (d.eylem < 1)], "yalnız eylem sıçraması": d[(d.eylem > 2) & (d.tehdit < 1)],
               "ikisi birden": d[(d.eylem > 2) & (d.tehdit > 2)], "sakin günler": d[(d.eylem < 1) & (d.tehdit < 1)]}
    for g, x in gruplar.items():
        out["sicrama"].append({"donem": per, "grup": g, "n": int(len(x)), "ort_r": round(float(x.r.mean()), 4) if len(x) else None,
                               "ort_abs_r": round(float(x.abs_r.mean()), 4) if len(x) else None})
# Türk medyası: yıllık oran (Türkçe kaynaklarda Hürmüz payı / küresel pay) ve ton farkı
H["y"] = H.index.year
yil = H.groupby("y")[["h_ad", "belge", "h_ad_tr", "tr_belge"]].sum()
out["tr_oran"] = {int(k): round(float((v.h_ad_tr / v.tr_belge) / (v.h_ad / v.belge)), 2) for k, v in yil.iterrows()}
(A.SONUC / "a6_gdelt_testi.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps(out, ensure_ascii=False)[:1500])
