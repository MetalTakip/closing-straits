"""Dört kritik analizin ortak veri yükleyicisi ve yardımcıları."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
XLSX = HERE.parent / "hurmuz_babulmendep_analiz_verisi.xlsx"
VERI = HERE / "veri"
SONUC = HERE / "sonuclar"
SONUC.mkdir(exist_ok=True)

SAVAS_ONCESI = pd.Timestamp("2026-02-27")

_REN = {
    "Brent (Dated, spot)": "dated", "Brent (ICE en yakın vadeli)": "ice_m1",
    "WTI (Cushing, spot)": "wti_spot", "WTI (NYMEX en yakın vadeli)": "wti_m1",
    "Altın (LBMA PM)": "altin", "Gümüş (LBMA)": "gumus",
    "ABD 10 yıllık reel faiz (TIPS)": "reel_faiz", "ABD dolar endeksi (DXY)": "dxy",
    "ABD 5 yıllık başabaş enflasyon": "basabas", "MSCI ACWI endeks fonu": "acwi",
    "Alüminyum (LME nakit)": "al_nakit", "Alüminyum (LME 3 ay)": "al_3ay",
    "Bakır (LME nakit)": "cu_nakit", "Bakır (LME 3 ay)": "cu_3ay", "Nikel (LME 3 ay)": "ni_3ay",
    "Petrol örtük oynaklık endeksi (OVX)": "ovx",
}


def excel() -> dict[str, pd.DataFrame]:
    return pd.read_excel(XLSX, sheet_name=None)


def fiyatlar(S: dict[str, pd.DataFrame]) -> pd.DataFrame:
    P = S["Fiyatlar"].set_index("Tarih").rename(columns=_REN)
    return P[[c for c in _REN.values() if c in P.columns]]


def getiriler(S: dict[str, pd.DataFrame]) -> pd.DataFrame:
    return S["Getiriler"].set_index("Tarih")


KURLAR = ["usdtry", "eurusd", "usdinr", "usdcny"]


def _londra_gunu(s: pd.Series) -> pd.Series:
    """Yahoo döviz çubukları Londra gün başında (00.00) damgalı. veri_cek.py'nin eski sürümü damgayı UTC tarihine çeviriyordu; İngiltere yaz
    saatinde 00.00 Londra = 23.00 UTC olduğu için her kur bir gün erken (pazar-perşembe) damgalanmıştı. Damga Londra takvim gününe taşınır.
    Yalnız eski damgalı dosyada çalışır (yazın pazar damgası varsa); düzeltilmiş indirmede seri olduğu gibi kalır."""
    s = s.dropna()
    if (s.index.dayofweek == 6).sum() < 10:
        return s
    ertesi = s.index + pd.Timedelta(days=1)
    yaz = ertesi.tz_localize("Europe/London").map(lambda t: t.utcoffset() != pd.Timedelta(0)).to_numpy(dtype=bool)
    s = s.copy()
    s.index = s.index.where(~yaz, ertesi)
    return s[~s.index.duplicated(keep="last")].sort_index()


def yahoo() -> pd.DataFrame:
    Y = pd.read_csv(VERI / "yahoo_vadeli_prim_kur.csv", index_col=0, parse_dates=True)
    sira = list(Y.columns)
    kur = {c: _londra_gunu(Y[c]) for c in KURLAR if c in Y}
    Y = Y.drop(columns=list(kur)).join(pd.DataFrame(kur), how="outer")
    return Y[sira].dropna(how="all").sort_index()


def lme() -> pd.DataFrame:
    return pd.read_csv(VERI / "lme_fiyat_stok_westmetall.csv", parse_dates=["tarih"])


def cot() -> pd.DataFrame:
    return pd.read_csv(VERI / "cftc_cot_disaggregated.csv", parse_dates=["tarih"])


def olaylar(S: dict[str, pd.DataFrame]) -> pd.DataFrame:
    EV = S["Olaylar"].copy()
    EV.columns = ["kod", "bogaz", "tur", "yon", "ana", "tarih", "t0", "aciklama", "gpr"]
    EV["t0"] = pd.to_datetime(EV["t0"])
    return EV


# Web doğrulamasında bulunan, rapordaki listede olmayan 2026 olayları.
# t0: fiyatların tepki verebildiği ilk işlem günü (hafta sonu olayları pazartesiye kayar).
EK_OLAYLAR = pd.DataFrame([
    ("E01", "Hürmüz", "Tırmanma", "2026-04-09", "İran boğazı yeniden kapattı"),
    ("E02", "Hürmüz", "Azalma", "2026-04-17", "Boğaz yeniden açıldı"),
    ("E03", "Hürmüz", "Tırmanma", "2026-04-20", "ABD ablukası sonrası yeniden kapanma (18 Nisan, Cumartesi)"),
    ("E04", "Hürmüz", "Azalma", "2026-06-17", "ABD-İran mutabakatı imzalandı"),
    ("E05", "Hürmüz", "İlan", "2026-06-22", "İran'ın kapanma ilanı (20 Haziran, Cumartesi)"),
    ("E06", "Hürmüz", "Tırmanma", "2026-07-13", "İran boğazı 'yeni bir duyuruya kadar' kapattı"),
    ("E07", "Bâbülmendep", "Tırmanma", "2026-09-11", "Doğu-Batı hattının kapatılması; Perim ve Hanish"),
], columns=["kod", "bogaz", "yon", "t0", "aciklama"])
EK_OLAYLAR["t0"] = pd.to_datetime(EK_OLAYLAR["t0"])


def car_table(R: pd.DataFrame, series: str, events: pd.DataFrame, all_t0: list[pd.Timestamp],
              windows=((-5, -1), (0, 0), (0, 1), (0, 5), (0, 10)), est_len=120, gap=11, excl=(-2, 10)):
    """Raporun yöntemi: sabit ortalama, 120 günlük tahmin penceresi (t0-11'de biter),
    diğer olayların [-2,+10] aralığı dışlanır. all_t0: dışlamada kullanılacak bütün olay günleri."""
    r = R[series].dropna()
    idx = r.index
    pos_all = [idx.searchsorted(t) for t in all_t0]
    out = []
    for _, ev in events.iterrows():
        p0 = idx.searchsorted(ev["t0"])
        if p0 >= len(idx):
            continue
        ex = set()
        for p2 in pos_all:
            if p2 != p0:
                ex.update(range(p2 + excl[0], p2 + excl[1] + 1))
        est, j = [], p0 - gap
        while len(est) < est_len and j >= 0:
            if j not in ex:
                est.append(j)
            j -= 1
        mu, sd = r.iloc[est].mean(), r.iloc[est].std(ddof=1)
        row = {"kod": ev["kod"], "t0": idx[p0].date(), "sd": sd}
        for a, b in windows:
            if p0 + b >= len(idx):
                row[f"[{a},{b}]"] = row[f"z[{a},{b}]"] = np.nan
                continue
            car = (r.iloc[p0 + a:p0 + b + 1] - mu).sum()
            row[f"[{a},{b}]"] = car
            row[f"z[{a},{b}]"] = car / (sd * np.sqrt(b - a + 1))
        out.append(row)
    return pd.DataFrame(out)


def pct(a: float, b: float) -> float:
    return 100.0 * (b / a - 1.0)
