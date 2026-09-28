"""NASA FIRMS VIIRS aktif ateş tespitleri (375 m) Körfez, Irak, İran ve Umman için: gaz alevi analizinin (analiz15_alev.py) ham verisi.

Kaynak: https://firms.modaps.eosdis.nasa.gov (NASA LANCE ve FIRMS; veri serbest, atıf gerekir). Anahtar depo kökündeki .env dosyasında
(FIRMS_MAP_KEY) ya da aynı adlı ortam değişkeninde; ekrana ve dosyalara yazılmaz. API alan isteği en fazla 5 gün kapsar; dönem 5 günlük parçalarla çekilir.
Ana seri NOAA-21 NRT (API'de Eylül 2024'ten bugüne tek işleme akışı; daha eskisi boş dönüyor). Sağlamlık için Suomi NPP: Haziran 2026'ya kadar standart işleme (SP),
sonrası NRT.
Kullanım: python firms_cek.py [kaynak ...]   Çıktı: veri/firms/<kaynak>.csv.gz
"""
from __future__ import annotations

import gzip
import io
import sys
import time
import urllib.error
import urllib.request
from datetime import date, timedelta
from pathlib import Path

import pandas as pd

BURASI = Path(__file__).resolve().parent
CIKTI = BURASI / "veri" / "firms"
CIKTI.mkdir(parents=True, exist_ok=True)
KUTU = "43.5,17.0,60.0,37.5"  # batı, güney, doğu, kuzey: Irak, Kuveyt, İran'ın güneyi, Körfez ülkeleri, Umman
KAYNAKLAR = {
    "VIIRS_NOAA21_NRT": ("2024-09-01", "2026-09-24"),
    "VIIRS_SNPP_SP": ("2024-09-01", "2026-06-30"),
    "VIIRS_SNPP_NRT": ("2026-07-01", "2026-09-24"),
}


def anahtar() -> str:
    """Önce FIRMS_MAP_KEY ortam değişkeni, sonra betiğin bir üst ve kendi klasöründeki .env dosyası."""
    import os
    if os.environ.get("FIRMS_MAP_KEY"):
        return os.environ["FIRMS_MAP_KEY"].strip()
    for env in (BURASI.parent / ".env", BURASI / ".env"):
        if env.exists():
            for satir in env.read_text(encoding="utf-8").splitlines():
                if satir.startswith("FIRMS_MAP_KEY="):
                    return satir.split("=", 1)[1].strip()
    raise SystemExit("FIRMS_MAP_KEY bulunamadı (ortam değişkeni ya da .env)")


def parca(k: str, kaynak: str, bas: date, gun: int) -> pd.DataFrame:
    url = f"https://firms.modaps.eosdis.nasa.gov/api/area/csv/{k}/{kaynak}/{KUTU}/{gun}/{bas.isoformat()}"
    for deneme in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                metin = r.read().decode("utf-8")
            if metin.startswith("Invalid") or "Exceeding" in metin[:200]:
                raise RuntimeError(metin[:200])
            return pd.read_csv(io.StringIO(metin)) if metin.strip() else pd.DataFrame()
        except (urllib.error.URLError, RuntimeError, TimeoutError) as e:
            bekle = 30 * (deneme + 1)
            print(f"  {bas} hata ({str(e)[:80]}); {bekle} sn bekleniyor", flush=True)
            time.sleep(bekle)
    raise SystemExit(f"{kaynak} {bas}: 5 denemede alınamadı")


def cek(kaynak: str):
    k = anahtar()
    a, b = (date.fromisoformat(x) for x in KAYNAKLAR[kaynak])
    hedef = CIKTI / f"{kaynak}.csv.gz"
    parcalar, g = [], a
    while g <= b:
        gun = min(5, (b - g).days + 1)
        d = parca(k, kaynak, g, gun)
        parcalar.append(d)
        print(f"{kaynak} {g} +{gun} gün: {len(d)} tespit", flush=True)
        g += timedelta(days=gun)
        time.sleep(1)
    D = pd.concat([d for d in parcalar if len(d)], ignore_index=True).drop_duplicates()
    with gzip.open(hedef, "wt", encoding="utf-8") as f:
        D.to_csv(f, index=False)
    print(f"{kaynak}: {len(D)} tespit -> {hedef.name}")


if __name__ == "__main__":
    for kaynak in (sys.argv[1:] or list(KAYNAKLAR)):
        cek(kaynak)
