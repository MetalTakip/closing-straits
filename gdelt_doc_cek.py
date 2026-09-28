"""GDELT DOC 2.0 API'den boğaz haber yoğunluğu ve tonu (günlük, 2017-2026).

Her sorgu ve mod için tek istek atılır (API çok yıllık aralıkta da günlük çözünürlük döndürüyor);
yanıtlar veri/gdelt_doc/ altına ham JSON olarak yazılır.
value = sorguya uyan makale sayısı, norm = o gün GDELT'in izlediği toplam makale sayısı.
API kuralı: en fazla 5 saniyede bir istek. 429 alınırsa artan beklemeyle yeniden denenir.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

KLASOR = Path(__file__).parent / "veri" / "gdelt_doc"
SORGULAR = {
    "hurmuz": '"strait of hormuz"',
    "babulmendep": '("bab el-mandeb" OR "bab al-mandab" OR "bab el mandeb")',
    "kizildeniz_husi": '"red sea" (houthi OR houthis)',
    "hurmuz_olay": '"strait of hormuz" (tanker OR vessel OR ship) (seized OR attacked OR attack OR mine OR drone OR missile)',
}
MODLAR = ("timelinevolraw", "timelinetone")
BAS, BIT = "20170101000000", "20260923235959"


def cek(sorgu: str, mod: str, bas: str, bit: str) -> dict | None:
    u = "https://api.gdeltproject.org/api/v2/doc/doc?" + urllib.parse.urlencode(
        {"query": sorgu, "mode": mod, "startdatetime": bas, "enddatetime": bit, "format": "json"})
    bekle = 30
    for _ in range(7):
        time.sleep(20)
        try:
            ham = urllib.request.urlopen(u, timeout=120).read().decode("utf-8", "replace")
            return json.loads(ham)
        except urllib.error.HTTPError as e:
            print(f"   HTTP {e.code}; {bekle} sn bekleniyor", flush=True)
        except json.JSONDecodeError:
            print(f"   JSON değil: {ham[:80]!r}; {bekle} sn bekleniyor", flush=True)
        except Exception as e:  # ağ hatası
            print(f"   hata {e!r}; {bekle} sn bekleniyor", flush=True)
        time.sleep(bekle)
        bekle = min(bekle * 2, 300)
    return None


if __name__ == "__main__":
    KLASOR.mkdir(parents=True, exist_ok=True)
    for ad, sorgu in SORGULAR.items():
        for mod in MODLAR:
            f = KLASOR / f"{ad}_{mod}_2017_2026.json"
            if f.exists() and f.stat().st_size > 200:
                continue
            j = cek(sorgu, mod, BAS, BIT)
            if j is None:
                print(f"BAŞARISIZ {f.name}", flush=True)
                continue
            f.write_text(json.dumps(j), encoding="utf-8")
            n = len(j["timeline"][0]["data"]) if j.get("timeline") else 0
            print(f"{f.name}: {n} gün, çözünürlük {j.get('query_details', {}).get('date_resolution')}", flush=True)
    print("BİTTİ", flush=True)
