"""Kapsama tam çalıştırmasının granül künyesi (ikinci çapraz denetim, 4. bölüm, madde 1: granül kimliği, UTC zaman, sensör, işleme sürümü ve
dosya özeti saklansın).

Her granül çifti için: UTC başlangıç zamanı, sensör, iki ürünün dosya adı (koleksiyon sürümü ve üretim zamanı adın içinde), üretim zamanı,
PGE sürümü (CMR'de varsa), bayt boyutu ve dosya özeti. Bulut ürünü (CLDMSK_L2_VIIRS_NOAA21, LAADS) için MD5 sağlayıcının CMR kaydından
alınır. Yangın ürünü (VJ214IMG, LP DAAC) için CMR özet vermiyor: dosyalar (her biri yaklaşık 1 MB) yeniden indirilir, SHA-256 hesaplanır,
dosya silinir. Adres üretim zamanını içerdiği için aynı adres işlenen dosyanın kendisidir.
Çıktı: veri/firms/kapsama/granul_kunyesi.csv
Kullanım: KAPSAMA_GECICI=<geçici klasör> python kapsama_kunye.py
"""
from __future__ import annotations

import hashlib
import json
import os
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

import kapsama_pilot as KP
import kapsama_tam as T

CIKTI = T.KLASOR / "granul_kunyesi.csv"


def cmr(kol: str, K: tuple) -> dict[str, dict]:
    q = urllib.parse.urlencode({"collection_concept_id": kol, "bounding_box": ",".join(map(str, K)), "day_night_flag": "night",
                                "temporal": f"{T.ZAMAN[0]}T00:00:00Z,{T.ZAMAN[1]}T23:59:59Z", "page_size": 2000})
    out, sa = {}, None
    while True:
        req = urllib.request.Request("https://cmr.earthdata.nasa.gov/search/granules.umm_json?" + q, headers={"CMR-Search-After": sa} if sa else {})
        with urllib.request.urlopen(req, timeout=120) as r:
            j = json.loads(r.read())
            sa = r.headers.get("CMR-Search-After")
        for it in j["items"]:
            u = it["umm"]
            dg = u.get("DataGranule", {})
            ad = next((i["Identifier"] for i in dg.get("Identifiers", []) if i.get("IdentifierType") == "ProducerGranuleId"), u.get("GranuleUR"))
            m = __import__("re").search(r"\.A(\d{7})\.(\d{4})\.", ad)
            if not m:
                continue
            ad_ = (dg.get("ArchiveAndDistributionInformation") or [{}])[0]
            ck = ad_.get("Checksum") or {}
            out[f"A{m.group(1)}.{m.group(2)}"] = {
                "dosya": ad if ad.endswith(".nc") else ad + ".nc", "uretim": dg.get("ProductionDateTime"),
                "pge": (u.get("PGEVersionClass") or {}).get("PGEVersion"),
                "bayt": ad_.get("SizeInBytes"), "ozet": ck.get("Value"), "ozet_alg": ck.get("Algorithm")}
        if not j["items"] or not sa:
            return out


def sha256(url: str, klasor: Path) -> tuple[str, int]:
    p = KP.indir(url, klasor)
    try:
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        return h, p.stat().st_size
    finally:
        p.unlink(missing_ok=True)


def main():
    H = T.hucreler()
    K = T.kutu(H)
    G = pd.read_csv(T.KLASOR / "granuller.csv")
    a, b = cmr(KP.KOL["ates"], K), cmr(KP.KOL["bulut"], K)
    klasor = Path(os.environ.get("KAPSAMA_GECICI", T.KLASOR / "_gecici")) / "kunye"
    klasor.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(8) as havuz:
        ozetler = list(havuz.map(lambda u: sha256(u, klasor), G["url_ates"]))
    satir = []
    for (_, r), (h, n) in zip(G.iterrows(), ozetler):
        ka, kb = a.get(r["anahtar"], {}), b.get(r["anahtar"], {})
        satir.append({"anahtar": r["anahtar"], "utc_baslangic": pd.to_datetime(r["anahtar"][1:], format="%Y%j.%H%M").strftime("%Y-%m-%dT%H:%MZ"),
                      "sensor": "NOAA-21 VIIRS",
                      "ates_dosya": r["url_ates"].rsplit("/", 1)[1], "ates_uretim": ka.get("uretim"), "ates_pge": ka.get("pge"),
                      "ates_bayt": n, "ates_sha256": h,
                      "bulut_dosya": r["url_bulut"].rsplit("/", 1)[1], "bulut_uretim": kb.get("uretim"), "bulut_bayt": kb.get("bayt"),
                      "bulut_md5_saglayici": kb.get("ozet") if kb.get("ozet_alg") == "MD5" else None})
    S = pd.DataFrame(satir)
    S.to_csv(CIKTI, index=False)
    print(f"{len(S)} çift | yangın SHA-256 {S['ates_sha256'].notna().sum()} | bulut MD5 {S['bulut_md5_saglayici'].notna().sum()} | "
          f"CMR eşleşmesi yangın {sum(k in a for k in G['anahtar'])}, bulut {sum(k in b for k in G['anahtar'])}")
    ad_uyumu = sum(r["bulut_dosya"] == b.get(r["anahtar"], {}).get("dosya") for _, r in S.iterrows())
    print(f"bulut dosya adı CMR ile aynı: {ad_uyumu} / {len(S)}")


if __name__ == "__main__":
    main()
