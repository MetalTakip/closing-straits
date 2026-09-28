"""SOMO (Irak Devlet Petrol Pazarlama Kurumu) aylık ham petrol ihracatı, çıkış noktasına göre.

Neden: SOMO'nun etkileşimli ihracat grafiği (somooil.gov.iq/en/exports/chart) yalnız Ocak-Haziran 2026'yı gösteriyor ve kategorileri
PDF raporlarıyla her ay örtüşmüyor: Şubat 2026'da "North" diye gösterdiği 971.130 varil, aynı ayın PDF'ine göre Basra yakınındaki
Khor al-Zubair'den (Körfez içi) yüklenen Qaiyarah ham petrolü; grafiğin "KRG" dediği 5.551.610 varil ise PDF'te "Kirkuk Crude Oil
Exports via Ceyhan". Aylık PDF raporları (somooil.gov.iq/ar/researches) ayrıca Haziran 2025-Mart 2026 ve Ağustos 2026'yı kapsıyor.

Yöntem: her PDF'te tablo başlıkları ile "Total" satırları konumlarına göre eşlenir (her toplam, aynı sayfada hemen üstündeki başlığa
aittir). Ham petrol kalemleri çıkış noktasına göre sınıflanır: Basrah Medium ve Heavy -> Basra (Körfez); Kirkuk ve (KRG) Kirkuk "via
Ceyhan" -> Ceyhan; Qaiyarah/Qayara floater -> Khor al-Zubair (Körfez). Fuel oil, nafta ve benzeri ürün tabloları dışarıda kalır.
PDF'i olmayan aylar (Nisan-Haziran 2026) grafikten alınır: Mart 2026 PDF'i grafiğin "North" kalemini "Kirkuk via Ceyhan", "KRG"
kalemini "(KRG) Kirkuk via Ceyhan" olarak doğruluyor; bu aylar kaynak sütununda ayrıca işaretlenir.
Doğrulama: PDF sonuçları grafiğin Ocak-Mart 2026 Basra ve Basra dışı toplamlarıyla, Mart 2026'da da iki Ceyhan kalemiyle tek tek
karşılaştırılır (fark 2 varilden büyükse betik durur). Nisan-Haziran 2026 eşlemesi Mart'a dayanan bir varsayımdır (kaynak sütununda yazılı).
Yeniden üretim iki düzeydir: (1) analizler CSV'den çalışır (analiz15_alev.py); (2) CSV bu betikle SOMO kaynaklarından üretilir. Kod paketi
PDF'leri dağıtmaz; veri/rafineri/irak_somo/kaynaklar.csv adresleri, erişim tarihlerini ve SHA256 değerlerini tutar.
Çıktı: veri/rafineri/irak_somo_cikis_aylik.csv, veri/rafineri/irak_somo/kaynaklar.csv; PDF'ler veri/rafineri/irak_somo/pdf/ altında.
Kullanım: python somo_cek.py   (ağ yoksa yalnız saklanan PDF'lerle çalışır)
"""
from __future__ import annotations

import calendar
import re
import urllib.request
from pathlib import Path

import pandas as pd
import pymupdf

import analiz_ortak as A

KLASOR = A.VERI / "rafineri" / "irak_somo"
PDFK = KLASOR / "pdf"
GRAFIK = KLASOR / "somo_exports_chart_2026-09-24.html"
CIKTI = A.VERI / "rafineri" / "irak_somo_cikis_aylik.csv"
SITE = "https://www.somooil.gov.iq"
ARASTIRMA = range(100, 131)  # rapor sayfası kimlikleri; bulunmayanlar atlanır
AYLAR = {m: i for i, m in enumerate(calendar.month_name) if m}


def indir():
    PDFK.mkdir(parents=True, exist_ok=True)
    for n in ARASTIRMA:
        if list(PDFK.glob(f"r{n}_*.pdf")):
            continue
        try:
            req = urllib.request.Request(f"{SITE}/ar/researches/{n}", headers={"User-Agent": "Mozilla/5.0"})
            s = urllib.request.urlopen(req, timeout=40).read().decode("utf-8", "ignore")
        except Exception:  # noqa: BLE001 (sayfa yok ya da ağ yok)
            continue
        m = re.search(r"storage/[0-9a-z-]+\.pdf", s)
        if not m:
            continue
        req = urllib.request.Request(f"{SITE}/{m.group(0)}", headers={"User-Agent": "Mozilla/5.0"})
        (PDFK / f"r{n}_{Path(m.group(0)).name}").write_bytes(urllib.request.urlopen(req, timeout=60).read())
        print("indirildi:", n, m.group(0))


def ay_bul(metin: str) -> str | None:
    m = re.search(r"(?:Report for|During|Exports\s*[–-]\S*)\s*([A-Z][a-z]+)\s+(20\d\d)", metin)
    return f"{m.group(2)}-{AYLAR[m.group(1)]:02d}" if m else None


def sinifla(baslik: str):
    b = baslik.lower()
    if "fuel oil" in b or "naphtha" in b or "gasoil" in b or "summary" in b:
        return None
    if "basrah medium" in b:
        return "Basra", "Basrah Medium"
    if "basrah heavy" in b:
        return "Basra", "Basrah Heavy"
    if "qaiyarah" in b or "qayara" in b:
        return "Khor al-Zubair", "Qaiyarah"
    if "kirkuk" in b and "ceyhan" in b:
        return "Ceyhan", "(KRG) Kirkuk" if "krg" in b else "Kirkuk"
    return None


# Satır: miktar (binlik virgüllü ya da virgülsüz en az dört hane; ondalıklı ton değerleri dışarıda) ve ardından yükleme tarihi
# (2-Jun-2025 ya da 26/06/2025). Haziran 2025 PDF'inde bir satır virgülsüz ve eğik çizgili (üçüncü denetim, Z05).
TARIH = r"(?:\d{1,2}-[A-Za-z]{3}(?:-\d{2,4})?|\d{1,2}/\d{1,2}/\d{2,4})"
SATIR = re.compile(r"(?<![\d.,])(\d{1,3}(?:,\d{3})+|\d{4,})(?![\d.,])\s+" + TARIH)
TARIH_RE = re.compile(r"(?<![\d/-])" + TARIH + r"(?![\d/])")
TOPLAM = re.compile(r"Total\s*(?:\(BBL\)|BBl|BBL)?\s*([\d,]{5,})")


def baslik_mi(t: str) -> bool:
    return bool(re.search(r"Crude Oil|Floater Tank|Fuel Oil|Naphtha", t) and re.search(r"Export|Exoprt|Floater|Summary|During", t))


def pdf_oku(yol: Path) -> list[dict]:
    """Sayfalar boyunca okuma sırasıyla: başlık satırı yeni tabloyu açar; tablo blokları (sayfaya taşan devam blokları dahil) son
    başlığa eklenir. Miktar: basılı toplam, yoksa satırların toplamı."""
    d = pymupdf.open(yol)
    ay = ay_bul(" ".join(p.get_text() for p in d))
    ogeler = []
    for i, p in enumerate(d):
        for b in p.get_text("dict")["blocks"]:
            for l in b.get("lines", []):
                t = " ".join(sp["text"] for sp in l["spans"]).strip()
                if baslik_mi(t):
                    ogeler.append((i, l["bbox"][1], 0, t))
        for b in p.get_text("blocks"):
            t = re.sub(r"\s+", " ", b[4]).strip()
            if "Vessel Name" in t or (SATIR.search(t) and not baslik_mi(t)):
                ogeler.append((i, b[1], 1, t))
    tablolar, simdiki = [], None
    for sayfa, y, tur, t in sorted(ogeler):
        if tur == 0:
            simdiki = {"baslik": t, "satir": [], "toplam": None, "sayfa": sayfa + 1, "tarih": 0}
            tablolar.append(simdiki)
        elif simdiki is not None:
            simdiki["satir"] += [int(x.replace(",", "")) for x in SATIR.findall(t)]
            simdiki["tarih"] += len(TARIH_RE.findall(t))
            m = TOPLAM.search(t)
            if m:
                simdiki["toplam"] = int(m.group(1).replace(",", ""))
    out = {}
    for tb in tablolar:
        s = sinifla(tb["baslik"])
        if not s or not (tb["satir"] or tb["toplam"]):
            continue
        miktar = tb["toplam"] if tb["toplam"] is not None else sum(tb["satir"])
        if tb["satir"] and tb["tarih"] != len(tb["satir"]):  # tarihli her satırın miktarı okunmuş olmalı
            print(f"  not {yol.name}: {s[1]} tablosunda {tb['tarih']} tarih, {len(tb['satir'])} miktar okundu")
        if tb["toplam"] is not None and tb["satir"] and abs(sum(tb["satir"]) - tb["toplam"]) > 2:
            print(f"  not {yol.name}: {s[1]} satır toplamı {sum(tb['satir']):,}, basılı toplam {tb['toplam']:,} (basılı toplam kullanıldı)")
        k = (ay, s[1])
        if k in out:  # aynı ayda aynı kalemden ikinci tablo (ör. iki Qaiyarah tankı): topla
            out[k]["miktar_varil"] += miktar
            out[k]["pdf_basligi"] += " + " + tb["baslik"][:80]
        else:
            out[k] = {"ay": ay, "cikis": s[0], "kalem": s[1], "miktar_varil": miktar, "pdf_basligi": tb["baslik"][:160],
                      "kaynak": f"SOMO aylık rapor PDF {yol.name} (s. {tb['sayfa']})"}
    return list(out.values())


def grafik_oku() -> pd.DataFrame:
    if not GRAFIK.exists():  # kaydedilmiş sayfa yoksa indir (Nisan-Haziran 2026 yalnız bu grafikte)
        req = urllib.request.Request(f"{SITE}/en/exports/chart", headers={"User-Agent": "Mozilla/5.0"})
        GRAFIK.parent.mkdir(parents=True, exist_ok=True)
        GRAFIK.write_bytes(urllib.request.urlopen(req, timeout=60).read())
    s = GRAFIK.read_text(encoding="utf-8", errors="ignore").replace('\\"', '"')
    satir = []
    for m in re.finditer(r'\{"m":"(20\d\d-\d\d)","reportId":\d+,"title":"[^"]*","q":\d+,"rev":[\d.]+,"price":[\d.]+,"qByT":\{([^}]*)\}', s):
        for k, v in re.findall(r'"(\d)":(\d+)', m.group(2)):
            satir.append({"ay": m.group(1), "terminal_id": int(k), "miktar_varil": int(v)})
    return pd.DataFrame(satir)


def manifest(G_yolu: Path) -> pd.DataFrame:
    """Kaynak manifesti: her PDF ve grafik anlık görüntüsü için adres, erişim tarihi, bayt ve SHA256 (paket PDF'leri dağıtmaz; bu dosyayla
    aynı kaynaklar yeniden indirilip doğrulanabilir)."""
    import hashlib
    satir = []
    for f in sorted(PDFK.glob("r*_*.pdf")):
        n, ad = f.name.split("_", 1)
        satir.append({"dosya": f"pdf/{f.name}", "tur": "aylık rapor PDF", "ay": ay_bul(" ".join(p.get_text() for p in pymupdf.open(f))),
                      "rapor_sayfasi": f"{SITE}/ar/researches/{n[1:]}", "adres": f"{SITE}/storage/{ad}",
                      "erisim_tarihi": ERISIM.get(f.name, "2026-09-25"), "bayt": f.stat().st_size, "sha256": hashlib.sha256(f.read_bytes()).hexdigest()})
    if G_yolu.exists():
        satir.append({"dosya": G_yolu.name, "tur": "ihracat grafiği anlık görüntüsü (HTML)", "ay": "2026-01..2026-06", "rapor_sayfasi": "",
                      "adres": f"{SITE}/en/exports/chart", "erisim_tarihi": "2026-09-24", "bayt": G_yolu.stat().st_size,
                      "sha256": hashlib.sha256(G_yolu.read_bytes()).hexdigest()})
    return pd.DataFrame(satir)


ERISIM: dict[str, str] = {}  # dosya adı -> erişim tarihi; listede olmayanlar 25 Eylül 2026'da indirildi


def main():
    try:
        indir()
    except Exception as e:  # noqa: BLE001
        print("indirme atlandı:", repr(e)[:80])
    pdfler = sorted(PDFK.glob("*.pdf"))
    if not pdfler:
        raise SystemExit(f"SOMO PDF'i yok: {PDFK} boş ve indirme başarısız. Kaynak listesi: {KLASOR / 'kaynaklar.csv'} (adres ve SHA256).")
    pdf = pd.DataFrame([r for f in pdfler for r in pdf_oku(f)])
    G = grafik_oku()
    # doğrulama 1: PDF'i olan grafik aylarında Basra ve Basra dışı toplamlar tutmalı
    for ay in sorted(set(pdf.ay) & set(G.ay)):
        gp = pdf[pdf.ay == ay]
        g = G[G.ay == ay]
        basra_g, diger_g = g[g.terminal_id == 1].miktar_varil.sum(), g[g.terminal_id != 1].miktar_varil.sum()
        basra_p, diger_p = gp[gp.cikis == "Basra"].miktar_varil.sum(), gp[gp.cikis != "Basra"].miktar_varil.sum()
        assert abs(basra_g - basra_p) <= 2 and abs(diger_g - diger_p) <= 2, (ay, basra_g, basra_p, diger_g, diger_p)
        print(f"doğrulandı {ay}: Basra {basra_p:,} | Basra dışı {diger_p:,} (grafikle aynı)")
    # doğrulama 2: Nisan-Haziran'da kullanılan eşlemenin dayanağı. Mart 2026'da grafiğin "North" (2) kalemi PDF'in "Kirkuk via Ceyhan",
    # "KRG" (3) kalemi PDF'in "(KRG) Kirkuk via Ceyhan" kalemine birebir eşit olmalı. (Şubat'ta grafiğin "North" kalemi Khor al-Zubair'den
    # yüklenen Qaiyarah petrolü; Ocak ve Şubat'ta tek Ceyhan kalemi grafikte "KRG" görünüyor: grafik etiketleri her ay aynı değil.)
    m = pdf[pdf.ay == "2026-03"].set_index("kalem").miktar_varil
    gm = G[G.ay == "2026-03"].set_index("terminal_id").miktar_varil
    assert abs(m["Kirkuk"] - gm[2]) <= 2 and abs(m["(KRG) Kirkuk"] - gm[3]) <= 2, (m.to_dict(), gm.to_dict())
    print(f"doğrulandı 2026-03 çıkış noktası: grafik North = PDF Kirkuk via Ceyhan ({gm[2]:,}), grafik KRG = PDF (KRG) Kirkuk ({gm[3]:,})")
    # PDF'i olmayan grafik ayları (Nisan-Haziran 2026): Mart eşlemesi VARSAYIM olarak uygulanır ve kaynak sütununda yazılır
    eklenecek = []
    for ay in sorted(set(G.ay) - set(pdf.ay)):
        for _, r in G[G.ay == ay].iterrows():
            cikis, kalem = {1: ("Basra", "Basra (toplam)"), 2: ("Ceyhan", "Kirkuk"), 3: ("Ceyhan", "(KRG) Kirkuk")}[r.terminal_id]
            eklenecek.append({"ay": ay, "cikis": cikis, "kalem": kalem, "miktar_varil": int(r.miktar_varil), "pdf_basligi": "",
                              "kaynak": "SOMO ihracat grafiği (somooil.gov.iq/en/exports/chart, 24 Eylül 2026); bu ay için aylık PDF yok; "
                                        + ("çıkış noktası Mart 2026 PDF eşlemesine göre varsayıldı" if r.terminal_id != 1 else "Basra toplamı")})
    T = pd.concat([pdf, pd.DataFrame(eklenecek)], ignore_index=True).sort_values(["ay", "cikis", "kalem"])
    T["gun"] = [calendar.monthrange(int(a[:4]), int(a[5:]))[1] for a in T.ay]
    T["varil_gun"] = (T.miktar_varil / T.gun).round(0).astype(int)
    T[["ay", "cikis", "kalem", "miktar_varil", "gun", "varil_gun", "pdf_basligi", "kaynak"]].to_csv(CIKTI, index=False, encoding="utf-8")
    M = manifest(GRAFIK)
    M.to_csv(KLASOR / "kaynaklar.csv", index=False, encoding="utf-8")
    ozet = T.pivot_table(index="ay", columns="cikis", values="miktar_varil", aggfunc="sum")
    print((ozet / 1e6).round(2).to_string())
    print("yazıldı:", CIKTI, "ve", KLASOR / "kaynaklar.csv", f"({len(M)} kaynak)")


if __name__ == "__main__":
    main()
