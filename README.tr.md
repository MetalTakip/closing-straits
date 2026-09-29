# Kapanan Boğazlar: 2026 Hürmüz kapanmaları için uydu radarıyla gemi sayımları, gaz alevleri ve GDELT haber endeksleri

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23009094.svg)](https://doi.org/10.5281/zenodo.23009094)

**Kapanan Boğazlar: Hürmüz ve Bâbülmendep şokları petrol, gaz, metal ve altın piyasalarını nasıl yeniden fiyatladı?**
(MetalTakip Sektörel Araştırma, Eylül 2026) raporundaki uydu radarı, gaz alevi ve haber verisi analizlerinin verisi ve kodu. Rapor:
[metaltakip.com](https://www.metaltakip.com/arastirma/kapanan-bogazlar). English: [README.md](README.md)

Yazarlar: Yusuf Oskay Dağlı ([ORCID 0009-0008-5513-9217](https://orcid.org/0009-0008-5513-9217)) ve Ebru Özpolat
([ORCID 0009-0006-3864-142X](https://orcid.org/0009-0006-3864-142X)), MetalTakip.

## İçerik

- **Beş noktada 407 Sentinel-1 radar görüntüsünde gemi sayımı**: Hürmüz geçiş koridoru, Füceyre demirleme alanı, Ras Tanura,
  Yanbu ve Bâbülmendep; Temmuz 2025'ten 22 Eylül 2026'ya kadar, görüntü başına bir satır.
- **Sayımı üreten kod**, raporda kullanılan dönem özetleri ve sağlamlık regresyonu, raporda gösterilen iki Hürmüz radar görüntüsü.
- **VIIRS gece yangın tespitlerinden gaz alevleri** (NASA FIRMS, Eylül 2024-Eylül 2026): Irak, Kuveyt, İran, Körfez ülkeleri ve Umman'daki
  kalıcı alev noktaları; bölge serilerini ve dönem etkilerini üreten kod.
- **Alevler için bulut ve kapsama denetimi**: kalıcı noktalar üzerinden geçen bütün NOAA-21 gece geçişleri, NASA'nın Level-2 yangın ve bulut
  maskelerinden (3.111 granül çifti); hücre x gece x geçiş tabloları ve bunları indirip modelleyen kod.
- **İki kapanmanın karşılaştırması ve Basra testi** (`analiz17_iki_kapanma.py`): aynı modellerde iki kapanmanın farkı ve Basra açık deniz
  terminallerinin çevresindeki altıncı bir kutuda, iki yörünge yönünde gemi sayımı (182 görüntü).
- **2017-2026 günlük GDELT endeksleri**: Hürmüz ya da Bâbülmendep adını anan haberler (GKG) ve iki boğazın çevresindeki kutularda
  konumlanan tehdit ve eylem olayları (Events); endeksleri üreten BigQuery sorguları.
- Raporun olay çalışmasındaki **39 olay** (Ek A), Türkçe ve İngilizce.

## Temel sonuçlar

Radar görüntüsü başına gemi sayısı, dönem medyanı; yalnız Sentinel-1C ve 1D görüntüleri. Son sütun, bütün görüntülerle kurulan
Poisson sayım modeline (dönem ve uydu kuklaları) göre ikinci kapanmadaki değişim (heteroskedastisiteye dayanıklı standart hata). Kaynak: `sonuclar/a7_sar_ozet.json`.

| Konum | Yörünge, çözünürlük | Dönemlerdeki 1C/1D görüntü | Savaş öncesi | İlk kapanma | İkinci kapanma | İkinci kapanmada değişim (%95 GA) |
|---|---|---:|---:|---:|---:|---|
| Hürmüz koridoru | inen, 40 m | 44 | 40,5 | 10 | 9,5 | −%81 (−%88,1 ile −%70,1) |
| Füceyre | inen, 20 m | 40 | 177 | 171 | 254 | +%46 (+%36,7 ile +%56,9) |
| Füceyre | çıkan, 20 m | 44 | 157 | 165 | 231 | +%43 (+%30,7 ile +%55,4) |
| Yanbu | çıkan, 20 m | 37 | 23 | 29 | 33 | +%32 (+%3,2 ile +%69,9) |
| Bâbülmendep | çıkan, 40 m | 47 | 46 | 50 | 54,5 | +%22 (−%1,9 ile +%52) |
| Ras Tanura | inen, 40 m | 44 | 115 | 118,5 | 108,5 | −%7 (−%14,4 ile +%0,5) |

Dönemler: savaş öncesi 1 Temmuz 2025-27 Şubat 2026; ilk kapanma 2 Mart-13 Haziran 2026; açılış 18 Haziran-10 Temmuz 2026;
ikinci kapanma 13 Temmuz-22 Eylül 2026.

Sağlamlık testleri her yörüngenin `saglamlik` anahtarında: log(1 + gemi) modeli (`log`), yalnız Sentinel-1C ve 1D (`tutarli`), yalnız
Sentinel-1C (`yalniz_1c`; her dönemde gözlenen tek uydu), yalnız büyük yankılar (`buyuk`: 20 m'de en az 12, 40 m'de en az 3 piksel; alan
eşiği, ölçülmüş gemi boyu değil), dört gecikmeli Newey-West (`HAC`), aya göre kümelenmiş hata (`ay`; t dağılımıyla sınırlar `alt_t`,
`ust_t`) ve ikinci kapanmada görüntüleri tek tek çıkarma (`tek_cikar`). Beş konum x iki kapanma için Holm düzeltmeli p değerleri `p_holm`
alanında. Bu kurgularda Hürmüz koridorundaki düşüşün nokta tahmini %75 ile %85 arasında kalıyor. Holm düzeltmesinden sonra %5
düzeyinde yalnız Hürmüz'ün iki düşüşü ve Füceyre'nin ikinci kapanma artışı anlamlı kalıyor; Yanbu'nun (altı görüntü) ve Bâbülmendep'in
(+%22; %95 aralığının alt sınırı −%2) ikinci kapanma artışları ise anlamlı kalmıyor.

İlk kapanmada AIS'e dayanan geçiş sayıları (IMF PortWatch) %93 düştü; EIA'ya göre boğazdan geçen petrol 2026'nın ikinci çeyreğinde
2025'in son çeyreğine göre %77 azdı. Koridordaki radar sayımı ise %80 düştü. Radar gemiyi AIS yayınından bağımsız görür, ama sayım
bir geçiş sayımı değildir: kutuda bekleyen, yavaşlayan ya da askerî görevdeki gemiler de sayılır. Üç seri farklı şeyleri farklı dönemlerde
ölçtüğü için aralarındaki fark AIS kaybına ya da yük hacmine çevrilemez.

## Gaz alevleri

Ölçü, kalıcı alev noktalarında gece başına tespit edilen ışınım gücüdür (FRP). Bir geçişteki pikseller toplanır; bir noktayı birden çok
geçiş gördüyse piksellerinin tarama yönüne dik ortalama boyu en küçük olan geçiş seçilir. Bu boy bakış açısıyla artar; tarama yönündeki boy
ise VIIRS'in piksel birleştirmesi yüzünden bakış açısını düzenli izlemez. Savaş öncesine göre değişim, dönem ve takvim ayı kuklalarıyla
Poisson sözde en çok olabilirlik (PPML) tahmininden gelir (HAC standart hata). `log`: aynı kuklalarla log(1 + FRP) modeli. Kalın: düzeltilmemiş %95 güven aralığı
sıfırı içermiyor. Holm düzeltmeli p değerleri (aile: aşağıdaki on bölge x iki kapanma) `p_holm` alanında; açılış dönemi, öteki bölgeler ve
güney-kuzey karşılaştırması bu ailenin dışında. Kaynak: `sonuclar/a15_alev.json`.

| Bölge | Savaş öncesi, MW | İlk kapanma (NOAA-21 / Suomi NPP / log) | İkinci kapanma (NOAA-21 / Suomi NPP / log) |
|---|---:|---|---|
| Irak (güney) | 1.053 | **−%84** (Holm <0,001) / **−%86** / **−%81** | **−%37** (Holm <0,001) / **−%36** / **−%41** |
| Kuveyt | 24 | −%4 (Holm 1,000) / −%2 / −%1 | **+%71** (Holm 0,005) / **+%86** / **+%58** |
| Katar | 67 | **−%31** (Holm 0,141) / −%38 / **−%58** | −%8 (Holm 0,770) / **−%18** / **−%13** |
| Bahreyn | 8 | **−%78** (Holm <0,001) / **−%73** / **−%69** | **−%74** (Holm 0,006) / **−%78** / **−%72** |
| Batı İran (çoğu Huzistan) | 1.052 | −%15 (Holm 0,498) / −%7 / −%11 | **−%26** (Holm <0,001) / **−%28** / **−%27** |
| İran (Asaluyeh) | 206 | −%4 (Holm 1,000) / −%15 / −%27 | **−%21** (Holm <0,001) / **−%23** / **−%27** |
| Suudi Arabistan | 179 | +%10 (Holm 1,000) / +%20 / −%3 | +%14 (Holm 1,000) / +%21 / +%14 |
| BAE | 62 | **+%304** (Holm <0,001) / **+%348** / **+%205** | **+%201** (Holm <0,001) / **+%206** / **+%203** |
| Umman | 81 | −%5 (Holm 1,000) / −%10 / −%10 | **+%43** (Holm 0,036) / **+%41** / **+%53** |
| Irak (kuzey) | 273 | **+%23** (Holm 0,006) / +%14 / +%5 | **+%30** (Holm <0,001) / **+%22** / **+%29** |

Irak'ın güney sahalarında tespit edilen alev gücü, Basra'dan ihracat çökerken düştü (SOMO: Şubat'ta günde 3,33 milyon varil, Mayıs'ta 79 bin varil); boğazı
atlayan Kerkük-Ceyhan hattına bağlı kuzeyde düşmedi. Tespit edilen alev gücü üretimle bire bir orantılı değildir: gazını satamayan bir tesis
daha çok gaz yakabilir; BAE ve Kuveyt'teki artışı bu açıklayabilir (test edilmedi). Katar'daki daha küçük düşüş sağlam değil: Holm
düzeltmesinden sonra ve Suomi NPP'de anlamlı değil, nokta tanımına göre değişiyor ve ilk kapanmada Katar'da hiçbir kalıcı noktada tespit
olmayan gecelerin payı %13 oldu (savaş öncesinde %2); aşağıdaki Level-2 denetimine göre bu gecelerin her birinde Katar'daki noktaların çoğu bulut
altındaydı. İlgili çalışma: Zhizhin ve Bazilian (2026, Payne Enstitüsü) VIIRS
Nightfire verisiyle Mart 2026'nın ilk on gününde Basra çevresinde alevlerin sert düştüğünü buldu; BAE, Katar ve Suudi Arabistan'da düşüş daha
küçüktü, İran'da ise alevler arttı. Bizim serilerimiz o günler için aynı işaretleri veriyor (`sonuclar/a15_alev.json` içinde `erken_mart`). Çıkış noktasına göre
ihracat (Basra, Ceyhan, Khor al-Zubair) `somo_cek.py` ile SOMO'nun aylık PDF raporlarından okunuyor.

## Bulut ve kapsama denetimi

FIRMS noktaları bulut ve kapsama bilgisi taşımaz. Alev sonuçlarını sınamak için her kalıcı hücre, NOAA-21'in Level-2 dosyalarındaki her gece
geçişiyle eşlendi: VJ214IMG'nin 375 m yangın maskesi ve konum ile bakış açısını da taşıyan 750 m bulut maskesi CLDMSK_L2_VIIRS_NOAA21
(748 gecede 3.111 granül çifti; tespit olmayan geçişler dahil). Hücre-geçiş, 375 m piksellerinin en az yarısı su, açık kara ya da
yangınsa geçerli sayılır (yangın pikseli açık sayılır: ürün bulut gördüğü yerde yangın aramaz); en az yarısı bulutsa bulutludur. Hücrenin
geçerli geçişleri arasında bakış açısı en küçük olanı FRP'ye bakılmadan seçilir ve kalan bowtie tekrarları (QA bit 22) çıkarılır. Ana model
bölgenin hücrelerinin en az %80, %90 ya da %95'inin geçerli gözlendiği gecelerde FIRMS ve Level-2 ölçüleriyle yeniden kurulur; geçerli
hücre-gecelerde hücre sabit etkili, Driscoll-Kraay standart hatalı PPML de tahmin edilir. Bağımsız bulut maskesi seçim için kullanılmaz: bulutlu
dediği pay tespit olmayan geçerli hücre-geçişlerde %25, 50 MW üstündekilerde %94; bu maskeyle seçim yapılsaydı güçlü alevler
dışarıda kalırdı. Kaynak: `sonuclar/a16_kapsama.json`; tablolar `veri/firms/kapsama/` altında.

| Bölge | İlk kapanma: bütün geceler | İlk: %90 gözlenen, FIRMS / L2 | İlk: hücre modeli | İkinci kapanma: bütün geceler | İkinci: %90 gözlenen, FIRMS / L2 | İkinci: hücre modeli |
|---|---|---|---|---|---|---|
| Irak (güney) | −%84 | −%83 / −%85 | −%85 | −%37 | −%37 / −%36 | −%36 |
| Kuveyt | −%4 | −%6 / −%11 | −%7 | +%71 | +%79 / +%85 | +%77 |
| Katar | −%31 | −%15 / −%19 | −%27 | −%8 | −%8 / −%10 | −%11 |
| Bahreyn | −%78 | −%70 / −%69 | −%71 | −%74 | −%75 / −%71 | −%69 |
| Batı İran (çoğu Huzistan) | −%15 | −%18 / −%18 | −%14 | −%26 | −%23 / −%23 | −%25 |
| İran (Asaluyeh) | −%4 | +%7 / +%10 | +%9 | −%21 | −%19 / −%18 | −%17 |
| Suudi Arabistan | +%10 | +%28 / +%27 | +%25 | +%14 | +%16 / +%10 | +%9 |
| BAE | +%304 | +%325 / +%324 | +%327 | +%201 | +%199 / +%200 | +%204 |
| Umman | −%5 | +%3 / +%0 | −%2 | +%43 | +%32 / +%41 | +%40 |
| Irak (kuzey) | +%23 | +%14 / +%13 | +%17 | +%30 | +%30 / +%27 | +%25 |

İyi gözlenen gecelerde güney Irak'ın ilk kapanmadaki değişimi on iki kurguda −%85 ile −%81 arasında, güney-kuzey farkı
−%87 ile −%86 arasında kalıyor. Katar'ın ilk kapanmadaki değişimi −%27 ile −%15 arasında ve hiçbir kurguda Holm
düzeltmesinden sonra anlamlı değil; ilk kapanma günleri Körfez'in güney kıyısında 2025'in aynı günlerinden çok daha bulutluydu (Katar'da
bulutlu hücre-gece payı 2025'te %8, 2026'da %21). Kuzey Irak'ın ilk kapanmadaki artışı +%7 ile +%22 arasında.

## İki kapanma ve Basra testi

İki kapanma aynı modellerde karşılaştırılır: fark 100 x [exp(β2 − β1) − 1] olarak hesaplanır ve aynı modelin kovaryansıyla sınanır. Hürmüz
koridorunda ikinci kapanmanın ilkinden farkı −%4 (p = 0,91); aralık geniş, veri ne fark ne de eşitlik gösteriyor. Güney Irak'ta
tespit edilen alev gücü ikinci kapanmada ilkinden %295, iki kapanmanın ilk 30 gününde %136 yüksekti; bu karşılaştırma üretimin değil,
tespit edilen alev gücünün. Altıncı radar kutusu Basra açık deniz terminallerini ve demirleme alanını kapsıyor (48,70-49,20 D, 29,40-29,85 K;
`KONUMLAR` içinde `basra`); iki yörünge yönü ayrı sabit hedef maskeleriyle sayıldı, inen seri 35 ve 108 numaralı göreli yörüngeleri
(`veri/sar/goreli_yorunge_basra.csv`) göreli yörünge kuklasıyla birleştiriyor. Bütün ikinci kapanmada doğrusal eğilim çıkan yörüngede 30 günde
+2,0 gemi [−0,5, +4,5], inen yörüngede −1,1 [−6,2, +3,9]. Yalnız Ağustos'ta inen yörüngede
eğilim 30 günde +12,8 gemi (HC3 aralığı −6,0 ile +31,5); geçici bir birikme dışlanamıyor. `birikim`,
`birikim_agustos` ve `birikim_goreli` alanlarında OLS (HC1, HC3) ve Poisson eğilimleri var. Gemi sayısı yük stoku değildir. Ölçek için: SOMO'nun
Ağustos 2026 Basra yüklemesinin her 2 milyon varili ayrı bir VLCC olarak kutuda kalsaydı sayı 30 günde yaklaşık 34 artardı. Çalışma
notu, IMF PortWatch'un AIS'e dayalı yük tahminlerini, Basra petrol terminalinde (liman verisi) ve boğazda (`capacity_tanker`: ton cinsinden
tahmini yük, DWT değil), Basra'nın yüklemesiyle karşılaştırıyor ve Irak'ın Ağustos ihracatına ilişkin bağımsız haberlere atıf yapıyor. Bu karşılaştırma paketin dağıtmadığı PortWatch verisini kullandığı için `analiz17_iki_kapanma.py`,
raporun çalışma kitabı yoksa onu ve Dated Brent farkını atlar. Kaynak: `sonuclar/a17_iki_kapanma.json`.

## Yöntem

1. Her konum için bir kutu (`sar_gemi_sayimi.py` içinde `KONUMLAR`). Görüntüler Microsoft Planetary Computer'daki Sentinel-1 RTC
   koleksiyonundan; hesap gerekmez.
2. VH polarizasyonu 20 m (Füceyre, Yanbu) ya da 40 m (Hürmüz, Ras Tanura, Bâbülmendep) çözünürlükte okunur, dB'ye çevrilir.
3. Sabit parlak hedefler (kara, liman yapıları, platformlar, şamandıralar) maskelenir: savaş öncesi medyanı eşiği aşan pikseller ve
   Natural Earth karası (300 m tampon).
4. Eşik: −18 dB ile deniz medyanı + 9 dB'den büyük olanı. 8-komşulu kümelerden 20 m'de en az 3, 40 m'de en az 2 piksellik olanlar gemi;
   20 m'de en az 12, 40 m'de en az 3 piksellik olanlar büyük yankı sayılır.
5. Kutunun deniz alanının %90'ından azını kapsayan görüntüler atlanır.
6. Sentinel-1A görüntülerinde deniz yüzeyinin medyanı 1C ve 1D'dekinden yaklaşık 2,5 dB yüksek ve 1A aynı dönemlerde daha az gemi
   sayıyor; ikinci kapanmada da görüntü vermedi. Dönem medyanları yalnız 1C ve 1D ile; Poisson modeli bütün görüntüleri uydu kuklalarıyla
   kullanır. Farkın mekanizmasını (gürültü tabanı ya da işleme) iddia etmiyoruz.
7. Çıkan ve inen yörüngeler geliş açıları farklı olduğu için ayrı tutulur.

İlgili çalışma: Cao ve diğerleri (2026), *The Innovation* 7(6), 101367 (doi:10.1016/j.xinn.2026.101367), Sentinel-1 ve AIS ile kapanmanın
ilk haftasında boğazın ana geçiş kanalında gemi sayısının %97 düştüğünü ölçtü. Bu paket 14 ayı, iki kapanmayı ve dört noktayı daha kapsıyor.

Sınırlar: gemi tipi ayırt edilmez (tanker, konteyner gemisi, savaş gemisi aynı sayılır); yan yana demirli gemiler tek küme görünebilir;
rüzgârlı günlerde yanlış alarm artabilir. Yanlış ve kaçırılan tespitler etiketli görüntülerle henüz ölçülmedi. Her görüntüye aynı sayıda
eklenen yanlış tespitler basit oranları bire (yüzde değişimleri sıfıra) yaklaştırır; döneme, uyduya, bakış geometrisine ya da deniz durumuna göre değişen hatalar ise değişimleri iki
yönde de saptırabilir.

## Yeniden üretme

Python 3.12. Kod rapor için nasıl çalıştıysa öyle yayımlandı.

```
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt

python analiz7_sar.py            # dönem özetleri ve regresyon -> sonuclar/a7_sar_ozet.json (saniyeler)
python analiz9_sar_yogunluk.py   # Füceyre yoğunluk haritaları -> pdf/gorsel/yogunluk_*.png, sonuclar/a9_yogunluk.json
bash run_radar.sh                # bütün görüntüleri Planetary Computer'dan yeniden sayar (birkaç saatlik indirme)
```

`python analiz15_alev.py` alev serilerini ve etkileri pakete dahil FIRMS dosyalarından yeniden üretir (`sonuclar/a15_alev.json`, saniyeler).
`firms_cek.py` bu dosyaları yeniden indirir; `FIRMS_MAP_KEY` ortam değişkeninde ücretsiz bir NASA FIRMS anahtarı ister.

`python analiz16_kapsama.py` bulut ve kapsama denetimini pakete dahil Level-2 özütlerinden (`veri/firms/kapsama/l2_*`) yaklaşık bir dakikada
yeniden üretir (`sonuclar/a16_kapsama.json` ve `veri/firms/kapsama/` altındaki tablolar). `kapsama_tam.py` Level-2 dosyalarını yeniden indirip
işler (3.111 granül çifti, yaklaşık 160 GB akış; her dosya işlendikten sonra silinir) ve `.env` dosyasında ücretsiz bir NASA Earthdata
anahtarı ister (`EARTHDATA_TOKEN=...`). `kapsama_kunye.py` granül künyesini, `kapsama_pilot.py` 59 gecelik pilotu yeniden üretir.

`python analiz17_iki_kapanma.py` iki kapanmanın karşılaştırmasını ve Basra testini yaklaşık bir dakikada yeniden üretir
(`sonuclar/a17_iki_kapanma.json`). Basra sayımları `veri/sar/gemi_sayimi_basra_40m.csv` (çıkan) ve `veri/sar/gemi_sayimi_basra_40m_desc.csv`
(inen) dosyalarında; yeniden sayım komutları `run_radar.sh` içinde. Basra'nın iki sabit hedef maskesinin kurulum künyesi var
(`veri/sar/basra_40m_sabit_maske_kunye.json` ve `veri/sar/basra_40m_desc_sabit_maske_kunye.json`: ızgara, dönüşüm, kullanılan ve atlanan
sahneler, parametreler, SHA-256). Çıkan yörünge çalıştırması iki yönün sahnelerini okudu; maske ızgarası (EPSG 32638) yalnız 174 numaralı
göreli yörüngeye uyduğu için medyana bu yörüngenin 6 savaş öncesi sahnesi girdi (inen yörünge maskesinde 31). Kayıtlı
sahnelerden 27 Eylül 2026'da yeniden kurulan iki maske de diskteki maskelerle piksel piksel aynı çıktı; `run_radar.sh`'teki ortam
değişkenleriyle `python sar_gemi_sayimi.py --maske-kunye basra` bu denetimi tekrarlar.

Irak ihracatı iki düzeyde yeniden üretilir. Analizler `veri/rafineri/irak_somo_cikis_aylik.csv` dosyasını okur. Bu tablo `python somo_cek.py`
ile SOMO'nun aylık PDF raporlarından ve ihracat grafiğinden üretilir; paket bu kaynakları dağıtmaz: betik onları indirir,
`veri/rafineri/irak_somo/kaynaklar.csv` her kaynağın adresini, erişim tarihini ve SHA256 değerini tutar. Nisan-Haziran 2026 için PDF rapor
yok; bu ayların çıkış noktaları grafik kategorilerinin Mart 2026 raporundaki eşlemesine göre varsayıldı (`kaynak` sütununda yazılı).

`sar_gemi_sayimi.py` çıktı CSV'sinde bulunan görüntüleri atlar; baştan saymak için pakete dahil `veri/sar/gemi_sayimi_*.csv`
dosyalarını önce başka yere taşıyın. Sabit hedef maskeleri (`veri/sar/*_sabit_maske.npy`) yeniden sayımın yayımlanan sayıları
vermesi için pakette; silinirse savaş öncesi görüntülerden yeniden üretilir.

**GDELT.** `bq/haber_endeksi_v2.sql` ve `bq/olay_endeksi_v2.sql` Google BigQuery konsolunda (Standard SQL) çalıştırılıp CSV olarak
kaydedilir. 25 Eylül 2026'daki kuru çalıştırmaya göre sorgular yaklaşık 0,63 TB ve 0,04 TB veri tarıyor; BigQuery'nin ücretsiz katmanı
ayda 1 TB sorguyu karşılıyor. `gdelt_doc_cek.py` ücretsiz GDELT DOC 2.0 API'sini kullanır (anahtar gerekmez; sık istekleri sınırlar).
Tehdit ve eylem testi (`analiz6_gdelt.py`) günlük Brent getirisi de ister; raporda ICE Brent ön vade kullanıldı ve bu veri yeniden
dağıtılamadığı için pakette yok.

## Pakette olmayanlar

Piyasa fiyatları (ICE, LME, LBMA, CME ve diğerleri) lisanslı olduğu için yeniden dağıtılamaz; raporun fiyat analizleri (vade eğrisi,
altın ayrıştırması, alüminyum primleri, senaryo simülasyonu) bu pakete girmedi. MetalTakip'in Türkiye iç piyasa primi göstergeleri
şirkete aittir.

Sütun adlarının açıklaması [README.md](README.md) içindeki sözlükte.

## Lisans

- Kod: MIT ([LICENSE](LICENSE)).
- Bu depodaki veri, sonuç ve görseller: CC BY 4.0 ([LICENSE-DATA.md](LICENSE-DATA.md)).
- Kaynaklar: Contains modified Copernicus Sentinel data 2025-2026; Sentinel-1 RTC işlemesi Catalyst'e ait, dağıtım Microsoft Planetary
  Computer (CC BY 4.0). Haber ve olay verisi [GDELT Project](https://www.gdeltproject.org/). Kara sınırları Natural Earth (kamu malı).
  Yangın tespitleri NASA FIRMS'ten (NASA ESDIS). Irak ihracat rakamları SOMO'dan. Olay listesini MetalTakip açık kaynaklardan derledi.

## Atıf

[CITATION.cff](CITATION.cff). Önerilen biçim:

> Dağlı, Y. O. ve Özpolat, E. (2026). *Closing Straits: satellite radar ship counts, gas flaring and GDELT news indices for the 2026 Strait of Hormuz closures* (sürüm 1.4.1) [veri ve kod]. Zenodo.
> https://doi.org/10.5281/zenodo.23009094

Paketin kayıtlı adı İngilizcedir; atıfta bu adı kullanın. Veri ve kod için paketi, bulgular için raporu ya da çalışma notunu kaynak gösterin:

- Rapor: Dağlı, Y. O. ve Özpolat, E. (2026). *Kapanan Boğazlar: Hürmüz ve Bâbülmendep şokları petrol, gaz, metal ve altın piyasalarını
  nasıl yeniden fiyatladı?* MetalTakip Sektörel Araştırma, Eylül 2026. https://www.metaltakip.com/arastirma/kapanan-bogazlar
- Çalışma notu: Dağlı, Y. O. ve Özpolat, E. (2026). *Measuring a chokepoint closure from space: radar ship counts and gas flaring during the 2026 Strait of Hormuz closures.* Working paper, version 1.6, Eylül 2026. https://doi.org/10.5281/zenodo.23009748

Sürümler: 1.4.1 ilk açık sürümdür ve çalışma notunun 1.6 sürümüne karşılık gelir. 1.0-1.4.0 sürümleri yayımlanmadı; her biri çalışma notunun
bir taslağına karşılık gelir (1.4.0 ile 1.5, 1.3.0 ile 1.4, 1.2.0 ile 1.3, 1.1.0 ile 1.2; 1.0 ve 1.1 taslaklarıyla denetlenen durumlar 1.0
etiketini taşıyordu). Taslak 1.0-1.2 üzerindeki üç denetim turu radar ve alev analizlerini kapsadı: ChatGPT (OpenAI) bunları koddan ve veriden
yeniden üretti, Claude (Anthropic) her bulguyu veriyle sınadı, Ebru Özpolat denetimleri inceledi. Üçüncü turun önerdiği protokolü izleyen tam
bulut ve kapsama denetimi (1.3.0) bu turlardan sonra eklendi. 1.4.0 üzerindeki dördüncü, dar kapsamlı tur iki kapanmanın karşılaştırmasını ve
Basra testini kapsadı: ChatGPT tahminleri yeniden üretti, Claude bulguları veriyle sınadı. Düzeltmeleri 1.4.1'de: PortWatch'un tanker alanı DWT
değil tahmini yük olarak okunuyor, Basra eğilimi Ağustos için ve HC3 ile Poisson'la da veriliyor, Basra maskelerinin kurulum künyesi var.
1.4.1'in yayın öncesi derlemesi üzerindeki kapanış kontrolü bu düzeltmeleri, çalışma notu ve rapordaki Basra karşılaştırmasının yorumu dışında
uygulanmış buldu; o yorum yayından önce yeniden yazıldı, eğilim tahminleri de tabloların tek kez yuvarlanması için dört ondalıkla saklanıyor. Bu
tur ham radar sayımını bağımsız olarak yeniden üretmedi; PortWatch liman verisi ve Basra'ya ilişkin haberler (bu pakette yok) ondan sonra eklendi.
Analizleri başka bir ortamda çalıştırmak p değerlerinin son basamaklarını (10⁻¹³ düzeyinde) değiştirebilir; tahminler ve aralıklar değişmez.

Bilgilendirme amaçlıdır; yatırım tavsiyesi değildir.
