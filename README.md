# Closing Straits: satellite radar ship counts, gas flaring and GDELT news indices for the 2026 Strait of Hormuz closures

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23009094.svg)](https://doi.org/10.5281/zenodo.23009094)

Replication data and code for the satellite radar, gas-flare and news analyses behind **Closing Straits: How the Hormuz and Bab el-Mandeb
shocks repriced oil, gas, metals and gold** (MetalTakip Sector Research, September 2026; Turkish edition: *Kapanan Boğazlar*).
The report is at [metaltakip.com](https://www.metaltakip.com/en/research/closing-straits).

Authors: Yusuf Oskay Dağlı ([ORCID 0009-0008-5513-9217](https://orcid.org/0009-0008-5513-9217)) and Ebru Özpolat
([ORCID 0009-0006-3864-142X](https://orcid.org/0009-0006-3864-142X)), MetalTakip.

Türkçe açıklama: [README.tr.md](README.tr.md)

## What is here

- **Ship counts from 407 Sentinel-1 radar images** over five locations: the Hormuz transit corridor, the Fujairah anchorage,
  Ras Tanura, Yanbu and Bab el-Mandeb, from July 2025 to 22 September 2026. One row per image.
- **The code that produced them**, the period summaries and the robustness regression used in the report, and the two radar
  images of the Hormuz corridor shown in the report.
- **Gas flares from VIIRS night-time fire detections** (NASA FIRMS, September 2024 to September 2026) at persistent flare sites across Iraq,
  Kuwait, Iran, the Gulf states and Oman, with the code that turns them into regional series and period effects.
- **A cloud and coverage check for the flares**: every NOAA-21 night pass over the persistent sites from NASA's Level-2 fire and cloud
  masks (3,111 granule pairs), as cell × night × pass tables, with the code that downloads and models them.
- **The two closures compared and a test at Basra** (`analiz17_iki_kapanma.py`): differences between the closures within the same
  models, and ship counts in a sixth box around the Basra offshore terminals on both orbit directions (182 images).
- **Daily GDELT indices, 2017 to 2026**: news documents naming Hormuz or Bab el-Mandeb (GKG), and threat and action events located
  in boxes around the two straits (Events), with the BigQuery SQL that built them.
- **The 39 events** used in the report's event study (Appendix A), in Turkish and English.

## Key results

Ships per radar image, median by period, Sentinel-1C and 1D images only. The last column is the change in the second closure
against the pre-war period from a Poisson count model with period and satellite dummies using all images (heteroskedasticity-robust standard
errors). Source: `sonuclar/a7_sar_ozet.json`.

| Location | Orbit, resolution | 1C/1D images in periods | Pre-war | First closure | Second closure | Change, second closure (95% CI) |
|---|---|---:|---:|---:|---:|---|
| Hormuz corridor | descending, 40 m | 44 | 40.5 | 10 | 9.5 | −81% (−88.1% to −70.1%) |
| Fujairah | descending, 20 m | 40 | 177 | 171 | 254 | +46% (+36.7% to +56.9%) |
| Fujairah | ascending, 20 m | 44 | 157 | 165 | 231 | +43% (+30.7% to +55.4%) |
| Yanbu | ascending, 20 m | 37 | 23 | 29 | 33 | +32% (+3.2% to +69.9%) |
| Bab el-Mandeb | ascending, 40 m | 47 | 46 | 50 | 54.5 | +22% (−1.9% to +52%) |
| Ras Tanura | descending, 40 m | 44 | 115 | 118.5 | 108.5 | −7% (−14.4% to +0.5%) |

Periods: pre-war 1 July 2025 to 27 February 2026; first closure 2 March to 13 June 2026; reopening 18 June to 10 July 2026;
second closure 13 July to 22 September 2026.

Robustness checks are stored under the `saglamlik` key of each orbit: a log(1 + ships) model (`log`), Sentinel-1C and 1D images only
(`tutarli`), Sentinel-1C only (`yalniz_1c`, the only satellite observed in every period), large returns only (`buyuk`: at least 12 pixels
at 20 m or 3 at 40 m; an area threshold, not a measured ship length), Newey-West standard errors with four lags (`HAC`), standard errors
clustered by month (`ay`, with t-based limits `alt_t` and `ust_t`) and leave-one-image-out estimates for the second closure (`tek_cikar`).
Holm-adjusted p-values for five locations x two closures are stored as `p_holm`. Across the point estimates, the fall in the Hormuz corridor
lies between 75% and 85%. After the Holm correction only the two Hormuz declines and Fujairah's second-closure rise remain significant at
the 5% level; the second-closure rises at Yanbu (six images) and Bab el-Mandeb (+22%, lower 95% limit −2%) do not.

In the first closure, AIS-based transit counts (IMF PortWatch) fell 93%, and the EIA's estimate of oil flows through the strait was
77% lower in the second quarter of 2026 than in the fourth quarter of 2025; the radar count in the corridor fell 80%. Radar sees ships whether or not they broadcast
AIS, but the count is not a transit count: ships waiting, moving slowly or on military duty in the box are counted too. The three series measure
different things over different periods, so the gaps between them cannot be read as lost AIS signals or as volumes.

## Gas flares

The measure is the fire radiative power (FRP) detected at persistent flare sites per night. Pixels are summed within a pass; when several
passes saw a site, we keep the pass with the smallest mean along-track pixel size, which grows with the viewing angle (the along-scan size does
not track the angle steadily, because of VIIRS pixel aggregation). Changes against the pre-war period come from a Poisson
pseudo-maximum-likelihood (PPML) regression on period and calendar-month dummies, with HAC standard errors.
`log`: log(1 + FRP) OLS with the same dummies. Bold: the unadjusted 95% confidence interval excludes zero. Holm-adjusted p-values (family: the ten regions below x two closures) are stored as `p_holm`; the reopening period, other regions and the south-north contrast are outside that family. Source: `sonuclar/a15_alev.json`.

| Region | Export route | Pre-war, MW | First closure (NOAA-21 / Suomi NPP / log) | Second closure (NOAA-21 / Suomi NPP / log) |
|---|---|---:|---|---|
| Iraq (south) | Basra terminals | 1,053 | **−84%** (Holm <0.001) / **−86%** / **−81%** | **−37%** (Holm <0.001) / **−36%** / **−41%** |
| Kuwait | Mina al-Ahmadi | 24 | −4% (Holm 1.000) / −2% / −1% | **+71%** (Holm 0.005) / **+86%** / **+58%** |
| Qatar | Ras Laffan | 67 | **−31%** (Holm 0.141) / −38% / **−58%** | −8% (Holm 0.770) / **−18%** / **−13%** |
| Bahrain | Sitra | 8 | **−78%** (Holm <0.001) / **−73%** / **−69%** | **−74%** (Holm 0.006) / **−78%** / **−72%** |
| Western Iran (mostly Khuzestan) | Kharg | 1,052 | −15% (Holm 0.498) / −7% / −11% | **−26%** (Holm <0.001) / **−28%** / **−27%** |
| Iran (Asaluyeh) | Asaluyeh | 206 | −4% (Holm 1.000) / −15% / −27% | **−21%** (Holm <0.001) / **−23%** / **−27%** |
| Saudi Arabia | East-West pipeline, Yanbu | 179 | +10% (Holm 1.000) / +20% / −3% | +14% (Holm 1.000) / +21% / +14% |
| UAE | Habshan-Fujairah pipeline | 62 | **+304%** (Holm <0.001) / **+348%** / **+205%** | **+201%** (Holm <0.001) / **+206%** / **+203%** |
| Oman | Gulf of Oman | 81 | −5% (Holm 1.000) / −10% / −10% | **+43%** (Holm 0.036) / **+41%** / **+53%** |
| Iraq (north) | Kirkuk-Ceyhan pipeline | 273 | **+23%** (Holm 0.006) / +14% / +5% | **+30%** (Holm <0.001) / **+22%** / **+29%** |

At Iraq's southern fields, detected flaring fell while exports from Basra collapsed (SOMO: 3.33 million barrels a day in February, 79,000 in May); at the
northern fields, whose crude bypasses the strait through the Kirkuk-Ceyhan pipeline, it did not fall. Detected flare power is not proportional
to output: a plant that cannot sell
its gas may burn more of it, which may explain the rises in the UAE and Kuwait (not tested). The smaller fall in Qatar is not robust: it is
not significant after the Holm correction or on Suomi NPP, it varies with the site definition, and no persistent site was detected on
13% of first-closure nights (2% before the war); the Level-2 check below shows that most Qatari sites were under cloud on
each of those nights. Related work: with VIIRS Nightfire, Zhizhin and Bazilian (2026, Payne Institute) found that flaring around Basra fell sharply in the first ten
days of March 2026, with smaller falls in the UAE, Qatar and Saudi Arabia and a rise in Iran; our series give the same signs for those days
(`erken_mart` in `sonuclar/a15_alev.json`). Export figures by outlet (Basra, Ceyhan, Khor al-Zubair) come from `somo_cek.py`, which reads SOMO's monthly PDF reports.

## Cloud and coverage check

FIRMS points carry no cloud or coverage information. To check the flare results, every persistent cell was matched with every NOAA-21
night pass in the Level-2 files: the 375 m fire mask of VJ214IMG and the 750 m cloud mask CLDMSK_L2_VIIRS_NOAA21, which also carries
geolocation and viewing angle (3,111 granule pairs on 748 nights, passes without a detection included). A cell-pass is valid when
at least half of its 375 m pixels are water, clear land or fire (fire pixels count as clear: the product looks for fire only where it sees no
cloud), and cloudy when at least half are cloud. Among a cell's valid passes the one with the smallest viewing angle is used, without
reference to FRP, and residual bow-tie duplicates (QA bit 22) are dropped. The main model is then re-estimated on nights when at least 80%,
90% or 95% of a region's cells were validly observed, with the FIRMS measure and the Level-2 measure, and a PPML with cell fixed effects and
Driscoll-Kraay standard errors is fitted on valid cell-nights. The independent cloud mask is not used as a filter: it calls 25% of valid
cell-passes without a detection cloudy and 94% of those above 50 MW, so selecting on it would drop strong flares. Source:
`sonuclar/a16_kapsama.json`; tables in `veri/firms/kapsama/`.

| Region | First closure: all nights | First: 90% observed, FIRMS / L2 | First: cell level | Second closure: all nights | Second: 90% observed, FIRMS / L2 | Second: cell level |
|---|---|---|---|---|---|---|
| Iraq (south) | −84% | −83% / −85% | −85% | −37% | −37% / −36% | −36% |
| Kuwait | −4% | −6% / −11% | −7% | +71% | +79% / +85% | +77% |
| Qatar | −31% | −15% / −19% | −27% | −8% | −8% / −10% | −11% |
| Bahrain | −78% | −70% / −69% | −71% | −74% | −75% / −71% | −69% |
| Western Iran (mostly Khuzestan) | −15% | −18% / −18% | −14% | −26% | −23% / −23% | −25% |
| Iran (Asaluyeh) | −4% | +7% / +10% | +9% | −21% | −19% / −18% | −17% |
| Saudi Arabia | +10% | +28% / +27% | +25% | +14% | +16% / +10% | +9% |
| UAE | +304% | +325% / +324% | +327% | +201% | +199% / +200% | +204% |
| Oman | −5% | +3% / +0% | −2% | +43% | +32% / +41% | +40% |
| Iraq (north) | +23% | +14% / +13% | +17% | +30% | +30% / +27% | +25% |

On well-observed nights the first-closure change in southern Iraq ranges from −85% to −81% across twelve specifications, and
the south-north difference from −87% to −86%. Qatar's first-closure change is −27% to −15% and not significant after a
Holm correction in any specification; the first closure was cloudier on the southern Gulf coast than the same days of 2025 (8% of
Qatari cell-nights cloudy in 2025, 21% in 2026). The first-closure rise in northern Iraq ranges from +7% to +22%.

## Two closures and the Basra test

The closures are compared within the same models: 100 × [exp(β2 − β1) − 1], tested with the same covariance. In the Hormuz corridor the
second closure differs from the first by −4% (p = 0.91); the interval is wide, so the data show neither a difference nor equality.
Detected flare power in southern Iraq was 295% higher in the second closure than in the first, and 136% higher over the first 30 days
of each; this compares detected flaring, not output. A sixth radar box covers the Basra offshore terminals and their anchorage
(48.70-49.20°E, 29.40-29.85°N; `basra` in `KONUMLAR`), counted on both orbit directions with separate fixed-target masks; the descending series
combines relative orbits 35 and 108 (`veri/sar/goreli_yorunge_basra.csv`) with a relative-orbit dummy. Over the whole second closure the linear
trend is +2.0 ships per 30 days [−0.5, +4.5] on the ascending orbit and −1.1 [−6.2, +3.9] on the
descending orbit. Within August alone the descending trend was +12.8 ships per 30 days (HC3 interval −6.0 to +31.5),
so a temporary build-up is not ruled out; `birikim`, `birikim_agustos` and `birikim_goreli` hold the OLS (HC1, HC3) and Poisson trends. Ship
counts are not cargo inventories. For scale, had every 2 million barrels of SOMO's August 2026 Basra loadings stayed in the box on a separate
VLCC, the count would have risen by about 34 per 30 days. The working paper also sets IMF PortWatch's AIS-based cargo estimates, at the
Basrah Oil Terminal (port data) and in the strait (`capacity_tanker`, an estimated cargo in tonnes, not deadweight), against Basra's loadings
and cites independent reports of Iraq's August exports. That comparison uses
PortWatch data, which this package does not redistribute, so `analiz17_iki_kapanma.py` skips it, together with the Dated Brent spread, when
the report's workbook is absent. Source: `sonuclar/a17_iki_kapanma.json`.

## Method in brief

1. A box is defined for each location (`KONUMLAR` in `sar_gemi_sayimi.py`). Scenes come from the Sentinel-1 RTC collection on
   Microsoft Planetary Computer; no account is needed.
2. The cross-polarised VH channel is read at 20 m (Fujairah, Yanbu) or 40 m (Hormuz, Ras Tanura, Bab el-Mandeb) and converted to dB.
3. Fixed bright targets (land, port structures, platforms, buoys) are masked: pixels whose pre-war median exceeds the threshold,
   plus Natural Earth land with a 300 m buffer.
4. Threshold: the larger of −18 dB and the sea median + 9 dB. Eight-connected clusters of at least 3 pixels at 20 m (2 at 40 m) count
   as ships; at least 12 pixels at 20 m (3 at 40 m) as large returns.
5. Scenes covering less than 90% of the box's sea area are skipped.
6. In Sentinel-1A images the sea-surface median is about 2.5 dB higher than in 1C and 1D, and 1A counts fewer ships in the same
   periods; 1A provided no images in the second closure. Period medians therefore use 1C and 1D images only; the Poisson model uses
   all images with satellite dummies. We do not claim a mechanism for the difference (noise floor or processing).
7. Ascending and descending orbits are kept apart because their incidence angles differ.

Related work: Cao et al. (2026), *The Innovation* 7(6), 101367 (doi:10.1016/j.xinn.2026.101367), measured a 97% fall in ships in the
strait's main shipping channel in the first week of the closure with Sentinel-1 and AIS. This package covers 14 months, both closures and
four more locations.

Limits: ship type is not identified (tankers, container ships and warships count alike); ships anchored side by side can merge into
one cluster; windy days can add false alarms. False and missed detections have not yet been measured against labelled images. The same
number of false detections in every image would pull simple ratios towards one (and percentage changes towards zero), but errors that vary with the period, the satellite, the viewing
geometry or the sea state can bias the changes in either direction.

## Reproduce

Python 3.12. Folder and variable names are Turkish because the code is published exactly as it ran for the report (see the
glossary below).

```
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt

python analiz7_sar.py            # period summaries and regression -> sonuclar/a7_sar_ozet.json (seconds)
python analiz9_sar_yogunluk.py   # Fujairah density maps -> pdf/gorsel/yogunluk_*.png, sonuclar/a9_yogunluk.json
bash run_radar.sh                # recount every image from Planetary Computer (several hours of downloads)
```

`python analiz15_alev.py` rebuilds the flare series and effects from the included FIRMS files (`sonuclar/a15_alev.json`, seconds).
`firms_cek.py` downloads them again; it needs a free NASA FIRMS MAP_KEY in the `FIRMS_MAP_KEY` environment variable.

`python analiz16_kapsama.py` rebuilds the cloud and coverage check from the included Level-2 extracts (`veri/firms/kapsama/l2_*`) in about a
minute (`sonuclar/a16_kapsama.json` and the tables in `veri/firms/kapsama/`). `kapsama_tam.py` downloads and processes the Level-2 files again
(3,111 granule pairs, about 160 GB streamed; each file is deleted after use) and needs a free NASA Earthdata token in `.env`
(`EARTHDATA_TOKEN=...`). `kapsama_kunye.py` rebuilds the granule register and `kapsama_pilot.py` the 59-night pilot.

`python analiz17_iki_kapanma.py` rebuilds the comparison of the two closures and the Basra test (`sonuclar/a17_iki_kapanma.json`, about a
minute). The Basra counts are in `veri/sar/gemi_sayimi_basra_40m.csv` (ascending) and `veri/sar/gemi_sayimi_basra_40m_desc.csv` (descending);
`run_radar.sh` lists the commands that recount them. Each Basra fixed-target mask has a build record (`veri/sar/basra_40m_sabit_maske_kunye.json`
and `veri/sar/basra_40m_desc_sabit_maske_kunye.json`: grid, transform, scenes used and skipped, parameters, SHA-256). The ascending run read
the scenes of both directions; its mask grid (EPSG 32638) fits relative orbit 174 only, so 6 pre-war scenes of that orbit entered the
median, against 31 for the descending mask. Rebuilt from the recorded scenes on 27 September 2026, both masks matched the stored ones
pixel for pixel; `python sar_gemi_sayimi.py --maske-kunye basra`, with the environment variables of `run_radar.sh`, repeats the check.

Iraqi exports are reproducible at two levels. The analyses read `veri/rafineri/irak_somo_cikis_aylik.csv`. That table is built by
`python somo_cek.py` from SOMO's monthly PDF reports and export chart, which this package does not redistribute: the script downloads them
and `veri/rafineri/irak_somo/kaynaklar.csv` lists each source with its address, access date and SHA256. April to June 2026 have no PDF
report; their outlets follow the March 2026 report's mapping of the chart categories (an assumption, stated in the `kaynak` column).

`sar_gemi_sayimi.py` skips images already in its output CSV, so move the included `veri/sar/gemi_sayimi_*.csv` files away before
recounting from scratch. The fixed-target masks (`veri/sar/*_sabit_maske.npy`) are included so that a recount reproduces the
published counts; delete them to rebuild the masks from the pre-war images.

**GDELT.** Run `bq/haber_endeksi_v2.sql` and `bq/olay_endeksi_v2.sql` in the Google BigQuery console (Standard SQL) and save the
results as CSV. A dry run on 25 September 2026 estimated about 0.63 TB and 0.04 TB of data processed; BigQuery's free tier
covers 1 TB of queries a month. `gdelt_doc_cek.py` uses the free GDELT DOC 2.0 API instead (no key; it rate-limits frequent
requests). `analiz6_gdelt.py`, the threat versus action test, also needs daily Brent returns, which are not included: the
report used ICE Brent front-month prices, which cannot be redistributed.

## What is not here

Market prices (ICE, LME, LBMA, CME and others) are licensed and cannot be redistributed, so the price analyses in the report
(forward curve, gold decomposition, aluminium premiums, scenario simulation) are not part of this package. MetalTakip's
Turkish market premium indicators are proprietary.

## Glossary

| Name | Meaning |
|---|---|
| `veri`, `sonuclar`, `pdf/gorsel` | data, results, figures |
| `gemi_sayimi_<location>.csv` | ship counts: `konum` location, `tarih` date (UTC), `saat` time (UTC), `sahne` scene ID, `yorunge` orbit direction, `kapsama` share of the box's sea area covered, `gemi` ships, `buyuk_gemi` large returns (Method in brief, step 4) |
| `tespit_konum_fuceyre.csv` | position of each detection at Fujairah: `boylam` longitude, `enlem` latitude, `alan_piksel` cluster size |
| `*_sabit_maske.npy` | fixed-target mask (land and permanent bright targets) |
| `goreli_yorunge_basra.csv` | relative orbit (`goreli_yorunge`) of each Sentinel-1 scene over the Basra box, from the STAC metadata |
| `haber_endeksi_v2.csv` | GKG news index: `gun` day, `belge` all documents, `tr_belge` Turkish-source documents, `h_ad` / `b_ad` documents naming Hormuz / Bab el-Mandeb, `h_url` / `b_url` strait name in the URL, `_tr` Turkish-source subset, `h_ton` / `b_ton` mean tone |
| `olay_endeksi_v2.csv` | Events index: `olay` events, `makale` articles, `h_` / `b_` Hormuz / Bab el-Mandeb box, `tehdit` threats (CAMEO 13, 15), `eylem` actions (CAMEO 17-20), `ton` tone, `goldstein` Goldstein scale; counts are article-weighted |
| `olaylar.csv` | event list: `kod` code, `bogaz` / `strait`, `tur` type, `yon` / `direction`, `olay_tarihi` event date, `t0` first trading day on which prices could react, `aciklama` / `description`, `ana_orneklem` in the main sample, `kaynak` added after web verification |
| `veri/firms/*.csv.gz` | NASA FIRMS VIIRS 375 m active fire detections as downloaded (NOAA-21 NRT; Suomi NPP SP and NRT) |
| `kalici_hucreler.csv` | persistent flare cells (0.01 degree): `hx`, `hy` cell index, `gece` pre-war nights with a detection, `bolge` region |
| `irak_somo_cikis_aylik.csv` | Iraq's monthly crude exports by outlet (SOMO monthly reports and export chart): `ay` month, `cikis` outlet (Basra, Ceyhan, Khor al-Zubair), `kalem` grade, `miktar_varil` barrels, `kaynak` source |
| `l2_hucre_granul.csv.gz` | Level-2 extract, one row per persistent cell and granule: `anahtar` cell key (hx × 100000 + hy), `n_m` 750 m pixels whose centre falls in the cell, `fm_0` to `fm_9` 375 m pixels by fire-mask class (0 not processed, 1 bow-tie deletion, 2 sun glint, 3 water, 4 cloud, 5 clear land, 6 unclassified, 7 to 9 fire), `cm_y1`, `cm_0` to `cm_3` 750 m pixels by cloud-mask class (−1 no result, 0 cloudy, 1 probably cloudy, 2 probably clear, 3 confident clear), `zenit` median viewing angle, `tarih` date (UTC), `gecis` granule |
| `l2_yangin_pikselleri.csv.gz` | fire pixels in persistent cells: `frp` MW, `vza` viewing angle, `guven` confidence, `bowtie` QA bit 22 |
| `hucre_gece_gecis.csv.gz` | cell × night × pass: `durum` state (`geçerli` valid, `bulutlu` cloudy, `belirsiz` uncertain, `durum yok` no 750 m pixel centre in the cell), `gecerli_pay` / `bulut_pay` valid / cloud share, `frp` FRP without bow-tie duplicates, `frp_tum` all FRP |
| `hucre_gece.csv.gz` | cell-night: `durum` (`granül yok` no granule, `kapsanmadı` not covered, `bulutlu`, `belirsiz`, `geçerli, tespit yok` / `var` valid without / with a detection), `gecis` chosen pass, `frp` Level-2 measure, `frp_firms` FIRMS measure |
| `bolge_gece.csv`, `granul_kunyesi.csv` | region-night table (`kapsama` share of validly observed cells, `l2`, `firms`); granule register (file names, production times, PGE version, SHA-256 of the fire file, provider MD5 of the cloud file) |
| `once`, `kapanma`, `acilis` | pre-war, closure, reopening |

## Licences

- Code: MIT, see [LICENSE](LICENSE).
- Data, results and figures in this repository: CC BY 4.0, see [LICENSE-DATA.md](LICENSE-DATA.md).
- Upstream sources: contains modified Copernicus Sentinel data 2025-2026; Sentinel-1 RTC processed by Catalyst and distributed by
  Microsoft Planetary Computer under CC BY 4.0. News and event data from the [GDELT Project](https://www.gdeltproject.org/).
  Land boundaries from Natural Earth (public domain). Fire detections from NASA's Fire Information for Resource Management System (FIRMS),
  part of NASA's Earth Science Data and Information System (ESDIS). Iraqi export figures from SOMO. The event list was compiled by MetalTakip
  from public sources.

## Citation

See [CITATION.cff](CITATION.cff). Suggested form:

> Dağlı, Y. O. and Özpolat, E. (2026). *Closing Straits: satellite radar ship counts, gas flaring and GDELT news indices for the 2026 Strait of Hormuz closures* (version 1.4.1) [data and code]. Zenodo.
> https://doi.org/10.5281/zenodo.23009094

Please cite the package for the data and code, and the report or the working paper for their findings:

- Report: Dağlı, Y. O. and Özpolat, E. (2026). *Closing Straits: How the Hormuz and Bab el-Mandeb shocks repriced oil, gas, metals and gold.*
  MetalTakip Sector Research, September 2026. https://www.metaltakip.com/en/research/closing-straits
- Working paper: Dağlı, Y. O. and Özpolat, E. (2026). *Measuring a chokepoint closure from space: radar ship counts and gas flaring during the 2026 Strait of Hormuz closures.* Working paper, version 1.6, September 2026. https://doi.org/10.5281/zenodo.23009748

Versions: 1.4.1 is the first public release and accompanies version 1.6 of the working paper. Versions 1.0 to 1.4.0 were not published;
each accompanied a draft of the working paper (1.4.0 with 1.5, 1.3.0 with 1.4, 1.2.0 with 1.3, 1.1.0 with 1.2; states audited with
drafts 1.0 and 1.1 were labelled 1.0). Three audit rounds, on drafts 1.0 to 1.2, covered the radar and flare analyses: ChatGPT (OpenAI)
reproduced them from the code and data, Claude (Anthropic) checked each finding against the data, and Ebru Özpolat reviewed the audits.
The full cloud and coverage check (1.3.0), which follows a protocol proposed in the third round, was added after these rounds. A fourth,
narrower round on 1.4.0 covered the comparison of the two closures and the Basra test: ChatGPT reproduced the estimates and Claude checked
the findings against the data. Its corrections are in 1.4.1: PortWatch's tanker field is read as a cargo estimate rather than deadweight,
the Basra trend is also reported for August and with HC3 and Poisson, and the Basra masks have build records. A closing check on a
pre-release build of 1.4.1 found these corrections applied except in the interpretation of the Basra comparison in the working paper and
the report, which was rewritten before release; the trend estimates are now stored with four decimals so that tables round only once. That
round did not re-run the raw radar counts independently, and the PortWatch port data and news reports on Basra, which are not in this
package, were added after it. Rerunning the analyses in
another environment can change p-values in the last digits (about 10⁻¹³); estimates and intervals are unaffected.

For information only; not investment advice.
