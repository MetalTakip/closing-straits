#!/usr/bin/env bash
# Recount ships in every Sentinel-1 RTC image with the settings used in the report, then rebuild the summaries.
# Downloads from Microsoft Planetary Computer (no account needed); expect several hours.
# sar_gemi_sayimi.py skips images already in its output CSV: move veri/sar/gemi_sayimi_*.csv away first for a full recount.
set -euo pipefail
cd "$(dirname "$0")"

BAS=2025-07-01
BIT=2026-09-22

python sar_gemi_sayimi.py fuceyre "$BAS" "$BIT"                                                       # 20 m, both orbits
SAR_RES=40 SAR_YORUNGE=descending SAR_EK=_40m python sar_gemi_sayimi.py hurmuz "$BAS" "$BIT"          # 40 m, descending
SAR_RES=40 SAR_YORUNGE=descending SAR_EK=_40m python sar_gemi_sayimi.py rastanura "$BAS" "$BIT"       # 40 m, descending
python sar_gemi_sayimi.py yanbu "$BAS" "$BIT"                                                         # 20 m (only ascending images cover the box)
SAR_RES=40 SAR_YORUNGE=ascending SAR_EK=_40m_asc python sar_gemi_sayimi.py babulmendep "$BAS" "$BIT"  # 40 m, ascending
SAR_RES=40 SAR_EK=_40m python sar_gemi_sayimi.py basra "$BAS" "$BIT"                                  # 40 m; reads both directions, only relative orbit 174 fits the mask grid
SAR_RES=40 SAR_YORUNGE=descending SAR_EK=_40m_desc python sar_gemi_sayimi.py basra "$BAS" "$BIT"      # 40 m, descending (relative orbits 35 and 108)
python sar_gemi_sayimi.py --goreli basra                                                              # relative orbit of each Basra scene
SAR_RES=40 SAR_EK=_40m python sar_gemi_sayimi.py --maske-kunye basra                                  # rebuild the Basra masks from the
SAR_RES=40 SAR_YORUNGE=descending SAR_EK=_40m_desc python sar_gemi_sayimi.py --maske-kunye basra      #   recorded scenes and compare

python sar_tespit_konum.py                                                                            # Fujairah detection positions
SAR_RES=40 SAR_YORUNGE=descending SAR_EK=_40m python sar_hurmuz_gorsel.py                             # Hormuz image pair
python analiz7_sar.py                                                                                 # period summaries, regression
python analiz9_sar_yogunluk.py                                                                        # Fujairah density maps
python analiz17_iki_kapanma.py                                                                        # two closures compared, Basra test
