# Data licence

The data, results and figures in this repository (the `veri/`, `sonuclar/` and `pdf/gorsel/` folders) are licensed under the
Creative Commons Attribution 4.0 International licence (CC BY 4.0): https://creativecommons.org/licenses/by/4.0/

Attribution: Dağlı, Y. O. and Özpolat, E. (2026). *Closing Straits: satellite radar ship counts, gas flaring and GDELT news indices for the 2026 Strait of Hormuz closures* (version 1.4.1) [data and code]. MetalTakip. https://doi.org/10.5281/zenodo.23009094

They are derived from the following sources, whose own terms also apply:

| Source | Used for | Terms |
|---|---|---|
| Copernicus Sentinel-1, radiometrically terrain corrected by Catalyst, via Microsoft Planetary Computer | ship counts, detection positions, masks, radar images | CC BY 4.0. Contains modified Copernicus Sentinel data 2025-2026. |
| The GDELT Project (GKG 2.0 and Events 2.0 on Google BigQuery; DOC 2.0 API) | news and event indices | Free for any use and redistribution with a citation to the GDELT Project and a link to https://www.gdeltproject.org/ |
| Natural Earth, Admin 0 countries, 1:10m | land mask, flare regions | Public domain |
| NASA FIRMS, VIIRS 375 m active fires (NOAA-21 NRT; Suomi NPP SP and NRT) | gas flare series (`veri/firms/`) | NASA data are free and open; acknowledgement: "We acknowledge the use of data from NASA's FIRMS, part of NASA's Earth Science Data and Information System (ESDIS)." |
| SOMO (State Organization for Marketing of Oil, Iraq), monthly crude exports | `veri/rafineri/irak_somo_cikis_aylik.csv` | Public figures reproduced from somooil.gov.iq with attribution; not relicensed under CC BY |
| Public reporting compiled by MetalTakip | event list | CC BY 4.0 as above |

The code is licensed separately under the MIT licence (see `LICENSE`).
