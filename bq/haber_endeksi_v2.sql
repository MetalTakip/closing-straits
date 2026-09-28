-- GDELT GKG 2.0: boğaz adı geçen belgeler (AllNames: metinde tanınan özel adlar; DocumentIdentifier: haber bağlantısı).
-- Türkçe kaynak: TranslationInfo 'srclc:tur' ya da alan adı .tr ile biter. Ton: V2Tone ilk alanı.
WITH g AS (
  SELECT
    DATE(_PARTITIONTIME) AS gun,
    (TranslationInfo LIKE '%srclc:tur%' OR REGEXP_CONTAINS(SourceCommonName, r'\.tr$')) AS tr,
    REGEXP_CONTAINS(AllNames, r'(?i)hormuz|hürmüz|hurmuz') AS h_ad,
    REGEXP_CONTAINS(AllNames, r'(?i)mand[ae]b|babülmendep|bâbülmendep|babulmendep') AS b_ad,
    REGEXP_CONTAINS(LOWER(DocumentIdentifier), r'hormuz|hurmuz|h%c3%bcrm%c3%bcz') AS h_url,
    REGEXP_CONTAINS(LOWER(DocumentIdentifier), r'mandeb|mandab|babulmendep|b%c3%a2b%c3%bclmendep') AS b_url,
    SAFE_CAST(SPLIT(V2Tone, ',')[SAFE_OFFSET(0)] AS FLOAT64) AS ton
  FROM `gdelt-bq.gdeltv2.gkg_partitioned`
  WHERE _PARTITIONTIME >= TIMESTAMP('2017-01-01') AND _PARTITIONTIME < TIMESTAMP('2026-09-24')
)
SELECT
  gun,
  COUNT(*) AS belge, COUNTIF(tr) AS tr_belge,
  COUNTIF(h_ad) AS h_ad, COUNTIF(b_ad) AS b_ad, COUNTIF(h_url) AS h_url, COUNTIF(b_url) AS b_url,
  COUNTIF(h_ad AND tr) AS h_ad_tr, COUNTIF(b_ad AND tr) AS b_ad_tr, COUNTIF(h_url AND tr) AS h_url_tr,
  AVG(IF(h_ad, ton, NULL)) AS h_ton, AVG(IF(b_ad, ton, NULL)) AS b_ton,
  AVG(IF(h_ad AND tr, ton, NULL)) AS h_ton_tr, AVG(IF(tr, ton, NULL)) AS tr_ton_genel
FROM g
GROUP BY gun
ORDER BY gun
