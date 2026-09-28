-- GDELT 2.0 Events, düzeltilmiş kutular. Bâbülmendep kutusu Aden (45,0 D) ve Taiz'i (13,58 K) dışarıda bırakır.
-- Hürmüz: 25,0-27,4 K, 55,0-57,9 D. Bâbülmendep: 11,8-13,45 K, 42,6-44,2 D. Tehdit: CAMEO 13, 15. Eylem: 17-20.
WITH e AS (
  SELECT
    DATE(_PARTITIONTIME) AS gun, EventRootCode AS kok, NumArticles AS n, AvgTone AS ton, GoldsteinScale AS gs,
    (ActionGeo_Lat BETWEEN 25.0 AND 27.4 AND ActionGeo_Long BETWEEN 55.0 AND 57.9 AND ActionGeo_Type IN (2, 3, 4, 5)) AS h,
    (ActionGeo_Lat BETWEEN 11.8 AND 13.45 AND ActionGeo_Long BETWEEN 42.6 AND 44.2 AND ActionGeo_Type IN (2, 3, 4, 5)) AS b
  FROM `gdelt-bq.gdeltv2.events_partitioned`
  WHERE _PARTITIONTIME >= TIMESTAMP('2017-01-01') AND _PARTITIONTIME < TIMESTAMP('2026-09-24')
)
SELECT
  gun, COUNT(*) AS olay, SUM(n) AS makale,
  SUM(IF(h, n, 0)) AS h_makale, SUM(IF(h AND kok IN ('13', '15'), n, 0)) AS h_tehdit, SUM(IF(h AND kok IN ('17', '18', '19', '20'), n, 0)) AS h_eylem,
  AVG(IF(h, ton, NULL)) AS h_ton, AVG(IF(h, gs, NULL)) AS h_goldstein,
  SUM(IF(b, n, 0)) AS b_makale, SUM(IF(b AND kok IN ('13', '15'), n, 0)) AS b_tehdit, SUM(IF(b AND kok IN ('17', '18', '19', '20'), n, 0)) AS b_eylem,
  AVG(IF(b, ton, NULL)) AS b_ton, AVG(IF(b, gs, NULL)) AS b_goldstein
FROM e
GROUP BY gun
ORDER BY gun
