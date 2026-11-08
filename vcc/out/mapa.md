# Mapa completo — local y S3, archivo por archivo

Solo lectura. Cero cambios: no se movió, copió, borró ni renombró nada; no se tocó
`systemd`, `S3`, ni `X` de ningún archivo en esta tarea.

## 1 · Los tres archivos de mapa, completos

| archivo | comando | líneas |
|---|---|---|
| `~/vcc/out/mapa_vcc.txt` | `find /home/ubuntu/vcc -printf '%y\t%p\t%s\t%TY-%Tm-%Td %TH:%TM\t%l\n' \| sort -k2` | 525 |
| `~/vcc/out/mapa_datos.txt` | mismo comando sobre `/datos` | 77 |
| `~/vcc/out/mapa_s3.txt` | `aws s3 ls s3://vc-challenge-data/ --recursive --human-readable` | 968 |

`mapa_s3.txt` tiene **968 líneas = 968 objetos** (confirmado también vía
`aws s3api list-objects-v2`, que da el mismo total exacto).

## 2 · Árbol ENTERO de `~/vcc` (525 líneas — tipo, ruta, tamaño en bytes, fecha)

```
d	/home/ubuntu/vcc	4096	2026-10-28 05:57	
d	/home/ubuntu/vcc/.claude	4096	2026-10-28 07:26	
f	/home/ubuntu/vcc/PLAN.md	9637	2026-10-27 06:07	
d	/home/ubuntu/vcc/corpus	4096	2026-10-28 21:12	
f	/home/ubuntu/vcc/corpus/inventario.parquet	819766	2026-10-28 21:07	
f	/home/ubuntu/vcc/corpus/particion.md	6503	2026-10-28 21:12	
f	/home/ubuntu/vcc/mf.sh	128128881	2026-10-27 01:58	
d	/home/ubuntu/vcc/out	4096	2026-10-28 21:43	
f	/home/ubuntu/vcc/out/b_panel_tf.md	8504	2026-10-28 05:55	
f	/home/ubuntu/vcc/out/b_panel_tf_raw.json	59667	2026-10-28 05:54	
f	/home/ubuntu/vcc/out/b_techo_i6_ruido.md	11639	2026-10-28 07:19	
f	/home/ubuntu/vcc/out/degenerado_check3_table.csv	2916	2026-10-27 20:13	
f	/home/ubuntu/vcc/out/estructura.md	13322	2026-10-28 21:27	
f	/home/ubuntu/vcc/out/geo_zenodo.md	16685	2026-10-27 08:50	
f	/home/ubuntu/vcc/out/gwps_fiabilidad.md	12417	2026-10-28 20:43	
f	/home/ubuntu/vcc/out/gwps_fiabilidad_raw.json	230761	2026-10-28 20:41	
f	/home/ubuntu/vcc/out/gwps_resolucion.md	10515	2026-10-28 20:28	
f	/home/ubuntu/vcc/out/gwps_resolucion_raw.json	33562	2026-10-28 20:26	
f	/home/ubuntu/vcc/out/mapa_vcc.txt	0	2026-10-28 21:43	
f	/home/ubuntu/vcc/out/nreal_marcador.md	8986	2026-10-28 19:53	
f	/home/ubuntu/vcc/out/nreal_marcador_basal_raw.json	15550	2026-10-28 19:52	
f	/home/ubuntu/vcc/out/nreal_predictores.md	13236	2026-10-28 19:44	
f	/home/ubuntu/vcc/out/nreal_predictores_raw.json	27972	2026-10-28 19:42	
f	/home/ubuntu/vcc/out/nreal_transfiere.md	8136	2026-10-28 19:58	
f	/home/ubuntu/vcc/out/nreal_transfiere_raw.json	24575	2026-10-28 19:57	
f	/home/ubuntu/vcc/out/panel_tf.md	3202	2026-10-28 05:02	
f	/home/ubuntu/vcc/out/panel_tf_raw.json	50035	2026-10-28 05:01	
f	/home/ubuntu/vcc/out/paso0_inventario.md	5548	2026-10-27 06:01	
f	/home/ubuntu/vcc/out/paso0_obs.md	13966	2026-10-27 06:15	
f	/home/ubuntu/vcc/out/rejilla_table.csv	14643	2026-10-27 19:58	
f	/home/ubuntu/vcc/out/replogle_addendum_k562_essential.json	3222	2026-10-27 20:37	
f	/home/ubuntu/vcc/out/replogle_addendum_k562_gwps.json	3233	2026-10-27 20:47	
f	/home/ubuntu/vcc/out/replogle_addendum_rpe1.json	3203	2026-10-27 20:39	
f	/home/ubuntu/vcc/out/replogle_bulk_parte1.json	28080	2026-10-27 20:26	
f	/home/ubuntu/vcc/out/replogle_bulk_parte3.json	9245	2026-10-27 20:27	
f	/home/ubuntu/vcc/out/replogle_contaminacion.md	5064	2026-10-28 04:38	
f	/home/ubuntu/vcc/out/replogle_contaminacion_raw.json	5556	2026-10-28 04:37	
f	/home/ubuntu/vcc/out/replogle_contaminacion_v2.md	5235	2026-10-28 04:48	
f	/home/ubuntu/vcc/out/replogle_contaminacion_v2_raw.json	8151	2026-10-28 04:48	
f	/home/ubuntu/vcc/out/replogle_controlexpr.md	7399	2026-10-28 05:46	
f	/home/ubuntu/vcc/out/replogle_controlexpr_k562_essential.json	1942	2026-10-28 05:28	
f	/home/ubuntu/vcc/out/replogle_controlexpr_k562_gwps.json	1933	2026-10-28 05:37	
f	/home/ubuntu/vcc/out/replogle_controlexpr_porgrupo.md	11280	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/replogle_controlexpr_porgrupo_k562_essential.json	15170	2026-10-28 07:07	
f	/home/ubuntu/vcc/out/replogle_controlexpr_porgrupo_k562_gwps.json	15293	2026-10-28 07:16	
f	/home/ubuntu/vcc/out/replogle_controlexpr_porgrupo_rpe1.json	15111	2026-10-28 07:06	
f	/home/ubuntu/vcc/out/replogle_controlexpr_rpe1.json	1930	2026-10-28 05:29	
f	/home/ubuntu/vcc/out/replogle_controlexpr_sample20_k562_essential.json	4736	2026-10-28 05:38	
f	/home/ubuntu/vcc/out/replogle_controlexpr_sample20_k562_gwps.json	4772	2026-10-28 05:45	
f	/home/ubuntu/vcc/out/replogle_controlexpr_sample20_rpe1.json	4703	2026-10-28 05:38	
f	/home/ubuntu/vcc/out/replogle_controlexpr_unfiltered.md	9740	2026-10-28 06:06	
f	/home/ubuntu/vcc/out/replogle_controlexpr_unfiltered_k562_essential.json	12968	2026-10-28 05:53	
f	/home/ubuntu/vcc/out/replogle_controlexpr_unfiltered_k562_gwps.json	13060	2026-10-28 06:05	
f	/home/ubuntu/vcc/out/replogle_controlexpr_unfiltered_rpe1.json	12883	2026-10-28 05:54	
f	/home/ubuntu/vcc/out/replogle_criba.md	28650	2026-10-27 20:49	
f	/home/ubuntu/vcc/out/replogle_dosetapas_detail_k562_essential.csv	217035	2026-10-28 05:07	
f	/home/ubuntu/vcc/out/replogle_dosetapas_detail_k562_gwps.csv	913246	2026-10-28 05:16	
f	/home/ubuntu/vcc/out/replogle_dosetapas_detail_rpe1.csv	237244	2026-10-28 05:09	
f	/home/ubuntu/vcc/out/replogle_dosetapas_k562_essential.json	2177	2026-10-28 05:07	
f	/home/ubuntu/vcc/out/replogle_dosetapas_k562_gwps.json	2154	2026-10-28 05:16	
f	/home/ubuntu/vcc/out/replogle_dosetapas_rpe1.json	2127	2026-10-28 05:09	
f	/home/ubuntu/vcc/out/replogle_escalon2.md	8598	2026-10-27 20:58	
f	/home/ubuntu/vcc/out/replogle_escalon2_partea.json	5332	2026-10-27 20:50	
f	/home/ubuntu/vcc/out/replogle_escalon2_parteb.json	1646	2026-10-27 20:56	
f	/home/ubuntu/vcc/out/replogle_foldexpr_descomposicion.md	4204	2026-10-28 04:52	
f	/home/ubuntu/vcc/out/replogle_foldexpr_descomposicion_raw.json	2536	2026-10-28 04:51	
f	/home/ubuntu/vcc/out/replogle_foldexpr_detail_k562_essential.csv	145628	2026-10-27 20:37	
f	/home/ubuntu/vcc/out/replogle_foldexpr_detail_k562_gwps.csv	610843	2026-10-27 20:47	
f	/home/ubuntu/vcc/out/replogle_foldexpr_detail_rpe1.csv	159364	2026-10-27 20:39	
f	/home/ubuntu/vcc/out/replogle_foldexpr_dosetapas.md	5414	2026-10-28 05:17	
f	/home/ubuntu/vcc/out/replogle_foldexpr_unfiltered.md	5833	2026-10-28 05:24	
f	/home/ubuntu/vcc/out/replogle_foldexpr_unfiltered_raw.json	7024	2026-10-28 05:23	
f	/home/ubuntu/vcc/out/replogle_tests_estimando.md	7134	2026-10-28 18:56	
f	/home/ubuntu/vcc/out/replogle_tests_estimando_mw_sample_k562_essential.json	49390	2026-10-28 18:54	
f	/home/ubuntu/vcc/out/replogle_tests_estimando_mw_sample_run.log	1996	2026-10-28 18:54	
f	/home/ubuntu/vcc/out/replogle_tests_estimando_sample_run.log	434	2026-10-28 17:54	
f	/home/ubuntu/vcc/out/respaldo.md	12985	2026-10-28 21:38	
f	/home/ubuntu/vcc/out/restore.DONE	34	2026-10-27 18:56	
f	/home/ubuntu/vcc/out/restore.log	6930	2026-10-27 18:56	
f	/home/ubuntu/vcc/out/s3_convertido.DONE	34	2026-10-27 08:55	
f	/home/ubuntu/vcc/out/s3_convertido.log	783	2026-10-27 08:55	
f	/home/ubuntu/vcc/out/s3_datos.DONE	34	2026-10-27 07:25	
f	/home/ubuntu/vcc/out/s3_datos.log	5325	2026-10-27 07:25	
f	/home/ubuntu/vcc/out/t13_comprobaciones.md	5134	2026-10-27 06:23	
f	/home/ubuntu/vcc/out/t13_pools.md	16374	2026-10-27 06:31	
f	/home/ubuntu/vcc/out/t13_replicas.md	7502	2026-10-27 06:46	
f	/home/ubuntu/vcc/out/t13_replicas.parquet	491113	2026-10-27 06:44	
f	/home/ubuntu/vcc/out/t3_inventario.md	17644	2026-10-28 21:14	
f	/home/ubuntu/vcc/out/t3_inventario_run.log	1689	2026-10-28 21:07	
f	/home/ubuntu/vcc/out/t3_q6_q7_raw.json	2493	2026-10-28 21:10	
f	/home/ubuntu/vcc/out/t3_queries_raw.json	23455	2026-10-28 21:10	
f	/home/ubuntu/vcc/out/techo_estado.md	7100	2026-10-28 04:55	
f	/home/ubuntu/vcc/out/techo_metricas.md	16617	2026-10-28 04:56	
d	/home/ubuntu/vcc/out/techo_metricas_work	4096	2026-10-28 05:56	
d	/home/ubuntu/vcc/out/techo_metricas_work/baseline	4096	2026-10-27 23:09	
f	/home/ubuntu/vcc/out/techo_metricas_work/baseline/baseline_agg.csv	1229	2026-10-27 23:09	
f	/home/ubuntu/vcc/out/techo_metricas_work/baseline/baseline_meta.json	12532	2026-10-27 23:09	
f	/home/ubuntu/vcc/out/techo_metricas_work/baseline/baseline_results.csv	41762	2026-10-27 23:09	
f	/home/ubuntu/vcc/out/techo_metricas_work/baseline/run_params.yaml	991	2026-10-27 23:09	
f	/home/ubuntu/vcc/out/techo_metricas_work/baseline_build.log	27013	2026-10-27 23:09	
d	/home/ubuntu/vcc/out/techo_metricas_work/bundle	4096	2026-10-28 00:55	
f	/home/ubuntu/vcc/out/techo_metricas_work/bundle/anchor_agg.parquet	3772	2026-10-28 00:55	
f	/home/ubuntu/vcc/out/techo_metricas_work/bundle/anchor_meta.json	2283	2026-10-28 00:55	
f	/home/ubuntu/vcc/out/techo_metricas_work/bundle/anchor_splits.parquet	2747	2026-10-28 00:55	
f	/home/ubuntu/vcc/out/techo_metricas_work/bundle/baseline_agg.csv	1229	2026-10-28 00:55	
f	/home/ubuntu/vcc/out/techo_metricas_work/bundle/baseline_meta.json	2184	2026-10-28 00:55	
f	/home/ubuntu/vcc/out/techo_metricas_work/bundle/config.yaml	990	2026-10-28 00:55	
f	/home/ubuntu/vcc/out/techo_metricas_work/bundle/manifest.json	2437	2026-10-28 00:55	
f	/home/ubuntu/vcc/out/techo_metricas_work/bundle_build.log	293158	2026-10-28 00:55	
f	/home/ubuntu/vcc/out/techo_metricas_work/driver_i1_i5.log	2244	2026-10-28 04:14	
d	/home/ubuntu/vcc/out/techo_metricas_work/i5	4096	2026-10-28 04:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i5/de_pred_scaled_0p0.parquet	14076994	2026-10-28 04:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i5/de_pred_scaled_0p25.parquet	30059967	2026-10-28 04:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/i5/de_pred_scaled_0p5.parquet	30057772	2026-10-28 04:05	
f	/home/ubuntu/vcc/out/techo_metricas_work/i5/de_pred_scaled_0p75.parquet	30053455	2026-10-28 04:02	
f	/home/ubuntu/vcc/out/techo_metricas_work/i5/de_pred_scaled_1p0.parquet	30056989	2026-10-28 03:59	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6	4096	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/b_run_summary.txt	888	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/b_techo_i6_gen_meta.json	22952	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p25_f0p25.parquet	27625041	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p25_f0p4.parquet	27613327	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p25_f0p5.parquet	27622452	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p25_f0p6.parquet	27624328	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p25_f0p75.parquet	27620155	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p25_f0p9.parquet	27623143	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p25_f1p0.parquet	27619665	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p5_f0p25.parquet	26406428	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p5_f0p4.parquet	26399428	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p5_f0p5.parquet	26402406	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p5_f0p6.parquet	26407209	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p5_f0p75.parquet	26402630	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p5_f0p9.parquet	26407324	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s0p5_f1p0.parquet	26402808	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s1p0_f0p25.parquet	25244566	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s1p0_f0p4.parquet	25237920	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s1p0_f0p5.parquet	25242244	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s1p0_f0p6.parquet	25244907	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s1p0_f0p75.parquet	25243504	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s1p0_f0p9.parquet	25247118	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/de_pred_s1p0_f1p0.parquet	25244735	2026-10-28 05:56	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p25_f0p25.log	592	2026-10-28 06:33	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p25_f0p4.log	589	2026-10-28 06:30	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p25_f0p5.log	589	2026-10-28 06:27	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p25_f0p6.log	589	2026-10-28 06:24	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p25_f0p75.log	592	2026-10-28 06:21	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p25_f0p9.log	594	2026-10-28 06:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p25_f1p0.log	589	2026-10-28 06:15	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p5_f0p25.log	591	2026-10-28 06:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p5_f0p4.log	590	2026-10-28 06:51	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p5_f0p5.log	593	2026-10-28 06:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p5_f0p6.log	590	2026-10-28 06:45	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p5_f0p75.log	591	2026-10-28 06:42	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p5_f0p9.log	590	2026-10-28 06:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s0p5_f1p0.log	588	2026-10-28 06:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s1p0_f0p25.log	594	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s1p0_f0p4.log	588	2026-10-28 07:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s1p0_f0p5.log	593	2026-10-28 07:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s1p0_f0p6.log	593	2026-10-28 07:07	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s1p0_f0p75.log	591	2026-10-28 07:04	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s1p0_f0p9.log	588	2026-10-28 07:01	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/run_s1p0_f1p0.log	588	2026-10-28 06:57	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p25_f0p25.log	915	2026-10-28 06:33	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p25_f0p4.log	914	2026-10-28 06:30	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p25_f0p5.log	914	2026-10-28 06:27	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p25_f0p6.log	914	2026-10-28 06:24	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p25_f0p75.log	915	2026-10-28 06:21	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p25_f0p9.log	914	2026-10-28 06:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p25_f1p0.log	914	2026-10-28 06:15	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p5_f0p25.log	914	2026-10-28 06:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p5_f0p4.log	913	2026-10-28 06:51	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p5_f0p5.log	913	2026-10-28 06:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p5_f0p6.log	913	2026-10-28 06:45	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p5_f0p75.log	914	2026-10-28 06:42	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p5_f0p9.log	913	2026-10-28 06:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s0p5_f1p0.log	913	2026-10-28 06:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s1p0_f0p25.log	914	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s1p0_f0p4.log	913	2026-10-28 07:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s1p0_f0p5.log	913	2026-10-28 07:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s1p0_f0p6.log	913	2026-10-28 07:07	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s1p0_f0p75.log	914	2026-10-28 07:04	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s1p0_f0p9.log	913	2026-10-28 07:01	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/score_s1p0_f1p0.log	913	2026-10-28 06:58	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p25_f0p25_baseonly.csv	323	2026-10-28 06:33	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p25_f0p4_baseonly.csv	322	2026-10-28 06:30	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p25_f0p5_baseonly.csv	321	2026-10-28 06:27	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p25_f0p6_baseonly.csv	323	2026-10-28 06:24	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p25_f0p75_baseonly.csv	321	2026-10-28 06:21	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p25_f0p9_baseonly.csv	322	2026-10-28 06:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p25_f1p0_baseonly.csv	322	2026-10-28 06:15	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p5_f0p25_baseonly.csv	323	2026-10-28 06:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p5_f0p4_baseonly.csv	322	2026-10-28 06:51	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p5_f0p5_baseonly.csv	322	2026-10-28 06:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p5_f0p6_baseonly.csv	322	2026-10-28 06:45	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p5_f0p75_baseonly.csv	322	2026-10-28 06:42	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p5_f0p9_baseonly.csv	322	2026-10-28 06:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s0p5_f1p0_baseonly.csv	320	2026-10-28 06:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s1p0_f0p25_baseonly.csv	323	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s1p0_f0p4_baseonly.csv	322	2026-10-28 07:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s1p0_f0p5_baseonly.csv	323	2026-10-28 07:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s1p0_f0p6_baseonly.csv	322	2026-10-28 07:07	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s1p0_f0p75_baseonly.csv	322	2026-10-28 07:04	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s1p0_f0p9_baseonly.csv	323	2026-10-28 07:01	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/scored_i6_s1p0_f1p0_baseonly.csv	322	2026-10-28 06:58	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p25	4096	2026-10-28 06:33	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p25/agg_results.csv	1209	2026-10-28 06:33	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p25/evaluation_results.json	35226	2026-10-28 06:33	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p25/metric_aggregation.csv	523	2026-10-28 06:33	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p25/results.csv	39683	2026-10-28 06:33	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p25/run_meta.json	2264	2026-10-28 06:33	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p25/run_params.yaml	1006	2026-10-28 06:33	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p4	4096	2026-10-28 06:30	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p4/agg_results.csv	1209	2026-10-28 06:30	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p4/evaluation_results.json	35225	2026-10-28 06:30	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p4/metric_aggregation.csv	523	2026-10-28 06:30	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p4/results.csv	39682	2026-10-28 06:30	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p4/run_meta.json	2264	2026-10-28 06:30	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p4/run_params.yaml	1005	2026-10-28 06:30	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p5	4096	2026-10-28 06:27	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p5/agg_results.csv	1210	2026-10-28 06:27	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p5/evaluation_results.json	35223	2026-10-28 06:27	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p5/metric_aggregation.csv	523	2026-10-28 06:27	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p5/results.csv	39679	2026-10-28 06:27	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p5/run_meta.json	2264	2026-10-28 06:27	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p5/run_params.yaml	1005	2026-10-28 06:27	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p6	4096	2026-10-28 06:24	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p6/agg_results.csv	1210	2026-10-28 06:24	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p6/evaluation_results.json	35255	2026-10-28 06:24	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p6/metric_aggregation.csv	523	2026-10-28 06:24	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p6/results.csv	39711	2026-10-28 06:24	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p6/run_meta.json	2264	2026-10-28 06:24	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p6/run_params.yaml	1005	2026-10-28 06:24	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p75	4096	2026-10-28 06:21	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p75/agg_results.csv	1211	2026-10-28 06:21	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p75/evaluation_results.json	35274	2026-10-28 06:21	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p75/metric_aggregation.csv	523	2026-10-28 06:21	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p75/results.csv	39729	2026-10-28 06:21	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p75/run_meta.json	2264	2026-10-28 06:21	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p75/run_params.yaml	1006	2026-10-28 06:21	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p9	4096	2026-10-28 06:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p9/agg_results.csv	1212	2026-10-28 06:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p9/evaluation_results.json	35280	2026-10-28 06:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p9/metric_aggregation.csv	523	2026-10-28 06:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p9/results.csv	39734	2026-10-28 06:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p9/run_meta.json	2264	2026-10-28 06:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f0p9/run_params.yaml	1005	2026-10-28 06:18	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f1p0	4096	2026-10-28 06:15	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f1p0/agg_results.csv	1212	2026-10-28 06:15	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f1p0/evaluation_results.json	35283	2026-10-28 06:15	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f1p0/metric_aggregation.csv	523	2026-10-28 06:15	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f1p0/results.csv	39737	2026-10-28 06:15	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f1p0/run_meta.json	2264	2026-10-28 06:15	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p25_f1p0/run_params.yaml	1005	2026-10-28 06:15	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p25	4096	2026-10-28 06:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p25/agg_results.csv	1210	2026-10-28 06:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p25/evaluation_results.json	35228	2026-10-28 06:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p25/metric_aggregation.csv	523	2026-10-28 06:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p25/results.csv	39684	2026-10-28 06:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p25/run_meta.json	2264	2026-10-28 06:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p25/run_params.yaml	1005	2026-10-28 06:54	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p4	4096	2026-10-28 06:51	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p4/agg_results.csv	1211	2026-10-28 06:51	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p4/evaluation_results.json	35230	2026-10-28 06:51	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p4/metric_aggregation.csv	523	2026-10-28 06:51	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p4/results.csv	39685	2026-10-28 06:51	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p4/run_meta.json	2264	2026-10-28 06:51	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p4/run_params.yaml	1004	2026-10-28 06:51	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p5	4096	2026-10-28 06:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p5/agg_results.csv	1208	2026-10-28 06:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p5/evaluation_results.json	35223	2026-10-28 06:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p5/metric_aggregation.csv	523	2026-10-28 06:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p5/results.csv	39681	2026-10-28 06:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p5/run_meta.json	2264	2026-10-28 06:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p5/run_params.yaml	1004	2026-10-28 06:48	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p6	4096	2026-10-28 06:45	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p6/agg_results.csv	1210	2026-10-28 06:45	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p6/evaluation_results.json	35250	2026-10-28 06:45	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p6/metric_aggregation.csv	523	2026-10-28 06:45	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p6/results.csv	39706	2026-10-28 06:45	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p6/run_meta.json	2264	2026-10-28 06:45	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p6/run_params.yaml	1004	2026-10-28 06:45	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p75	4096	2026-10-28 06:42	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p75/agg_results.csv	1212	2026-10-28 06:42	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p75/evaluation_results.json	35263	2026-10-28 06:42	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p75/metric_aggregation.csv	523	2026-10-28 06:42	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p75/results.csv	39717	2026-10-28 06:42	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p75/run_meta.json	2264	2026-10-28 06:42	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p75/run_params.yaml	1005	2026-10-28 06:42	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p9	4096	2026-10-28 06:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p9/agg_results.csv	1211	2026-10-28 06:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p9/evaluation_results.json	35266	2026-10-28 06:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p9/metric_aggregation.csv	523	2026-10-28 06:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p9/results.csv	39721	2026-10-28 06:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p9/run_meta.json	2264	2026-10-28 06:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f0p9/run_params.yaml	1004	2026-10-28 06:39	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f1p0	4096	2026-10-28 06:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f1p0/agg_results.csv	1211	2026-10-28 06:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f1p0/evaluation_results.json	35263	2026-10-28 06:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f1p0/metric_aggregation.csv	523	2026-10-28 06:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f1p0/results.csv	39718	2026-10-28 06:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f1p0/run_meta.json	2264	2026-10-28 06:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s0p5_f1p0/run_params.yaml	1004	2026-10-28 06:36	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p25	4096	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p25/agg_results.csv	1211	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p25/evaluation_results.json	35229	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p25/metric_aggregation.csv	523	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p25/results.csv	39684	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p25/run_meta.json	2264	2026-10-28 07:18	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p25/run_params.yaml	1005	2026-10-28 07:18	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p4	4096	2026-10-28 07:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p4/agg_results.csv	1208	2026-10-28 07:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p4/evaluation_results.json	35220	2026-10-28 07:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p4/metric_aggregation.csv	523	2026-10-28 07:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p4/results.csv	39678	2026-10-28 07:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p4/run_meta.json	2264	2026-10-28 07:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p4/run_params.yaml	1004	2026-10-28 07:14	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p5	4096	2026-10-28 07:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p5/agg_results.csv	1211	2026-10-28 07:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p5/evaluation_results.json	35228	2026-10-28 07:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p5/metric_aggregation.csv	523	2026-10-28 07:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p5/results.csv	39683	2026-10-28 07:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p5/run_meta.json	2264	2026-10-28 07:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p5/run_params.yaml	1004	2026-10-28 07:11	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p6	4096	2026-10-28 07:07	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p6/agg_results.csv	1211	2026-10-28 07:07	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p6/evaluation_results.json	35242	2026-10-28 07:07	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p6/metric_aggregation.csv	523	2026-10-28 07:07	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p6/results.csv	39697	2026-10-28 07:07	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p6/run_meta.json	2264	2026-10-28 07:07	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p6/run_params.yaml	1004	2026-10-28 07:07	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p75	4096	2026-10-28 07:04	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p75/agg_results.csv	1211	2026-10-28 07:04	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p75/evaluation_results.json	35244	2026-10-28 07:04	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p75/metric_aggregation.csv	523	2026-10-28 07:04	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p75/results.csv	39699	2026-10-28 07:04	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p75/run_meta.json	2264	2026-10-28 07:04	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p75/run_params.yaml	1005	2026-10-28 07:04	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p9	4096	2026-10-28 07:01	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p9/agg_results.csv	1210	2026-10-28 07:01	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p9/evaluation_results.json	35251	2026-10-28 07:01	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p9/metric_aggregation.csv	523	2026-10-28 07:01	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p9/results.csv	39707	2026-10-28 07:01	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p9/run_meta.json	2264	2026-10-28 07:01	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f0p9/run_params.yaml	1004	2026-10-28 07:01	
d	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f1p0	4096	2026-10-28 06:57	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f1p0/agg_results.csv	1208	2026-10-28 06:57	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f1p0/evaluation_results.json	35242	2026-10-28 06:57	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f1p0/metric_aggregation.csv	523	2026-10-28 06:57	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f1p0/results.csv	39700	2026-10-28 06:57	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f1p0/run_meta.json	2264	2026-10-28 06:57	
f	/home/ubuntu/vcc/out/techo_metricas_work/i6/user_i6_s1p0_f1p0/run_params.yaml	1004	2026-10-28 06:57	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_base.log	26046	2026-10-28 01:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_i1_cpm1e6.log	26052	2026-10-28 03:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_i2_invert_agg.log	26057	2026-10-28 03:28	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_i3_log1p.log	26052	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_i4_log.log	16923	2026-10-28 03:59	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_i4_log_retry.log	26052	2026-10-28 04:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_i5_0p0.log	579	2026-10-28 04:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_i5_0p25.log	579	2026-10-28 04:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_i5_0p5.log	584	2026-10-28 04:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_i5_0p75.log	580	2026-10-28 04:05	
f	/home/ubuntu/vcc/out/techo_metricas_work/run_i5_1p0.log	581	2026-10-28 04:02	
f	/home/ubuntu/vcc/out/techo_metricas_work/score_i1_cpm1e6.log	1209	2026-10-28 03:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/score_i2_invert_agg.log	1208	2026-10-28 03:28	
f	/home/ubuntu/vcc/out/techo_metricas_work/score_i3_log1p.log	1209	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/score_i4_log.log	222	2026-10-28 03:59	
f	/home/ubuntu/vcc/out/techo_metricas_work/score_i5_0p0.log	1373	2026-10-28 04:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/score_i5_0p25.log	1373	2026-10-28 04:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/score_i5_0p5.log	1373	2026-10-28 04:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/score_i5_0p75.log	1373	2026-10-28 04:05	
f	/home/ubuntu/vcc/out/techo_metricas_work/score_i5_1p0.log	1373	2026-10-28 04:02	
f	/home/ubuntu/vcc/out/techo_metricas_work/scored_base.csv	1540	2026-10-28 02:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/scored_i1_cpm1e6_baseonly.csv	324	2026-10-28 03:09	
f	/home/ubuntu/vcc/out/techo_metricas_work/scored_i2_invert_agg_baseonly.csv	323	2026-10-28 03:28	
f	/home/ubuntu/vcc/out/techo_metricas_work/scored_i3_log1p_baseonly.csv	316	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/scored_i4_log_baseonly.csv	313	2026-10-28 04:55	
f	/home/ubuntu/vcc/out/techo_metricas_work/scored_i5_0p0_baseonly.csv	328	2026-10-28 04:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/scored_i5_0p25_baseonly.csv	322	2026-10-28 04:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/scored_i5_0p5_baseonly.csv	322	2026-10-28 04:36	
f	/home/ubuntu/vcc/out/techo_metricas_work/scored_i5_0p75_baseonly.csv	322	2026-10-28 04:35	
f	/home/ubuntu/vcc/out/techo_metricas_work/scored_i5_1p0_baseonly.csv	320	2026-10-28 04:35	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_base	4096	2026-10-28 01:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_base/agg_results.csv	1214	2026-10-28 01:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_base/de_pred.parquet	30056989	2026-10-28 01:38	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_base/de_real.parquet	30103831	2026-10-28 01:38	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_base/evaluation_results.json	35410	2026-10-28 01:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_base/metric_aggregation.csv	523	2026-10-28 01:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_base/results.csv	39739	2026-10-28 01:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_base/run_meta.json	2219	2026-10-28 01:39	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_base/run_params.yaml	993	2026-10-28 01:39	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i1_cpm1e6	4096	2026-10-28 03:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i1_cpm1e6/agg_results.csv	1214	2026-10-28 03:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i1_cpm1e6/evaluation_results.json	35324	2026-10-28 03:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i1_cpm1e6/metric_aggregation.csv	523	2026-10-28 03:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i1_cpm1e6/results.csv	39821	2026-10-28 03:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i1_cpm1e6/run_meta.json	2219	2026-10-28 03:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i1_cpm1e6/run_params.yaml	998	2026-10-28 03:08	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i1_cpm1e6_default_attempt	4096	2026-10-28 02:49	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i1_cpm1e6_default_attempt.log	2377	2026-10-28 02:49	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i2_invert_agg	4096	2026-10-28 03:28	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i2_invert_agg/agg_results.csv	1213	2026-10-28 03:28	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i2_invert_agg/evaluation_results.json	35318	2026-10-28 03:28	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i2_invert_agg/metric_aggregation.csv	523	2026-10-28 03:28	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i2_invert_agg/results.csv	39816	2026-10-28 03:28	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i2_invert_agg/run_meta.json	2219	2026-10-28 03:28	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i2_invert_agg/run_params.yaml	1002	2026-10-28 03:28	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i2_invert_agg_default_attempt	4096	2026-10-28 03:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i2_invert_agg_default_attempt.log	2376	2026-10-28 03:09	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i3_log1p	4096	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i3_log1p/agg_results.csv	1214	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i3_log1p/evaluation_results.json	37359	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i3_log1p/metric_aggregation.csv	523	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i3_log1p/results.csv	41704	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i3_log1p/run_meta.json	2219	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i3_log1p/run_params.yaml	997	2026-10-28 03:48	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i3_log1p_default_attempt	4096	2026-10-28 03:28	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i3_log1p_default_attempt.log	2377	2026-10-28 03:29	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i4_log	4096	2026-10-28 04:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i4_log/agg_results.csv	1234	2026-10-28 04:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i4_log/evaluation_results.json	37153	2026-10-28 04:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i4_log/metric_aggregation.csv	523	2026-10-28 04:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i4_log/results.csv	41478	2026-10-28 04:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i4_log/run_meta.json	2219	2026-10-28 04:54	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i4_log/run_params.yaml	995	2026-10-28 04:54	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i4_log_default_attempt	4096	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i4_log_default_attempt.log	2377	2026-10-28 03:48	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p0	4096	2026-10-28 04:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p0/agg_results.csv	1025	2026-10-28 04:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p0/evaluation_results.json	32478	2026-10-28 04:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p0/metric_aggregation.csv	523	2026-10-28 04:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p0/results.csv	37119	2026-10-28 04:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p0/run_meta.json	2264	2026-10-28 04:14	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p0/run_params.yaml	995	2026-10-28 04:14	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p25	4096	2026-10-28 04:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p25/agg_results.csv	1211	2026-10-28 04:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p25/evaluation_results.json	35230	2026-10-28 04:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p25/metric_aggregation.csv	523	2026-10-28 04:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p25/results.csv	39685	2026-10-28 04:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p25/run_meta.json	2264	2026-10-28 04:11	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p25/run_params.yaml	996	2026-10-28 04:11	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p5	4096	2026-10-28 04:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p5/agg_results.csv	1209	2026-10-28 04:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p5/evaluation_results.json	35234	2026-10-28 04:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p5/metric_aggregation.csv	523	2026-10-28 04:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p5/results.csv	39691	2026-10-28 04:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p5/run_meta.json	2264	2026-10-28 04:08	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p5/run_params.yaml	995	2026-10-28 04:08	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p75	4096	2026-10-28 04:05	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p75/agg_results.csv	1212	2026-10-28 04:05	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p75/evaluation_results.json	35266	2026-10-28 04:05	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p75/metric_aggregation.csv	523	2026-10-28 04:05	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p75/results.csv	39720	2026-10-28 04:05	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p75/run_meta.json	2264	2026-10-28 04:05	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_0p75/run_params.yaml	996	2026-10-28 04:05	
d	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_1p0	4096	2026-10-28 04:02	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_1p0/agg_results.csv	1213	2026-10-28 04:02	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_1p0/evaluation_results.json	35286	2026-10-28 04:02	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_1p0/metric_aggregation.csv	523	2026-10-28 04:02	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_1p0/results.csv	39739	2026-10-28 04:02	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_1p0/run_meta.json	2264	2026-10-28 04:02	
f	/home/ubuntu/vcc/out/techo_metricas_work/user_i5_1p0/run_params.yaml	995	2026-10-28 04:02	
f	/home/ubuntu/vcc/out/techo_metricas_work/validation_summary.txt	780	2026-10-28 03:48	
f	/home/ubuntu/vcc/out/vcc2025.json	6134	2026-10-27 15:50	
f	/home/ubuntu/vcc/out/vcc2025_analisis.md	21494	2026-10-27 19:30	
f	/home/ubuntu/vcc/out/vcc2025_analysis_raw.json	13815	2026-10-27 16:12	
f	/home/ubuntu/vcc/out/vcc2025_degenerado.md	11624	2026-10-27 20:15	
f	/home/ubuntu/vcc/out/vcc2025_degenerado_raw.json	26655	2026-10-27 20:13	
f	/home/ubuntu/vcc/out/vcc2025_knockdown_grid.json	4057	2026-10-27 19:26	
f	/home/ubuntu/vcc/out/vcc2025_knockdown_raw_arrays.npz	1884816	2026-10-27 19:26	
f	/home/ubuntu/vcc/out/vcc2025_percell_umi.json	696	2026-10-27 19:28	
f	/home/ubuntu/vcc/out/vcc2025_rejilla.md	14460	2026-10-27 20:06	
f	/home/ubuntu/vcc/out/vcc2025_rejilla_grid.json	71690	2026-10-27 19:58	
f	/home/ubuntu/vcc/out/vcc2025_rejilla_raw_arrays_test_train.npz	12113673	2026-10-27 19:57	
f	/home/ubuntu/vcc/out/vcc2025_rejilla_raw_arrays_validation.npz	2503056	2026-10-27 19:56	
f	/home/ubuntu/vcc/out/vcc2025_rejilla_scan_summary_test_train.json	201	2026-10-27 19:57	
f	/home/ubuntu/vcc/out/vcc2025_rejilla_scan_summary_validation.json	105	2026-10-27 19:56	
f	/home/ubuntu/vcc/out/vcc2025_structure_full.json	13569	2026-10-27 19:21	
d	/home/ubuntu/vcc/registro	4096	2026-10-28 21:34	
f	/home/ubuntu/vcc/registro/CONSTANTES_A_VERIFICAR.md	7699	2026-10-27 05:53	
f	/home/ubuntu/vcc/registro/limitaciones.md	9313	2026-10-27 07:43	
f	/home/ubuntu/vcc/registro/syncout.service.bak	299	2026-10-28 21:34	
d	/home/ubuntu/vcc/src	4096	2026-10-28 21:03	
d	/home/ubuntu/vcc/src/__pycache__	4096	2026-10-27 02:00	
f	/home/ubuntu/vcc/src/__pycache__/cpm_stats.cpython-311.pyc	5775	2026-10-27 02:00	
f	/home/ubuntu/vcc/src/__pycache__/pseudobulk.cpython-311.pyc	2555	2026-10-27 02:00	
f	/home/ubuntu/vcc/src/analyze_vcc2025.py	8167	2026-10-27 16:12	
f	/home/ubuntu/vcc/src/b_panel_tf.py	5465	2026-10-28 05:54	
f	/home/ubuntu/vcc/src/b_techo_i6_ruido_gen.py	4766	2026-10-28 05:56	
f	/home/ubuntu/vcc/src/b_techo_i6_run.sh	1313	2026-10-28 06:05	
f	/home/ubuntu/vcc/src/cpm_stats.py	2961	2026-10-27 01:56	
f	/home/ubuntu/vcc/src/degenerado_check.py	8412	2026-10-27 20:13	
f	/home/ubuntu/vcc/src/explore_structure.py	8074	2026-10-27 19:21	
f	/home/ubuntu/vcc/src/knockdown_grid.py	9536	2026-10-27 19:25	
f	/home/ubuntu/vcc/src/knockdown_grid_v2_compute.py	8994	2026-10-27 19:58	
f	/home/ubuntu/vcc/src/knockdown_grid_v2_scan.py	6185	2026-10-27 19:56	
f	/home/ubuntu/vcc/src/nreal_predictores.py	8402	2026-10-28 19:42	
f	/home/ubuntu/vcc/src/panel_tf.py	2965	2026-10-28 05:01	
f	/home/ubuntu/vcc/src/percell_umi.py	2061	2026-10-27 19:28	
f	/home/ubuntu/vcc/src/pseudobulk.py	1033	2026-10-27 01:56	
f	/home/ubuntu/vcc/src/replogle_contaminacion.py	4606	2026-10-28 04:37	
f	/home/ubuntu/vcc/src/replogle_contaminacion_v2.py	5693	2026-10-28 04:47	
f	/home/ubuntu/vcc/src/replogle_controlexpr.py	11765	2026-10-28 05:28	
f	/home/ubuntu/vcc/src/replogle_controlexpr_porgrupo.py	10269	2026-10-28 07:04	
f	/home/ubuntu/vcc/src/replogle_controlexpr_sample20.py	3111	2026-10-28 05:38	
f	/home/ubuntu/vcc/src/replogle_controlexpr_unfiltered.py	8962	2026-10-28 05:52	
f	/home/ubuntu/vcc/src/replogle_criba_addendum.py	12414	2026-10-27 20:36	
f	/home/ubuntu/vcc/src/replogle_criba_parte1.py	3920	2026-10-27 20:26	
f	/home/ubuntu/vcc/src/replogle_criba_parte3.py	5668	2026-10-27 20:27	
f	/home/ubuntu/vcc/src/replogle_escalon2_partea.py	3655	2026-10-27 20:50	
f	/home/ubuntu/vcc/src/replogle_escalon2_parteb.py	5786	2026-10-27 20:51	
f	/home/ubuntu/vcc/src/replogle_foldexpr_descomposicion.py	3690	2026-10-28 04:51	
f	/home/ubuntu/vcc/src/replogle_foldexpr_dosetapas.py	10405	2026-10-28 05:06	
f	/home/ubuntu/vcc/src/replogle_foldexpr_unfiltered.py	5279	2026-10-28 05:23	
f	/home/ubuntu/vcc/src/replogle_tests_estimando_mw_sample.py	6950	2026-10-28 17:56	
f	/home/ubuntu/vcc/src/replogle_tests_estimando_sample.py	8944	2026-10-28 17:48	
f	/home/ubuntu/vcc/src/t3_inventario.py	20634	2026-10-28 21:03	
f	/home/ubuntu/vcc/src/techo_metricas_run.sh	3800	2026-10-27 22:42	
f	/home/ubuntu/vcc/src/techo_metricas_run_i1_i4_i5.sh	3134	2026-10-28 02:48	
f	/home/ubuntu/vcc/src/techo_metricas_split.py	4155	2026-10-27 22:32	
f	/home/ubuntu/vcc/src/techo_metricas_transforms.py	3280	2026-10-27 22:35	
d	/home/ubuntu/vcc/workflow	4096	2026-10-27 02:00	
d	/home/ubuntu/vcc/workflow/vigilantes	4096	2026-10-27 20:03	
f	/home/ubuntu/vcc/workflow/vigilantes/spot_watch.sh	1220	2026-10-27 02:00	
f	/home/ubuntu/vcc/workflow/vigilantes/sync_out.log	58	2026-10-27 20:03	
f	/home/ubuntu/vcc/workflow/vigilantes/sync_out.pid	10	2026-10-27 20:03	
f	/home/ubuntu/vcc/workflow/vigilantes/sync_out.sh	235	2026-10-27 02:00	
```

## 3 · Árbol ENTERO de `/datos` (77 líneas)

```
d	/datos	4096	2026-10-27 22:31	
d	/datos/_scratch_techo_metricas	4096	2026-10-28 06:02	
f	/datos/_scratch_techo_metricas/predA_base.h5ad	5973295832	2026-10-28 06:02	
f	/datos/_scratch_techo_metricas/realB.h5ad	5971459718	2026-10-28 06:04	
f	/datos/_scratch_techo_metricas/split_log.json	8227	2026-10-28 06:04	
d	/datos/lost+found	16384	2026-10-27 17:28	
d	/datos/ref	4096	2026-10-27 18:30	
d	/datos/ref/convertido	4096	2026-10-27 18:32	
f	/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad	11816834569	2026-10-27 08:50	
f	/datos/ref/convertido/Seurat_object_IFNG_Perturb_seq.h5ad	8035282509	2026-10-27 08:50	
f	/datos/ref/convertido/Seurat_object_INS_Perturb_seq.h5ad	15420720085	2026-10-27 08:50	
f	/datos/ref/convertido/Seurat_object_TGFB_Perturb_seq.h5ad	7298609725	2026-10-27 08:50	
f	/datos/ref/convertido/Seurat_object_TNFA_Perturb_seq.h5ad	12849409777	2026-10-27 08:50	
f	/datos/ref/dl.log	1388	2026-10-27 07:11	
d	/datos/ref/feng	4096	2026-10-27 18:32	
f	/datos/ref/feng/GenomeWideScreen_LFC_byGene.tsv.gz	794610300	2026-10-27 07:11	
f	/datos/ref/feng/TargetedScreen_LFC_byGene-perLine.tsv.gz	1919661029	2026-10-27 07:11	
f	/datos/ref/feng/TargetedScreen_LFC_byGene.tsv.gz	77563731	2026-10-27 07:11	
d	/datos/ref/nadig	4096	2026-10-27 18:34	
f	/datos/ref/nadig/hepg2_raw_sc.h5ad	5614460941	2026-10-27 07:11	
f	/datos/ref/nadig/jurkat_raw_sc.h5ad	9366490264	2026-10-27 07:11	
d	/datos/ref/orion	4096	2026-10-27 18:35	
f	/datos/ref/orion/hct116.parquet	3395023058	2026-10-27 07:11	
f	/datos/ref/orion/hek293t.parquet	1820568723	2026-10-27 07:11	
d	/datos/ref/otros	4096	2026-10-27 18:40	
f	/datos/ref/otros/A_readme.txt	1539	2026-10-27 07:11	
f	/datos/ref/otros/Bulk_RNAseq_Seurat_object_IFNG_and_TGFB_stim.rds	3369918	2026-10-27 07:11	
f	/datos/ref/otros/DE_results_all_pathway.zip	324112587	2026-10-27 07:11	
f	/datos/ref/otros/GSE139944_RAW.tar	9902397440	2026-10-27 07:11	
f	/datos/ref/otros/GSE225775_RAW.tar	13392046080	2026-10-27 07:11	
f	/datos/ref/otros/HClust_Pathway_celltype_specific_genelist.rds	69497	2026-10-27 07:11	
f	/datos/ref/otros/Parse_Guide_Capture_Protocol.pdf	303899	2026-10-27 07:11	
f	/datos/ref/otros/Pathway_Exclusive_genelist.rds	9038	2026-10-27 07:11	
f	/datos/ref/otros/Pathway_genelist.rds	20666	2026-10-27 07:11	
f	/datos/ref/otros/Seurat_object_IFNB_Perturb_seq.rds	4326548669	2026-10-27 07:11	
f	/datos/ref/otros/Seurat_object_IFNG_Perturb_seq.rds	2915636149	2026-10-27 07:12	
f	/datos/ref/otros/Seurat_object_INS_Perturb_seq.rds	5601176410	2026-10-27 07:12	
f	/datos/ref/otros/Seurat_object_TGFB_Perturb_seq.rds	2642041433	2026-10-27 07:13	
f	/datos/ref/otros/Seurat_object_TNFA_Perturb_seq.rds	4656209976	2026-10-27 07:13	
d	/datos/ref/replogle	4096	2026-10-27 18:52	
f	/datos/ref/replogle/k562_essential_raw_bulk.h5ad	79766954	2026-10-27 07:13	
f	/datos/ref/replogle/k562_essential_raw_sc.h5ad	10661879995	2026-10-27 07:13	
f	/datos/ref/replogle/k562_gwps_raw_bulk.h5ad	374587922	2026-10-27 07:13	
f	/datos/ref/replogle/k562_gwps_raw_sc.h5ad	65830941948	2026-10-27 07:13	
f	/datos/ref/replogle/rpe1_raw_bulk.h5ad	95350546	2026-10-27 07:13	
f	/datos/ref/replogle/rpe1_raw_sc.h5ad	8700873216	2026-10-27 07:13	
d	/datos/ref/tfact	4096	2026-10-27 18:48	
f	/datos/ref/tfact/TFactprofiler_all.csv	85747741	2026-10-27 07:13	
d	/datos/reto	4096	2026-10-27 18:25	
d	/datos/reto/final	4096	2026-10-27 18:49	
f	/datos/reto/final/context_D.h5ad	217126334	2026-10-27 07:15	
f	/datos/reto/final/context_E.h5ad	213808470	2026-10-27 07:16	
f	/datos/reto/final/context_F.h5ad	212299193	2026-10-27 07:16	
f	/datos/reto/final/gene_names.csv	119295	2026-10-27 07:16	
f	/datos/reto/final/manifest.json	717	2026-10-27 07:16	
f	/datos/reto/final/pert_counts.csv	1913	2026-10-27 07:16	
f	/datos/reto/final/vcc_2026_controls_final.zip	643356568	2026-10-27 07:16	
d	/datos/reto/val	4096	2026-10-27 18:50	
f	/datos/reto/val/context_A.h5ad	224973316	2026-10-27 07:16	
f	/datos/reto/val/context_B.h5ad	210966111	2026-10-27 07:16	
f	/datos/reto/val/context_C.h5ad	226056691	2026-10-27 07:16	
f	/datos/reto/val/gene_names.csv	119295	2026-10-27 07:17	
f	/datos/reto/val/manifest.json	715	2026-10-27 07:17	
f	/datos/reto/val/pert_counts.csv	1906	2026-10-27 07:17	
f	/datos/reto/val/vcc_2026_controls.zip	662118680	2026-10-27 07:17	
d	/datos/vcc2025	4096	2026-10-27 18:50	
f	/datos/vcc2025/gene_names.csv	116023	2026-10-27 16:47	
d	/datos/vcc2025/test	4096	2026-10-27 18:53	
f	/datos/vcc2025/test/adata_Test.h5ad	11950739168	2026-10-27 16:47	
f	/datos/vcc2025/test/pert_counts_Test.csv	1894	2026-10-27 16:47	
d	/datos/vcc2025/train	4096	2026-10-27 18:55	
f	/datos/vcc2025/train/adata_Training.h5ad	15482497461	2026-10-27 16:47	
f	/datos/vcc2025/train/pert_counts_Training.csv	2824	2026-10-27 16:47	
d	/datos/vcc2025/validation	4096	2026-10-27 18:56	
f	/datos/vcc2025/validation/adata_Validation.h5ad	6928967541	2026-10-27 16:47	
f	/datos/vcc2025/validation/pert_counts_Validation.csv	978	2026-10-27 16:47	
bfs: error: /datos/lost+found: Permission denied.
```

## 4 · S3 por prefijo y subprefijo (dos niveles) — resumen

Los ~968 objetos uno a uno están en `~/vcc/out/mapa_s3.txt` (968 líneas), no aquí. Lo
que sigue es el árbol a dos niveles, con objetos y tamaño total por cada (prefijo,
subprefijo):

| prefijo (nivel 1) | subprefijo (nivel 2) | objetos | tamaño total |
|---|---|---|---|
| `datos/` | `ref/` | 34 | 193.7GiB |
| `datos/` | `reto/` | 14 | 2.4GiB |
| `datos/` | `vcc2025/` | 7 | 32.0GiB |
| `out/` | `(archivos sueltos en out/)` | 101 | 19.6MiB |
| `out/` | `techo_metricas_work/` | 328 | 717.6MiB |
| `scratch_techo/` | `(archivos sueltos en scratch_techo/)` | 3 | 11.1GiB |
| `vcc/` | `(archivos sueltos en vcc/)` | 2 | 122.2MiB |
| `vcc/` | `.claude/` | 1 | 113B |
| `vcc/` | `corpus/` | 2 | 806.9KiB |
| `vcc/` | `out/` | 427 | 737.2MiB |
| `vcc/` | `registro/` | 3 | 16.9KiB |
| `vcc/` | `src/` | 37 | 227.3KiB |
| `vcc/` | `workflow/` | 4 | 1.5KiB |
| `virtual cell/` | `(archivos sueltos en virtual cell/)` | 2 | 10.1KiB |
| `virtual cell/` | `registro/` | 1 | 5.2KiB |
| `virtual cell/` | `src/` | 2 | 3.9KiB |

**Totales por prefijo de nivel 1:**

| prefijo | objetos | tamaño total |
|---|---|---|
| `datos/` | 55 | 228.2GiB |
| `out/` | 429 | 737.2MiB |
| `scratch_techo/` | 3 | 11.1GiB |
| `vcc/` | 476 | 860.4MiB |
| `virtual cell/` | 5 | 19.3KiB |

**Total general: 968 objetos.**

## 5 · Dos columnas por ruta de primer y segundo nivel — volumen y respaldo

### `~/vcc` -- primer y segundo nivel

| ruta | tipo | volumen | respaldo (prefijo S3 vigente) |
|---|---|---|---|
| `/home/ubuntu/vcc/.claude` | dir | nvme0n1p1 (raiz) | `vcc/.claude/` |
| `/home/ubuntu/vcc/PLAN.md` | file | nvme0n1p1 (raiz) | `vcc/PLAN.md` |
| `/home/ubuntu/vcc/corpus` | dir | nvme0n1p1 (raiz) | `vcc/corpus/` |
| `/home/ubuntu/vcc/corpus/inventario.parquet` | file | nvme0n1p1 (raiz) | `vcc/corpus/inventario.parquet` |
| `/home/ubuntu/vcc/corpus/particion.md` | file | nvme0n1p1 (raiz) | `vcc/corpus/particion.md` |
| `/home/ubuntu/vcc/mf.sh` | file | nvme0n1p1 (raiz) | `vcc/mf.sh` |
| `/home/ubuntu/vcc/out` | dir | nvme0n1p1 (raiz) | `vcc/out/` |
| `/home/ubuntu/vcc/out/b_panel_tf.md` | file | nvme0n1p1 (raiz) | `vcc/out/b_panel_tf.md` |
| `/home/ubuntu/vcc/out/b_panel_tf_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/b_panel_tf_raw.json` |
| `/home/ubuntu/vcc/out/b_techo_i6_ruido.md` | file | nvme0n1p1 (raiz) | `vcc/out/b_techo_i6_ruido.md` |
| `/home/ubuntu/vcc/out/degenerado_check3_table.csv` | file | nvme0n1p1 (raiz) | `vcc/out/degenerado_check3_table.csv` |
| `/home/ubuntu/vcc/out/estructura.md` | file | nvme0n1p1 (raiz) | `vcc/out/estructura.md` |
| `/home/ubuntu/vcc/out/geo_zenodo.md` | file | nvme0n1p1 (raiz) | `vcc/out/geo_zenodo.md` |
| `/home/ubuntu/vcc/out/gwps_fiabilidad.md` | file | nvme0n1p1 (raiz) | `vcc/out/gwps_fiabilidad.md` |
| `/home/ubuntu/vcc/out/gwps_fiabilidad_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/gwps_fiabilidad_raw.json` |
| `/home/ubuntu/vcc/out/gwps_resolucion.md` | file | nvme0n1p1 (raiz) | `vcc/out/gwps_resolucion.md` |
| `/home/ubuntu/vcc/out/gwps_resolucion_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/gwps_resolucion_raw.json` |
| `/home/ubuntu/vcc/out/mapa_vcc.txt` | file | nvme0n1p1 (raiz) | `vcc/out/mapa_vcc.txt` |
| `/home/ubuntu/vcc/out/nreal_marcador.md` | file | nvme0n1p1 (raiz) | `vcc/out/nreal_marcador.md` |
| `/home/ubuntu/vcc/out/nreal_marcador_basal_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/nreal_marcador_basal_raw.json` |
| `/home/ubuntu/vcc/out/nreal_predictores.md` | file | nvme0n1p1 (raiz) | `vcc/out/nreal_predictores.md` |
| `/home/ubuntu/vcc/out/nreal_predictores_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/nreal_predictores_raw.json` |
| `/home/ubuntu/vcc/out/nreal_transfiere.md` | file | nvme0n1p1 (raiz) | `vcc/out/nreal_transfiere.md` |
| `/home/ubuntu/vcc/out/nreal_transfiere_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/nreal_transfiere_raw.json` |
| `/home/ubuntu/vcc/out/panel_tf.md` | file | nvme0n1p1 (raiz) | `vcc/out/panel_tf.md` |
| `/home/ubuntu/vcc/out/panel_tf_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/panel_tf_raw.json` |
| `/home/ubuntu/vcc/out/paso0_inventario.md` | file | nvme0n1p1 (raiz) | `vcc/out/paso0_inventario.md` |
| `/home/ubuntu/vcc/out/paso0_obs.md` | file | nvme0n1p1 (raiz) | `vcc/out/paso0_obs.md` |
| `/home/ubuntu/vcc/out/rejilla_table.csv` | file | nvme0n1p1 (raiz) | `vcc/out/rejilla_table.csv` |
| `/home/ubuntu/vcc/out/replogle_addendum_k562_essential.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_addendum_k562_essential.json` |
| `/home/ubuntu/vcc/out/replogle_addendum_k562_gwps.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_addendum_k562_gwps.json` |
| `/home/ubuntu/vcc/out/replogle_addendum_rpe1.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_addendum_rpe1.json` |
| `/home/ubuntu/vcc/out/replogle_bulk_parte1.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_bulk_parte1.json` |
| `/home/ubuntu/vcc/out/replogle_bulk_parte3.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_bulk_parte3.json` |
| `/home/ubuntu/vcc/out/replogle_contaminacion.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_contaminacion.md` |
| `/home/ubuntu/vcc/out/replogle_contaminacion_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_contaminacion_raw.json` |
| `/home/ubuntu/vcc/out/replogle_contaminacion_v2.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_contaminacion_v2.md` |
| `/home/ubuntu/vcc/out/replogle_contaminacion_v2_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_contaminacion_v2_raw.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr.md` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_k562_essential.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_k562_essential.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_k562_gwps.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_k562_gwps.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_porgrupo.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_porgrupo.md` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_porgrupo_k562_essential.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_porgrupo_k562_essential.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_porgrupo_k562_gwps.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_porgrupo_k562_gwps.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_porgrupo_rpe1.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_porgrupo_rpe1.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_rpe1.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_rpe1.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_sample20_k562_essential.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_sample20_k562_essential.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_sample20_k562_gwps.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_sample20_k562_gwps.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_sample20_rpe1.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_sample20_rpe1.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_unfiltered.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_unfiltered.md` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_unfiltered_k562_essential.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_unfiltered_k562_essential.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_unfiltered_k562_gwps.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_unfiltered_k562_gwps.json` |
| `/home/ubuntu/vcc/out/replogle_controlexpr_unfiltered_rpe1.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_controlexpr_unfiltered_rpe1.json` |
| `/home/ubuntu/vcc/out/replogle_criba.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_criba.md` |
| `/home/ubuntu/vcc/out/replogle_dosetapas_detail_k562_essential.csv` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_dosetapas_detail_k562_essential.csv` |
| `/home/ubuntu/vcc/out/replogle_dosetapas_detail_k562_gwps.csv` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_dosetapas_detail_k562_gwps.csv` |
| `/home/ubuntu/vcc/out/replogle_dosetapas_detail_rpe1.csv` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_dosetapas_detail_rpe1.csv` |
| `/home/ubuntu/vcc/out/replogle_dosetapas_k562_essential.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_dosetapas_k562_essential.json` |
| `/home/ubuntu/vcc/out/replogle_dosetapas_k562_gwps.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_dosetapas_k562_gwps.json` |
| `/home/ubuntu/vcc/out/replogle_dosetapas_rpe1.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_dosetapas_rpe1.json` |
| `/home/ubuntu/vcc/out/replogle_escalon2.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_escalon2.md` |
| `/home/ubuntu/vcc/out/replogle_escalon2_partea.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_escalon2_partea.json` |
| `/home/ubuntu/vcc/out/replogle_escalon2_parteb.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_escalon2_parteb.json` |
| `/home/ubuntu/vcc/out/replogle_foldexpr_descomposicion.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_foldexpr_descomposicion.md` |
| `/home/ubuntu/vcc/out/replogle_foldexpr_descomposicion_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_foldexpr_descomposicion_raw.json` |
| `/home/ubuntu/vcc/out/replogle_foldexpr_detail_k562_essential.csv` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_foldexpr_detail_k562_essential.csv` |
| `/home/ubuntu/vcc/out/replogle_foldexpr_detail_k562_gwps.csv` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_foldexpr_detail_k562_gwps.csv` |
| `/home/ubuntu/vcc/out/replogle_foldexpr_detail_rpe1.csv` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_foldexpr_detail_rpe1.csv` |
| `/home/ubuntu/vcc/out/replogle_foldexpr_dosetapas.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_foldexpr_dosetapas.md` |
| `/home/ubuntu/vcc/out/replogle_foldexpr_unfiltered.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_foldexpr_unfiltered.md` |
| `/home/ubuntu/vcc/out/replogle_foldexpr_unfiltered_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_foldexpr_unfiltered_raw.json` |
| `/home/ubuntu/vcc/out/replogle_tests_estimando.md` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_tests_estimando.md` |
| `/home/ubuntu/vcc/out/replogle_tests_estimando_mw_sample_k562_essential.json` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_tests_estimando_mw_sample_k562_essential.json` |
| `/home/ubuntu/vcc/out/replogle_tests_estimando_mw_sample_run.log` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_tests_estimando_mw_sample_run.log` |
| `/home/ubuntu/vcc/out/replogle_tests_estimando_sample_run.log` | file | nvme0n1p1 (raiz) | `vcc/out/replogle_tests_estimando_sample_run.log` |
| `/home/ubuntu/vcc/out/respaldo.md` | file | nvme0n1p1 (raiz) | `vcc/out/respaldo.md` |
| `/home/ubuntu/vcc/out/restore.DONE` | file | nvme0n1p1 (raiz) | `vcc/out/restore.DONE` |
| `/home/ubuntu/vcc/out/restore.log` | file | nvme0n1p1 (raiz) | `vcc/out/restore.log` |
| `/home/ubuntu/vcc/out/s3_convertido.DONE` | file | nvme0n1p1 (raiz) | `vcc/out/s3_convertido.DONE` |
| `/home/ubuntu/vcc/out/s3_convertido.log` | file | nvme0n1p1 (raiz) | `vcc/out/s3_convertido.log` |
| `/home/ubuntu/vcc/out/s3_datos.DONE` | file | nvme0n1p1 (raiz) | `vcc/out/s3_datos.DONE` |
| `/home/ubuntu/vcc/out/s3_datos.log` | file | nvme0n1p1 (raiz) | `vcc/out/s3_datos.log` |
| `/home/ubuntu/vcc/out/t13_comprobaciones.md` | file | nvme0n1p1 (raiz) | `vcc/out/t13_comprobaciones.md` |
| `/home/ubuntu/vcc/out/t13_pools.md` | file | nvme0n1p1 (raiz) | `vcc/out/t13_pools.md` |
| `/home/ubuntu/vcc/out/t13_replicas.md` | file | nvme0n1p1 (raiz) | `vcc/out/t13_replicas.md` |
| `/home/ubuntu/vcc/out/t13_replicas.parquet` | file | nvme0n1p1 (raiz) | `vcc/out/t13_replicas.parquet` |
| `/home/ubuntu/vcc/out/t3_inventario.md` | file | nvme0n1p1 (raiz) | `vcc/out/t3_inventario.md` |
| `/home/ubuntu/vcc/out/t3_inventario_run.log` | file | nvme0n1p1 (raiz) | `vcc/out/t3_inventario_run.log` |
| `/home/ubuntu/vcc/out/t3_q6_q7_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/t3_q6_q7_raw.json` |
| `/home/ubuntu/vcc/out/t3_queries_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/t3_queries_raw.json` |
| `/home/ubuntu/vcc/out/techo_estado.md` | file | nvme0n1p1 (raiz) | `vcc/out/techo_estado.md` |
| `/home/ubuntu/vcc/out/techo_metricas.md` | file | nvme0n1p1 (raiz) | `vcc/out/techo_metricas.md` |
| `/home/ubuntu/vcc/out/techo_metricas_work` | dir | nvme0n1p1 (raiz) | `vcc/out/techo_metricas_work/` |
| `/home/ubuntu/vcc/out/vcc2025.json` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025.json` |
| `/home/ubuntu/vcc/out/vcc2025_analisis.md` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_analisis.md` |
| `/home/ubuntu/vcc/out/vcc2025_analysis_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_analysis_raw.json` |
| `/home/ubuntu/vcc/out/vcc2025_degenerado.md` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_degenerado.md` |
| `/home/ubuntu/vcc/out/vcc2025_degenerado_raw.json` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_degenerado_raw.json` |
| `/home/ubuntu/vcc/out/vcc2025_knockdown_grid.json` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_knockdown_grid.json` |
| `/home/ubuntu/vcc/out/vcc2025_knockdown_raw_arrays.npz` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_knockdown_raw_arrays.npz` |
| `/home/ubuntu/vcc/out/vcc2025_percell_umi.json` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_percell_umi.json` |
| `/home/ubuntu/vcc/out/vcc2025_rejilla.md` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_rejilla.md` |
| `/home/ubuntu/vcc/out/vcc2025_rejilla_grid.json` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_rejilla_grid.json` |
| `/home/ubuntu/vcc/out/vcc2025_rejilla_raw_arrays_test_train.npz` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_rejilla_raw_arrays_test_train.npz` |
| `/home/ubuntu/vcc/out/vcc2025_rejilla_raw_arrays_validation.npz` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_rejilla_raw_arrays_validation.npz` |
| `/home/ubuntu/vcc/out/vcc2025_rejilla_scan_summary_test_train.json` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_rejilla_scan_summary_test_train.json` |
| `/home/ubuntu/vcc/out/vcc2025_rejilla_scan_summary_validation.json` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_rejilla_scan_summary_validation.json` |
| `/home/ubuntu/vcc/out/vcc2025_structure_full.json` | file | nvme0n1p1 (raiz) | `vcc/out/vcc2025_structure_full.json` |
| `/home/ubuntu/vcc/registro` | dir | nvme0n1p1 (raiz) | `vcc/registro/` |
| `/home/ubuntu/vcc/registro/CONSTANTES_A_VERIFICAR.md` | file | nvme0n1p1 (raiz) | `vcc/registro/CONSTANTES_A_VERIFICAR.md` |
| `/home/ubuntu/vcc/registro/limitaciones.md` | file | nvme0n1p1 (raiz) | `vcc/registro/limitaciones.md` |
| `/home/ubuntu/vcc/registro/syncout.service.bak` | file | nvme0n1p1 (raiz) | `vcc/registro/syncout.service.bak` |
| `/home/ubuntu/vcc/src` | dir | nvme0n1p1 (raiz) | `vcc/src/` |
| `/home/ubuntu/vcc/src/__pycache__` | dir | nvme0n1p1 (raiz) | `vcc/src/__pycache__/` |
| `/home/ubuntu/vcc/src/analyze_vcc2025.py` | file | nvme0n1p1 (raiz) | `vcc/src/analyze_vcc2025.py` |
| `/home/ubuntu/vcc/src/b_panel_tf.py` | file | nvme0n1p1 (raiz) | `vcc/src/b_panel_tf.py` |
| `/home/ubuntu/vcc/src/b_techo_i6_ruido_gen.py` | file | nvme0n1p1 (raiz) | `vcc/src/b_techo_i6_ruido_gen.py` |
| `/home/ubuntu/vcc/src/b_techo_i6_run.sh` | file | nvme0n1p1 (raiz) | `vcc/src/b_techo_i6_run.sh` |
| `/home/ubuntu/vcc/src/cpm_stats.py` | file | nvme0n1p1 (raiz) | `vcc/src/cpm_stats.py` |
| `/home/ubuntu/vcc/src/degenerado_check.py` | file | nvme0n1p1 (raiz) | `vcc/src/degenerado_check.py` |
| `/home/ubuntu/vcc/src/explore_structure.py` | file | nvme0n1p1 (raiz) | `vcc/src/explore_structure.py` |
| `/home/ubuntu/vcc/src/knockdown_grid.py` | file | nvme0n1p1 (raiz) | `vcc/src/knockdown_grid.py` |
| `/home/ubuntu/vcc/src/knockdown_grid_v2_compute.py` | file | nvme0n1p1 (raiz) | `vcc/src/knockdown_grid_v2_compute.py` |
| `/home/ubuntu/vcc/src/knockdown_grid_v2_scan.py` | file | nvme0n1p1 (raiz) | `vcc/src/knockdown_grid_v2_scan.py` |
| `/home/ubuntu/vcc/src/nreal_predictores.py` | file | nvme0n1p1 (raiz) | `vcc/src/nreal_predictores.py` |
| `/home/ubuntu/vcc/src/panel_tf.py` | file | nvme0n1p1 (raiz) | `vcc/src/panel_tf.py` |
| `/home/ubuntu/vcc/src/percell_umi.py` | file | nvme0n1p1 (raiz) | `vcc/src/percell_umi.py` |
| `/home/ubuntu/vcc/src/pseudobulk.py` | file | nvme0n1p1 (raiz) | `vcc/src/pseudobulk.py` |
| `/home/ubuntu/vcc/src/replogle_contaminacion.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_contaminacion.py` |
| `/home/ubuntu/vcc/src/replogle_contaminacion_v2.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_contaminacion_v2.py` |
| `/home/ubuntu/vcc/src/replogle_controlexpr.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_controlexpr.py` |
| `/home/ubuntu/vcc/src/replogle_controlexpr_porgrupo.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_controlexpr_porgrupo.py` |
| `/home/ubuntu/vcc/src/replogle_controlexpr_sample20.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_controlexpr_sample20.py` |
| `/home/ubuntu/vcc/src/replogle_controlexpr_unfiltered.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_controlexpr_unfiltered.py` |
| `/home/ubuntu/vcc/src/replogle_criba_addendum.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_criba_addendum.py` |
| `/home/ubuntu/vcc/src/replogle_criba_parte1.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_criba_parte1.py` |
| `/home/ubuntu/vcc/src/replogle_criba_parte3.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_criba_parte3.py` |
| `/home/ubuntu/vcc/src/replogle_escalon2_partea.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_escalon2_partea.py` |
| `/home/ubuntu/vcc/src/replogle_escalon2_parteb.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_escalon2_parteb.py` |
| `/home/ubuntu/vcc/src/replogle_foldexpr_descomposicion.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_foldexpr_descomposicion.py` |
| `/home/ubuntu/vcc/src/replogle_foldexpr_dosetapas.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_foldexpr_dosetapas.py` |
| `/home/ubuntu/vcc/src/replogle_foldexpr_unfiltered.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_foldexpr_unfiltered.py` |
| `/home/ubuntu/vcc/src/replogle_tests_estimando_mw_sample.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_tests_estimando_mw_sample.py` |
| `/home/ubuntu/vcc/src/replogle_tests_estimando_sample.py` | file | nvme0n1p1 (raiz) | `vcc/src/replogle_tests_estimando_sample.py` |
| `/home/ubuntu/vcc/src/t3_inventario.py` | file | nvme0n1p1 (raiz) | `vcc/src/t3_inventario.py` |
| `/home/ubuntu/vcc/src/techo_metricas_run.sh` | file | nvme0n1p1 (raiz) | `vcc/src/techo_metricas_run.sh` |
| `/home/ubuntu/vcc/src/techo_metricas_run_i1_i4_i5.sh` | file | nvme0n1p1 (raiz) | `vcc/src/techo_metricas_run_i1_i4_i5.sh` |
| `/home/ubuntu/vcc/src/techo_metricas_split.py` | file | nvme0n1p1 (raiz) | `vcc/src/techo_metricas_split.py` |
| `/home/ubuntu/vcc/src/techo_metricas_transforms.py` | file | nvme0n1p1 (raiz) | `vcc/src/techo_metricas_transforms.py` |
| `/home/ubuntu/vcc/workflow` | dir | nvme0n1p1 (raiz) | `vcc/workflow/` |
| `/home/ubuntu/vcc/workflow/vigilantes` | dir | nvme0n1p1 (raiz) | `vcc/workflow/vigilantes/` |

### `/datos` -- primer y segundo nivel

| ruta | tipo | volumen | respaldo (prefijo S3 vigente) |
|---|---|---|---|
| `/datos/_scratch_techo_metricas` | dir | nvme1n1 (/datos) | `scratch_techo/` (sincronizado por syncout.service, continuo) |
| `/datos/_scratch_techo_metricas/predA_base.h5ad` | file | nvme1n1 (/datos) | `scratch_techo/` (sincronizado por syncout.service, continuo) |
| `/datos/_scratch_techo_metricas/realB.h5ad` | file | nvme1n1 (/datos) | `scratch_techo/` (sincronizado por syncout.service, continuo) |
| `/datos/_scratch_techo_metricas/split_log.json` | file | nvme1n1 (/datos) | `scratch_techo/` (sincronizado por syncout.service, continuo) |
| `/datos/lost+found` | dir | nvme1n1 (/datos) | NINGUNO (no accesible, permiso denegado) |
| `/datos/ref` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/ref/convertido` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/ref/dl.log` | file | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/ref/feng` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/ref/nadig` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/ref/orion` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/ref/otros` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/ref/replogle` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/ref/tfact` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/reto` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/reto/final` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/reto/val` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/vcc2025` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/vcc2025/gene_names.csv` | file | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/vcc2025/test` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/vcc2025/train` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |
| `/datos/vcc2025/validation` | dir | nvme1n1 (/datos) | `datos/` -- **solo la carga inicial del 27-oct, no se sincroniza mas** |

## 6 · Solapes — qué rutas locales aparecen bajo más de un prefijo de S3

### `out/` (legado, raíz del bucket) vs. `vcc/out/` (vigente)

- `out/` (raíz, legado): **429** objetos.
- `vcc/out/` (vigente): **427** objetos.
- Mismo nombre relativo bajo ambos: **427**, de los cuales **0** tienen tamaño distinto entre copias (la copia vieja quedó más chica -- se escribió antes de que el archivo creciera).
- Solo en `out/` (nunca llegaron a `vcc/out/` en esta ejecución, o fueron borrados de la fuente): **2**.
- Solo en `vcc/out/` (nunca llegaron a `out/`): **0**.

**Único objeto presente en `out/` que NO tiene ninguna ruta local actual**:
`out/techo_metricas_DRAFT.md` (153 bytes, `LastModified` 2026-10-27T23:22:08Z). No
existe en `~/vcc` bajo ningún nombre (`find ~/vcc -iname "*techo_metricas_DRAFT*"` no
devuelve nada) — es el borrador previo a `techo_metricas.md`, sincronizado una vez por
el mecanismo viejo y nunca limpiado de S3 cuando el archivo local se renombró/completó.
Contenido completo (153 bytes, cabe entero):

```
# Techo de las metricas (cell-eval2, perfil vcc2026) — BORRADOR EN CONSTRUCCION

(placeholder; se completa cuando terminen las corridas de cell-eval2)
```

El otro objeto exclusivo de `out/` es `LEGADO.txt` — intencional, añadido en la tarea
anterior, no es una duplicación accidental.

### `virtual cell/` (legado) vs. `vcc/` (vigente) — 4 archivos, 2 de ellos YA divergen en contenido

| ruta relativa | prefijos que la contienen | tamaños por copia |
|---|---|---|
| `LEGADO.txt` | `out/`, `virtual cell/` | `out/`=1417B, `virtual cell/`=753B |
| `PLAN.md` | `vcc/`, `virtual cell/` | `vcc/`=9637B, `virtual cell/`=9637B |
| `registro/CONSTANTES_A_VERIFICAR.md` | `vcc/`, `virtual cell/` | `vcc/`=7699B, `virtual cell/`=5333B |
| `src/cpm_stats.py` | `vcc/`, `virtual cell/` | `vcc/`=2961B, `virtual cell/`=2961B |
| `src/pseudobulk.py` | `vcc/`, `virtual cell/` | `vcc/`=1033B, `virtual cell/`=1033B |

**Nota de lectura de esta tabla**: `LEGADO.txt` aparece en ambos prefijos por
coincidencia de nombre, no por ser el mismo archivo — son los dos marcadores subidos a
propósito en la tarea anterior, con contenido distinto a propósito (cada uno dice cuál
es su propio prefijo vigente). Los otros 3 sí son la misma ruta local real:
`src/cpm_stats.py` y `src/pseudobulk.py` son copias byte-idénticas (mismo tamaño);
**`PLAN.md` también es idéntico** (9,637B en ambos); pero
**`registro/CONSTANTES_A_VERIFICAR.md` ya DIVERGE: 7,699B en `vcc/` (vigente) contra
5,333B en `virtual cell/` (legado)** — el archivo creció después del 27 de octubre y el
prefijo abandonado quedó con la versión vieja, más chica. Es exactamente el tipo de
"copia congelada pareciendo válida" que motivó el `LEGADO.txt`.

## 7 · S3 prefijos sin ninguna ruta local actual

Ninguno de los 5 prefijos de nivel 1 (`datos/`, `out/`, `scratch_techo/`, `vcc/`,
`virtual cell/`) carece por completo de correspondencia local — los 5 tienen origen en
algo que existió o existe en `~/vcc` o `/datos`. Dentro de ellos, el único objeto
individual sin ruta local actual es `out/techo_metricas_DRAFT.md` (§6).

## 8 · Estado de lo que quedó a medias — tal cual, sin arreglar

**`syncout.service`, textual, ahora mismo:**

```ini
# /etc/systemd/system/syncout.service
[Unit]
Description=Sync vcc out a S3
[Service]
User=ubuntu
ExecStart=/bin/bash -c "while true; do aws s3 sync /home/ubuntu/vcc/ s3://vc-challenge-data/vcc/ --no-progress >/dev/null 2>&1; aws s3 sync /datos/_scratch_techo_metricas/ s3://vc-challenge-data/scratch_techo/ --no-progress >/dev/null 2>&1; sleep 120; done"
Restart=always
[Install]
WantedBy=multi-user.target
```

**Prefijos en la raíz del bucket, ahora mismo:**

```
$ aws s3 ls s3://vc-challenge-data/
                           PRE datos/
                           PRE out/
                           PRE scratch_techo/
                           PRE vcc/
                           PRE virtual cell/
```

`scratch_techo/` ya existe (se creó en la tarea anterior) y tiene **3 objetos**:

```
2026-10-28 21:34:45 5973295832 scratch_techo/predA_base.h5ad
2026-10-28 21:34:45 5971459718 scratch_techo/realB.h5ad
2026-10-28 21:34:45       8227 scratch_techo/split_log.json
```

Los tres tamaños coinciden byte a byte con los archivos locales en
`/datos/_scratch_techo_metricas/` (ver §3 del árbol de `/datos` arriba). Nada quedó a
medias en lo que tocó la tarea anterior: la unidad corrigida está activa, el contenido
subió completo, y el escritor del prefijo `out/` legado sigue sin identificar — igual
que se dejó entonces.
