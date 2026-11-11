## v8 means

```
                       slide  runs  median_pcc_all  median_pcc_svg20  median_ssim_svg50  cmd_all  mean_rmse  median_rmse_pos
method                                                                                                                      
GHIST+-Soft         breast5k     9          0.0881            0.4354             0.2775   0.5770     0.2016           0.6383
Niche-Soft          breast5k     9          0.0881            0.4354             0.2775   0.5770     0.2016           0.6383
RRC-soft            breast5k     9          0.1117            0.4878             0.6584   0.4258     0.0981           0.7069
Niche-RRC           breast5k     9          0.1117            0.4878             0.6584   0.4258     0.0981           0.7069
Prog-RRC            breast5k     9          0.1215            0.4931             0.6191   0.5203     0.0978           0.7089
Niche-Prog-RRC      breast5k     9          0.1215            0.4931             0.6191   0.5203     0.0978           0.7089
GHIST+-Soft            atera     9          0.1445            0.5949             0.6349   0.5041     0.2841           0.7334
Niche-Soft             atera     9          0.1445            0.5949             0.6349   0.5041     0.2841           0.7334
RRC-soft               atera     9          0.1860            0.6921             0.7879   0.2955     0.2205           0.7411
Niche-RRC              atera     9          0.1860            0.6921             0.7879   0.2955     0.2205           0.7411
Prog-RRC               atera     9          0.1957            0.6932             0.7888   0.3439     0.2198           0.7417
Niche-Prog-RRC         atera     9          0.1957            0.6932             0.7888   0.3439     0.2198           0.7417
GHIST+-Soft     breast5k_pan     9          0.0715            0.4459             0.3007   0.6026     0.1716           0.6406
Niche-Soft      breast5k_pan     9          0.0715            0.4459             0.3007   0.6026     0.1716           0.6406
RRC-soft        breast5k_pan     9          0.0822            0.5239             0.5908   0.4671     0.0930           0.6965
Niche-RRC       breast5k_pan     9          0.0822            0.5239             0.5908   0.4671     0.0930           0.6965
Prog-RRC        breast5k_pan     9          0.1015            0.5342             0.5317   0.5512     0.0914           0.7001
Niche-Prog-RRC  breast5k_pan     9          0.1015            0.5342             0.5317   0.5512     0.0914           0.7001
```

## v8 vs RRC-soft

```
           slide          method        vs  runs  CI>0  CI<0    dPCC   dSSIM%   dRMSE%   dCMD%
0       breast5k     GHIST+-Soft  RRC-soft     9     0     9 -0.0236 -57.8532 105.5349 35.5154
1       breast5k      Niche-Soft  RRC-soft     9     0     9 -0.0236 -57.8532 105.5349 35.5154
2       breast5k       Niche-RRC  RRC-soft     9     0     0  0.0000   0.0000   0.0000  0.0000
3       breast5k        Prog-RRC  RRC-soft     9     9     0  0.0099  -5.9595  -0.2871 22.2119
4       breast5k  Niche-Prog-RRC  RRC-soft     9     9     0  0.0099  -5.9595  -0.2871 22.2119
5          atera     GHIST+-Soft  RRC-soft     9     0     9 -0.0416 -19.4146  28.8569 70.6155
6          atera      Niche-Soft  RRC-soft     9     0     9 -0.0416 -19.4146  28.8569 70.6155
7          atera       Niche-RRC  RRC-soft     9     0     0  0.0000   0.0000   0.0000  0.0000
8          atera        Prog-RRC  RRC-soft     9     9     0  0.0097   0.1226  -0.2847 16.3792
9          atera  Niche-Prog-RRC  RRC-soft     9     9     0  0.0097   0.1226  -0.2847 16.3792
10  breast5k_pan     GHIST+-Soft  RRC-soft     9     0     9 -0.0107 -49.1028  84.4878 29.0273
11  breast5k_pan      Niche-Soft  RRC-soft     9     0     9 -0.0107 -49.1028  84.4878 29.0273
12  breast5k_pan       Niche-RRC  RRC-soft     9     0     0  0.0000   0.0000   0.0000  0.0000
13  breast5k_pan        Prog-RRC  RRC-soft     9     9     0  0.0192  -9.9922  -1.7424 18.0084
14  breast5k_pan  Niche-Prog-RRC  RRC-soft     9     9     0  0.0192  -9.9922  -1.7424 18.0084
```

## v9 means

```
                       slide  runs  median_pcc_all  median_pcc_svg20  median_ssim_svg50  cmd_all  mean_rmse  median_rmse_pos
method                                                                                                                      
GHIST+-Soft         breast5k     9          0.0881            0.4354             0.2775   0.5770     0.2016           0.6383
RRC-soft            breast5k     9          0.1117            0.4878             0.6584   0.4258     0.0981           0.7069
HAC                 breast5k     9          0.1053            0.4884             0.6987   0.4025     0.0976           0.7120
BioHAN              breast5k     9          0.1176            0.4862             0.7757   0.5013     0.0976           0.7125
BioHAN-noniche      breast5k     9          0.1194            0.4857             0.7757   0.4978     0.0976           0.7127
GHIST+-Soft            atera     9          0.1445            0.5949             0.6349   0.5041     0.2841           0.7334
RRC-soft               atera     9          0.1860            0.6921             0.7879   0.2955     0.2205           0.7411
HAC                    atera     9          0.1897            0.6893             0.7799   0.3185     0.2210           0.7446
BioHAN                 atera     9          0.1918            0.6958             0.8073   0.3137     0.2200           0.7441
BioHAN-noniche         atera     9          0.1888            0.6926             0.8038   0.3129     0.2205           0.7441
GHIST+-Soft     breast5k_pan     9          0.0715            0.4459             0.3007   0.6026     0.1716           0.6406
RRC-soft        breast5k_pan     9          0.0822            0.5239             0.5908   0.4671     0.0930           0.6965
HAC             breast5k_pan     9          0.0735            0.5322             0.7223   0.3794     0.0876           0.7077
BioHAN          breast5k_pan     9          0.0923            0.5253             0.7854   0.4992     0.0874           0.7095
BioHAN-noniche  breast5k_pan     9          0.0919            0.5259             0.7899   0.4913     0.0873           0.7097
```

## v9 vs RRC-soft

```
           slide          method        vs  runs  CI>0  CI<0    dPCC   dSSIM%   dRMSE%    dCMD%
0       breast5k     GHIST+-Soft  RRC-soft     9     0     9 -0.0236 -57.8532 105.5349  35.5154
1       breast5k             HAC  RRC-soft     9     0     7 -0.0064   6.1277  -0.5430  -5.4514
2       breast5k          BioHAN  RRC-soft     9     9     0  0.0060  17.8285  -0.4990  17.7386
3       breast5k  BioHAN-noniche  RRC-soft     9     8     0  0.0077  17.8235  -0.5437  16.9270
4          atera     GHIST+-Soft  RRC-soft     9     0     9 -0.0416 -19.4146  28.8569  70.6155
5          atera             HAC  RRC-soft     9     3     0  0.0037  -1.0178   0.2293   7.7814
6          atera          BioHAN  RRC-soft     9     5     0  0.0057   2.4594  -0.2232   6.1784
7          atera  BioHAN-noniche  RRC-soft     9     2     0  0.0028   2.0187   0.0067   5.8820
8   breast5k_pan     GHIST+-Soft  RRC-soft     9     0     9 -0.0107 -49.1028  84.4878  29.0273
9   breast5k_pan             HAC  RRC-soft     9     0     9 -0.0088  22.2702  -5.8074 -18.7642
10  breast5k_pan          BioHAN  RRC-soft     9     9     0  0.0101  32.9405  -6.0095   6.8705
11  breast5k_pan  BioHAN-noniche  RRC-soft     9     9     0  0.0097  33.7161  -6.0851   5.1815
```

## v9_ablate means

```
                            slide  runs  median_pcc_all  median_pcc_svg20  median_ssim_svg50  cmd_all  mean_rmse  median_rmse_pos
method                                                                                                                           
RRC-soft                 breast5k     6          0.1105            0.4888             0.6535   0.4253     0.0982           0.7071
BioHAN                   breast5k     6          0.1173            0.4863             0.7720   0.4963     0.0977           0.7130
BioHAN-nobottleneck      breast5k     6          0.1176            0.4857             0.7697   0.5004     0.0977           0.7135
BioHAN-noweight          breast5k     6          0.1179            0.4867             0.7761   0.5007     0.0977           0.7128
BioHAN-noniche           breast5k     6          0.1187            0.4859             0.7720   0.4951     0.0977           0.7128
RRC-soft                    atera     6          0.1865            0.6909             0.7862   0.2957     0.2204           0.7415
BioHAN                      atera     6          0.1923            0.6935             0.8031   0.3127     0.2200           0.7449
BioHAN-nobottleneck         atera     6          0.1896            0.6911             0.7972   0.3107     0.2202           0.7444
BioHAN-noweight             atera     6          0.1928            0.6920             0.8041   0.3161     0.2201           0.7440
BioHAN-noniche              atera     6          0.1897            0.6903             0.8011   0.3134     0.2206           0.7450
RRC-soft             breast5k_pan     6          0.0825            0.5212             0.5660   0.4866     0.0944           0.6944
BioHAN               breast5k_pan     6          0.0922            0.5264             0.7890   0.4988     0.0874           0.7101
BioHAN-nobottleneck  breast5k_pan     6          0.0920            0.5268             0.7867   0.5024     0.0875           0.7090
BioHAN-noweight      breast5k_pan     6          0.0919            0.5267             0.7985   0.5033     0.0873           0.7102
BioHAN-noniche       breast5k_pan     6          0.0915            0.5268             0.7916   0.4963     0.0874           0.7099
```

## v9_ablate vs BioHAN

```
           slide               method      vs  runs  CI>0  CI<0    dPCC   dSSIM%  dRMSE%    dCMD%
0       breast5k             RRC-soft  BioHAN     6     0     6 -0.0068 -15.3441  0.5012 -14.3121
1       breast5k  BioHAN-nobottleneck  BioHAN     6     0     1  0.0003  -0.2995 -0.0007   0.8284
2       breast5k      BioHAN-noweight  BioHAN     6     0     0  0.0006   0.5299 -0.0776   0.8800
3       breast5k       BioHAN-noniche  BioHAN     6     1     0  0.0014  -0.0008 -0.0250  -0.2307
4          atera             RRC-soft  BioHAN     6     0     4 -0.0058  -2.1073  0.1971  -5.4224
5          atera  BioHAN-nobottleneck  BioHAN     6     0     1 -0.0028  -0.7339  0.1008  -0.6241
6          atera      BioHAN-noweight  BioHAN     6     0     0  0.0004   0.1219  0.0436   1.0846
7          atera       BioHAN-noniche  BioHAN     6     0     1 -0.0027  -0.2433  0.2491   0.2193
8   breast5k_pan             RRC-soft  BioHAN     6     0     6 -0.0097 -28.2661  8.0288  -2.4371
9   breast5k_pan  BioHAN-nobottleneck  BioHAN     6     0     0 -0.0003  -0.2868  0.1129   0.7372
10  breast5k_pan      BioHAN-noweight  BioHAN     6     1     0 -0.0003   1.2023 -0.0623   0.9066
11  breast5k_pan       BioHAN-noniche  BioHAN     6     0     0 -0.0007   0.3381 -0.0135  -0.4893
```
