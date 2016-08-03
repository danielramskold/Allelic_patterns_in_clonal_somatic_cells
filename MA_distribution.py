# from python ERCCsum_vs_monoallelic_nonDia_minrpkm_v8.py ../bjorn_reinius_complete_set/rpkmforgenes/star_merged_cast1_mm9/ensembl/rpkms_counts_rmnameoverlap.txt nonDia_singlecell_withERCC_faliedorduplicatesamplesremoved.txt nonDia_newS20allelehits_clonalprchrblanking.txt ../snp-validation/chromosome_check/autosomal.genelist :
data = '''clone_3_A_x1_CxB fibroblast 23.3168605413 0.913385826772
clone_3_C_x3_CxB fibroblast 23.4867555325 0.910026462805
clone_3_D_x4_CxB fibroblast 23.4032349294 0.934722222222
clone4_A_1_CxB fibroblast 23.0863763073 0.918564920273
clone4_B_2_CxB fibroblast 21.6107752303 0.94030261348
clone6_A_5_CxB fibroblast 23.2042812329 0.956101601325
clone6_B_6_CxB fibroblast 22.296188261 0.939742536291
clone6_C_7_CxB fibroblast 23.1085237528 0.954831645223
clone6_D_8_CxB fibroblast 23.425457534 0.96075740944
clone8_Div_11_CxB fibroblast 23.6934829535 0.960033076075
clone9_Div_12_CxB fibroblast 23.2487042801 0.963700873362
clone10_Div_13_CxB fibroblast 23.1383945014 0.932328767123
MAF_CxB_clone_A_1 fibroblast 23.9492187914 0.857908847185
MAF_CxB_clone_A_2 fibroblast 24.1259986084 0.833753943218
MAF_CxB_clone_A_3 fibroblast 24.8080988954 0.861990305104
MAF_CxB_clone_A_4 fibroblast 23.9557304076 0.810042946812
MAF_CxB_clone_A_6 fibroblast 24.0597704434 0.807140580172
MAF_CxB_clone_A_7 fibroblast 24.5903473787 0.884739782417
MAF_CxB_clone_A_8 fibroblast 23.9727324307 0.904288554388
MAF_CxB_clone_A_9 fibroblast 23.4111390356 0.904507042254
MAF_CxB_clone_B_21 fibroblast 22.1644943459 0.838912732475
MAF_CxB_clone_B_22 fibroblast 22.5930836183 0.868045848476
MAF_CxB_clone_B_23 fibroblast 21.7160035547 0.816065192084
MAF_CxB_clone_B_25 fibroblast 24.3496922034 0.889710062536
MAF_CxB_clone_B_26 fibroblast 24.5123457881 0.928170594837
MAF_CxB_clone_B_27 fibroblast 24.1524312103 0.91678500426
MAF_CxB_clone_B_28 fibroblast 22.1592752514 0.843812233286
MAF_CxB_clone_B_30 fibroblast 23.9672930406 0.945803224013
MAF_CxB_clone_B_32 fibroblast 23.7091974724 0.896743668244
MAF_CxB_clone_B_33 fibroblast 24.3106780049 0.934044344653
MAF_CxB_clone_B_35 fibroblast 22.6039219447 0.863458310017
MAF_CxB_clone_B_36 fibroblast 23.6011942178 0.913418708241
MAF_CxB_clone_B_37 fibroblast 24.1266595561 0.935457291953
MAF_CxB_clone_B_38 fibroblast 24.1808270477 0.891922005571
MAF_CxB_clone_B_39 fibroblast 23.9129563306 0.933612273361
MAF_CxB_clone_B_40 fibroblast 24.6844810368 0.950988582568
MAF_CxB_clone_B_42 fibroblast 22.5509321515 0.870976770221
MAF_CxB_clone_B_43 fibroblast 22.0298042481 0.842195052602
MAF_CxB_clone_B_44 fibroblast 22.0164032435 0.831671888351
MAF_CxB_clone_B_45 fibroblast 23.5316174568 0.890486725664
MAF_CxB_clone_B_47 fibroblast 22.3063601593 0.872640180332
MAF_CxB_clone_B_48 fibroblast 23.713757681 0.915488976823
MAF_CxB_clone_B_49 fibroblast 23.6048164309 0.91628830874
MAF_CxB_clone_B_51 fibroblast 23.661212355 0.891176470588
MAF_CxB_clone_B_52 fibroblast 22.2939027914 0.855649717514
MAF_CxB_clone_B_54 fibroblast 22.679640434 0.866046249295
MAF_CxB_clone_B_55 fibroblast 22.3393915773 0.84
MAF_CxB_clone_B_56 fibroblast 22.6699587998 0.907417893545
MAF_CxB_clone_B_57 fibroblast 23.1369491759 0.928291316527
MAF_CxB_clone_B_59 fibroblast 21.6667758198 0.787905346188
MAF_CxB_clone_B_60 fibroblast 21.6027507631 0.772923076923
MAF_CxB_clone_B_62 fibroblast 22.3103811913 0.794144556267
MAF_CxB_clone_B_63 fibroblast 21.9814778007 0.827990841442
MAF_CxB_clone_B_64 fibroblast 23.2528879496 0.842538190364
MAF_CxB_clone_B_66 fibroblast 22.3722164651 0.856901408451
MAF_CxB_clone_B_67 fibroblast 22.1950131259 0.848647125141
MAF_CxB_clone_B_70 fibroblast 21.7898002461 0.815893665158
MAF_CxB_clone_B_71 fibroblast 22.1963703732 0.851768033946
MAF_CxB_clone_B_72 fibroblast 22.185365455 0.86012998022
MAF_CxB_clone_B_73 fibroblast 21.4658096031 0.796264367816
MAF_CxB_clone_B_74 fibroblast 22.6265244973 0.872876557191
MAF_CxB_clone_B_78 fibroblast 23.7200657777 0.914430951709
MAF_CxB_clone_B_79 fibroblast 22.1185636963 0.854143337066
MAF_CxB_clone_B_80 fibroblast 22.7266157648 0.878144214645
MAF_CxB_clone_B_81 fibroblast 22.0118454672 0.846457804121
MAF_CxB_clone_B_82-A fibroblast 22.806079303 0.882862734135
MAF_CxB_clone_B_83 fibroblast 21.9675148285 0.84338028169
MAF_CxB_clone_B_84 fibroblast 21.7599719656 0.830580865604
MAF_CxB_clone_B_85 fibroblast 22.2661819757 0.839569160998
MAF_CxB_clone_B_89 fibroblast 21.6750455955 0.825092883681
MAF_CxB_clone_B_90 fibroblast 23.1403776154 0.887894883981
MAF_CxB_clone_B_92 fibroblast 22.3015371309 0.858796950014
MAF_CxB_clone_B_93 fibroblast 22.8213066283 0.872691662003
MAF_CxB_clone_B_94 fibroblast 22.1329040327 0.833851370444
MAF_CxB_clone_B_96 fibroblast 22.5376346759 0.864577504197
MAF_CxB_clone_B_97 fibroblast 22.0742412759 0.842179559571
MAF_CxB_clone_B_98 fibroblast 22.2330013496 0.841216216216
MAF_CxB_clone_B_99 fibroblast 22.535288833 0.854940034266
MAF_CxB_RE_clone_C_105-A fibroblast 24.7982099508 0.911903358869
MEF_E14_Clone_BxC_H_1A fibroblast 22.7561301678 0.90316091954
MEF_E14_Clone_BxC_H_2A fibroblast 22.7854798246 0.844045080719
MEF_E14_Clone_BxC_H_3A fibroblast 22.3802146774 0.822959483265
MEF_E14_Clone_BxC_H_6 fibroblast 22.6165403208 0.86419383787
MEF_E14_Clone_BxC_H_7 fibroblast 22.0705935203 0.853896103896
MEF_E14_Clone_BxC_H_9 fibroblast 22.6586336035 0.802197802198
MEF_E14_Clone_BxC_H_11 fibroblast 21.7096169011 0.756869772999
MEF_E14_Clone_BxC_H_15 fibroblast 22.4765634221 0.890642615558
MEF_E14_Clone_BxC_H_16 fibroblast 22.4412010301 0.916341327384
MEF_E14_Clone_BxC_H_17 fibroblast 22.6234796922 0.916407794408
MEF_E14_Clone_BxC_H_19 fibroblast 22.6440107174 0.817495807714
MEF_E14_Clone_BxC_H_20 fibroblast 21.2863647974 0.85131617009
MEF_E14_Clone_BxC_H_21 fibroblast 22.732741632 0.909963307931
MEF_E14_Clone_BxC_H_23 fibroblast 21.2344438632 0.733195143619
MEF_E14_Clone_BxC_H_24 fibroblast 22.0653948416 0.901616104338
MEF_E14_Clone_BxC_H_25 fibroblast 22.4926385373 0.916364650125
MEF_E14_Clone_BxC_H_26 fibroblast 20.7080881176 0.790759270159
MEF_E14_Clone_BxC_H_27 fibroblast 21.5903481366 0.88403313339
MEF_E14_Clone_BxC_H_28 fibroblast 21.8805658269 0.881693299692
MEF_E14_Clone_BxC_H_29 fibroblast 21.0669992968 0.782088681447
MEF_E14_Clone_BxC_H_30 fibroblast 22.3540083128 0.912315550511
MEF_E14_Clone_BxC_H_31 fibroblast 23.2512817785 0.941747572816
MEF_E14_Clone_BxC_H_32 fibroblast 21.4135082121 0.820498139135
MEF_E14_Clone_BxC_H_33 fibroblast 22.4464556509 0.915900863269
MEF_E14_Clone_BxC_H_34 fibroblast 21.0627186727 0.805202312139
MEF_E14_Clone_BxC_H_35 fibroblast 20.9752173487 0.796031514444
MEF_E14_Clone_BxC_H_36 fibroblast 22.2438128506 0.906689151209
MEF_E14_Clone_BxC_H_37 fibroblast 20.2059918195 0.741729893778
MEF_E14_Clone_BxC_H_38 fibroblast 21.618738578 0.814782858786
MEF_E14_Clone_BxC_H_39 fibroblast 21.7167006262 0.670624281884
MEF_E14_Clone_BxC_H_40 fibroblast 21.8952258402 0.804947398351
MEF_E14_Clone_BxC_H_41 fibroblast 22.4261412217 0.916130835896
MEF_E14_Clone_BxC_H_42 fibroblast 22.0589882745 0.894692502106
MEF_E14_Clone_BxC_H_43 fibroblast 22.3745603991 0.875844594595
MEF_E14_Clone_BxC_H_44 fibroblast 22.1934249371 0.84093637455
MEF_E14_Clone_BxC_H_45 fibroblast 21.8278983563 0.823034522773
MEF_E14_Clone_BxC_H_47 fibroblast 23.11877347 0.937603993344
MEF_E14_Clone_BxC_H_48 fibroblast 21.1914070025 0.820289855072
BQx1_indG_EmbryoMEF_BxC fibroblast 24.1875597174 0.933106575964
BQx2_indG_EmbryoMEF_BxC fibroblast 23.1038497252 0.88082152155
BQx3_indG_EmbryoMEF_BxC fibroblast 21.9434611818 0.830867420946
BQx4_indG_EmbryoMEF_BxC fibroblast 21.8053934292 0.761620603015
BQx5_indG_EmbryoMEF_BxC fibroblast 23.0172824506 0.899381675098
BQx6_indG_EmbryoMEF_BxC fibroblast 22.2806888394 0.853113983549
BQx7_indG_EmbryoMEF_BxC fibroblast 23.4148322442 0.902192729371
BQx8_indG_EmbryoMEF_BxC fibroblast 22.6738839202 0.876603020804
BQx33_indG_EmbryoMEF_BxC fibroblast 23.0561405038 0.741183328476
BQx34_indG_EmbryoMEF_BxC fibroblast 23.1558312903 0.91109868605
BQx35_indG_EmbryoMEF_BxC fibroblast 24.747801154 0.936312849162
BQx36_indG_EmbryoMEF_BxC fibroblast 21.9792724877 0.836953373878
BQx37_indG_EmbryoMEF_BxC fibroblast 22.3477365492 0.86197021764
BQx38_indG_EmbryoMEF_BxC fibroblast 23.6330485299 0.940783986656
BQx39_indG_EmbryoMEF_BxC fibroblast 22.818767871 0.903216783217
BQx40_indG_EmbryoMEF_BxC fibroblast 24.2365361122 0.943137254902
BQx58_indH_EmbryoMEF_BxC fibroblast 22.4318540006 0.855783097779
BQx59_indH_EmbryoMEF_BxC fibroblast 22.0616997267 0.825014679977
BQx60_indH_EmbryoMEF_BxC fibroblast 22.21273463 0.805671641791
BQx61_indH_EmbryoMEF_BxC fibroblast 22.1460959017 0.800509878904
BQx69_indH_EmbryoMEF_BxC fibroblast 23.313708642 0.923400673401
BQx70_indH_EmbryoMEF_BxC fibroblast 22.8221328151 0.918101054969
BQx71_indH_EmbryoMEF_BxC fibroblast 22.6056947522 0.87027818448
BQx72_indH_EmbryoMEF_BxC fibroblast 22.1982835228 0.882553191489
BQx41_indD_EmbryoMEF_BxC fibroblast 23.5160256485 0.877144438296
BQx42_indD_EmbryoMEF_BxC fibroblast 22.7789569981 0.900754822477
BQx43_indD_EmbryoMEF_BxC fibroblast 21.8798790266 0.853428571429
BQx44_indD_EmbryoMEF_BxC fibroblast 22.8289110195 0.911403018446
BQx45_indD_EmbryoMEF_BxC fibroblast 22.219371057 0.838431838432
BQx47_indD_EmbryoMEF_BxC fibroblast 23.3711405119 0.945779765232
BQx48_indD_EmbryoMEF_BxC fibroblast 22.1131691783 0.87074829932
BQx49_indD_EmbryoMEF_BxC fibroblast 22.2083075334 0.874894574079
BQx50_indD_EmbryoMEF_BxC fibroblast 22.143993423 0.887630908576
BQx51_indD_EmbryoMEF_BxC fibroblast 22.9705183508 0.915445321308
BQx52_indD_EmbryoMEF_BxC fibroblast 23.383621892 0.925034965035
BQx53_indD_EmbryoMEF_BxC fibroblast 22.3565501372 0.824416619237
BQx54_indD_EmbryoMEF_BxC fibroblast 22.3577414266 0.858554488989
BQx55_indD_EmbryoMEF_BxC fibroblast 22.2172218845 0.867531003382
BQx57_indH_EmbryoMEF_BxC fibroblast 22.2973429111 0.79813302217
BQx62_indH_EmbryoMEF_BxC fibroblast 23.0141227094 0.836659064994
BQx63_indH_EmbryoMEF_BxC fibroblast 22.5718164289 0.913201228021
BQx64_indH_EmbryoMEF_BxC fibroblast 22.0842151961 0.863506567676
BQx65_indH_EmbryoMEF_BxC fibroblast 21.8981375725 0.84350816853
BQx66_indH_EmbryoMEF_BxC fibroblast 22.9019708657 0.776422764228
BQx67_indH_EmbryoMEF_BxC fibroblast 22.7568401612 0.902762119504
BQx68_indH_EmbryoMEF_BxC fibroblast 22.8272651544 0.900828808231'''

import pylab, dr_tools

if '__main__' == __name__:
	values = [(1-float(line.split()[-1]))*100 for line in data.split('\n')]
	pylab.hist(values, 8)
	pylab.xlabel('number of cells')
	pylab.xlabel('% monoallelic >=20 rpkm (mean across all)')
	pylab.savefig('MA_distribution.pdf')
