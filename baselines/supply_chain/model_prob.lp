\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj1:
Subject To
 c1:  - 6 mu_t0_f0_0 - 9 mu_t0_f0_1 - 7.95 mu_t0_f0_2 - 7.05 mu_t0_f0_3
      - 1.95 mu_t0_f0_4 - 13.05 mu_t0_f0_5 - 4.05 mu_t0_f0_6 - 10.95 mu_t0_f0_7
      - 1.95 mu_t0_f0_8 - 13.05 mu_t0_f0_9 - 12 mu_t0_f0_10 - 3 mu_t0_f0_11
      - 1.95 mu_t0_f0_12 - 13.05 mu_t0_f0_13 - 10.95 mu_t0_f0_14
      - 4.05 mu_t0_f0_15 + x68 + 2 x69 + 4 x70 + 8 x71 <= -1
 c2:  - 13.05 mu_t0_f1_0 - 1.95 mu_t0_f1_1 - 9 mu_t0_f1_2 - 6 mu_t0_f1_3
      - 15 mu_t0_f1_4 - 13.95 mu_t0_f1_6 - 1.05 mu_t0_f1_7 - 7.05 mu_t0_f1_8
      - 7.95 mu_t0_f1_9 - 4.05 mu_t0_f1_10 - 10.95 mu_t0_f1_11
      - 1.05 mu_t0_f1_12 - 13.95 mu_t0_f1_13 - 6 mu_t0_f1_14 - 9 mu_t0_f1_15
      - 1.05 mu_t0_f1_16 - 13.95 mu_t0_f1_17 - 7.95 mu_t0_f1_18
      - 7.05 mu_t0_f1_19 - 13.05 mu_t0_f1_20 - 1.95 mu_t0_f1_21
      - 7.95 mu_t0_f1_22 - 7.05 mu_t0_f1_23 - 7.05 mu_t0_f1_24
      - 7.95 mu_t0_f1_25 - 1.95 mu_t0_f1_26 - 13.05 mu_t0_f1_27
      - 10.05 mu_t0_f1_28 - 4.95 mu_t0_f1_29 - 10.05 mu_t0_f1_30
      - 4.95 mu_t0_f1_31 + x72 + 2 x73 + 4 x74 + 8 x75 <= -1
 c3:  - 9 mu_t0_f2_0 - 6 mu_t0_f2_1 + x76 + 2 x77 + 4 x78 + 8 x79 <= -1
 c4:  - 15 mu_t0_f3_1 - 1.05 mu_t0_f3_2 - 13.95 mu_t0_f3_3 - 3 mu_t0_f3_4
      - 12 mu_t0_f3_5 - 9 mu_t0_f3_6 - 6 mu_t0_f3_7 + x80 + 2 x81 + 4 x82
      + 8 x83 <= -1
 c5:  - 6 mu_t0_f4_0 - 9 mu_t0_f4_1 - 1.05 mu_t0_f4_2 - 13.95 mu_t0_f4_3 + x84
      + 2 x85 + 4 x86 + 8 x87 <= -1
 c6:  - de_t0_12_23 + mu_t0_f0_8 + mu_t0_f0_9 + mu_t0_f0_10 + mu_t0_f0_11
      + mu_t0_f0_12 + mu_t0_f0_13 + mu_t0_f0_14 + mu_t0_f0_15  = 0
 c7:  de_t0_12_23 + mu_t0_f0_0 + mu_t0_f0_1 + mu_t0_f0_2 + mu_t0_f0_3
      + mu_t0_f0_4 + mu_t0_f0_5 + mu_t0_f0_6 + mu_t0_f0_7  = 1
 c8:  - de_t0_8_10 + mu_t0_f0_4 + mu_t0_f0_5 + mu_t0_f0_6 + mu_t0_f0_7
      + mu_t0_f0_12 + mu_t0_f0_13 + mu_t0_f0_14 + mu_t0_f0_15  = 0
 c9:  de_t0_8_10 + mu_t0_f0_0 + mu_t0_f0_1 + mu_t0_f0_2 + mu_t0_f0_3
      + mu_t0_f0_8 + mu_t0_f0_9 + mu_t0_f0_10 + mu_t0_f0_11  = 1
 c10: - de_t0_19_21 + mu_t0_f0_2 + mu_t0_f0_3 + mu_t0_f0_6 + mu_t0_f0_7
      + mu_t0_f0_10 + mu_t0_f0_11 + mu_t0_f0_14 + mu_t0_f0_15  = 0
 c11: de_t0_19_21 + mu_t0_f0_0 + mu_t0_f0_1 + mu_t0_f0_4 + mu_t0_f0_5
      + mu_t0_f0_8 + mu_t0_f0_9 + mu_t0_f0_12 + mu_t0_f0_13  = 1
 c12: - de_t0_12_26 + mu_t0_f0_1 + mu_t0_f0_3 + mu_t0_f0_5 + mu_t0_f0_7
      + mu_t0_f0_9 + mu_t0_f0_11 + mu_t0_f0_13 + mu_t0_f0_15  = 0
 c13: de_t0_12_26 + mu_t0_f0_0 + mu_t0_f0_2 + mu_t0_f0_4 + mu_t0_f0_6
      + mu_t0_f0_8 + mu_t0_f0_10 + mu_t0_f0_12 + mu_t0_f0_14  = 1
 c14: - de_t0_8_10 + mu_t0_f1_16 + mu_t0_f1_17 + mu_t0_f1_18 + mu_t0_f1_19
      + mu_t0_f1_20 + mu_t0_f1_21 + mu_t0_f1_22 + mu_t0_f1_23 + mu_t0_f1_24
      + mu_t0_f1_25 + mu_t0_f1_26 + mu_t0_f1_27 + mu_t0_f1_28 + mu_t0_f1_29
      + mu_t0_f1_30 + mu_t0_f1_31  = 0
 c15: de_t0_8_10 + mu_t0_f1_0 + mu_t0_f1_1 + mu_t0_f1_2 + mu_t0_f1_3
      + mu_t0_f1_4 + mu_t0_f1_5 + mu_t0_f1_6 + mu_t0_f1_7 + mu_t0_f1_8
      + mu_t0_f1_9 + mu_t0_f1_10 + mu_t0_f1_11 + mu_t0_f1_12 + mu_t0_f1_13
      + mu_t0_f1_14 + mu_t0_f1_15  = 1
 c16: - de_t0_19_21 + mu_t0_f1_8 + mu_t0_f1_9 + mu_t0_f1_10 + mu_t0_f1_11
      + mu_t0_f1_12 + mu_t0_f1_13 + mu_t0_f1_14 + mu_t0_f1_15 + mu_t0_f1_24
      + mu_t0_f1_25 + mu_t0_f1_26 + mu_t0_f1_27 + mu_t0_f1_28 + mu_t0_f1_29
      + mu_t0_f1_30 + mu_t0_f1_31  = 0
 c17: de_t0_19_21 + mu_t0_f1_0 + mu_t0_f1_1 + mu_t0_f1_2 + mu_t0_f1_3
      + mu_t0_f1_4 + mu_t0_f1_5 + mu_t0_f1_6 + mu_t0_f1_7 + mu_t0_f1_16
      + mu_t0_f1_17 + mu_t0_f1_18 + mu_t0_f1_19 + mu_t0_f1_20 + mu_t0_f1_21
      + mu_t0_f1_22 + mu_t0_f1_23  = 1
 c18: - de_t0_12_23 + mu_t0_f1_4 + mu_t0_f1_5 + mu_t0_f1_6 + mu_t0_f1_7
      + mu_t0_f1_12 + mu_t0_f1_13 + mu_t0_f1_14 + mu_t0_f1_15 + mu_t0_f1_20
      + mu_t0_f1_21 + mu_t0_f1_22 + mu_t0_f1_23 + mu_t0_f1_28 + mu_t0_f1_29
      + mu_t0_f1_30 + mu_t0_f1_31  = 0
 c19: de_t0_12_23 + mu_t0_f1_0 + mu_t0_f1_1 + mu_t0_f1_2 + mu_t0_f1_3
      + mu_t0_f1_8 + mu_t0_f1_9 + mu_t0_f1_10 + mu_t0_f1_11 + mu_t0_f1_16
      + mu_t0_f1_17 + mu_t0_f1_18 + mu_t0_f1_19 + mu_t0_f1_24 + mu_t0_f1_25
      + mu_t0_f1_26 + mu_t0_f1_27  = 1
 c20: - de_t0_12_26 + mu_t0_f1_2 + mu_t0_f1_3 + mu_t0_f1_6 + mu_t0_f1_7
      + mu_t0_f1_10 + mu_t0_f1_11 + mu_t0_f1_14 + mu_t0_f1_15 + mu_t0_f1_18
      + mu_t0_f1_19 + mu_t0_f1_22 + mu_t0_f1_23 + mu_t0_f1_26 + mu_t0_f1_27
      + mu_t0_f1_30 + mu_t0_f1_31  = 0
 c21: de_t0_12_26 + mu_t0_f1_0 + mu_t0_f1_1 + mu_t0_f1_4 + mu_t0_f1_5
      + mu_t0_f1_8 + mu_t0_f1_9 + mu_t0_f1_12 + mu_t0_f1_13 + mu_t0_f1_16
      + mu_t0_f1_17 + mu_t0_f1_20 + mu_t0_f1_21 + mu_t0_f1_24 + mu_t0_f1_25
      + mu_t0_f1_28 + mu_t0_f1_29  = 1
 c22: - de_t0_11_24 + mu_t0_f1_1 + mu_t0_f1_3 + mu_t0_f1_5 + mu_t0_f1_7
      + mu_t0_f1_9 + mu_t0_f1_11 + mu_t0_f1_13 + mu_t0_f1_15 + mu_t0_f1_17
      + mu_t0_f1_19 + mu_t0_f1_21 + mu_t0_f1_23 + mu_t0_f1_25 + mu_t0_f1_27
      + mu_t0_f1_29 + mu_t0_f1_31  = 0
 c23: de_t0_11_24 + mu_t0_f1_0 + mu_t0_f1_2 + mu_t0_f1_4 + mu_t0_f1_6
      + mu_t0_f1_8 + mu_t0_f1_10 + mu_t0_f1_12 + mu_t0_f1_14 + mu_t0_f1_16
      + mu_t0_f1_18 + mu_t0_f1_20 + mu_t0_f1_22 + mu_t0_f1_24 + mu_t0_f1_26
      + mu_t0_f1_28 + mu_t0_f1_30  = 1
 c24: - de_t0_19_21 + mu_t0_f2_1  = 0
 c25: de_t0_19_21 + mu_t0_f2_0  = 1
 c26: - de_t0_8_10 + mu_t0_f3_4 + mu_t0_f3_5 + mu_t0_f3_6 + mu_t0_f3_7  = 0
 c27: de_t0_8_10 + mu_t0_f3_0 + mu_t0_f3_1 + mu_t0_f3_2 + mu_t0_f3_3  = 1
 c28: - de_t0_19_21 + mu_t0_f3_2 + mu_t0_f3_3 + mu_t0_f3_6 + mu_t0_f3_7  = 0
 c29: de_t0_19_21 + mu_t0_f3_0 + mu_t0_f3_1 + mu_t0_f3_4 + mu_t0_f3_5  = 1
 c30: - de_t0_12_23 + mu_t0_f3_1 + mu_t0_f3_3 + mu_t0_f3_5 + mu_t0_f3_7  = 0
 c31: de_t0_12_23 + mu_t0_f3_0 + mu_t0_f3_2 + mu_t0_f3_4 + mu_t0_f3_6  = 1
 c32: - de_t0_19_21 + mu_t0_f4_2 + mu_t0_f4_3  = 0
 c33: de_t0_19_21 + mu_t0_f4_0 + mu_t0_f4_1  = 1
 c34: - de_t0_8_10 + mu_t0_f4_1 + mu_t0_f4_3  = 0
 c35: de_t0_8_10 + mu_t0_f4_0 + mu_t0_f4_2  = 1
Bounds
 0 <= de_t0_12_26 <= 1
 0 <= de_t0_11_24 <= 1
 0 <= de_t0_19_21 <= 1
 0 <= de_t0_12_23 <= 1
 0 <= de_t0_8_10 <= 1
 0 <= mu_t0_f0_0 <= 1
 0 <= mu_t0_f0_1 <= 1
 0 <= mu_t0_f0_2 <= 1
 0 <= mu_t0_f0_3 <= 1
 0 <= mu_t0_f0_4 <= 1
 0 <= mu_t0_f0_5 <= 1
 0 <= mu_t0_f0_6 <= 1
 0 <= mu_t0_f0_7 <= 1
 0 <= mu_t0_f0_8 <= 1
 0 <= mu_t0_f0_9 <= 1
 0 <= mu_t0_f0_10 <= 1
 0 <= mu_t0_f0_11 <= 1
 0 <= mu_t0_f0_12 <= 1
 0 <= mu_t0_f0_13 <= 1
 0 <= mu_t0_f0_14 <= 1
 0 <= mu_t0_f0_15 <= 1
 0 <= mu_t0_f1_0 <= 1
 0 <= mu_t0_f1_1 <= 1
 0 <= mu_t0_f1_2 <= 1
 0 <= mu_t0_f1_3 <= 1
 0 <= mu_t0_f1_4 <= 1
 0 <= mu_t0_f1_5 <= 1
 0 <= mu_t0_f1_6 <= 1
 0 <= mu_t0_f1_7 <= 1
 0 <= mu_t0_f1_8 <= 1
 0 <= mu_t0_f1_9 <= 1
 0 <= mu_t0_f1_10 <= 1
 0 <= mu_t0_f1_11 <= 1
 0 <= mu_t0_f1_12 <= 1
 0 <= mu_t0_f1_13 <= 1
 0 <= mu_t0_f1_14 <= 1
 0 <= mu_t0_f1_15 <= 1
 0 <= mu_t0_f1_16 <= 1
 0 <= mu_t0_f1_17 <= 1
 0 <= mu_t0_f1_18 <= 1
 0 <= mu_t0_f1_19 <= 1
 0 <= mu_t0_f1_20 <= 1
 0 <= mu_t0_f1_21 <= 1
 0 <= mu_t0_f1_22 <= 1
 0 <= mu_t0_f1_23 <= 1
 0 <= mu_t0_f1_24 <= 1
 0 <= mu_t0_f1_25 <= 1
 0 <= mu_t0_f1_26 <= 1
 0 <= mu_t0_f1_27 <= 1
 0 <= mu_t0_f1_28 <= 1
 0 <= mu_t0_f1_29 <= 1
 0 <= mu_t0_f1_30 <= 1
 0 <= mu_t0_f1_31 <= 1
 0 <= mu_t0_f2_0 <= 1
 0 <= mu_t0_f2_1 <= 1
 0 <= mu_t0_f3_0 <= 1
 0 <= mu_t0_f3_1 <= 1
 0 <= mu_t0_f3_2 <= 1
 0 <= mu_t0_f3_3 <= 1
 0 <= mu_t0_f3_4 <= 1
 0 <= mu_t0_f3_5 <= 1
 0 <= mu_t0_f3_6 <= 1
 0 <= mu_t0_f3_7 <= 1
 0 <= mu_t0_f4_0 <= 1
 0 <= mu_t0_f4_1 <= 1
 0 <= mu_t0_f4_2 <= 1
 0 <= mu_t0_f4_3 <= 1
 0 <= x68 <= 1
 0 <= x69 <= 1
 0 <= x70 <= 1
 0 <= x71 <= 1
 0 <= x72 <= 1
 0 <= x73 <= 1
 0 <= x74 <= 1
 0 <= x75 <= 1
 0 <= x76 <= 1
 0 <= x77 <= 1
 0 <= x78 <= 1
 0 <= x79 <= 1
 0 <= x80 <= 1
 0 <= x81 <= 1
 0 <= x82 <= 1
 0 <= x83 <= 1
 0 <= x84 <= 1
 0 <= x85 <= 1
 0 <= x86 <= 1
 0 <= x87 <= 1
Binaries
 de_t0_12_26  de_t0_11_24  de_t0_19_21  de_t0_12_23  de_t0_8_10  mu_t0_f0_0 
 mu_t0_f0_1  mu_t0_f0_2  mu_t0_f0_3  mu_t0_f0_4  mu_t0_f0_5  mu_t0_f0_6 
 mu_t0_f0_7  mu_t0_f0_8  mu_t0_f0_9  mu_t0_f0_10  mu_t0_f0_11  mu_t0_f0_12 
 mu_t0_f0_13  mu_t0_f0_14  mu_t0_f0_15  mu_t0_f1_0  mu_t0_f1_1  mu_t0_f1_2 
 mu_t0_f1_3  mu_t0_f1_4  mu_t0_f1_5  mu_t0_f1_6  mu_t0_f1_7  mu_t0_f1_8 
 mu_t0_f1_9  mu_t0_f1_10  mu_t0_f1_11  mu_t0_f1_12  mu_t0_f1_13  mu_t0_f1_14 
 mu_t0_f1_15  mu_t0_f1_16  mu_t0_f1_17  mu_t0_f1_18  mu_t0_f1_19  mu_t0_f1_20 
 mu_t0_f1_21  mu_t0_f1_22  mu_t0_f1_23  mu_t0_f1_24  mu_t0_f1_25  mu_t0_f1_26 
 mu_t0_f1_27  mu_t0_f1_28  mu_t0_f1_29  mu_t0_f1_30  mu_t0_f1_31  mu_t0_f2_0 
 mu_t0_f2_1  mu_t0_f3_0  mu_t0_f3_1  mu_t0_f3_2  mu_t0_f3_3  mu_t0_f3_4 
 mu_t0_f3_5  mu_t0_f3_6  mu_t0_f3_7  mu_t0_f4_0  mu_t0_f4_1  mu_t0_f4_2 
 mu_t0_f4_3  x68  x69  x70  x71  x72  x73  x74  x75  x76  x77  x78  x79  x80 
 x81  x82  x83  x84  x85  x86  x87 
End
