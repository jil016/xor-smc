\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj1:
Subject To
 c1:   dummy  = 1
 c2:   alpha_0_0 + alpha_0_2 + alpha_0_4 + alpha_0_6 + alpha_0_8 + alpha_0_10
       + alpha_0_12  = 1
 c3:   alpha_1_0 + alpha_1_2 + alpha_1_4 + alpha_1_6 + alpha_1_8  = 1
 c4:   alpha_2_0 + alpha_2_2 + alpha_2_4 + alpha_2_6 + alpha_2_8 + alpha_2_10
       + alpha_2_12 + alpha_2_14 + alpha_2_16  = 1
 c5:   alpha_3_0 + alpha_3_2 + alpha_3_4 + alpha_3_6 + alpha_3_8 + alpha_3_10
       + alpha_3_12 + alpha_3_14  = 1
 c6:   - alpha_0_0 + zeta_0_0_0 <= 0
 c7:   - alpha_0_2 + zeta_0_0_2 <= 0
 c8:   - alpha_0_4 + zeta_0_0_4 <= 0
 c9:   - alpha_0_6 + zeta_0_0_6 <= 0
 c10:  - alpha_0_8 + zeta_0_0_8 <= 0
 c11:  - alpha_0_10 + zeta_0_0_10 <= 0
 c12:  - alpha_0_12 + zeta_0_0_12 <= 0
 c13:  - de_t0_12_26 + zeta_0_0_0 + zeta_0_0_2 + zeta_0_0_4 + zeta_0_0_6
       + zeta_0_0_8 + zeta_0_0_10 + zeta_0_0_12  = 0
 c14:  - alpha_0_0 + zeta_1_0_0 <= 0
 c15:  - alpha_0_2 + zeta_1_0_2 <= 0
 c16:  - alpha_0_4 + zeta_1_0_4 <= 0
 c17:  - alpha_0_6 + zeta_1_0_6 <= 0
 c18:  - alpha_0_8 + zeta_1_0_8 <= 0
 c19:  - alpha_0_10 + zeta_1_0_10 <= 0
 c20:  - alpha_0_12 + zeta_1_0_12 <= 0
 c21:  - de_t0_11_24 + zeta_1_0_0 + zeta_1_0_2 + zeta_1_0_4 + zeta_1_0_6
       + zeta_1_0_8 + zeta_1_0_10 + zeta_1_0_12  = 0
 c22:  - alpha_1_0 + zeta_1_1_0 <= 0
 c23:  - alpha_1_2 + zeta_1_1_2 <= 0
 c24:  - alpha_1_4 + zeta_1_1_4 <= 0
 c25:  - alpha_1_6 + zeta_1_1_6 <= 0
 c26:  - alpha_1_8 + zeta_1_1_8 <= 0
 c27:  - de_t0_11_24 + zeta_1_1_0 + zeta_1_1_2 + zeta_1_1_4 + zeta_1_1_6
       + zeta_1_1_8  = 0
 c28:  - alpha_2_0 + zeta_1_2_0 <= 0
 c29:  - alpha_2_2 + zeta_1_2_2 <= 0
 c30:  - alpha_2_4 + zeta_1_2_4 <= 0
 c31:  - alpha_2_6 + zeta_1_2_6 <= 0
 c32:  - alpha_2_8 + zeta_1_2_8 <= 0
 c33:  - alpha_2_10 + zeta_1_2_10 <= 0
 c34:  - alpha_2_12 + zeta_1_2_12 <= 0
 c35:  - alpha_2_14 + zeta_1_2_14 <= 0
 c36:  - alpha_2_16 + zeta_1_2_16 <= 0
 c37:  - de_t0_11_24 + zeta_1_2_0 + zeta_1_2_2 + zeta_1_2_4 + zeta_1_2_6
       + zeta_1_2_8 + zeta_1_2_10 + zeta_1_2_12 + zeta_1_2_14 + zeta_1_2_16
        = 0
 c38:  - alpha_0_0 + zeta_2_0_0 <= 0
 c39:  - alpha_0_2 + zeta_2_0_2 <= 0
 c40:  - alpha_0_4 + zeta_2_0_4 <= 0
 c41:  - alpha_0_6 + zeta_2_0_6 <= 0
 c42:  - alpha_0_8 + zeta_2_0_8 <= 0
 c43:  - alpha_0_10 + zeta_2_0_10 <= 0
 c44:  - alpha_0_12 + zeta_2_0_12 <= 0
 c45:  - de_t0_19_21 + zeta_2_0_0 + zeta_2_0_2 + zeta_2_0_4 + zeta_2_0_6
       + zeta_2_0_8 + zeta_2_0_10 + zeta_2_0_12  = 0
 c46:  - alpha_2_0 + zeta_2_2_0 <= 0
 c47:  - alpha_2_2 + zeta_2_2_2 <= 0
 c48:  - alpha_2_4 + zeta_2_2_4 <= 0
 c49:  - alpha_2_6 + zeta_2_2_6 <= 0
 c50:  - alpha_2_8 + zeta_2_2_8 <= 0
 c51:  - alpha_2_10 + zeta_2_2_10 <= 0
 c52:  - alpha_2_12 + zeta_2_2_12 <= 0
 c53:  - alpha_2_14 + zeta_2_2_14 <= 0
 c54:  - alpha_2_16 + zeta_2_2_16 <= 0
 c55:  - de_t0_19_21 + zeta_2_2_0 + zeta_2_2_2 + zeta_2_2_4 + zeta_2_2_6
       + zeta_2_2_8 + zeta_2_2_10 + zeta_2_2_12 + zeta_2_2_14 + zeta_2_2_16
        = 0
 c56:  - alpha_3_0 + zeta_2_3_0 <= 0
 c57:  - alpha_3_2 + zeta_2_3_2 <= 0
 c58:  - alpha_3_4 + zeta_2_3_4 <= 0
 c59:  - alpha_3_6 + zeta_2_3_6 <= 0
 c60:  - alpha_3_8 + zeta_2_3_8 <= 0
 c61:  - alpha_3_10 + zeta_2_3_10 <= 0
 c62:  - alpha_3_12 + zeta_2_3_12 <= 0
 c63:  - alpha_3_14 + zeta_2_3_14 <= 0
 c64:  - de_t0_19_21 + zeta_2_3_0 + zeta_2_3_2 + zeta_2_3_4 + zeta_2_3_6
       + zeta_2_3_8 + zeta_2_3_10 + zeta_2_3_12 + zeta_2_3_14  = 0
 c65:  - alpha_0_0 + zeta_3_0_0 <= 0
 c66:  - alpha_0_2 + zeta_3_0_2 <= 0
 c67:  - alpha_0_4 + zeta_3_0_4 <= 0
 c68:  - alpha_0_6 + zeta_3_0_6 <= 0
 c69:  - alpha_0_8 + zeta_3_0_8 <= 0
 c70:  - alpha_0_10 + zeta_3_0_10 <= 0
 c71:  - alpha_0_12 + zeta_3_0_12 <= 0
 c72:  - de_t0_12_23 + zeta_3_0_0 + zeta_3_0_2 + zeta_3_0_4 + zeta_3_0_6
       + zeta_3_0_8 + zeta_3_0_10 + zeta_3_0_12  = 0
 c73:  - alpha_2_0 + zeta_3_2_0 <= 0
 c74:  - alpha_2_2 + zeta_3_2_2 <= 0
 c75:  - alpha_2_4 + zeta_3_2_4 <= 0
 c76:  - alpha_2_6 + zeta_3_2_6 <= 0
 c77:  - alpha_2_8 + zeta_3_2_8 <= 0
 c78:  - alpha_2_10 + zeta_3_2_10 <= 0
 c79:  - alpha_2_12 + zeta_3_2_12 <= 0
 c80:  - alpha_2_14 + zeta_3_2_14 <= 0
 c81:  - alpha_2_16 + zeta_3_2_16 <= 0
 c82:  - de_t0_12_23 + zeta_3_2_0 + zeta_3_2_2 + zeta_3_2_4 + zeta_3_2_6
       + zeta_3_2_8 + zeta_3_2_10 + zeta_3_2_12 + zeta_3_2_14 + zeta_3_2_16
        = 0
 c83:  - alpha_1_0 + zeta_4_1_0 <= 0
 c84:  - alpha_1_2 + zeta_4_1_2 <= 0
 c85:  - alpha_1_4 + zeta_4_1_4 <= 0
 c86:  - alpha_1_6 + zeta_4_1_6 <= 0
 c87:  - alpha_1_8 + zeta_4_1_8 <= 0
 c88:  - de_t0_8_10 + zeta_4_1_0 + zeta_4_1_2 + zeta_4_1_4 + zeta_4_1_6
       + zeta_4_1_8  = 0
 c89:  - alpha_2_0 + zeta_4_2_0 <= 0
 c90:  - alpha_2_2 + zeta_4_2_2 <= 0
 c91:  - alpha_2_4 + zeta_4_2_4 <= 0
 c92:  - alpha_2_6 + zeta_4_2_6 <= 0
 c93:  - alpha_2_8 + zeta_4_2_8 <= 0
 c94:  - alpha_2_10 + zeta_4_2_10 <= 0
 c95:  - alpha_2_12 + zeta_4_2_12 <= 0
 c96:  - alpha_2_14 + zeta_4_2_14 <= 0
 c97:  - alpha_2_16 + zeta_4_2_16 <= 0
 c98:  - de_t0_8_10 + zeta_4_2_0 + zeta_4_2_2 + zeta_4_2_4 + zeta_4_2_6
       + zeta_4_2_8 + zeta_4_2_10 + zeta_4_2_12 + zeta_4_2_14 + zeta_4_2_16
        = 0
 c99:  - alpha_1_0 + zeta_5_1_0 <= 0
 c100: - alpha_1_2 + zeta_5_1_2 <= 0
 c101: - alpha_1_4 + zeta_5_1_4 <= 0
 c102: - alpha_1_6 + zeta_5_1_6 <= 0
 c103: - alpha_1_8 + zeta_5_1_8 <= 0
 c104: - x6 + zeta_5_1_0 + zeta_5_1_2 + zeta_5_1_4 + zeta_5_1_6 + zeta_5_1_8
        = 0
 c105: - alpha_3_0 + zeta_5_3_0 <= 0
 c106: - alpha_3_2 + zeta_5_3_2 <= 0
 c107: - alpha_3_4 + zeta_5_3_4 <= 0
 c108: - alpha_3_6 + zeta_5_3_6 <= 0
 c109: - alpha_3_8 + zeta_5_3_8 <= 0
 c110: - alpha_3_10 + zeta_5_3_10 <= 0
 c111: - alpha_3_12 + zeta_5_3_12 <= 0
 c112: - alpha_3_14 + zeta_5_3_14 <= 0
 c113: - x6 + zeta_5_3_0 + zeta_5_3_2 + zeta_5_3_4 + zeta_5_3_6 + zeta_5_3_8
       + zeta_5_3_10 + zeta_5_3_12 + zeta_5_3_14  = 0
 c114: - alpha_3_0 + zeta_6_3_0 <= 0
 c115: - alpha_3_2 + zeta_6_3_2 <= 0
 c116: - alpha_3_4 + zeta_6_3_4 <= 0
 c117: - alpha_3_6 + zeta_6_3_6 <= 0
 c118: - alpha_3_8 + zeta_6_3_8 <= 0
 c119: - alpha_3_10 + zeta_6_3_10 <= 0
 c120: - alpha_3_12 + zeta_6_3_12 <= 0
 c121: - alpha_3_14 + zeta_6_3_14 <= 0
 c122: - x7 + zeta_6_3_0 + zeta_6_3_2 + zeta_6_3_4 + zeta_6_3_6 + zeta_6_3_8
       + zeta_6_3_10 + zeta_6_3_12 + zeta_6_3_14  = 0
 c123: - alpha_2_0 + zeta_7_2_0 <= 0
 c124: - alpha_2_2 + zeta_7_2_2 <= 0
 c125: - alpha_2_4 + zeta_7_2_4 <= 0
 c126: - alpha_2_6 + zeta_7_2_6 <= 0
 c127: - alpha_2_8 + zeta_7_2_8 <= 0
 c128: - alpha_2_10 + zeta_7_2_10 <= 0
 c129: - alpha_2_12 + zeta_7_2_12 <= 0
 c130: - alpha_2_14 + zeta_7_2_14 <= 0
 c131: - alpha_2_16 + zeta_7_2_16 <= 0
 c132: - x8 + zeta_7_2_0 + zeta_7_2_2 + zeta_7_2_4 + zeta_7_2_6 + zeta_7_2_8
       + zeta_7_2_10 + zeta_7_2_12 + zeta_7_2_14 + zeta_7_2_16  = 0
 c133: - alpha_3_0 + zeta_7_3_0 <= 0
 c134: - alpha_3_2 + zeta_7_3_2 <= 0
 c135: - alpha_3_4 + zeta_7_3_4 <= 0
 c136: - alpha_3_6 + zeta_7_3_6 <= 0
 c137: - alpha_3_8 + zeta_7_3_8 <= 0
 c138: - alpha_3_10 + zeta_7_3_10 <= 0
 c139: - alpha_3_12 + zeta_7_3_12 <= 0
 c140: - alpha_3_14 + zeta_7_3_14 <= 0
 c141: - x8 + zeta_7_3_0 + zeta_7_3_2 + zeta_7_3_4 + zeta_7_3_6 + zeta_7_3_8
       + zeta_7_3_10 + zeta_7_3_12 + zeta_7_3_14  = 0
 c142: - alpha_0_0 + zeta_8_0_0 <= 0
 c143: - alpha_0_2 + zeta_8_0_2 <= 0
 c144: - alpha_0_4 + zeta_8_0_4 <= 0
 c145: - alpha_0_6 + zeta_8_0_6 <= 0
 c146: - alpha_0_8 + zeta_8_0_8 <= 0
 c147: - alpha_0_10 + zeta_8_0_10 <= 0
 c148: - alpha_0_12 + zeta_8_0_12 <= 0
 c149: - x9 + zeta_8_0_0 + zeta_8_0_2 + zeta_8_0_4 + zeta_8_0_6 + zeta_8_0_8
       + zeta_8_0_10 + zeta_8_0_12  = 0
 c150: - alpha_2_0 + zeta_8_2_0 <= 0
 c151: - alpha_2_2 + zeta_8_2_2 <= 0
 c152: - alpha_2_4 + zeta_8_2_4 <= 0
 c153: - alpha_2_6 + zeta_8_2_6 <= 0
 c154: - alpha_2_8 + zeta_8_2_8 <= 0
 c155: - alpha_2_10 + zeta_8_2_10 <= 0
 c156: - alpha_2_12 + zeta_8_2_12 <= 0
 c157: - alpha_2_14 + zeta_8_2_14 <= 0
 c158: - alpha_2_16 + zeta_8_2_16 <= 0
 c159: - x9 + zeta_8_2_0 + zeta_8_2_2 + zeta_8_2_4 + zeta_8_2_6 + zeta_8_2_8
       + zeta_8_2_10 + zeta_8_2_12 + zeta_8_2_14 + zeta_8_2_16  = 0
 c160: - alpha_3_0 + zeta_8_3_0 <= 0
 c161: - alpha_3_2 + zeta_8_3_2 <= 0
 c162: - alpha_3_4 + zeta_8_3_4 <= 0
 c163: - alpha_3_6 + zeta_8_3_6 <= 0
 c164: - alpha_3_8 + zeta_8_3_8 <= 0
 c165: - alpha_3_10 + zeta_8_3_10 <= 0
 c166: - alpha_3_12 + zeta_8_3_12 <= 0
 c167: - alpha_3_14 + zeta_8_3_14 <= 0
 c168: - x9 + zeta_8_3_0 + zeta_8_3_2 + zeta_8_3_4 + zeta_8_3_6 + zeta_8_3_8
       + zeta_8_3_10 + zeta_8_3_12 + zeta_8_3_14  = 0
 c169: - alpha_0_0 + zeta_9_0_0 <= 0
 c170: - alpha_0_2 + zeta_9_0_2 <= 0
 c171: - alpha_0_4 + zeta_9_0_4 <= 0
 c172: - alpha_0_6 + zeta_9_0_6 <= 0
 c173: - alpha_0_8 + zeta_9_0_8 <= 0
 c174: - alpha_0_10 + zeta_9_0_10 <= 0
 c175: - alpha_0_12 + zeta_9_0_12 <= 0
 c176: - x10 + zeta_9_0_0 + zeta_9_0_2 + zeta_9_0_4 + zeta_9_0_6 + zeta_9_0_8
       + zeta_9_0_10 + zeta_9_0_12  = 0
 c177: - alpha_2_0 + zeta_9_2_0 <= 0
 c178: - alpha_2_2 + zeta_9_2_2 <= 0
 c179: - alpha_2_4 + zeta_9_2_4 <= 0
 c180: - alpha_2_6 + zeta_9_2_6 <= 0
 c181: - alpha_2_8 + zeta_9_2_8 <= 0
 c182: - alpha_2_10 + zeta_9_2_10 <= 0
 c183: - alpha_2_12 + zeta_9_2_12 <= 0
 c184: - alpha_2_14 + zeta_9_2_14 <= 0
 c185: - alpha_2_16 + zeta_9_2_16 <= 0
 c186: - x10 + zeta_9_2_0 + zeta_9_2_2 + zeta_9_2_4 + zeta_9_2_6 + zeta_9_2_8
       + zeta_9_2_10 + zeta_9_2_12 + zeta_9_2_14 + zeta_9_2_16  = 0
 c187: - alpha_1_0 + zeta_10_1_0 <= 0
 c188: - alpha_1_2 + zeta_10_1_2 <= 0
 c189: - alpha_1_4 + zeta_10_1_4 <= 0
 c190: - alpha_1_6 + zeta_10_1_6 <= 0
 c191: - alpha_1_8 + zeta_10_1_8 <= 0
 c192: - x11 + zeta_10_1_0 + zeta_10_1_2 + zeta_10_1_4 + zeta_10_1_6
       + zeta_10_1_8  = 0
 c193: - alpha_2_0 + zeta_10_2_0 <= 0
 c194: - alpha_2_2 + zeta_10_2_2 <= 0
 c195: - alpha_2_4 + zeta_10_2_4 <= 0
 c196: - alpha_2_6 + zeta_10_2_6 <= 0
 c197: - alpha_2_8 + zeta_10_2_8 <= 0
 c198: - alpha_2_10 + zeta_10_2_10 <= 0
 c199: - alpha_2_12 + zeta_10_2_12 <= 0
 c200: - alpha_2_14 + zeta_10_2_14 <= 0
 c201: - alpha_2_16 + zeta_10_2_16 <= 0
 c202: - x11 + zeta_10_2_0 + zeta_10_2_2 + zeta_10_2_4 + zeta_10_2_6
       + zeta_10_2_8 + zeta_10_2_10 + zeta_10_2_12 + zeta_10_2_14
       + zeta_10_2_16  = 0
 c203: - alpha_0_0 + zeta_11_0_0 <= 0
 c204: - alpha_0_2 + zeta_11_0_2 <= 0
 c205: - alpha_0_4 + zeta_11_0_4 <= 0
 c206: - alpha_0_6 + zeta_11_0_6 <= 0
 c207: - alpha_0_8 + zeta_11_0_8 <= 0
 c208: - alpha_0_10 + zeta_11_0_10 <= 0
 c209: - alpha_0_12 + zeta_11_0_12 <= 0
 c210: - x12 + zeta_11_0_0 + zeta_11_0_2 + zeta_11_0_4 + zeta_11_0_6
       + zeta_11_0_8 + zeta_11_0_10 + zeta_11_0_12  = 0
 c211: - alpha_2_0 + zeta_11_2_0 <= 0
 c212: - alpha_2_2 + zeta_11_2_2 <= 0
 c213: - alpha_2_4 + zeta_11_2_4 <= 0
 c214: - alpha_2_6 + zeta_11_2_6 <= 0
 c215: - alpha_2_8 + zeta_11_2_8 <= 0
 c216: - alpha_2_10 + zeta_11_2_10 <= 0
 c217: - alpha_2_12 + zeta_11_2_12 <= 0
 c218: - alpha_2_14 + zeta_11_2_14 <= 0
 c219: - alpha_2_16 + zeta_11_2_16 <= 0
 c220: - x12 + zeta_11_2_0 + zeta_11_2_2 + zeta_11_2_4 + zeta_11_2_6
       + zeta_11_2_8 + zeta_11_2_10 + zeta_11_2_12 + zeta_11_2_14
       + zeta_11_2_16  = 0
 c221: - alpha_0_0 + zeta_12_0_0 <= 0
 c222: - alpha_0_2 + zeta_12_0_2 <= 0
 c223: - alpha_0_4 + zeta_12_0_4 <= 0
 c224: - alpha_0_6 + zeta_12_0_6 <= 0
 c225: - alpha_0_8 + zeta_12_0_8 <= 0
 c226: - alpha_0_10 + zeta_12_0_10 <= 0
 c227: - alpha_0_12 + zeta_12_0_12 <= 0
 c228: - x13 + zeta_12_0_0 + zeta_12_0_2 + zeta_12_0_4 + zeta_12_0_6
       + zeta_12_0_8 + zeta_12_0_10 + zeta_12_0_12  = 0
 c229: - alpha_3_0 + zeta_12_3_0 <= 0
 c230: - alpha_3_2 + zeta_12_3_2 <= 0
 c231: - alpha_3_4 + zeta_12_3_4 <= 0
 c232: - alpha_3_6 + zeta_12_3_6 <= 0
 c233: - alpha_3_8 + zeta_12_3_8 <= 0
 c234: - alpha_3_10 + zeta_12_3_10 <= 0
 c235: - alpha_3_12 + zeta_12_3_12 <= 0
 c236: - alpha_3_14 + zeta_12_3_14 <= 0
 c237: - x13 + zeta_12_3_0 + zeta_12_3_2 + zeta_12_3_4 + zeta_12_3_6
       + zeta_12_3_8 + zeta_12_3_10 + zeta_12_3_12 + zeta_12_3_14  = 0
 c238: - alpha_2_0 + zeta_13_2_0 <= 0
 c239: - alpha_2_2 + zeta_13_2_2 <= 0
 c240: - alpha_2_4 + zeta_13_2_4 <= 0
 c241: - alpha_2_6 + zeta_13_2_6 <= 0
 c242: - alpha_2_8 + zeta_13_2_8 <= 0
 c243: - alpha_2_10 + zeta_13_2_10 <= 0
 c244: - alpha_2_12 + zeta_13_2_12 <= 0
 c245: - alpha_2_14 + zeta_13_2_14 <= 0
 c246: - alpha_2_16 + zeta_13_2_16 <= 0
 c247: - x14 + zeta_13_2_0 + zeta_13_2_2 + zeta_13_2_4 + zeta_13_2_6
       + zeta_13_2_8 + zeta_13_2_10 + zeta_13_2_12 + zeta_13_2_14
       + zeta_13_2_16  = 0
 c248: - alpha_1_0 + zeta_14_1_0 <= 0
 c249: - alpha_1_2 + zeta_14_1_2 <= 0
 c250: - alpha_1_4 + zeta_14_1_4 <= 0
 c251: - alpha_1_6 + zeta_14_1_6 <= 0
 c252: - alpha_1_8 + zeta_14_1_8 <= 0
 c253: - x15 + zeta_14_1_0 + zeta_14_1_2 + zeta_14_1_4 + zeta_14_1_6
       + zeta_14_1_8  = 0
 c254: - alpha_2_0 + zeta_14_2_0 <= 0
 c255: - alpha_2_2 + zeta_14_2_2 <= 0
 c256: - alpha_2_4 + zeta_14_2_4 <= 0
 c257: - alpha_2_6 + zeta_14_2_6 <= 0
 c258: - alpha_2_8 + zeta_14_2_8 <= 0
 c259: - alpha_2_10 + zeta_14_2_10 <= 0
 c260: - alpha_2_12 + zeta_14_2_12 <= 0
 c261: - alpha_2_14 + zeta_14_2_14 <= 0
 c262: - alpha_2_16 + zeta_14_2_16 <= 0
 c263: - x15 + zeta_14_2_0 + zeta_14_2_2 + zeta_14_2_4 + zeta_14_2_6
       + zeta_14_2_8 + zeta_14_2_10 + zeta_14_2_12 + zeta_14_2_14
       + zeta_14_2_16  = 0
 c264: - alpha_3_0 + zeta_14_3_0 <= 0
 c265: - alpha_3_2 + zeta_14_3_2 <= 0
 c266: - alpha_3_4 + zeta_14_3_4 <= 0
 c267: - alpha_3_6 + zeta_14_3_6 <= 0
 c268: - alpha_3_8 + zeta_14_3_8 <= 0
 c269: - alpha_3_10 + zeta_14_3_10 <= 0
 c270: - alpha_3_12 + zeta_14_3_12 <= 0
 c271: - alpha_3_14 + zeta_14_3_14 <= 0
 c272: - x15 + zeta_14_3_0 + zeta_14_3_2 + zeta_14_3_4 + zeta_14_3_6
       + zeta_14_3_8 + zeta_14_3_10 + zeta_14_3_12 + zeta_14_3_14  = 0
 c273: - alpha_0_0 + zeta_15_0_0 <= 0
 c274: - alpha_0_2 + zeta_15_0_2 <= 0
 c275: - alpha_0_4 + zeta_15_0_4 <= 0
 c276: - alpha_0_6 + zeta_15_0_6 <= 0
 c277: - alpha_0_8 + zeta_15_0_8 <= 0
 c278: - alpha_0_10 + zeta_15_0_10 <= 0
 c279: - alpha_0_12 + zeta_15_0_12 <= 0
 c280: - x16 + zeta_15_0_0 + zeta_15_0_2 + zeta_15_0_4 + zeta_15_0_6
       + zeta_15_0_8 + zeta_15_0_10 + zeta_15_0_12  = 0
 c281: - alpha_3_0 + zeta_15_3_0 <= 0
 c282: - alpha_3_2 + zeta_15_3_2 <= 0
 c283: - alpha_3_4 + zeta_15_3_4 <= 0
 c284: - alpha_3_6 + zeta_15_3_6 <= 0
 c285: - alpha_3_8 + zeta_15_3_8 <= 0
 c286: - alpha_3_10 + zeta_15_3_10 <= 0
 c287: - alpha_3_12 + zeta_15_3_12 <= 0
 c288: - alpha_3_14 + zeta_15_3_14 <= 0
 c289: - x16 + zeta_15_3_0 + zeta_15_3_2 + zeta_15_3_4 + zeta_15_3_6
       + zeta_15_3_8 + zeta_15_3_10 + zeta_15_3_12 + zeta_15_3_14  = 0
 c290: - alpha_2_0 + zeta_16_2_0 <= 0
 c291: - alpha_2_2 + zeta_16_2_2 <= 0
 c292: - alpha_2_4 + zeta_16_2_4 <= 0
 c293: - alpha_2_6 + zeta_16_2_6 <= 0
 c294: - alpha_2_8 + zeta_16_2_8 <= 0
 c295: - alpha_2_10 + zeta_16_2_10 <= 0
 c296: - alpha_2_12 + zeta_16_2_12 <= 0
 c297: - alpha_2_14 + zeta_16_2_14 <= 0
 c298: - alpha_2_16 + zeta_16_2_16 <= 0
 c299: - x17 + zeta_16_2_0 + zeta_16_2_2 + zeta_16_2_4 + zeta_16_2_6
       + zeta_16_2_8 + zeta_16_2_10 + zeta_16_2_12 + zeta_16_2_14
       + zeta_16_2_16  = 0
 c300: - alpha_0_0 + zeta_17_0_0 <= 0
 c301: - alpha_0_2 + zeta_17_0_2 <= 0
 c302: - alpha_0_4 + zeta_17_0_4 <= 0
 c303: - alpha_0_6 + zeta_17_0_6 <= 0
 c304: - alpha_0_8 + zeta_17_0_8 <= 0
 c305: - alpha_0_10 + zeta_17_0_10 <= 0
 c306: - alpha_0_12 + zeta_17_0_12 <= 0
 c307: - x18 + zeta_17_0_0 + zeta_17_0_2 + zeta_17_0_4 + zeta_17_0_6
       + zeta_17_0_8 + zeta_17_0_10 + zeta_17_0_12  = 0
 c308: - alpha_3_0 + zeta_18_3_0 <= 0
 c309: - alpha_3_2 + zeta_18_3_2 <= 0
 c310: - alpha_3_4 + zeta_18_3_4 <= 0
 c311: - alpha_3_6 + zeta_18_3_6 <= 0
 c312: - alpha_3_8 + zeta_18_3_8 <= 0
 c313: - alpha_3_10 + zeta_18_3_10 <= 0
 c314: - alpha_3_12 + zeta_18_3_12 <= 0
 c315: - alpha_3_14 + zeta_18_3_14 <= 0
 c316: - x19 + zeta_18_3_0 + zeta_18_3_2 + zeta_18_3_4 + zeta_18_3_6
       + zeta_18_3_8 + zeta_18_3_10 + zeta_18_3_12 + zeta_18_3_14  = 0
 c317: - alpha_2_0 + zeta_19_2_0 <= 0
 c318: - alpha_2_2 + zeta_19_2_2 <= 0
 c319: - alpha_2_4 + zeta_19_2_4 <= 0
 c320: - alpha_2_6 + zeta_19_2_6 <= 0
 c321: - alpha_2_8 + zeta_19_2_8 <= 0
 c322: - alpha_2_10 + zeta_19_2_10 <= 0
 c323: - alpha_2_12 + zeta_19_2_12 <= 0
 c324: - alpha_2_14 + zeta_19_2_14 <= 0
 c325: - alpha_2_16 + zeta_19_2_16 <= 0
 c326: - x20 + zeta_19_2_0 + zeta_19_2_2 + zeta_19_2_4 + zeta_19_2_6
       + zeta_19_2_8 + zeta_19_2_10 + zeta_19_2_12 + zeta_19_2_14
       + zeta_19_2_16  = 0
 c327: - alpha_3_0 + zeta_19_3_0 <= 0
 c328: - alpha_3_2 + zeta_19_3_2 <= 0
 c329: - alpha_3_4 + zeta_19_3_4 <= 0
 c330: - alpha_3_6 + zeta_19_3_6 <= 0
 c331: - alpha_3_8 + zeta_19_3_8 <= 0
 c332: - alpha_3_10 + zeta_19_3_10 <= 0
 c333: - alpha_3_12 + zeta_19_3_12 <= 0
 c334: - alpha_3_14 + zeta_19_3_14 <= 0
 c335: - x20 + zeta_19_3_0 + zeta_19_3_2 + zeta_19_3_4 + zeta_19_3_6
       + zeta_19_3_8 + zeta_19_3_10 + zeta_19_3_12 + zeta_19_3_14  = 0
 c336: - alpha_1_0 + zeta_20_1_0 <= 0
 c337: - alpha_1_2 + zeta_20_1_2 <= 0
 c338: - alpha_1_4 + zeta_20_1_4 <= 0
 c339: - alpha_1_6 + zeta_20_1_6 <= 0
 c340: - alpha_1_8 + zeta_20_1_8 <= 0
 c341: - x21 + zeta_20_1_0 + zeta_20_1_2 + zeta_20_1_4 + zeta_20_1_6
       + zeta_20_1_8  = 0
 c342: - alpha_2_0 + zeta_21_2_0 <= 0
 c343: - alpha_2_2 + zeta_21_2_2 <= 0
 c344: - alpha_2_4 + zeta_21_2_4 <= 0
 c345: - alpha_2_6 + zeta_21_2_6 <= 0
 c346: - alpha_2_8 + zeta_21_2_8 <= 0
 c347: - alpha_2_10 + zeta_21_2_10 <= 0
 c348: - alpha_2_12 + zeta_21_2_12 <= 0
 c349: - alpha_2_14 + zeta_21_2_14 <= 0
 c350: - alpha_2_16 + zeta_21_2_16 <= 0
 c351: - x22 + zeta_21_2_0 + zeta_21_2_2 + zeta_21_2_4 + zeta_21_2_6
       + zeta_21_2_8 + zeta_21_2_10 + zeta_21_2_12 + zeta_21_2_14
       + zeta_21_2_16  = 0
 c352: - alpha_0_0 + zeta_22_0_0 <= 0
 c353: - alpha_0_2 + zeta_22_0_2 <= 0
 c354: - alpha_0_4 + zeta_22_0_4 <= 0
 c355: - alpha_0_6 + zeta_22_0_6 <= 0
 c356: - alpha_0_8 + zeta_22_0_8 <= 0
 c357: - alpha_0_10 + zeta_22_0_10 <= 0
 c358: - alpha_0_12 + zeta_22_0_12 <= 0
 c359: - x23 + zeta_22_0_0 + zeta_22_0_2 + zeta_22_0_4 + zeta_22_0_6
       + zeta_22_0_8 + zeta_22_0_10 + zeta_22_0_12  = 0
 c360: - alpha_3_0 + zeta_23_3_0 <= 0
 c361: - alpha_3_2 + zeta_23_3_2 <= 0
 c362: - alpha_3_4 + zeta_23_3_4 <= 0
 c363: - alpha_3_6 + zeta_23_3_6 <= 0
 c364: - alpha_3_8 + zeta_23_3_8 <= 0
 c365: - alpha_3_10 + zeta_23_3_10 <= 0
 c366: - alpha_3_12 + zeta_23_3_12 <= 0
 c367: - alpha_3_14 + zeta_23_3_14 <= 0
 c368: - x24 + zeta_23_3_0 + zeta_23_3_2 + zeta_23_3_4 + zeta_23_3_6
       + zeta_23_3_8 + zeta_23_3_10 + zeta_23_3_12 + zeta_23_3_14  = 0
 c369: - alpha_2_0 + zeta_24_2_0 <= 0
 c370: - alpha_2_2 + zeta_24_2_2 <= 0
 c371: - alpha_2_4 + zeta_24_2_4 <= 0
 c372: - alpha_2_6 + zeta_24_2_6 <= 0
 c373: - alpha_2_8 + zeta_24_2_8 <= 0
 c374: - alpha_2_10 + zeta_24_2_10 <= 0
 c375: - alpha_2_12 + zeta_24_2_12 <= 0
 c376: - alpha_2_14 + zeta_24_2_14 <= 0
 c377: - alpha_2_16 + zeta_24_2_16 <= 0
 c378: - x25 + zeta_24_2_0 + zeta_24_2_2 + zeta_24_2_4 + zeta_24_2_6
       + zeta_24_2_8 + zeta_24_2_10 + zeta_24_2_12 + zeta_24_2_14
       + zeta_24_2_16  = 0
 c379: - alpha_3_0 + zeta_24_3_0 <= 0
 c380: - alpha_3_2 + zeta_24_3_2 <= 0
 c381: - alpha_3_4 + zeta_24_3_4 <= 0
 c382: - alpha_3_6 + zeta_24_3_6 <= 0
 c383: - alpha_3_8 + zeta_24_3_8 <= 0
 c384: - alpha_3_10 + zeta_24_3_10 <= 0
 c385: - alpha_3_12 + zeta_24_3_12 <= 0
 c386: - alpha_3_14 + zeta_24_3_14 <= 0
 c387: - x25 + zeta_24_3_0 + zeta_24_3_2 + zeta_24_3_4 + zeta_24_3_6
       + zeta_24_3_8 + zeta_24_3_10 + zeta_24_3_12 + zeta_24_3_14  = 0
 c388: - alpha_1_0 + zeta_25_1_0 <= 0
 c389: - alpha_1_2 + zeta_25_1_2 <= 0
 c390: - alpha_1_4 + zeta_25_1_4 <= 0
 c391: - alpha_1_6 + zeta_25_1_6 <= 0
 c392: - alpha_1_8 + zeta_25_1_8 <= 0
 c393: - x26 + zeta_25_1_0 + zeta_25_1_2 + zeta_25_1_4 + zeta_25_1_6
       + zeta_25_1_8  = 0
 c394: - alpha_2_0 + zeta_26_2_0 <= 0
 c395: - alpha_2_2 + zeta_26_2_2 <= 0
 c396: - alpha_2_4 + zeta_26_2_4 <= 0
 c397: - alpha_2_6 + zeta_26_2_6 <= 0
 c398: - alpha_2_8 + zeta_26_2_8 <= 0
 c399: - alpha_2_10 + zeta_26_2_10 <= 0
 c400: - alpha_2_12 + zeta_26_2_12 <= 0
 c401: - alpha_2_14 + zeta_26_2_14 <= 0
 c402: - alpha_2_16 + zeta_26_2_16 <= 0
 c403: - x27 + zeta_26_2_0 + zeta_26_2_2 + zeta_26_2_4 + zeta_26_2_6
       + zeta_26_2_8 + zeta_26_2_10 + zeta_26_2_12 + zeta_26_2_14
       + zeta_26_2_16  = 0
 c404: - alpha_0_0 + zeta_27_0_0 <= 0
 c405: - alpha_0_2 + zeta_27_0_2 <= 0
 c406: - alpha_0_4 + zeta_27_0_4 <= 0
 c407: - alpha_0_6 + zeta_27_0_6 <= 0
 c408: - alpha_0_8 + zeta_27_0_8 <= 0
 c409: - alpha_0_10 + zeta_27_0_10 <= 0
 c410: - alpha_0_12 + zeta_27_0_12 <= 0
 c411: - x28 + zeta_27_0_0 + zeta_27_0_2 + zeta_27_0_4 + zeta_27_0_6
       + zeta_27_0_8 + zeta_27_0_10 + zeta_27_0_12  = 0
 c412: - alpha_3_0 + zeta_28_3_0 <= 0
 c413: - alpha_3_2 + zeta_28_3_2 <= 0
 c414: - alpha_3_4 + zeta_28_3_4 <= 0
 c415: - alpha_3_6 + zeta_28_3_6 <= 0
 c416: - alpha_3_8 + zeta_28_3_8 <= 0
 c417: - alpha_3_10 + zeta_28_3_10 <= 0
 c418: - alpha_3_12 + zeta_28_3_12 <= 0
 c419: - alpha_3_14 + zeta_28_3_14 <= 0
 c420: - x29 + zeta_28_3_0 + zeta_28_3_2 + zeta_28_3_4 + zeta_28_3_6
       + zeta_28_3_8 + zeta_28_3_10 + zeta_28_3_12 + zeta_28_3_14  = 0
 c421: - alpha_0_0 + zeta_29_0_0 <= 0
 c422: - alpha_0_2 + zeta_29_0_2 <= 0
 c423: - alpha_0_4 + zeta_29_0_4 <= 0
 c424: - alpha_0_6 + zeta_29_0_6 <= 0
 c425: - alpha_0_8 + zeta_29_0_8 <= 0
 c426: - alpha_0_10 + zeta_29_0_10 <= 0
 c427: - alpha_0_12 + zeta_29_0_12 <= 0
 c428: - x30 + zeta_29_0_0 + zeta_29_0_2 + zeta_29_0_4 + zeta_29_0_6
       + zeta_29_0_8 + zeta_29_0_10 + zeta_29_0_12  = 0
 c429: - alpha_1_0 + zeta_29_1_0 <= 0
 c430: - alpha_1_2 + zeta_29_1_2 <= 0
 c431: - alpha_1_4 + zeta_29_1_4 <= 0
 c432: - alpha_1_6 + zeta_29_1_6 <= 0
 c433: - alpha_1_8 + zeta_29_1_8 <= 0
 c434: - x30 + zeta_29_1_0 + zeta_29_1_2 + zeta_29_1_4 + zeta_29_1_6
       + zeta_29_1_8  = 0
 c435: - alpha_2_0 + zeta_30_2_0 <= 0
 c436: - alpha_2_2 + zeta_30_2_2 <= 0
 c437: - alpha_2_4 + zeta_30_2_4 <= 0
 c438: - alpha_2_6 + zeta_30_2_6 <= 0
 c439: - alpha_2_8 + zeta_30_2_8 <= 0
 c440: - alpha_2_10 + zeta_30_2_10 <= 0
 c441: - alpha_2_12 + zeta_30_2_12 <= 0
 c442: - alpha_2_14 + zeta_30_2_14 <= 0
 c443: - alpha_2_16 + zeta_30_2_16 <= 0
 c444: - x31 + zeta_30_2_0 + zeta_30_2_2 + zeta_30_2_4 + zeta_30_2_6
       + zeta_30_2_8 + zeta_30_2_10 + zeta_30_2_12 + zeta_30_2_14
       + zeta_30_2_16  = 0
 c445: - alpha_3_0 + zeta_30_3_0 <= 0
 c446: - alpha_3_2 + zeta_30_3_2 <= 0
 c447: - alpha_3_4 + zeta_30_3_4 <= 0
 c448: - alpha_3_6 + zeta_30_3_6 <= 0
 c449: - alpha_3_8 + zeta_30_3_8 <= 0
 c450: - alpha_3_10 + zeta_30_3_10 <= 0
 c451: - alpha_3_12 + zeta_30_3_12 <= 0
 c452: - alpha_3_14 + zeta_30_3_14 <= 0
 c453: - x31 + zeta_30_3_0 + zeta_30_3_2 + zeta_30_3_4 + zeta_30_3_6
       + zeta_30_3_8 + zeta_30_3_10 + zeta_30_3_12 + zeta_30_3_14  = 0
 c454: - alpha_3_0 + zeta_31_3_0 <= 0
 c455: - alpha_3_2 + zeta_31_3_2 <= 0
 c456: - alpha_3_4 + zeta_31_3_4 <= 0
 c457: - alpha_3_6 + zeta_31_3_6 <= 0
 c458: - alpha_3_8 + zeta_31_3_8 <= 0
 c459: - alpha_3_10 + zeta_31_3_10 <= 0
 c460: - alpha_3_12 + zeta_31_3_12 <= 0
 c461: - alpha_3_14 + zeta_31_3_14 <= 0
 c462: zeta_31_3_0 + zeta_31_3_2 + zeta_31_3_4 + zeta_31_3_6 + zeta_31_3_8
       + zeta_31_3_10 + zeta_31_3_12 + zeta_31_3_14 - dummy  = 0
 c463: zeta_0_0_0 + zeta_1_0_0 + zeta_2_0_0 + zeta_3_0_0 + zeta_8_0_0
       + zeta_9_0_0 + zeta_11_0_0 + zeta_12_0_0 + zeta_15_0_0 + zeta_17_0_0
       + zeta_22_0_0 + zeta_27_0_0 + zeta_29_0_0  = 0
 c464: - 2 alpha_0_2 + zeta_0_0_2 + zeta_1_0_2 + zeta_2_0_2 + zeta_3_0_2
       + zeta_8_0_2 + zeta_9_0_2 + zeta_11_0_2 + zeta_12_0_2 + zeta_15_0_2
       + zeta_17_0_2 + zeta_22_0_2 + zeta_27_0_2 + zeta_29_0_2  = 0
 c465: - 4 alpha_0_4 + zeta_0_0_4 + zeta_1_0_4 + zeta_2_0_4 + zeta_3_0_4
       + zeta_8_0_4 + zeta_9_0_4 + zeta_11_0_4 + zeta_12_0_4 + zeta_15_0_4
       + zeta_17_0_4 + zeta_22_0_4 + zeta_27_0_4 + zeta_29_0_4  = 0
 c466: - 6 alpha_0_6 + zeta_0_0_6 + zeta_1_0_6 + zeta_2_0_6 + zeta_3_0_6
       + zeta_8_0_6 + zeta_9_0_6 + zeta_11_0_6 + zeta_12_0_6 + zeta_15_0_6
       + zeta_17_0_6 + zeta_22_0_6 + zeta_27_0_6 + zeta_29_0_6  = 0
 c467: - 8 alpha_0_8 + zeta_0_0_8 + zeta_1_0_8 + zeta_2_0_8 + zeta_3_0_8
       + zeta_8_0_8 + zeta_9_0_8 + zeta_11_0_8 + zeta_12_0_8 + zeta_15_0_8
       + zeta_17_0_8 + zeta_22_0_8 + zeta_27_0_8 + zeta_29_0_8  = 0
 c468: - 10 alpha_0_10 + zeta_0_0_10 + zeta_1_0_10 + zeta_2_0_10 + zeta_3_0_10
       + zeta_8_0_10 + zeta_9_0_10 + zeta_11_0_10 + zeta_12_0_10 + zeta_15_0_10
       + zeta_17_0_10 + zeta_22_0_10 + zeta_27_0_10 + zeta_29_0_10  = 0
 c469: - 12 alpha_0_12 + zeta_0_0_12 + zeta_1_0_12 + zeta_2_0_12 + zeta_3_0_12
       + zeta_8_0_12 + zeta_9_0_12 + zeta_11_0_12 + zeta_12_0_12 + zeta_15_0_12
       + zeta_17_0_12 + zeta_22_0_12 + zeta_27_0_12 + zeta_29_0_12  = 0
 c470: zeta_1_1_0 + zeta_4_1_0 + zeta_5_1_0 + zeta_10_1_0 + zeta_14_1_0
       + zeta_20_1_0 + zeta_25_1_0 + zeta_29_1_0  = 0
 c471: - 2 alpha_1_2 + zeta_1_1_2 + zeta_4_1_2 + zeta_5_1_2 + zeta_10_1_2
       + zeta_14_1_2 + zeta_20_1_2 + zeta_25_1_2 + zeta_29_1_2  = 0
 c472: - 4 alpha_1_4 + zeta_1_1_4 + zeta_4_1_4 + zeta_5_1_4 + zeta_10_1_4
       + zeta_14_1_4 + zeta_20_1_4 + zeta_25_1_4 + zeta_29_1_4  = 0
 c473: - 6 alpha_1_6 + zeta_1_1_6 + zeta_4_1_6 + zeta_5_1_6 + zeta_10_1_6
       + zeta_14_1_6 + zeta_20_1_6 + zeta_25_1_6 + zeta_29_1_6  = 0
 c474: - 8 alpha_1_8 + zeta_1_1_8 + zeta_4_1_8 + zeta_5_1_8 + zeta_10_1_8
       + zeta_14_1_8 + zeta_20_1_8 + zeta_25_1_8 + zeta_29_1_8  = 0
 c475: zeta_1_2_0 + zeta_2_2_0 + zeta_3_2_0 + zeta_4_2_0 + zeta_7_2_0
       + zeta_8_2_0 + zeta_9_2_0 + zeta_10_2_0 + zeta_11_2_0 + zeta_13_2_0
       + zeta_14_2_0 + zeta_16_2_0 + zeta_19_2_0 + zeta_21_2_0 + zeta_24_2_0
       + zeta_26_2_0 + zeta_30_2_0  = 0
 c476: - 2 alpha_2_2 + zeta_1_2_2 + zeta_2_2_2 + zeta_3_2_2 + zeta_4_2_2
       + zeta_7_2_2 + zeta_8_2_2 + zeta_9_2_2 + zeta_10_2_2 + zeta_11_2_2
       + zeta_13_2_2 + zeta_14_2_2 + zeta_16_2_2 + zeta_19_2_2 + zeta_21_2_2
       + zeta_24_2_2 + zeta_26_2_2 + zeta_30_2_2  = 0
 c477: - 4 alpha_2_4 + zeta_1_2_4 + zeta_2_2_4 + zeta_3_2_4 + zeta_4_2_4
       + zeta_7_2_4 + zeta_8_2_4 + zeta_9_2_4 + zeta_10_2_4 + zeta_11_2_4
       + zeta_13_2_4 + zeta_14_2_4 + zeta_16_2_4 + zeta_19_2_4 + zeta_21_2_4
       + zeta_24_2_4 + zeta_26_2_4 + zeta_30_2_4  = 0
 c478: - 6 alpha_2_6 + zeta_1_2_6 + zeta_2_2_6 + zeta_3_2_6 + zeta_4_2_6
       + zeta_7_2_6 + zeta_8_2_6 + zeta_9_2_6 + zeta_10_2_6 + zeta_11_2_6
       + zeta_13_2_6 + zeta_14_2_6 + zeta_16_2_6 + zeta_19_2_6 + zeta_21_2_6
       + zeta_24_2_6 + zeta_26_2_6 + zeta_30_2_6  = 0
 c479: - 8 alpha_2_8 + zeta_1_2_8 + zeta_2_2_8 + zeta_3_2_8 + zeta_4_2_8
       + zeta_7_2_8 + zeta_8_2_8 + zeta_9_2_8 + zeta_10_2_8 + zeta_11_2_8
       + zeta_13_2_8 + zeta_14_2_8 + zeta_16_2_8 + zeta_19_2_8 + zeta_21_2_8
       + zeta_24_2_8 + zeta_26_2_8 + zeta_30_2_8  = 0
 c480: - 10 alpha_2_10 + zeta_1_2_10 + zeta_2_2_10 + zeta_3_2_10 + zeta_4_2_10
       + zeta_7_2_10 + zeta_8_2_10 + zeta_9_2_10 + zeta_10_2_10 + zeta_11_2_10
       + zeta_13_2_10 + zeta_14_2_10 + zeta_16_2_10 + zeta_19_2_10
       + zeta_21_2_10 + zeta_24_2_10 + zeta_26_2_10 + zeta_30_2_10  = 0
 c481: - 12 alpha_2_12 + zeta_1_2_12 + zeta_2_2_12 + zeta_3_2_12 + zeta_4_2_12
       + zeta_7_2_12 + zeta_8_2_12 + zeta_9_2_12 + zeta_10_2_12 + zeta_11_2_12
       + zeta_13_2_12 + zeta_14_2_12 + zeta_16_2_12 + zeta_19_2_12
       + zeta_21_2_12 + zeta_24_2_12 + zeta_26_2_12 + zeta_30_2_12  = 0
 c482: - 14 alpha_2_14 + zeta_1_2_14 + zeta_2_2_14 + zeta_3_2_14 + zeta_4_2_14
       + zeta_7_2_14 + zeta_8_2_14 + zeta_9_2_14 + zeta_10_2_14 + zeta_11_2_14
       + zeta_13_2_14 + zeta_14_2_14 + zeta_16_2_14 + zeta_19_2_14
       + zeta_21_2_14 + zeta_24_2_14 + zeta_26_2_14 + zeta_30_2_14  = 0
 c483: - 16 alpha_2_16 + zeta_1_2_16 + zeta_2_2_16 + zeta_3_2_16 + zeta_4_2_16
       + zeta_7_2_16 + zeta_8_2_16 + zeta_9_2_16 + zeta_10_2_16 + zeta_11_2_16
       + zeta_13_2_16 + zeta_14_2_16 + zeta_16_2_16 + zeta_19_2_16
       + zeta_21_2_16 + zeta_24_2_16 + zeta_26_2_16 + zeta_30_2_16  = 0
 c484: zeta_2_3_0 + zeta_5_3_0 + zeta_6_3_0 + zeta_7_3_0 + zeta_8_3_0
       + zeta_12_3_0 + zeta_14_3_0 + zeta_15_3_0 + zeta_18_3_0 + zeta_19_3_0
       + zeta_23_3_0 + zeta_24_3_0 + zeta_28_3_0 + zeta_30_3_0 + zeta_31_3_0
        = 0
 c485: - 2 alpha_3_2 + zeta_2_3_2 + zeta_5_3_2 + zeta_6_3_2 + zeta_7_3_2
       + zeta_8_3_2 + zeta_12_3_2 + zeta_14_3_2 + zeta_15_3_2 + zeta_18_3_2
       + zeta_19_3_2 + zeta_23_3_2 + zeta_24_3_2 + zeta_28_3_2 + zeta_30_3_2
       + zeta_31_3_2  = 0
 c486: - 4 alpha_3_4 + zeta_2_3_4 + zeta_5_3_4 + zeta_6_3_4 + zeta_7_3_4
       + zeta_8_3_4 + zeta_12_3_4 + zeta_14_3_4 + zeta_15_3_4 + zeta_18_3_4
       + zeta_19_3_4 + zeta_23_3_4 + zeta_24_3_4 + zeta_28_3_4 + zeta_30_3_4
       + zeta_31_3_4  = 0
 c487: - 6 alpha_3_6 + zeta_2_3_6 + zeta_5_3_6 + zeta_6_3_6 + zeta_7_3_6
       + zeta_8_3_6 + zeta_12_3_6 + zeta_14_3_6 + zeta_15_3_6 + zeta_18_3_6
       + zeta_19_3_6 + zeta_23_3_6 + zeta_24_3_6 + zeta_28_3_6 + zeta_30_3_6
       + zeta_31_3_6  = 0
 c488: - 8 alpha_3_8 + zeta_2_3_8 + zeta_5_3_8 + zeta_6_3_8 + zeta_7_3_8
       + zeta_8_3_8 + zeta_12_3_8 + zeta_14_3_8 + zeta_15_3_8 + zeta_18_3_8
       + zeta_19_3_8 + zeta_23_3_8 + zeta_24_3_8 + zeta_28_3_8 + zeta_30_3_8
       + zeta_31_3_8  = 0
 c489: - 10 alpha_3_10 + zeta_2_3_10 + zeta_5_3_10 + zeta_6_3_10 + zeta_7_3_10
       + zeta_8_3_10 + zeta_12_3_10 + zeta_14_3_10 + zeta_15_3_10
       + zeta_18_3_10 + zeta_19_3_10 + zeta_23_3_10 + zeta_24_3_10
       + zeta_28_3_10 + zeta_30_3_10 + zeta_31_3_10  = 0
 c490: - 12 alpha_3_12 + zeta_2_3_12 + zeta_5_3_12 + zeta_6_3_12 + zeta_7_3_12
       + zeta_8_3_12 + zeta_12_3_12 + zeta_14_3_12 + zeta_15_3_12
       + zeta_18_3_12 + zeta_19_3_12 + zeta_23_3_12 + zeta_24_3_12
       + zeta_28_3_12 + zeta_30_3_12 + zeta_31_3_12  = 0
 c491: - 14 alpha_3_14 + zeta_2_3_14 + zeta_5_3_14 + zeta_6_3_14 + zeta_7_3_14
       + zeta_8_3_14 + zeta_12_3_14 + zeta_14_3_14 + zeta_15_3_14
       + zeta_18_3_14 + zeta_19_3_14 + zeta_23_3_14 + zeta_24_3_14
       + zeta_28_3_14 + zeta_30_3_14 + zeta_31_3_14  = 0
Bounds
 0 <= de_t0_12_26 <= 1
 0 <= de_t0_11_24 <= 1
 0 <= de_t0_19_21 <= 1
 0 <= de_t0_12_23 <= 1
 0 <= de_t0_8_10 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1
 0 <= x8 <= 1
 0 <= x9 <= 1
 0 <= x10 <= 1
 0 <= x11 <= 1
 0 <= x12 <= 1
 0 <= x13 <= 1
 0 <= x14 <= 1
 0 <= x15 <= 1
 0 <= x16 <= 1
 0 <= x17 <= 1
 0 <= x18 <= 1
 0 <= x19 <= 1
 0 <= x20 <= 1
 0 <= x21 <= 1
 0 <= x22 <= 1
 0 <= x23 <= 1
 0 <= x24 <= 1
 0 <= x25 <= 1
 0 <= x26 <= 1
 0 <= x27 <= 1
 0 <= x28 <= 1
 0 <= x29 <= 1
 0 <= x30 <= 1
 0 <= x31 <= 1
 0 <= alpha_0_0 <= 1
 0 <= alpha_0_2 <= 1
 0 <= alpha_0_4 <= 1
 0 <= alpha_0_6 <= 1
 0 <= alpha_0_8 <= 1
 0 <= alpha_0_10 <= 1
 0 <= alpha_0_12 <= 1
 0 <= alpha_1_0 <= 1
 0 <= alpha_1_2 <= 1
 0 <= alpha_1_4 <= 1
 0 <= alpha_1_6 <= 1
 0 <= alpha_1_8 <= 1
 0 <= alpha_2_0 <= 1
 0 <= alpha_2_2 <= 1
 0 <= alpha_2_4 <= 1
 0 <= alpha_2_6 <= 1
 0 <= alpha_2_8 <= 1
 0 <= alpha_2_10 <= 1
 0 <= alpha_2_12 <= 1
 0 <= alpha_2_14 <= 1
 0 <= alpha_2_16 <= 1
 0 <= alpha_3_0 <= 1
 0 <= alpha_3_2 <= 1
 0 <= alpha_3_4 <= 1
 0 <= alpha_3_6 <= 1
 0 <= alpha_3_8 <= 1
 0 <= alpha_3_10 <= 1
 0 <= alpha_3_12 <= 1
 0 <= alpha_3_14 <= 1
 0 <= zeta_0_0_0 <= 1
 0 <= zeta_0_0_2 <= 1
 0 <= zeta_0_0_4 <= 1
 0 <= zeta_0_0_6 <= 1
 0 <= zeta_0_0_8 <= 1
 0 <= zeta_0_0_10 <= 1
 0 <= zeta_0_0_12 <= 1
 0 <= zeta_1_0_0 <= 1
 0 <= zeta_1_0_2 <= 1
 0 <= zeta_1_0_4 <= 1
 0 <= zeta_1_0_6 <= 1
 0 <= zeta_1_0_8 <= 1
 0 <= zeta_1_0_10 <= 1
 0 <= zeta_1_0_12 <= 1
 0 <= zeta_1_1_0 <= 1
 0 <= zeta_1_1_2 <= 1
 0 <= zeta_1_1_4 <= 1
 0 <= zeta_1_1_6 <= 1
 0 <= zeta_1_1_8 <= 1
 0 <= zeta_1_2_0 <= 1
 0 <= zeta_1_2_2 <= 1
 0 <= zeta_1_2_4 <= 1
 0 <= zeta_1_2_6 <= 1
 0 <= zeta_1_2_8 <= 1
 0 <= zeta_1_2_10 <= 1
 0 <= zeta_1_2_12 <= 1
 0 <= zeta_1_2_14 <= 1
 0 <= zeta_1_2_16 <= 1
 0 <= zeta_2_0_0 <= 1
 0 <= zeta_2_0_2 <= 1
 0 <= zeta_2_0_4 <= 1
 0 <= zeta_2_0_6 <= 1
 0 <= zeta_2_0_8 <= 1
 0 <= zeta_2_0_10 <= 1
 0 <= zeta_2_0_12 <= 1
 0 <= zeta_2_2_0 <= 1
 0 <= zeta_2_2_2 <= 1
 0 <= zeta_2_2_4 <= 1
 0 <= zeta_2_2_6 <= 1
 0 <= zeta_2_2_8 <= 1
 0 <= zeta_2_2_10 <= 1
 0 <= zeta_2_2_12 <= 1
 0 <= zeta_2_2_14 <= 1
 0 <= zeta_2_2_16 <= 1
 0 <= zeta_2_3_0 <= 1
 0 <= zeta_2_3_2 <= 1
 0 <= zeta_2_3_4 <= 1
 0 <= zeta_2_3_6 <= 1
 0 <= zeta_2_3_8 <= 1
 0 <= zeta_2_3_10 <= 1
 0 <= zeta_2_3_12 <= 1
 0 <= zeta_2_3_14 <= 1
 0 <= zeta_3_0_0 <= 1
 0 <= zeta_3_0_2 <= 1
 0 <= zeta_3_0_4 <= 1
 0 <= zeta_3_0_6 <= 1
 0 <= zeta_3_0_8 <= 1
 0 <= zeta_3_0_10 <= 1
 0 <= zeta_3_0_12 <= 1
 0 <= zeta_3_2_0 <= 1
 0 <= zeta_3_2_2 <= 1
 0 <= zeta_3_2_4 <= 1
 0 <= zeta_3_2_6 <= 1
 0 <= zeta_3_2_8 <= 1
 0 <= zeta_3_2_10 <= 1
 0 <= zeta_3_2_12 <= 1
 0 <= zeta_3_2_14 <= 1
 0 <= zeta_3_2_16 <= 1
 0 <= zeta_4_1_0 <= 1
 0 <= zeta_4_1_2 <= 1
 0 <= zeta_4_1_4 <= 1
 0 <= zeta_4_1_6 <= 1
 0 <= zeta_4_1_8 <= 1
 0 <= zeta_4_2_0 <= 1
 0 <= zeta_4_2_2 <= 1
 0 <= zeta_4_2_4 <= 1
 0 <= zeta_4_2_6 <= 1
 0 <= zeta_4_2_8 <= 1
 0 <= zeta_4_2_10 <= 1
 0 <= zeta_4_2_12 <= 1
 0 <= zeta_4_2_14 <= 1
 0 <= zeta_4_2_16 <= 1
 0 <= zeta_5_1_0 <= 1
 0 <= zeta_5_1_2 <= 1
 0 <= zeta_5_1_4 <= 1
 0 <= zeta_5_1_6 <= 1
 0 <= zeta_5_1_8 <= 1
 0 <= zeta_5_3_0 <= 1
 0 <= zeta_5_3_2 <= 1
 0 <= zeta_5_3_4 <= 1
 0 <= zeta_5_3_6 <= 1
 0 <= zeta_5_3_8 <= 1
 0 <= zeta_5_3_10 <= 1
 0 <= zeta_5_3_12 <= 1
 0 <= zeta_5_3_14 <= 1
 0 <= zeta_6_3_0 <= 1
 0 <= zeta_6_3_2 <= 1
 0 <= zeta_6_3_4 <= 1
 0 <= zeta_6_3_6 <= 1
 0 <= zeta_6_3_8 <= 1
 0 <= zeta_6_3_10 <= 1
 0 <= zeta_6_3_12 <= 1
 0 <= zeta_6_3_14 <= 1
 0 <= zeta_7_2_0 <= 1
 0 <= zeta_7_2_2 <= 1
 0 <= zeta_7_2_4 <= 1
 0 <= zeta_7_2_6 <= 1
 0 <= zeta_7_2_8 <= 1
 0 <= zeta_7_2_10 <= 1
 0 <= zeta_7_2_12 <= 1
 0 <= zeta_7_2_14 <= 1
 0 <= zeta_7_2_16 <= 1
 0 <= zeta_7_3_0 <= 1
 0 <= zeta_7_3_2 <= 1
 0 <= zeta_7_3_4 <= 1
 0 <= zeta_7_3_6 <= 1
 0 <= zeta_7_3_8 <= 1
 0 <= zeta_7_3_10 <= 1
 0 <= zeta_7_3_12 <= 1
 0 <= zeta_7_3_14 <= 1
 0 <= zeta_8_0_0 <= 1
 0 <= zeta_8_0_2 <= 1
 0 <= zeta_8_0_4 <= 1
 0 <= zeta_8_0_6 <= 1
 0 <= zeta_8_0_8 <= 1
 0 <= zeta_8_0_10 <= 1
 0 <= zeta_8_0_12 <= 1
 0 <= zeta_8_2_0 <= 1
 0 <= zeta_8_2_2 <= 1
 0 <= zeta_8_2_4 <= 1
 0 <= zeta_8_2_6 <= 1
 0 <= zeta_8_2_8 <= 1
 0 <= zeta_8_2_10 <= 1
 0 <= zeta_8_2_12 <= 1
 0 <= zeta_8_2_14 <= 1
 0 <= zeta_8_2_16 <= 1
 0 <= zeta_8_3_0 <= 1
 0 <= zeta_8_3_2 <= 1
 0 <= zeta_8_3_4 <= 1
 0 <= zeta_8_3_6 <= 1
 0 <= zeta_8_3_8 <= 1
 0 <= zeta_8_3_10 <= 1
 0 <= zeta_8_3_12 <= 1
 0 <= zeta_8_3_14 <= 1
 0 <= zeta_9_0_0 <= 1
 0 <= zeta_9_0_2 <= 1
 0 <= zeta_9_0_4 <= 1
 0 <= zeta_9_0_6 <= 1
 0 <= zeta_9_0_8 <= 1
 0 <= zeta_9_0_10 <= 1
 0 <= zeta_9_0_12 <= 1
 0 <= zeta_9_2_0 <= 1
 0 <= zeta_9_2_2 <= 1
 0 <= zeta_9_2_4 <= 1
 0 <= zeta_9_2_6 <= 1
 0 <= zeta_9_2_8 <= 1
 0 <= zeta_9_2_10 <= 1
 0 <= zeta_9_2_12 <= 1
 0 <= zeta_9_2_14 <= 1
 0 <= zeta_9_2_16 <= 1
 0 <= zeta_10_1_0 <= 1
 0 <= zeta_10_1_2 <= 1
 0 <= zeta_10_1_4 <= 1
 0 <= zeta_10_1_6 <= 1
 0 <= zeta_10_1_8 <= 1
 0 <= zeta_10_2_0 <= 1
 0 <= zeta_10_2_2 <= 1
 0 <= zeta_10_2_4 <= 1
 0 <= zeta_10_2_6 <= 1
 0 <= zeta_10_2_8 <= 1
 0 <= zeta_10_2_10 <= 1
 0 <= zeta_10_2_12 <= 1
 0 <= zeta_10_2_14 <= 1
 0 <= zeta_10_2_16 <= 1
 0 <= zeta_11_0_0 <= 1
 0 <= zeta_11_0_2 <= 1
 0 <= zeta_11_0_4 <= 1
 0 <= zeta_11_0_6 <= 1
 0 <= zeta_11_0_8 <= 1
 0 <= zeta_11_0_10 <= 1
 0 <= zeta_11_0_12 <= 1
 0 <= zeta_11_2_0 <= 1
 0 <= zeta_11_2_2 <= 1
 0 <= zeta_11_2_4 <= 1
 0 <= zeta_11_2_6 <= 1
 0 <= zeta_11_2_8 <= 1
 0 <= zeta_11_2_10 <= 1
 0 <= zeta_11_2_12 <= 1
 0 <= zeta_11_2_14 <= 1
 0 <= zeta_11_2_16 <= 1
 0 <= zeta_12_0_0 <= 1
 0 <= zeta_12_0_2 <= 1
 0 <= zeta_12_0_4 <= 1
 0 <= zeta_12_0_6 <= 1
 0 <= zeta_12_0_8 <= 1
 0 <= zeta_12_0_10 <= 1
 0 <= zeta_12_0_12 <= 1
 0 <= zeta_12_3_0 <= 1
 0 <= zeta_12_3_2 <= 1
 0 <= zeta_12_3_4 <= 1
 0 <= zeta_12_3_6 <= 1
 0 <= zeta_12_3_8 <= 1
 0 <= zeta_12_3_10 <= 1
 0 <= zeta_12_3_12 <= 1
 0 <= zeta_12_3_14 <= 1
 0 <= zeta_13_2_0 <= 1
 0 <= zeta_13_2_2 <= 1
 0 <= zeta_13_2_4 <= 1
 0 <= zeta_13_2_6 <= 1
 0 <= zeta_13_2_8 <= 1
 0 <= zeta_13_2_10 <= 1
 0 <= zeta_13_2_12 <= 1
 0 <= zeta_13_2_14 <= 1
 0 <= zeta_13_2_16 <= 1
 0 <= zeta_14_1_0 <= 1
 0 <= zeta_14_1_2 <= 1
 0 <= zeta_14_1_4 <= 1
 0 <= zeta_14_1_6 <= 1
 0 <= zeta_14_1_8 <= 1
 0 <= zeta_14_2_0 <= 1
 0 <= zeta_14_2_2 <= 1
 0 <= zeta_14_2_4 <= 1
 0 <= zeta_14_2_6 <= 1
 0 <= zeta_14_2_8 <= 1
 0 <= zeta_14_2_10 <= 1
 0 <= zeta_14_2_12 <= 1
 0 <= zeta_14_2_14 <= 1
 0 <= zeta_14_2_16 <= 1
 0 <= zeta_14_3_0 <= 1
 0 <= zeta_14_3_2 <= 1
 0 <= zeta_14_3_4 <= 1
 0 <= zeta_14_3_6 <= 1
 0 <= zeta_14_3_8 <= 1
 0 <= zeta_14_3_10 <= 1
 0 <= zeta_14_3_12 <= 1
 0 <= zeta_14_3_14 <= 1
 0 <= zeta_15_0_0 <= 1
 0 <= zeta_15_0_2 <= 1
 0 <= zeta_15_0_4 <= 1
 0 <= zeta_15_0_6 <= 1
 0 <= zeta_15_0_8 <= 1
 0 <= zeta_15_0_10 <= 1
 0 <= zeta_15_0_12 <= 1
 0 <= zeta_15_3_0 <= 1
 0 <= zeta_15_3_2 <= 1
 0 <= zeta_15_3_4 <= 1
 0 <= zeta_15_3_6 <= 1
 0 <= zeta_15_3_8 <= 1
 0 <= zeta_15_3_10 <= 1
 0 <= zeta_15_3_12 <= 1
 0 <= zeta_15_3_14 <= 1
 0 <= zeta_16_2_0 <= 1
 0 <= zeta_16_2_2 <= 1
 0 <= zeta_16_2_4 <= 1
 0 <= zeta_16_2_6 <= 1
 0 <= zeta_16_2_8 <= 1
 0 <= zeta_16_2_10 <= 1
 0 <= zeta_16_2_12 <= 1
 0 <= zeta_16_2_14 <= 1
 0 <= zeta_16_2_16 <= 1
 0 <= zeta_17_0_0 <= 1
 0 <= zeta_17_0_2 <= 1
 0 <= zeta_17_0_4 <= 1
 0 <= zeta_17_0_6 <= 1
 0 <= zeta_17_0_8 <= 1
 0 <= zeta_17_0_10 <= 1
 0 <= zeta_17_0_12 <= 1
 0 <= zeta_18_3_0 <= 1
 0 <= zeta_18_3_2 <= 1
 0 <= zeta_18_3_4 <= 1
 0 <= zeta_18_3_6 <= 1
 0 <= zeta_18_3_8 <= 1
 0 <= zeta_18_3_10 <= 1
 0 <= zeta_18_3_12 <= 1
 0 <= zeta_18_3_14 <= 1
 0 <= zeta_19_2_0 <= 1
 0 <= zeta_19_2_2 <= 1
 0 <= zeta_19_2_4 <= 1
 0 <= zeta_19_2_6 <= 1
 0 <= zeta_19_2_8 <= 1
 0 <= zeta_19_2_10 <= 1
 0 <= zeta_19_2_12 <= 1
 0 <= zeta_19_2_14 <= 1
 0 <= zeta_19_2_16 <= 1
 0 <= zeta_19_3_0 <= 1
 0 <= zeta_19_3_2 <= 1
 0 <= zeta_19_3_4 <= 1
 0 <= zeta_19_3_6 <= 1
 0 <= zeta_19_3_8 <= 1
 0 <= zeta_19_3_10 <= 1
 0 <= zeta_19_3_12 <= 1
 0 <= zeta_19_3_14 <= 1
 0 <= zeta_20_1_0 <= 1
 0 <= zeta_20_1_2 <= 1
 0 <= zeta_20_1_4 <= 1
 0 <= zeta_20_1_6 <= 1
 0 <= zeta_20_1_8 <= 1
 0 <= zeta_21_2_0 <= 1
 0 <= zeta_21_2_2 <= 1
 0 <= zeta_21_2_4 <= 1
 0 <= zeta_21_2_6 <= 1
 0 <= zeta_21_2_8 <= 1
 0 <= zeta_21_2_10 <= 1
 0 <= zeta_21_2_12 <= 1
 0 <= zeta_21_2_14 <= 1
 0 <= zeta_21_2_16 <= 1
 0 <= zeta_22_0_0 <= 1
 0 <= zeta_22_0_2 <= 1
 0 <= zeta_22_0_4 <= 1
 0 <= zeta_22_0_6 <= 1
 0 <= zeta_22_0_8 <= 1
 0 <= zeta_22_0_10 <= 1
 0 <= zeta_22_0_12 <= 1
 0 <= zeta_23_3_0 <= 1
 0 <= zeta_23_3_2 <= 1
 0 <= zeta_23_3_4 <= 1
 0 <= zeta_23_3_6 <= 1
 0 <= zeta_23_3_8 <= 1
 0 <= zeta_23_3_10 <= 1
 0 <= zeta_23_3_12 <= 1
 0 <= zeta_23_3_14 <= 1
 0 <= zeta_24_2_0 <= 1
 0 <= zeta_24_2_2 <= 1
 0 <= zeta_24_2_4 <= 1
 0 <= zeta_24_2_6 <= 1
 0 <= zeta_24_2_8 <= 1
 0 <= zeta_24_2_10 <= 1
 0 <= zeta_24_2_12 <= 1
 0 <= zeta_24_2_14 <= 1
 0 <= zeta_24_2_16 <= 1
 0 <= zeta_24_3_0 <= 1
 0 <= zeta_24_3_2 <= 1
 0 <= zeta_24_3_4 <= 1
 0 <= zeta_24_3_6 <= 1
 0 <= zeta_24_3_8 <= 1
 0 <= zeta_24_3_10 <= 1
 0 <= zeta_24_3_12 <= 1
 0 <= zeta_24_3_14 <= 1
 0 <= zeta_25_1_0 <= 1
 0 <= zeta_25_1_2 <= 1
 0 <= zeta_25_1_4 <= 1
 0 <= zeta_25_1_6 <= 1
 0 <= zeta_25_1_8 <= 1
 0 <= zeta_26_2_0 <= 1
 0 <= zeta_26_2_2 <= 1
 0 <= zeta_26_2_4 <= 1
 0 <= zeta_26_2_6 <= 1
 0 <= zeta_26_2_8 <= 1
 0 <= zeta_26_2_10 <= 1
 0 <= zeta_26_2_12 <= 1
 0 <= zeta_26_2_14 <= 1
 0 <= zeta_26_2_16 <= 1
 0 <= zeta_27_0_0 <= 1
 0 <= zeta_27_0_2 <= 1
 0 <= zeta_27_0_4 <= 1
 0 <= zeta_27_0_6 <= 1
 0 <= zeta_27_0_8 <= 1
 0 <= zeta_27_0_10 <= 1
 0 <= zeta_27_0_12 <= 1
 0 <= zeta_28_3_0 <= 1
 0 <= zeta_28_3_2 <= 1
 0 <= zeta_28_3_4 <= 1
 0 <= zeta_28_3_6 <= 1
 0 <= zeta_28_3_8 <= 1
 0 <= zeta_28_3_10 <= 1
 0 <= zeta_28_3_12 <= 1
 0 <= zeta_28_3_14 <= 1
 0 <= zeta_29_0_0 <= 1
 0 <= zeta_29_0_2 <= 1
 0 <= zeta_29_0_4 <= 1
 0 <= zeta_29_0_6 <= 1
 0 <= zeta_29_0_8 <= 1
 0 <= zeta_29_0_10 <= 1
 0 <= zeta_29_0_12 <= 1
 0 <= zeta_29_1_0 <= 1
 0 <= zeta_29_1_2 <= 1
 0 <= zeta_29_1_4 <= 1
 0 <= zeta_29_1_6 <= 1
 0 <= zeta_29_1_8 <= 1
 0 <= zeta_30_2_0 <= 1
 0 <= zeta_30_2_2 <= 1
 0 <= zeta_30_2_4 <= 1
 0 <= zeta_30_2_6 <= 1
 0 <= zeta_30_2_8 <= 1
 0 <= zeta_30_2_10 <= 1
 0 <= zeta_30_2_12 <= 1
 0 <= zeta_30_2_14 <= 1
 0 <= zeta_30_2_16 <= 1
 0 <= zeta_30_3_0 <= 1
 0 <= zeta_30_3_2 <= 1
 0 <= zeta_30_3_4 <= 1
 0 <= zeta_30_3_6 <= 1
 0 <= zeta_30_3_8 <= 1
 0 <= zeta_30_3_10 <= 1
 0 <= zeta_30_3_12 <= 1
 0 <= zeta_30_3_14 <= 1
 0 <= zeta_31_3_0 <= 1
 0 <= zeta_31_3_2 <= 1
 0 <= zeta_31_3_4 <= 1
 0 <= zeta_31_3_6 <= 1
 0 <= zeta_31_3_8 <= 1
 0 <= zeta_31_3_10 <= 1
 0 <= zeta_31_3_12 <= 1
 0 <= zeta_31_3_14 <= 1
 0 <= dummy <= 1
Binaries
 de_t0_12_26  de_t0_11_24  de_t0_19_21  de_t0_12_23  de_t0_8_10  x6  x7  x8 
 x9  x10  x11  x12  x13  x14  x15  x16  x17  x18  x19  x20  x21  x22  x23  x24 
 x25  x26  x27  x28  x29  x30  x31  alpha_0_0  alpha_0_2  alpha_0_4  alpha_0_6 
 alpha_0_8  alpha_0_10  alpha_0_12  alpha_1_0  alpha_1_2  alpha_1_4  alpha_1_6 
 alpha_1_8  alpha_2_0  alpha_2_2  alpha_2_4  alpha_2_6  alpha_2_8  alpha_2_10 
 alpha_2_12  alpha_2_14  alpha_2_16  alpha_3_0  alpha_3_2  alpha_3_4 
 alpha_3_6  alpha_3_8  alpha_3_10  alpha_3_12  alpha_3_14  zeta_0_0_0 
 zeta_0_0_2  zeta_0_0_4  zeta_0_0_6  zeta_0_0_8  zeta_0_0_10  zeta_0_0_12 
 zeta_1_0_0  zeta_1_0_2  zeta_1_0_4  zeta_1_0_6  zeta_1_0_8  zeta_1_0_10 
 zeta_1_0_12  zeta_1_1_0  zeta_1_1_2  zeta_1_1_4  zeta_1_1_6  zeta_1_1_8 
 zeta_1_2_0  zeta_1_2_2  zeta_1_2_4  zeta_1_2_6  zeta_1_2_8  zeta_1_2_10 
 zeta_1_2_12  zeta_1_2_14  zeta_1_2_16  zeta_2_0_0  zeta_2_0_2  zeta_2_0_4 
 zeta_2_0_6  zeta_2_0_8  zeta_2_0_10  zeta_2_0_12  zeta_2_2_0  zeta_2_2_2 
 zeta_2_2_4  zeta_2_2_6  zeta_2_2_8  zeta_2_2_10  zeta_2_2_12  zeta_2_2_14 
 zeta_2_2_16  zeta_2_3_0  zeta_2_3_2  zeta_2_3_4  zeta_2_3_6  zeta_2_3_8 
 zeta_2_3_10  zeta_2_3_12  zeta_2_3_14  zeta_3_0_0  zeta_3_0_2  zeta_3_0_4 
 zeta_3_0_6  zeta_3_0_8  zeta_3_0_10  zeta_3_0_12  zeta_3_2_0  zeta_3_2_2 
 zeta_3_2_4  zeta_3_2_6  zeta_3_2_8  zeta_3_2_10  zeta_3_2_12  zeta_3_2_14 
 zeta_3_2_16  zeta_4_1_0  zeta_4_1_2  zeta_4_1_4  zeta_4_1_6  zeta_4_1_8 
 zeta_4_2_0  zeta_4_2_2  zeta_4_2_4  zeta_4_2_6  zeta_4_2_8  zeta_4_2_10 
 zeta_4_2_12  zeta_4_2_14  zeta_4_2_16  zeta_5_1_0  zeta_5_1_2  zeta_5_1_4 
 zeta_5_1_6  zeta_5_1_8  zeta_5_3_0  zeta_5_3_2  zeta_5_3_4  zeta_5_3_6 
 zeta_5_3_8  zeta_5_3_10  zeta_5_3_12  zeta_5_3_14  zeta_6_3_0  zeta_6_3_2 
 zeta_6_3_4  zeta_6_3_6  zeta_6_3_8  zeta_6_3_10  zeta_6_3_12  zeta_6_3_14 
 zeta_7_2_0  zeta_7_2_2  zeta_7_2_4  zeta_7_2_6  zeta_7_2_8  zeta_7_2_10 
 zeta_7_2_12  zeta_7_2_14  zeta_7_2_16  zeta_7_3_0  zeta_7_3_2  zeta_7_3_4 
 zeta_7_3_6  zeta_7_3_8  zeta_7_3_10  zeta_7_3_12  zeta_7_3_14  zeta_8_0_0 
 zeta_8_0_2  zeta_8_0_4  zeta_8_0_6  zeta_8_0_8  zeta_8_0_10  zeta_8_0_12 
 zeta_8_2_0  zeta_8_2_2  zeta_8_2_4  zeta_8_2_6  zeta_8_2_8  zeta_8_2_10 
 zeta_8_2_12  zeta_8_2_14  zeta_8_2_16  zeta_8_3_0  zeta_8_3_2  zeta_8_3_4 
 zeta_8_3_6  zeta_8_3_8  zeta_8_3_10  zeta_8_3_12  zeta_8_3_14  zeta_9_0_0 
 zeta_9_0_2  zeta_9_0_4  zeta_9_0_6  zeta_9_0_8  zeta_9_0_10  zeta_9_0_12 
 zeta_9_2_0  zeta_9_2_2  zeta_9_2_4  zeta_9_2_6  zeta_9_2_8  zeta_9_2_10 
 zeta_9_2_12  zeta_9_2_14  zeta_9_2_16  zeta_10_1_0  zeta_10_1_2  zeta_10_1_4 
 zeta_10_1_6  zeta_10_1_8  zeta_10_2_0  zeta_10_2_2  zeta_10_2_4  zeta_10_2_6 
 zeta_10_2_8  zeta_10_2_10  zeta_10_2_12  zeta_10_2_14  zeta_10_2_16 
 zeta_11_0_0  zeta_11_0_2  zeta_11_0_4  zeta_11_0_6  zeta_11_0_8  zeta_11_0_10 
 zeta_11_0_12  zeta_11_2_0  zeta_11_2_2  zeta_11_2_4  zeta_11_2_6  zeta_11_2_8 
 zeta_11_2_10  zeta_11_2_12  zeta_11_2_14  zeta_11_2_16  zeta_12_0_0 
 zeta_12_0_2  zeta_12_0_4  zeta_12_0_6  zeta_12_0_8  zeta_12_0_10 
 zeta_12_0_12  zeta_12_3_0  zeta_12_3_2  zeta_12_3_4  zeta_12_3_6  zeta_12_3_8 
 zeta_12_3_10  zeta_12_3_12  zeta_12_3_14  zeta_13_2_0  zeta_13_2_2 
 zeta_13_2_4  zeta_13_2_6  zeta_13_2_8  zeta_13_2_10  zeta_13_2_12 
 zeta_13_2_14  zeta_13_2_16  zeta_14_1_0  zeta_14_1_2  zeta_14_1_4 
 zeta_14_1_6  zeta_14_1_8  zeta_14_2_0  zeta_14_2_2  zeta_14_2_4  zeta_14_2_6 
 zeta_14_2_8  zeta_14_2_10  zeta_14_2_12  zeta_14_2_14  zeta_14_2_16 
 zeta_14_3_0  zeta_14_3_2  zeta_14_3_4  zeta_14_3_6  zeta_14_3_8  zeta_14_3_10 
 zeta_14_3_12  zeta_14_3_14  zeta_15_0_0  zeta_15_0_2  zeta_15_0_4 
 zeta_15_0_6  zeta_15_0_8  zeta_15_0_10  zeta_15_0_12  zeta_15_3_0 
 zeta_15_3_2  zeta_15_3_4  zeta_15_3_6  zeta_15_3_8  zeta_15_3_10 
 zeta_15_3_12  zeta_15_3_14  zeta_16_2_0  zeta_16_2_2  zeta_16_2_4 
 zeta_16_2_6  zeta_16_2_8  zeta_16_2_10  zeta_16_2_12  zeta_16_2_14 
 zeta_16_2_16  zeta_17_0_0  zeta_17_0_2  zeta_17_0_4  zeta_17_0_6  zeta_17_0_8 
 zeta_17_0_10  zeta_17_0_12  zeta_18_3_0  zeta_18_3_2  zeta_18_3_4 
 zeta_18_3_6  zeta_18_3_8  zeta_18_3_10  zeta_18_3_12  zeta_18_3_14 
 zeta_19_2_0  zeta_19_2_2  zeta_19_2_4  zeta_19_2_6  zeta_19_2_8  zeta_19_2_10 
 zeta_19_2_12  zeta_19_2_14  zeta_19_2_16  zeta_19_3_0  zeta_19_3_2 
 zeta_19_3_4  zeta_19_3_6  zeta_19_3_8  zeta_19_3_10  zeta_19_3_12 
 zeta_19_3_14  zeta_20_1_0  zeta_20_1_2  zeta_20_1_4  zeta_20_1_6  zeta_20_1_8 
 zeta_21_2_0  zeta_21_2_2  zeta_21_2_4  zeta_21_2_6  zeta_21_2_8  zeta_21_2_10 
 zeta_21_2_12  zeta_21_2_14  zeta_21_2_16  zeta_22_0_0  zeta_22_0_2 
 zeta_22_0_4  zeta_22_0_6  zeta_22_0_8  zeta_22_0_10  zeta_22_0_12 
 zeta_23_3_0  zeta_23_3_2  zeta_23_3_4  zeta_23_3_6  zeta_23_3_8  zeta_23_3_10 
 zeta_23_3_12  zeta_23_3_14  zeta_24_2_0  zeta_24_2_2  zeta_24_2_4 
 zeta_24_2_6  zeta_24_2_8  zeta_24_2_10  zeta_24_2_12  zeta_24_2_14 
 zeta_24_2_16  zeta_24_3_0  zeta_24_3_2  zeta_24_3_4  zeta_24_3_6  zeta_24_3_8 
 zeta_24_3_10  zeta_24_3_12  zeta_24_3_14  zeta_25_1_0  zeta_25_1_2 
 zeta_25_1_4  zeta_25_1_6  zeta_25_1_8  zeta_26_2_0  zeta_26_2_2  zeta_26_2_4 
 zeta_26_2_6  zeta_26_2_8  zeta_26_2_10  zeta_26_2_12  zeta_26_2_14 
 zeta_26_2_16  zeta_27_0_0  zeta_27_0_2  zeta_27_0_4  zeta_27_0_6  zeta_27_0_8 
 zeta_27_0_10  zeta_27_0_12  zeta_28_3_0  zeta_28_3_2  zeta_28_3_4 
 zeta_28_3_6  zeta_28_3_8  zeta_28_3_10  zeta_28_3_12  zeta_28_3_14 
 zeta_29_0_0  zeta_29_0_2  zeta_29_0_4  zeta_29_0_6  zeta_29_0_8  zeta_29_0_10 
 zeta_29_0_12  zeta_29_1_0  zeta_29_1_2  zeta_29_1_4  zeta_29_1_6  zeta_29_1_8 
 zeta_30_2_0  zeta_30_2_2  zeta_30_2_4  zeta_30_2_6  zeta_30_2_8  zeta_30_2_10 
 zeta_30_2_12  zeta_30_2_14  zeta_30_2_16  zeta_30_3_0  zeta_30_3_2 
 zeta_30_3_4  zeta_30_3_6  zeta_30_3_8  zeta_30_3_10  zeta_30_3_12 
 zeta_30_3_14  zeta_31_3_0  zeta_31_3_2  zeta_31_3_4  zeta_31_3_6  zeta_31_3_8 
 zeta_31_3_10  zeta_31_3_12  zeta_31_3_14  dummy 
End
