package main

var blosum62 map[string]float64 = map[string]float64{
	"AA": 4.0,
	"AR": -1.0,
	"AN": -2.0,
	"AD": -2.0,
	"AC": 0.0,
	"AQ": -1.0,
	"AE": -1.0,
	"AG": 0.0,
	"AH": -2.0,
	"AI": -1.0,
	"AL": -1.0,
	"AK": -1.0,
	"AM": -1.0,
	"AF": -2.0,
	"AP": -1.0,
	"AS": 1.0,
	"AT": 0.0,
	"AW": -3.0,
	"AY": -2.0,
	"AV": 0.0,
	"AB": -2.0,
	"AJ": -1.0,
	"AZ": -1.0,
	"AX": -1.0,
	"A*": -4.0,
	"RA": -1.0,
	"RR": 5.0,
	"RN": 0.0,
	"RD": -2.0,
	"RC": -3.0,
	"RQ": 1.0,
	"RE": 0.0,
	"RG": -2.0,
	"RH": 0.0,
	"RI": -3.0,
	"RL": -2.0,
	"RK": 2.0,
	"RM": -1.0,
	"RF": -3.0,
	"RP": -2.0,
	"RS": -1.0,
	"RT": -1.0,
	"RW": -3.0,
	"RY": -2.0,
	"RV": -3.0,
	"RB": -1.0,
	"RJ": -2.0,
	"RZ": 0.0,
	"RX": -1.0,
	"R*": -4.0,
	"NA": -2.0,
	"NR": 0.0,
	"NN": 6.0,
	"ND": 1.0,
	"NC": -3.0,
	"NQ": 0.0,
	"NE": 0.0,
	"NG": 0.0,
	"NH": 1.0,
	"NI": -3.0,
	"NL": -3.0,
	"NK": 0.0,
	"NM": -2.0,
	"NF": -3.0,
	"NP": -2.0,
	"NS": 1.0,
	"NT": 0.0,
	"NW": -4.0,
	"NY": -2.0,
	"NV": -3.0,
	"NB": 4.0,
	"NJ": -3.0,
	"NZ": 0.0,
	"NX": -1.0,
	"N*": -4.0,
	"DA": -2.0,
	"DR": -2.0,
	"DN": 1.0,
	"DD": 6.0,
	"DC": -3.0,
	"DQ": 0.0,
	"DE": 2.0,
	"DG": -1.0,
	"DH": -1.0,
	"DI": -3.0,
	"DL": -4.0,
	"DK": -1.0,
	"DM": -3.0,
	"DF": -3.0,
	"DP": -1.0,
	"DS": 0.0,
	"DT": -1.0,
	"DW": -4.0,
	"DY": -3.0,
	"DV": -3.0,
	"DB": 4.0,
	"DJ": -3.0,
	"DZ": 1.0,
	"DX": -1.0,
	"D*": -4.0,
	"CA": 0.0,
	"CR": -3.0,
	"CN": -3.0,
	"CD": -3.0,
	"CC": 9.0,
	"CQ": -3.0,
	"CE": -4.0,
	"CG": -3.0,
	"CH": -3.0,
	"CI": -1.0,
	"CL": -1.0,
	"CK": -3.0,
	"CM": -1.0,
	"CF": -2.0,
	"CP": -3.0,
	"CS": -1.0,
	"CT": -1.0,
	"CW": -2.0,
	"CY": -2.0,
	"CV": -1.0,
	"CB": -3.0,
	"CJ": -1.0,
	"CZ": -3.0,
	"CX": -1.0,
	"C*": -4.0,
	"QA": -1.0,
	"QR": 1.0,
	"QN": 0.0,
	"QD": 0.0,
	"QC": -3.0,
	"QQ": 5.0,
	"QE": 2.0,
	"QG": -2.0,
	"QH": 0.0,
	"QI": -3.0,
	"QL": -2.0,
	"QK": 1.0,
	"QM": 0.0,
	"QF": -3.0,
	"QP": -1.0,
	"QS": 0.0,
	"QT": -1.0,
	"QW": -2.0,
	"QY": -1.0,
	"QV": -2.0,
	"QB": 0.0,
	"QJ": -2.0,
	"QZ": 4.0,
	"QX": -1.0,
	"Q*": -4.0,
	"EA": -1.0,
	"ER": 0.0,
	"EN": 0.0,
	"ED": 2.0,
	"EC": -4.0,
	"EQ": 2.0,
	"EE": 5.0,
	"EG": -2.0,
	"EH": 0.0,
	"EI": -3.0,
	"EL": -3.0,
	"EK": 1.0,
	"EM": -2.0,
	"EF": -3.0,
	"EP": -1.0,
	"ES": 0.0,
	"ET": -1.0,
	"EW": -3.0,
	"EY": -2.0,
	"EV": -2.0,
	"EB": 1.0,
	"EJ": -3.0,
	"EZ": 4.0,
	"EX": -1.0,
	"E*": -4.0,
	"GA": 0.0,
	"GR": -2.0,
	"GN": 0.0,
	"GD": -1.0,
	"GC": -3.0,
	"GQ": -2.0,
	"GE": -2.0,
	"GG": 6.0,
	"GH": -2.0,
	"GI": -4.0,
	"GL": -4.0,
	"GK": -2.0,
	"GM": -3.0,
	"GF": -3.0,
	"GP": -2.0,
	"GS": 0.0,
	"GT": -2.0,
	"GW": -2.0,
	"GY": -3.0,
	"GV": -3.0,
	"GB": -1.0,
	"GJ": -4.0,
	"GZ": -2.0,
	"GX": -1.0,
	"G*": -4.0,
	"HA": -2.0,
	"HR": 0.0,
	"HN": 1.0,
	"HD": -1.0,
	"HC": -3.0,
	"HQ": 0.0,
	"HE": 0.0,
	"HG": -2.0,
	"HH": 8.0,
	"HI": -3.0,
	"HL": -3.0,
	"HK": -1.0,
	"HM": -2.0,
	"HF": -1.0,
	"HP": -2.0,
	"HS": -1.0,
	"HT": -2.0,
	"HW": -2.0,
	"HY": 2.0,
	"HV": -3.0,
	"HB": 0.0,
	"HJ": -3.0,
	"HZ": 0.0,
	"HX": -1.0,
	"H*": -4.0,
	"IA": -1.0,
	"IR": -3.0,
	"IN": -3.0,
	"ID": -3.0,
	"IC": -1.0,
	"IQ": -3.0,
	"IE": -3.0,
	"IG": -4.0,
	"IH": -3.0,
	"II": 4.0,
	"IL": 2.0,
	"IK": -3.0,
	"IM": 1.0,
	"IF": 0.0,
	"IP": -3.0,
	"IS": -2.0,
	"IT": -1.0,
	"IW": -3.0,
	"IY": -1.0,
	"IV": 3.0,
	"IB": -3.0,
	"IJ": 3.0,
	"IZ": -3.0,
	"IX": -1.0,
	"I*": -4.0,
	"LA": -1.0,
	"LR": -2.0,
	"LN": -3.0,
	"LD": -4.0,
	"LC": -1.0,
	"LQ": -2.0,
	"LE": -3.0,
	"LG": -4.0,
	"LH": -3.0,
	"LI": 2.0,
	"LL": 4.0,
	"LK": -2.0,
	"LM": 2.0,
	"LF": 0.0,
	"LP": -3.0,
	"LS": -2.0,
	"LT": -1.0,
	"LW": -2.0,
	"LY": -1.0,
	"LV": 1.0,
	"LB": -4.0,
	"LJ": 3.0,
	"LZ": -3.0,
	"LX": -1.0,
	"L*": -4.0,
	"KA": -1.0,
	"KR": 2.0,
	"KN": 0.0,
	"KD": -1.0,
	"KC": -3.0,
	"KQ": 1.0,
	"KE": 1.0,
	"KG": -2.0,
	"KH": -1.0,
	"KI": -3.0,
	"KL": -2.0,
	"KK": 5.0,
	"KM": -1.0,
	"KF": -3.0,
	"KP": -1.0,
	"KS": 0.0,
	"KT": -1.0,
	"KW": -3.0,
	"KY": -2.0,
	"KV": -2.0,
	"KB": 0.0,
	"KJ": -3.0,
	"KZ": 1.0,
	"KX": -1.0,
	"K*": -4.0,
	"MA": -1.0,
	"MR": -1.0,
	"MN": -2.0,
	"MD": -3.0,
	"MC": -1.0,
	"MQ": 0.0,
	"ME": -2.0,
	"MG": -3.0,
	"MH": -2.0,
	"MI": 1.0,
	"ML": 2.0,
	"MK": -1.0,
	"MM": 5.0,
	"MF": 0.0,
	"MP": -2.0,
	"MS": -1.0,
	"MT": -1.0,
	"MW": -1.0,
	"MY": -1.0,
	"MV": 1.0,
	"MB": -3.0,
	"MJ": 2.0,
	"MZ": -1.0,
	"MX": -1.0,
	"M*": -4.0,
	"FA": -2.0,
	"FR": -3.0,
	"FN": -3.0,
	"FD": -3.0,
	"FC": -2.0,
	"FQ": -3.0,
	"FE": -3.0,
	"FG": -3.0,
	"FH": -1.0,
	"FI": 0.0,
	"FL": 0.0,
	"FK": -3.0,
	"FM": 0.0,
	"FF": 6.0,
	"FP": -4.0,
	"FS": -2.0,
	"FT": -2.0,
	"FW": 1.0,
	"FY": 3.0,
	"FV": -1.0,
	"FB": -3.0,
	"FJ": 0.0,
	"FZ": -3.0,
	"FX": -1.0,
	"F*": -4.0,
	"PA": -1.0,
	"PR": -2.0,
	"PN": -2.0,
	"PD": -1.0,
	"PC": -3.0,
	"PQ": -1.0,
	"PE": -1.0,
	"PG": -2.0,
	"PH": -2.0,
	"PI": -3.0,
	"PL": -3.0,
	"PK": -1.0,
	"PM": -2.0,
	"PF": -4.0,
	"PP": 7.0,
	"PS": -1.0,
	"PT": -1.0,
	"PW": -4.0,
	"PY": -3.0,
	"PV": -2.0,
	"PB": -2.0,
	"PJ": -3.0,
	"PZ": -1.0,
	"PX": -1.0,
	"P*": -4.0,
	"SA": 1.0,
	"SR": -1.0,
	"SN": 1.0,
	"SD": 0.0,
	"SC": -1.0,
	"SQ": 0.0,
	"SE": 0.0,
	"SG": 0.0,
	"SH": -1.0,
	"SI": -2.0,
	"SL": -2.0,
	"SK": 0.0,
	"SM": -1.0,
	"SF": -2.0,
	"SP": -1.0,
	"SS": 4.0,
	"ST": 1.0,
	"SW": -3.0,
	"SY": -2.0,
	"SV": -2.0,
	"SB": 0.0,
	"SJ": -2.0,
	"SZ": 0.0,
	"SX": -1.0,
	"S*": -4.0,
	"TA": 0.0,
	"TR": -1.0,
	"TN": 0.0,
	"TD": -1.0,
	"TC": -1.0,
	"TQ": -1.0,
	"TE": -1.0,
	"TG": -2.0,
	"TH": -2.0,
	"TI": -1.0,
	"TL": -1.0,
	"TK": -1.0,
	"TM": -1.0,
	"TF": -2.0,
	"TP": -1.0,
	"TS": 1.0,
	"TT": 5.0,
	"TW": -2.0,
	"TY": -2.0,
	"TV": 0.0,
	"TB": -1.0,
	"TJ": -1.0,
	"TZ": -1.0,
	"TX": -1.0,
	"T*": -4.0,
	"WA": -3.0,
	"WR": -3.0,
	"WN": -4.0,
	"WD": -4.0,
	"WC": -2.0,
	"WQ": -2.0,
	"WE": -3.0,
	"WG": -2.0,
	"WH": -2.0,
	"WI": -3.0,
	"WL": -2.0,
	"WK": -3.0,
	"WM": -1.0,
	"WF": 1.0,
	"WP": -4.0,
	"WS": -3.0,
	"WT": -2.0,
	"WW": 11.0,
	"WY": 2.0,
	"WV": -3.0,
	"WB": -4.0,
	"WJ": -2.0,
	"WZ": -2.0,
	"WX": -1.0,
	"W*": -4.0,
	"YA": -2.0,
	"YR": -2.0,
	"YN": -2.0,
	"YD": -3.0,
	"YC": -2.0,
	"YQ": -1.0,
	"YE": -2.0,
	"YG": -3.0,
	"YH": 2.0,
	"YI": -1.0,
	"YL": -1.0,
	"YK": -2.0,
	"YM": -1.0,
	"YF": 3.0,
	"YP": -3.0,
	"YS": -2.0,
	"YT": -2.0,
	"YW": 2.0,
	"YY": 7.0,
	"YV": -1.0,
	"YB": -3.0,
	"YJ": -1.0,
	"YZ": -2.0,
	"YX": -1.0,
	"Y*": -4.0,
	"VA": 0.0,
	"VR": -3.0,
	"VN": -3.0,
	"VD": -3.0,
	"VC": -1.0,
	"VQ": -2.0,
	"VE": -2.0,
	"VG": -3.0,
	"VH": -3.0,
	"VI": 3.0,
	"VL": 1.0,
	"VK": -2.0,
	"VM": 1.0,
	"VF": -1.0,
	"VP": -2.0,
	"VS": -2.0,
	"VT": 0.0,
	"VW": -3.0,
	"VY": -1.0,
	"VV": 4.0,
	"VB": -3.0,
	"VJ": 2.0,
	"VZ": -2.0,
	"VX": -1.0,
	"V*": -4.0,
	"BA": -2.0,
	"BR": -1.0,
	"BN": 4.0,
	"BD": 4.0,
	"BC": -3.0,
	"BQ": 0.0,
	"BE": 1.0,
	"BG": -1.0,
	"BH": 0.0,
	"BI": -3.0,
	"BL": -4.0,
	"BK": 0.0,
	"BM": -3.0,
	"BF": -3.0,
	"BP": -2.0,
	"BS": 0.0,
	"BT": -1.0,
	"BW": -4.0,
	"BY": -3.0,
	"BV": -3.0,
	"BB": 4.0,
	"BJ": -3.0,
	"BZ": 0.0,
	"BX": -1.0,
	"B*": -4.0,
	"JA": -1.0,
	"JR": -2.0,
	"JN": -3.0,
	"JD": -3.0,
	"JC": -1.0,
	"JQ": -2.0,
	"JE": -3.0,
	"JG": -4.0,
	"JH": -3.0,
	"JI": 3.0,
	"JL": 3.0,
	"JK": -3.0,
	"JM": 2.0,
	"JF": 0.0,
	"JP": -3.0,
	"JS": -2.0,
	"JT": -1.0,
	"JW": -2.0,
	"JY": -1.0,
	"JV": 2.0,
	"JB": -3.0,
	"JJ": 3.0,
	"JZ": -3.0,
	"JX": -1.0,
	"J*": -4.0,
	"ZA": -1.0,
	"ZR": 0.0,
	"ZN": 0.0,
	"ZD": 1.0,
	"ZC": -3.0,
	"ZQ": 4.0,
	"ZE": 4.0,
	"ZG": -2.0,
	"ZH": 0.0,
	"ZI": -3.0,
	"ZL": -3.0,
	"ZK": 1.0,
	"ZM": -1.0,
	"ZF": -3.0,
	"ZP": -1.0,
	"ZS": 0.0,
	"ZT": -1.0,
	"ZW": -2.0,
	"ZY": -2.0,
	"ZV": -2.0,
	"ZB": 0.0,
	"ZJ": -3.0,
	"ZZ": 4.0,
	"ZX": -1.0,
	"Z*": -4.0,
	"XA": -1.0,
	"XR": -1.0,
	"XN": -1.0,
	"XD": -1.0,
	"XC": -1.0,
	"XQ": -1.0,
	"XE": -1.0,
	"XG": -1.0,
	"XH": -1.0,
	"XI": -1.0,
	"XL": -1.0,
	"XK": -1.0,
	"XM": -1.0,
	"XF": -1.0,
	"XP": -1.0,
	"XS": -1.0,
	"XT": -1.0,
	"XW": -1.0,
	"XY": -1.0,
	"XV": -1.0,
	"XB": -1.0,
	"XJ": -1.0,
	"XZ": -1.0,
	"XX": -1.0,
	"X*": -4.0,
	"*A": -4.0,
	"*R": -4.0,
	"*N": -4.0,
	"*D": -4.0,
	"*C": -4.0,
	"*Q": -4.0,
	"*E": -4.0,
	"*G": -4.0,
	"*H": -4.0,
	"*I": -4.0,
	"*L": -4.0,
	"*K": -4.0,
	"*M": -4.0,
	"*F": -4.0,
	"*P": -4.0,
	"*S": -4.0,
	"*T": -4.0,
	"*W": -4.0,
	"*Y": -4.0,
	"*V": -4.0,
	"*B": -4.0,
	"*J": -4.0,
	"*Z": -4.0,
	"*X": -4.0,
	"**": 1.0,
}