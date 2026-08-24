import cython
import numpy as np
cimport numpy as np

ctypedef np.double_t DTYPEDBL_t

if cython.compiled:
    print("Yep, I'm compiled.")
else:
    print("Just a lowly interpreted script.")

####
#Definition of Functor class:
####

cdef class stepClass:

    cdef public params

    def __init__(self, np.ndarray[DTYPEDBL_t, ndim=1] newParams):
        cdef int i
        cdef int n = len(newParams)
        self.params = np.zeros([n], dtype=np.double)
        for i in range(n):
            self.params[i] = newParams[i]


####
#Definition of reaction flux function:
####

    def calcFlux(self, float t, np.ndarray[DTYPEDBL_t, ndim=1] y):
        return self.calcFlux_c(t, y) 

    cdef np.ndarray[DTYPEDBL_t, ndim=1] calcFlux_c(self, float t, np.ndarray[DTYPEDBL_t, ndim=1] y):



        # (M_ACP_c) M_ACP_c
        cdef double M_ACP_c = y[0]

        # (M_ACP_R_c) M_ACP_R_c
        cdef double M_ACP_R_c = y[1]

        # (M_apoACP_c) M_apoACP_c
        cdef double M_apoACP_c = y[2]

        # (M_trdrd_c) M_trdrd_c
        cdef double M_trdrd_c = y[3]

        # (M_trdox_c) M_trdox_c
        cdef double M_trdox_c = y[4]

        # (M_dhlpl_PdhC_c) M_dhlpl_PdhC_c
        cdef double M_dhlpl_PdhC_c = y[5]

        # (M_acdhlpl_PdhC_c) M_acdhlpl_PdhC_c
        cdef double M_acdhlpl_PdhC_c = y[6]

        # (M_lpl_PdhC_c) M_lpl_PdhC_c
        cdef double M_lpl_PdhC_c = y[7]

        # (ptsi_P) ptsi_P
        cdef double ptsi_P = y[8]

        # (ptsi) ptsi
        cdef double ptsi = y[9]

        # (ptsh_P) ptsh_P
        cdef double ptsh_P = y[10]

        # (ptsh) ptsh
        cdef double ptsh = y[11]

        # (crr_P) crr_P
        cdef double crr_P = y[12]

        # (crr) crr
        cdef double crr = y[13]

        # (ptsg_P) ptsg_P
        cdef double ptsg_P = y[14]

        # (ptsg) ptsg
        cdef double ptsg = y[15]

        # (M_g6p_c) M_g6p_c
        cdef double M_g6p_c = y[16]

        # (M_f6p_c) M_f6p_c
        cdef double M_f6p_c = y[17]

        # (M_atp_c) M_atp_c
        cdef double M_atp_c = y[18]

        # (M_adp_c) M_adp_c
        cdef double M_adp_c = y[19]

        # (M_fdp_c) M_fdp_c
        cdef double M_fdp_c = y[20]

        # (M_dhap_c) M_dhap_c
        cdef double M_dhap_c = y[21]

        # (M_g3p_c) M_g3p_c
        cdef double M_g3p_c = y[22]

        # (M_nad_c) M_nad_c
        cdef double M_nad_c = y[23]

        # (M_pi_c) M_pi_c
        cdef double M_pi_c = y[24]

        # (M_13dpg_c) M_13dpg_c
        cdef double M_13dpg_c = y[25]

        # (M_nadh_c) M_nadh_c
        cdef double M_nadh_c = y[26]

        # (M_nadp_c) M_nadp_c
        cdef double M_nadp_c = y[27]

        # (M_3pg_c) M_3pg_c
        cdef double M_3pg_c = y[28]

        # (M_nadph_c) M_nadph_c
        cdef double M_nadph_c = y[29]

        # (M_2pg_c) M_2pg_c
        cdef double M_2pg_c = y[30]

        # (M_pep_c) M_pep_c
        cdef double M_pep_c = y[31]

        # (M_pyr_c) M_pyr_c
        cdef double M_pyr_c = y[32]

        # (M_lac__L_c) M_lac__L_c
        cdef double M_lac__L_c = y[33]

        # (M_acald_c) M_acald_c
        cdef double M_acald_c = y[34]

        # (M_coa_c) M_coa_c
        cdef double M_coa_c = y[35]

        # (M_accoa_c) M_accoa_c
        cdef double M_accoa_c = y[36]

        # (M_actp_c) M_actp_c
        cdef double M_actp_c = y[37]

        # (M_ac_c) M_ac_c
        cdef double M_ac_c = y[38]

        # (M_o2_c) M_o2_c
        cdef double M_o2_c = y[39]

        # (M_e4p_c) M_e4p_c
        cdef double M_e4p_c = y[40]

        # (M_s7p_c) M_s7p_c
        cdef double M_s7p_c = y[41]

        # (M_r5p_c) M_r5p_c
        cdef double M_r5p_c = y[42]

        # (M_xu5p__D_c) M_xu5p__D_c
        cdef double M_xu5p__D_c = y[43]

        # (M_ru5p__D_c) M_ru5p__D_c
        cdef double M_ru5p__D_c = y[44]

        # (M_amp_c) M_amp_c
        cdef double M_amp_c = y[45]

        # (M_prpp_c) M_prpp_c
        cdef double M_prpp_c = y[46]

        # (M_r1p_c) M_r1p_c
        cdef double M_r1p_c = y[47]

        # (M_2dr1p_c) M_2dr1p_c
        cdef double M_2dr1p_c = y[48]

        # (M_2dr5p_c) M_2dr5p_c
        cdef double M_2dr5p_c = y[49]

        # (M_ade_c) M_ade_c
        cdef double M_ade_c = y[50]

        # (M_ppi_c) M_ppi_c
        cdef double M_ppi_c = y[51]

        # (M_adn_c) M_adn_c
        cdef double M_adn_c = y[52]

        # (M_dadp_c) M_dadp_c
        cdef double M_dadp_c = y[53]

        # (M_dad_2_c) M_dad_2_c
        cdef double M_dad_2_c = y[54]

        # (M_damp_c) M_damp_c
        cdef double M_damp_c = y[55]

        # (M_datp_c) M_datp_c
        cdef double M_datp_c = y[56]

        # (M_gua_c) M_gua_c
        cdef double M_gua_c = y[57]

        # (M_gmp_c) M_gmp_c
        cdef double M_gmp_c = y[58]

        # (M_gsn_c) M_gsn_c
        cdef double M_gsn_c = y[59]

        # (M_gdp_c) M_gdp_c
        cdef double M_gdp_c = y[60]

        # (M_dgdp_c) M_dgdp_c
        cdef double M_dgdp_c = y[61]

        # (M_gtp_c) M_gtp_c
        cdef double M_gtp_c = y[62]

        # (M_dgsn_c) M_dgsn_c
        cdef double M_dgsn_c = y[63]

        # (M_dgmp_c) M_dgmp_c
        cdef double M_dgmp_c = y[64]

        # (M_dgtp_c) M_dgtp_c
        cdef double M_dgtp_c = y[65]

        # (M_ura_c) M_ura_c
        cdef double M_ura_c = y[66]

        # (M_ump_c) M_ump_c
        cdef double M_ump_c = y[67]

        # (M_udp_c) M_udp_c
        cdef double M_udp_c = y[68]

        # (M_utp_c) M_utp_c
        cdef double M_utp_c = y[69]

        # (M_dudp_c) M_dudp_c
        cdef double M_dudp_c = y[70]

        # (M_cmp_c) M_cmp_c
        cdef double M_cmp_c = y[71]

        # (M_cdp_c) M_cdp_c
        cdef double M_cdp_c = y[72]

        # (M_ctp_c) M_ctp_c
        cdef double M_ctp_c = y[73]

        # (M_dcdp_c) M_dcdp_c
        cdef double M_dcdp_c = y[74]

        # (M_dcyt_c) M_dcyt_c
        cdef double M_dcyt_c = y[75]

        # (M_dcmp_c) M_dcmp_c
        cdef double M_dcmp_c = y[76]

        # (M_dctp_c) M_dctp_c
        cdef double M_dctp_c = y[77]

        # (M_duri_c) M_duri_c
        cdef double M_duri_c = y[78]

        # (M_dump_c) M_dump_c
        cdef double M_dump_c = y[79]

        # (M_thymd_c) M_thymd_c
        cdef double M_thymd_c = y[80]

        # (M_dtmp_c) M_dtmp_c
        cdef double M_dtmp_c = y[81]

        # (M_dtdp_c) M_dtdp_c
        cdef double M_dtdp_c = y[82]

        # (M_dttp_c) M_dttp_c
        cdef double M_dttp_c = y[83]

        # (M_dutp_c) M_dutp_c
        cdef double M_dutp_c = y[84]

        # (M_gln__L_c) M_gln__L_c
        cdef double M_gln__L_c = y[85]

        # (M_glu__L_c) M_glu__L_c
        cdef double M_glu__L_c = y[86]

        # (M_nh3_c) M_nh3_c
        cdef double M_nh3_c = y[87]

        # (M_uri_c) M_uri_c
        cdef double M_uri_c = y[88]

        # (M_glyc_c) M_glyc_c
        cdef double M_glyc_c = y[89]

        # (M_glyc3p_c) M_glyc3p_c
        cdef double M_glyc3p_c = y[90]

        # (M_pap_c) M_pap_c
        cdef double M_pap_c = y[91]

        # (M_fa_c) M_fa_c
        cdef double M_fa_c = y[92]

        # (M_ap_c) M_ap_c
        cdef double M_ap_c = y[93]

        # (M_1ag3p_c) M_1ag3p_c
        cdef double M_1ag3p_c = y[94]

        # (M_pa_c) M_pa_c
        cdef double M_pa_c = y[95]

        # (M_cdpdag_c) M_cdpdag_c
        cdef double M_cdpdag_c = y[96]

        # (M_pg3p_c) M_pg3p_c
        cdef double M_pg3p_c = y[97]

        # (M_pg_c) M_pg_c
        cdef double M_pg_c = y[98]

        # (M_clpn_c) M_clpn_c
        cdef double M_clpn_c = y[99]

        # (M_12dgr_c) M_12dgr_c
        cdef double M_12dgr_c = y[100]

        # (M_g1p_c) M_g1p_c
        cdef double M_g1p_c = y[101]

        # (M_udpg_c) M_udpg_c
        cdef double M_udpg_c = y[102]

        # (M_udpgal_c) M_udpgal_c
        cdef double M_udpgal_c = y[103]

        # (M_udpgalfur_c) M_udpgalfur_c
        cdef double M_udpgalfur_c = y[104]

        # (M_galfur12dgr_c) M_galfur12dgr_c
        cdef double M_galfur12dgr_c = y[105]

        # (M_nac_c) M_nac_c
        cdef double M_nac_c = y[106]

        # (M_nicrnt_c) M_nicrnt_c
        cdef double M_nicrnt_c = y[107]

        # (M_dnad_c) M_dnad_c
        cdef double M_dnad_c = y[108]

        # (M_ribflv_c) M_ribflv_c
        cdef double M_ribflv_c = y[109]

        # (M_fmn_c) M_fmn_c
        cdef double M_fmn_c = y[110]

        # (M_fad_c) M_fad_c
        cdef double M_fad_c = y[111]

        # (M_5fthf_c) M_5fthf_c
        cdef double M_5fthf_c = y[112]

        # (M_5fthfglu3_c) M_5fthfglu3_c
        cdef double M_5fthfglu3_c = y[113]

        # (M_methfglu3_c) M_methfglu3_c
        cdef double M_methfglu3_c = y[114]

        # (M_10fthfglu3_c) M_10fthfglu3_c
        cdef double M_10fthfglu3_c = y[115]

        # (M_mettrna_c) M_mettrna_c
        cdef double M_mettrna_c = y[116]

        # (M_fmettrna_c) M_fmettrna_c
        cdef double M_fmettrna_c = y[117]

        # (M_thfglu3_c) M_thfglu3_c
        cdef double M_thfglu3_c = y[118]

        # (M_ser__L_c) M_ser__L_c
        cdef double M_ser__L_c = y[119]

        # (M_gly_c) M_gly_c
        cdef double M_gly_c = y[120]

        # (M_mlthfglu3_c) M_mlthfglu3_c
        cdef double M_mlthfglu3_c = y[121]

        # (M_pydx5p_c) M_pydx5p_c
        cdef double M_pydx5p_c = y[122]

        # (M_thmpp_c) M_thmpp_c
        cdef double M_thmpp_c = y[123]

        # (M_sprm_c) M_sprm_c
        cdef double M_sprm_c = y[124]

        # (M_na1_c) M_na1_c
        cdef double M_na1_c = y[125]

        # (M_k_c) M_k_c
        cdef double M_k_c = y[126]

        # (M_mg2_c) M_mg2_c
        cdef double M_mg2_c = y[127]

        # (M_ca2_c) M_ca2_c
        cdef double M_ca2_c = y[128]

        # (M_arg__L_c) M_arg__L_c
        cdef double M_arg__L_c = y[129]

        # (M_asp__L_c) M_asp__L_c
        cdef double M_asp__L_c = y[130]

        # (M_ile__L_c) M_ile__L_c
        cdef double M_ile__L_c = y[131]

        # (M_ala__L_c) M_ala__L_c
        cdef double M_ala__L_c = y[132]

        # (M_asn__L_c) M_asn__L_c
        cdef double M_asn__L_c = y[133]

        # (M_leu__L_c) M_leu__L_c
        cdef double M_leu__L_c = y[134]

        # (M_his__L_c) M_his__L_c
        cdef double M_his__L_c = y[135]

        # (M_lys__L_c) M_lys__L_c
        cdef double M_lys__L_c = y[136]

        # (M_pro__L_c) M_pro__L_c
        cdef double M_pro__L_c = y[137]

        # (M_phe__L_c) M_phe__L_c
        cdef double M_phe__L_c = y[138]

        # (M_thr__L_c) M_thr__L_c
        cdef double M_thr__L_c = y[139]

        # (M_trp__L_c) M_trp__L_c
        cdef double M_trp__L_c = y[140]

        # (M_tyr__L_c) M_tyr__L_c
        cdef double M_tyr__L_c = y[141]

        # (M_val__L_c) M_val__L_c
        cdef double M_val__L_c = y[142]

        # (M_met__L_c) M_met__L_c
        cdef double M_met__L_c = y[143]

        # (M_cys__L_c) M_cys__L_c
        cdef double M_cys__L_c = y[144]

        # (M_sm_c) M_sm_c
        cdef double M_sm_c = y[145]

        # (M_pc_c) M_pc_c
        cdef double M_pc_c = y[146]

        # (M_tag_c) M_tag_c
        cdef double M_tag_c = y[147]

        # (M_chsterol_c) M_chsterol_c
        cdef double M_chsterol_c = y[148]



        # (PGI) Reaction PGI
        cdef double V1 = self.params[1] * self.params[0] * ( ( 804.3384 * ( M_g6p_c / 0.28 ) - 650 * ( M_f6p_c / 0.15 ) ) / ( ( 1 + M_g6p_c / 0.28 ) + ( 1 + M_f6p_c / 0.15 ) - 1 ) )

        # (PFK) Reaction PFK
        cdef double V2 = self.params[3] * self.params[2] * ( ( 111 * ( M_atp_c / 0.117 ) * ( M_f6p_c / 0.2133 ) - 0.3438 * ( M_adp_c / 0.1062 ) * ( M_fdp_c / 0.1062 ) ) / ( ( 1 + M_atp_c / 0.117 ) * ( 1 + M_f6p_c / 0.2133 ) + ( 1 + M_adp_c / 0.1062 ) * ( 1 + M_fdp_c / 0.1062 ) - 1 ) )

        # (FBA) Reaction FBA
        cdef double V3 = self.params[5] * self.params[4] * ( ( 64.5 * ( M_fdp_c / 0.12 ) - 12.6 * ( M_dhap_c / 0.095 ) * ( M_g3p_c / 0.052 ) ) / ( ( 1 + M_fdp_c / 0.12 ) + ( 1 + M_dhap_c / 0.095 ) * ( 1 + M_g3p_c / 0.052 ) - 1 ) )

        # (TPI) Reaction TPI
        cdef double V4 = self.params[7] * self.params[6] * ( ( 759.6105 * ( M_dhap_c / 0.077 ) - 67 * ( M_g3p_c / 0.084 ) ) / ( ( 1 + M_dhap_c / 0.077 ) + ( 1 + M_g3p_c / 0.084 ) - 1 ) )

        # (GAPD) Reaction GAPD
        cdef double V5 = self.params[9] * self.params[8] * ( ( 442 * ( M_g3p_c / 0.0171 ) * ( M_nad_c / 1.3 ) * ( M_pi_c / 6.5826 ) - 73.6 * ( M_13dpg_c / 0.47 ) * ( M_nadh_c / 0.061 ) ) / ( ( 1 + M_g3p_c / 0.0171 ) * ( 1 + M_nad_c / 1.3 ) * ( 1 + M_pi_c / 6.5826 ) + ( 1 + M_13dpg_c / 0.47 ) * ( 1 + M_nadh_c / 0.061 ) - 1 ) )

        # (GAPDP) Reaction GAPDP
        cdef double V6 = self.params[11] * self.params[10] * ( ( 2.8 * ( M_g3p_c / 0.052 ) * ( M_nadp_c / 0.385 ) - 0 * ( M_3pg_c / 0.1924 ) * ( M_nadph_c / 0.202 ) ) / ( ( 1 + M_g3p_c / 0.052 ) * ( 1 + M_nadp_c / 0.385 ) + ( 1 + M_3pg_c / 0.1924 ) * ( 1 + M_nadph_c / 0.202 ) - 1 ) )

        # (PGK) Reaction PGK
        cdef double V7 = self.params[13] * self.params[12] * ( ( 220 * ( M_13dpg_c / 0.01 ) * ( M_adp_c / 0.5273 ) - 3.4 * ( M_3pg_c / 0.1 ) * ( M_atp_c / 0.368 ) ) / ( ( 1 + M_13dpg_c / 0.01 ) * ( 1 + M_adp_c / 0.5273 ) + ( 1 + M_3pg_c / 0.1 ) * ( 1 + M_atp_c / 0.368 ) - 1 ) )

        # (PGM) Reaction PGM
        cdef double V8 = self.params[15] * self.params[14] * ( ( 434 * ( M_3pg_c / 3.6 ) - 14 * ( M_2pg_c / 0.2 ) ) / ( ( 1 + M_3pg_c / 3.6 ) + ( 1 + M_2pg_c / 0.2 ) - 1 ) )

        # (ENO) Reaction ENO
        cdef double V9 = self.params[17] * self.params[16] * ( ( 62.1825 * ( M_2pg_c / 0.246 ) - 23.453 * ( M_pep_c / 0.2874 ) ) / ( ( 1 + M_2pg_c / 0.246 ) + ( 1 + M_pep_c / 0.2874 ) - 1 ) )

        # (PYK) Reaction PYK
        cdef double V10 = self.params[19] * self.params[18] * ( ( 3204 * ( M_adp_c / 0.2814 ) * ( M_pep_c / 0.075 ) - 0.0021 * ( M_atp_c / 0.121 ) * ( M_pyr_c / 0.0585 ) ) / ( ( 1 + M_adp_c / 0.2814 ) * ( 1 + M_pep_c / 0.075 ) + ( 1 + M_atp_c / 0.121 ) * ( 1 + M_pyr_c / 0.0585 ) - 1 ) )

        # (LDH_L) Reaction LDH_L
        cdef double V11 = self.params[21] * self.params[20] * ( ( 388.468 * ( M_nadh_c / 0.0381 ) * ( M_pyr_c / 0.9649 ) - 400.2106 * ( M_lac__L_c / 59.656 ) * ( M_nad_c / 32.972 ) ) / ( ( 1 + M_nadh_c / 0.0381 ) * ( 1 + M_pyr_c / 0.9649 ) + ( 1 + M_lac__L_c / 59.656 ) * ( 1 + M_nad_c / 32.972 ) - 1 ) )

        # (PDH_acald) Reaction PDH_acald
        cdef double V12 = self.params[23] * self.params[22] * ( ( 38.1501 * ( M_acald_c / 0.1 ) * ( M_coa_c / 0.1 ) * ( M_lpl_PdhC_c / 0.1 ) - 2.6212 * ( M_accoa_c / 0.1 ) * ( M_dhlpl_PdhC_c / 0.1 ) ) / ( ( 1 + M_acald_c / 0.1 ) * ( 1 + M_coa_c / 0.1 ) * ( 1 + M_lpl_PdhC_c / 0.1 ) + ( 1 + M_accoa_c / 0.1 ) * ( 1 + M_dhlpl_PdhC_c / 0.1 ) - 1 ) )

        # (PDH_E3) Reaction PDH_E3
        cdef double V13 = self.params[25] * self.params[24] * ( ( 110 * ( M_dhlpl_PdhC_c / 0.1003 ) * ( M_nad_c / 0.14 ) - 448 * ( M_lpl_PdhC_c / 0.0997 ) * ( M_nadh_c / 2.9919 ) ) / ( ( 1 + M_dhlpl_PdhC_c / 0.1003 ) * ( 1 + M_nad_c / 0.14 ) + ( 1 + M_lpl_PdhC_c / 0.0997 ) * ( 1 + M_nadh_c / 2.9919 ) - 1 ) )

        # (PTAr) Reaction PTAr
        cdef double V14 = self.params[27] * self.params[26] * ( ( 895.9342 * ( M_accoa_c / 0.5535 ) * ( M_pi_c / 5.291 ) - 1509.0778 * ( M_actp_c / 0.2373 ) * ( M_coa_c / 0.022 ) ) / ( ( 1 + M_accoa_c / 0.5535 ) * ( 1 + M_pi_c / 5.291 ) + ( 1 + M_actp_c / 0.2373 ) * ( 1 + M_coa_c / 0.022 ) - 1 ) )

        # (ACKr) Reaction ACKr
        cdef double V15 = self.params[29] * self.params[28] * ( ( 1452.2498 * ( M_actp_c / 0.9669 ) * ( M_adp_c / 0.6657 ) - 655.0835 * ( M_ac_c / 55.7644 ) * ( M_atp_c / 1.6092 ) ) / ( ( 1 + M_actp_c / 0.9669 ) * ( 1 + M_adp_c / 0.6657 ) + ( 1 + M_ac_c / 55.7644 ) * ( 1 + M_atp_c / 1.6092 ) - 1 ) )

        # (NOX) Reaction NOX
        cdef double V16 = self.params[31] * self.params[30] * ( ( 1520.2025 * ( M_nadh_c / 0.0268 ) * ( M_nadh_c / 0.0268 ) * ( M_o2_c / 0.0204 ) - 0.0693 * ( M_nad_c / 0.0993 ) * ( M_nad_c / 0.0993 ) ) / ( ( 1 + M_nadh_c / 0.0268 ) * ( 1 + M_nadh_c / 0.0268 ) * ( 1 + M_o2_c / 0.0204 ) + ( 1 + M_nad_c / 0.0993 ) * ( 1 + M_nad_c / 0.0993 ) - 1 ) )

        # (TALA) Reaction TALA
        cdef double V17 = self.params[33] * self.params[32] * ( ( 22.3 * ( M_e4p_c / 0.0401 ) * ( M_f6p_c / 0.6688 ) - 0.54 * ( M_g3p_c / 1.9 ) * ( M_s7p_c / 0.285 ) ) / ( ( 1 + M_e4p_c / 0.0401 ) * ( 1 + M_f6p_c / 0.6688 ) + ( 1 + M_g3p_c / 1.9 ) * ( 1 + M_s7p_c / 0.285 ) - 1 ) )

        # (TKT1) Reaction TKT1
        cdef double V18 = self.params[35] * self.params[34] * ( ( 20.58 * ( M_g3p_c / 0.743 ) * ( M_s7p_c / 3.7298 ) - 0.8 * ( M_r5p_c / 0.4717 ) * ( M_xu5p__D_c / 0.134 ) ) / ( ( 1 + M_g3p_c / 0.743 ) * ( 1 + M_s7p_c / 3.7298 ) + ( 1 + M_r5p_c / 0.4717 ) * ( 1 + M_xu5p__D_c / 0.134 ) - 1 ) )

        # (TKT2) Reaction TKT2
        cdef double V19 = self.params[37] * self.params[36] * ( ( 26.87 * ( M_f6p_c / 0.25 ) * ( M_g3p_c / 0.743 ) - 1.4 * ( M_e4p_c / 0.0227 ) * ( M_xu5p__D_c / 0.134 ) ) / ( ( 1 + M_f6p_c / 0.25 ) * ( 1 + M_g3p_c / 0.743 ) + ( 1 + M_e4p_c / 0.0227 ) * ( 1 + M_xu5p__D_c / 0.134 ) - 1 ) )

        # (RPE) Reaction RPE
        cdef double V20 = self.params[39] * self.params[38] * ( ( 1337.133 * ( M_xu5p__D_c / 2.9022 ) - 1301.6531 * ( M_ru5p__D_c / 1.9567 ) ) / ( ( 1 + M_xu5p__D_c / 2.9022 ) + ( 1 + M_ru5p__D_c / 1.9567 ) - 1 ) )

        # (RPI) Reaction RPI
        cdef double V21 = self.params[41] * self.params[40] * ( ( 10 * ( M_ru5p__D_c / 2.1524 ) - 1 * ( M_r5p_c / 3.3617 ) ) / ( ( 1 + M_ru5p__D_c / 2.1524 ) + ( 1 + M_r5p_c / 3.3617 ) - 1 ) )

        # (PRPPS) Reaction PRPPS
        cdef double V22 = self.params[43] * self.params[42] * ( ( 30.2215 * ( M_atp_c / 0.0749 ) * ( M_r5p_c / 0.0676 ) - 10.9581 * ( M_amp_c / 0.0453 ) * ( M_prpp_c / 0.0141 ) ) / ( ( 1 + M_atp_c / 0.0749 ) * ( 1 + M_r5p_c / 0.0676 ) + ( 1 + M_amp_c / 0.0453 ) * ( 1 + M_prpp_c / 0.0141 ) - 1 ) )

        # (PPM) Reaction PPM
        cdef double V23 = self.params[45] * self.params[44] * ( ( 7.8902 * ( M_r1p_c / 0.073 ) - 0.0822 * ( M_r5p_c / 0.0854 ) ) / ( ( 1 + M_r1p_c / 0.073 ) + ( 1 + M_r5p_c / 0.0854 ) - 1 ) )

        # (PPM2) Reaction PPM2
        cdef double V24 = self.params[47] * self.params[46] * ( ( 173 * ( M_2dr1p_c / 0.013 ) - 10.201 * ( M_2dr5p_c / 1.2 ) ) / ( ( 1 + M_2dr1p_c / 0.013 ) + ( 1 + M_2dr5p_c / 1.2 ) - 1 ) )

        # (DRPA) Reaction DRPA
        cdef double V25 = self.params[49] * self.params[48] * ( ( 521 * ( M_2dr5p_c / 0.2456 ) - 34 * ( M_acald_c / 0.267 ) * ( M_g3p_c / 0.2 ) ) / ( ( 1 + M_2dr5p_c / 0.2456 ) + ( 1 + M_acald_c / 0.267 ) * ( 1 + M_g3p_c / 0.2 ) - 1 ) )

        # (ADPT) Reaction ADPT
        cdef double V26 = self.params[51] * self.params[50] * ( ( 6.5832 * ( M_ade_c / 0.0044 ) * ( M_prpp_c / 0.0085 ) - 0.895 * ( M_amp_c / 0.0319 ) * ( M_ppi_c / 0.1586 ) ) / ( ( 1 + M_ade_c / 0.0044 ) * ( 1 + M_prpp_c / 0.0085 ) + ( 1 + M_amp_c / 0.0319 ) * ( 1 + M_ppi_c / 0.1586 ) - 1 ) )

        # (PUNP1) Reaction PUNP1
        cdef double V27 = self.params[53] * self.params[52] * ( ( 26.2772 * ( M_adn_c / 0.7652 ) * ( M_pi_c / 2.1387 ) - 109.8453 * ( M_ade_c / 4.0132 ) * ( M_r1p_c / 0.0142 ) ) / ( ( 1 + M_adn_c / 0.7652 ) * ( 1 + M_pi_c / 2.1387 ) + ( 1 + M_ade_c / 4.0132 ) * ( 1 + M_r1p_c / 0.0142 ) - 1 ) )

        # (ADK1) Reaction ADK1
        cdef double V28 = self.params[55] * self.params[54] * ( ( 319.1556 * ( M_amp_c / 0.096 ) * ( M_atp_c / 0.1426 ) - 783.2866 * ( M_adp_c / 0.2669 ) * ( M_adp_c / 0.2669 ) ) / ( ( 1 + M_amp_c / 0.096 ) * ( 1 + M_atp_c / 0.1426 ) + ( 1 + M_adp_c / 0.2669 ) * ( 1 + M_adp_c / 0.2669 ) - 1 ) )

        # (RNDR1) Reaction RNDR1
        cdef double V29 = self.params[57] * self.params[56] * ( ( 4.3833 * ( M_adp_c / 0.0206 ) * ( M_trdrd_c / 0.0389 ) - 0 * ( M_dadp_c / 0.2569 ) * ( M_trdox_c / 0.2569 ) ) / ( ( 1 + M_adp_c / 0.0206 ) * ( 1 + M_trdrd_c / 0.0389 ) + ( 1 + M_dadp_c / 0.2569 ) * ( 1 + M_trdox_c / 0.2569 ) - 1 ) )

        # (TRDR) Reaction TRDR
        cdef double V30 = self.params[59] * self.params[58] * ( ( 19.0309 * ( M_nadph_c / 0.0468 ) * ( M_trdox_c / 0.0129 ) - 139.1897 * ( M_nadp_c / 0.7737 ) * ( M_trdrd_c / 0.7737 ) ) / ( ( 1 + M_nadph_c / 0.0468 ) * ( 1 + M_trdox_c / 0.0129 ) + ( 1 + M_nadp_c / 0.7737 ) * ( 1 + M_trdrd_c / 0.7737 ) - 1 ) )

        # (PPA) Reaction PPA
        cdef double V31 = self.params[61] * self.params[60] * ( ( 646.7271 * ( M_ppi_c / 0.0024 ) - 0.1874 * ( M_pi_c / 0.0976 ) * ( M_pi_c / 0.0976 ) ) / ( ( 1 + M_ppi_c / 0.0024 ) + ( 1 + M_pi_c / 0.0976 ) * ( 1 + M_pi_c / 0.0976 ) - 1 ) )

        # (PUNP2) Reaction PUNP2
        cdef double V32 = self.params[63] * self.params[62] * ( ( 3.3 * ( M_dad_2_c / 0.7005 ) * ( M_pi_c / 2.1382 ) - 109.8439 * ( M_2dr1p_c / 0.0143 ) * ( M_ade_c / 0.41 ) ) / ( ( 1 + M_dad_2_c / 0.7005 ) * ( 1 + M_pi_c / 2.1382 ) + ( 1 + M_2dr1p_c / 0.0143 ) * ( 1 + M_ade_c / 0.41 ) - 1 ) )

        # (DADNK) Reaction DADNK
        cdef double V33 = self.params[65] * self.params[64] * ( ( 3.7 * ( M_atp_c / 0.15 ) * ( M_dad_2_c / 0.012 ) - 0.4533 * ( M_adp_c / 0.1435 ) * ( M_damp_c / 0.1435 ) ) / ( ( 1 + M_atp_c / 0.15 ) * ( 1 + M_dad_2_c / 0.012 ) + ( 1 + M_adp_c / 0.1435 ) * ( 1 + M_damp_c / 0.1435 ) - 1 ) )

        # (DADK) Reaction DADK
        cdef double V34 = self.params[67] * self.params[66] * ( ( 2.9129 * ( M_atp_c / 0.1 ) * ( M_damp_c / 0.1 ) - 3.43299 * ( M_adp_c / 1 ) * ( M_dadp_c / 0.1 ) ) / ( ( 1 + M_atp_c / 0.1 ) * ( 1 + M_damp_c / 0.1 ) + ( 1 + M_adp_c / 1 ) * ( 1 + M_dadp_c / 0.1 ) - 1 ) )

        # (PYK5) Reaction PYK5
        cdef double V35 = self.params[69] * self.params[68] * ( ( 576.72 * ( M_dadp_c / 0.2814 ) * ( M_pep_c / 0.075 ) - 0.02 * ( M_datp_c / 0.1168 ) * ( M_pyr_c / 0.0553 ) ) / ( ( 1 + M_dadp_c / 0.2814 ) * ( 1 + M_pep_c / 0.075 ) + ( 1 + M_datp_c / 0.1168 ) * ( 1 + M_pyr_c / 0.0553 ) - 1 ) )

        # (PGK2) Reaction PGK2
        cdef double V36 = self.params[71] * self.params[70] * ( ( 37.4 * ( M_13dpg_c / 0.0057 ) * ( M_dadp_c / 0.5357 ) - 0.07699573353053 * ( M_3pg_c / 0.7375 ) * ( M_datp_c / 0.0902 ) ) / ( ( 1 + M_13dpg_c / 0.0057 ) * ( 1 + M_dadp_c / 0.5357 ) + ( 1 + M_3pg_c / 0.7375 ) * ( 1 + M_datp_c / 0.0902 ) - 1 ) )

        # (GUAPRT) Reaction GUAPRT
        cdef double V37 = self.params[73] * self.params[72] * ( ( 14.4958 * ( M_gua_c / 0.0005 ) * ( M_prpp_c / 0.009 ) - 5.2464 * ( M_gmp_c / 0.3931 ) * ( M_ppi_c / 0.4352 ) ) / ( ( 1 + M_gua_c / 0.0005 ) * ( 1 + M_prpp_c / 0.009 ) + ( 1 + M_gmp_c / 0.3931 ) * ( 1 + M_ppi_c / 0.4352 ) - 1 ) )

        # (PUNP3) Reaction PUNP3
        cdef double V38 = self.params[75] * self.params[74] * ( ( 26.2373 * ( M_gsn_c / 1.2701 ) * ( M_pi_c / 2.1432 ) - 109.8576 * ( M_gua_c / 1.8086 ) * ( M_r1p_c / 0.0132 ) ) / ( ( 1 + M_gsn_c / 1.2701 ) * ( 1 + M_pi_c / 2.1432 ) + ( 1 + M_gua_c / 1.8086 ) * ( 1 + M_r1p_c / 0.0132 ) - 1 ) )

        # (GK1) Reaction GK1
        cdef double V39 = self.params[77] * self.params[76] * ( ( 410.2268 * ( M_atp_c / 0.1501 ) * ( M_gmp_c / 0.0259 ) - 127.921 * ( M_adp_c / 0.0676 ) * ( M_gdp_c / 0.0236 ) ) / ( ( 1 + M_atp_c / 0.1501 ) * ( 1 + M_gmp_c / 0.0259 ) + ( 1 + M_adp_c / 0.0676 ) * ( 1 + M_gdp_c / 0.0236 ) - 1 ) )

        # (RNDR2) Reaction RNDR2
        cdef double V40 = self.params[79] * self.params[78] * ( ( 4.1518 * ( M_gdp_c / 0.24 ) * ( M_trdrd_c / 0.0429 ) - 0 * ( M_dgdp_c / 0.2331 ) * ( M_trdox_c / 0.2331 ) ) / ( ( 1 + M_gdp_c / 0.24 ) * ( 1 + M_trdrd_c / 0.0429 ) + ( 1 + M_dgdp_c / 0.2331 ) * ( 1 + M_trdox_c / 0.2331 ) - 1 ) )

        # (PYK3) Reaction PYK3
        cdef double V41 = self.params[81] * self.params[80] * ( ( 672.84 * ( M_gdp_c / 0.2814 ) * ( M_pep_c / 0.075 ) - 0.0014 * ( M_gtp_c / 0.1376 ) * ( M_pyr_c / 0.0589 ) ) / ( ( 1 + M_gdp_c / 0.2814 ) * ( 1 + M_pep_c / 0.075 ) + ( 1 + M_gtp_c / 0.1376 ) * ( 1 + M_pyr_c / 0.0589 ) - 1 ) )

        # (PGK3) Reaction PGK3
        cdef double V42 = self.params[83] * self.params[82] * ( ( 140.8 * ( M_13dpg_c / 0.0032 ) * ( M_gdp_c / 0.5357 ) - 0.030214276516037 * ( M_3pg_c / 0.7944 ) * ( M_gtp_c / 0.2033 ) ) / ( ( 1 + M_13dpg_c / 0.0032 ) * ( 1 + M_gdp_c / 0.5357 ) + ( 1 + M_3pg_c / 0.7944 ) * ( 1 + M_gtp_c / 0.2033 ) - 1 ) )

        # (PUNP4) Reaction PUNP4
        cdef double V43 = self.params[85] * self.params[84] * ( ( 5.1 * ( M_dgsn_c / 0.183 ) * ( M_pi_c / 2.2063 ) - 110.0264 * ( M_2dr1p_c / 0.0049 ) * ( M_gua_c / 0.8 ) ) / ( ( 1 + M_dgsn_c / 0.183 ) * ( 1 + M_pi_c / 2.2063 ) + ( 1 + M_2dr1p_c / 0.0049 ) * ( 1 + M_gua_c / 0.8 ) - 1 ) )

        # (DGSNK) Reaction DGSNK
        cdef double V44 = self.params[87] * self.params[86] * ( ( 2.25 * ( M_atp_c / 0.0698 ) * ( M_dgsn_c / 0.0061 ) - 0.4437 * ( M_adp_c / 0.1484 ) * ( M_dgmp_c / 0.1484 ) ) / ( ( 1 + M_atp_c / 0.0698 ) * ( 1 + M_dgsn_c / 0.0061 ) + ( 1 + M_adp_c / 0.1484 ) * ( 1 + M_dgmp_c / 0.1484 ) - 1 ) )

        # (DGK1) Reaction DGK1
        cdef double V45 = self.params[89] * self.params[88] * ( ( 410.1815 * ( M_atp_c / 0.1549 ) * ( M_dgmp_c / 0.0441 ) - 831.9981 * ( M_adp_c / 0.0601 ) * ( M_dgdp_c / 0.0601 ) ) / ( ( 1 + M_atp_c / 0.1549 ) * ( 1 + M_dgmp_c / 0.0441 ) + ( 1 + M_adp_c / 0.0601 ) * ( 1 + M_dgdp_c / 0.0601 ) - 1 ) )

        # (PYK6) Reaction PYK6
        cdef double V46 = self.params[91] * self.params[90] * ( ( 833.04 * ( M_dgdp_c / 0.2814 ) * ( M_pep_c / 0.075 ) - 0.02 * ( M_dgtp_c / 0.1168 ) * ( M_pyr_c / 0.0553 ) ) / ( ( 1 + M_dgdp_c / 0.2814 ) * ( 1 + M_pep_c / 0.075 ) + ( 1 + M_dgtp_c / 0.1168 ) * ( 1 + M_pyr_c / 0.0553 ) - 1 ) )

        # (PGK4) Reaction PGK4
        cdef double V47 = self.params[93] * self.params[92] * ( ( 35.2 * ( M_13dpg_c / 0.0057 ) * ( M_dgdp_c / 0.5357 ) - 0.072466572734617 * ( M_3pg_c / 0.7375 ) * ( M_dgtp_c / 0.0902 ) ) / ( ( 1 + M_13dpg_c / 0.0057 ) * ( 1 + M_dgdp_c / 0.5357 ) + ( 1 + M_3pg_c / 0.7375 ) * ( 1 + M_dgtp_c / 0.0902 ) - 1 ) )

        # (UPPRT) Reaction UPPRT
        cdef double V48 = self.params[95] * self.params[94] * ( ( 5.1 * ( M_prpp_c / 0.084 ) * ( M_ura_c / 0.00407 ) - 0.0035 * ( M_ppi_c / 1 ) * ( M_ump_c / 0.0138 ) ) / ( ( 1 + M_prpp_c / 0.084 ) * ( 1 + M_ura_c / 0.00407 ) + ( 1 + M_ppi_c / 1 ) * ( 1 + M_ump_c / 0.0138 ) - 1 ) )

        # (UMPK) Reaction UMPK
        cdef double V49 = self.params[97] * self.params[96] * ( ( 263.7845 * ( M_atp_c / 2.3428 ) * ( M_ump_c / 0.2951 ) - 0.7058 * ( M_adp_c / 0.0962 ) * ( M_udp_c / 0.0962 ) ) / ( ( 1 + M_atp_c / 2.3428 ) * ( 1 + M_ump_c / 0.2951 ) + ( 1 + M_adp_c / 0.0962 ) * ( 1 + M_udp_c / 0.0962 ) - 1 ) )

        # (PYK2) Reaction PYK2
        cdef double V50 = self.params[99] * self.params[98] * ( ( 352.44 * ( M_pep_c / 0.075 ) * ( M_udp_c / 0.2814 ) - 0.0004 * ( M_pyr_c / 0.0608 ) * ( M_utp_c / 0.1491 ) ) / ( ( 1 + M_pep_c / 0.075 ) * ( 1 + M_udp_c / 0.2814 ) + ( 1 + M_pyr_c / 0.0608 ) * ( 1 + M_utp_c / 0.1491 ) - 1 ) )

        # (RNDR4) Reaction RNDR4
        cdef double V51 = self.params[101] * self.params[100] * ( ( 0.29 * ( M_trdrd_c / 0.0395 ) * ( M_udp_c / 1.2 ) - 0 * ( M_dudp_c / 0.253 ) * ( M_trdox_c / 0.253 ) ) / ( ( 1 + M_trdrd_c / 0.0395 ) * ( 1 + M_udp_c / 1.2 ) + ( 1 + M_dudp_c / 0.253 ) * ( 1 + M_trdox_c / 0.253 ) - 1 ) )

        # (CYTK1) Reaction CYTK1
        cdef double V52 = self.params[103] * self.params[102] * ( ( 18.464 * ( M_atp_c / 0.1 ) * ( M_cmp_c / 0.1 ) - 5.416 * ( M_adp_c / 0.1 ) * ( M_cdp_c / 0.1 ) ) / ( ( 1 + M_atp_c / 0.1 ) * ( 1 + M_cmp_c / 0.1 ) + ( 1 + M_adp_c / 0.1 ) * ( 1 + M_cdp_c / 0.1 ) - 1 ) )

        # (PYK4) Reaction PYK4
        cdef double V53 = self.params[105] * self.params[104] * ( ( 544.68 * ( M_cdp_c / 0.2814 ) * ( M_pep_c / 0.075 ) - 0.0001 * ( M_ctp_c / 0.1611 ) * ( M_pyr_c / 0.0627 ) ) / ( ( 1 + M_cdp_c / 0.2814 ) * ( 1 + M_pep_c / 0.075 ) + ( 1 + M_ctp_c / 0.1611 ) * ( 1 + M_pyr_c / 0.0627 ) - 1 ) )

        # (RNDR3) Reaction RNDR3
        cdef double V54 = self.params[107] * self.params[106] * ( ( 4.1273 * ( M_cdp_c / 0.057 ) * ( M_trdrd_c / 0.0433 ) - 0 * ( M_dcdp_c / 0.2307 ) * ( M_trdox_c / 0.2307 ) ) / ( ( 1 + M_cdp_c / 0.057 ) * ( 1 + M_trdrd_c / 0.0433 ) + ( 1 + M_dcdp_c / 0.2307 ) * ( 1 + M_trdox_c / 0.2307 ) - 1 ) )

        # (DCYTK) Reaction DCYTK
        cdef double V55 = self.params[109] * self.params[108] * ( ( 1.8 * ( M_atp_c / 0.0136 ) * ( M_dcyt_c / 0.0023 ) - 0.0058 * ( M_adp_c / 0.2097 ) * ( M_dcmp_c / 0.2097 ) ) / ( ( 1 + M_atp_c / 0.0136 ) * ( 1 + M_dcyt_c / 0.0023 ) + ( 1 + M_adp_c / 0.2097 ) * ( 1 + M_dcmp_c / 0.2097 ) - 1 ) )

        # (CYTK2) Reaction CYTK2
        cdef double V56 = self.params[111] * self.params[110] * ( ( 18.566 * ( M_atp_c / 0.1 ) * ( M_dcmp_c / 0.1 ) - 5.3862 * ( M_adp_c / 0.1 ) * ( M_dcdp_c / 0.1 ) ) / ( ( 1 + M_atp_c / 0.1 ) * ( 1 + M_dcmp_c / 0.1 ) + ( 1 + M_adp_c / 0.1 ) * ( 1 + M_dcdp_c / 0.1 ) - 1 ) )

        # (PYK7) Reaction PYK7
        cdef double V57 = self.params[113] * self.params[112] * ( ( 416.52 * ( M_dcdp_c / 0.2814 ) * ( M_pep_c / 0.075 ) - 0.0308 * ( M_dctp_c / 0.1137 ) * ( M_pyr_c / 0.0547 ) ) / ( ( 1 + M_dcdp_c / 0.2814 ) * ( 1 + M_pep_c / 0.075 ) + ( 1 + M_dctp_c / 0.1137 ) * ( 1 + M_pyr_c / 0.0547 ) - 1 ) )

        # (DURIK1) Reaction DURIK1
        cdef double V58 = self.params[115] * self.params[114] * ( ( 19.3256 * ( M_atp_c / 0.4221 ) * ( M_duri_c / 0.0283 ) - 0.1482 * ( M_adp_c / 0.1249 ) * ( M_dump_c / 0.1249 ) ) / ( ( 1 + M_atp_c / 0.4221 ) * ( 1 + M_duri_c / 0.0283 ) + ( 1 + M_adp_c / 0.1249 ) * ( 1 + M_dump_c / 0.1249 ) - 1 ) )

        # (TMDK1) Reaction TMDK1
        cdef double V59 = self.params[117] * self.params[116] * ( ( 19.2637 * ( M_atp_c / 0.4246 ) * ( M_thymd_c / 0.0061 ) - 0.2901 * ( M_adp_c / 0.1198 ) * ( M_dtmp_c / 0.1198 ) ) / ( ( 1 + M_atp_c / 0.4246 ) * ( 1 + M_thymd_c / 0.0061 ) + ( 1 + M_adp_c / 0.1198 ) * ( 1 + M_dtmp_c / 0.1198 ) - 1 ) )

        # (TMPK) Reaction TMPK
        cdef double V60 = self.params[119] * self.params[118] * ( ( 4.7572 * ( M_atp_c / 0.0524 ) * ( M_dtmp_c / 0.0168 ) - 41.5553 * ( M_adp_c / 0.0958 ) * ( M_dtdp_c / 0.0958 ) ) / ( ( 1 + M_atp_c / 0.0524 ) * ( 1 + M_dtmp_c / 0.0168 ) + ( 1 + M_adp_c / 0.0958 ) * ( 1 + M_dtdp_c / 0.0958 ) - 1 ) )

        # (PYK8) Reaction PYK8
        cdef double V61 = self.params[121] * self.params[120] * ( ( 160.2 * ( M_dtdp_c / 0.2814 ) * ( M_pep_c / 0.075 ) - 0.0184 * ( M_dttp_c / 0.1174 ) * ( M_pyr_c / 0.0554 ) ) / ( ( 1 + M_dtdp_c / 0.2814 ) * ( 1 + M_pep_c / 0.075 ) + ( 1 + M_dttp_c / 0.1174 ) * ( 1 + M_pyr_c / 0.0554 ) - 1 ) )

        # (PYK9) Reaction PYK9
        cdef double V62 = self.params[123] * self.params[122] * ( ( 32.04 * ( M_dudp_c / 0.2814 ) * ( M_pep_c / 0.075 ) - 0.0126 * ( M_dutp_c / 0.1202 ) * ( M_pyr_c / 0.0559 ) ) / ( ( 1 + M_dudp_c / 0.2814 ) * ( 1 + M_pep_c / 0.075 ) + ( 1 + M_dutp_c / 0.1202 ) * ( 1 + M_pyr_c / 0.0559 ) - 1 ) )

        # (CTPS2) Reaction CTPS2
        cdef double V63 = self.params[125] * self.params[124] * ( ( 2.2 * ( M_atp_c / 0.0164 ) * ( M_gln__L_c / 0.354 ) * ( M_utp_c / 0.071 ) - 0.027698870056497 * ( M_adp_c / 0.1326 ) * ( M_ctp_c / 0.1326 ) * ( M_glu__L_c / 0.1326 ) * ( M_pi_c / 0.1326 ) ) / ( ( 1 + M_atp_c / 0.0164 ) * ( 1 + M_gln__L_c / 0.354 ) * ( 1 + M_utp_c / 0.071 ) + ( 1 + M_adp_c / 0.1326 ) * ( 1 + M_ctp_c / 0.1326 ) * ( 1 + M_glu__L_c / 0.1326 ) * ( 1 + M_pi_c / 0.1326 ) - 1 ) )

        # (CTPSDUMP) Reaction CTPSDUMP
        cdef double V64 = self.params[127] * self.params[126] * ( ( 4.9027 * ( M_atp_c / 0.0164 ) * ( M_dump_c / 0.0351 ) * ( M_gln__L_c / 0.0667 ) - 0.027698870056497 * ( M_adp_c / 0.1326 ) * ( M_dcmp_c / 0.1326 ) * ( M_glu__L_c / 0.1326 ) * ( M_pi_c / 0.1326 ) ) / ( ( 1 + M_atp_c / 0.0164 ) * ( 1 + M_dump_c / 0.0351 ) * ( 1 + M_gln__L_c / 0.0667 ) + ( 1 + M_adp_c / 0.1326 ) * ( 1 + M_dcmp_c / 0.1326 ) * ( 1 + M_glu__L_c / 0.1326 ) * ( 1 + M_pi_c / 0.1326 ) - 1 ) )

        # (DCMPDA) Reaction DCMPDA
        cdef double V65 = self.params[129] * self.params[128] * ( ( 69 * ( M_dcmp_c / 1 ) - 0.0041 * ( M_dump_c / 0.1259 ) * ( M_nh3_c / 0.1259 ) ) / ( ( 1 + M_dcmp_c / 1 ) + ( 1 + M_dump_c / 0.1259 ) * ( 1 + M_nh3_c / 0.1259 ) - 1 ) )

        # (DUTPDP) Reaction DUTPDP
        cdef double V66 = self.params[131] * self.params[130] * ( ( 12.4805 * ( M_dutp_c / 0.0055 ) - 0 * ( M_dump_c / 0.2224 ) * ( M_ppi_c / 0.2224 ) ) / ( ( 1 + M_dutp_c / 0.0055 ) + ( 1 + M_dump_c / 0.2224 ) * ( 1 + M_ppi_c / 0.2224 ) - 1 ) )

        # (NTD5) Reaction NTD5
        cdef double V67 = self.params[133] * self.params[132] * ( ( 1.1 * ( M_dtmp_c / 0.23 ) - 0.2985 * ( M_pi_c / 0.1 ) * ( M_thymd_c / 0.1 ) ) / ( ( 1 + M_dtmp_c / 0.23 ) + ( 1 + M_pi_c / 0.1 ) * ( 1 + M_thymd_c / 0.1 ) - 1 ) )

        # (NTD1) Reaction NTD1
        cdef double V68 = self.params[135] * self.params[134] * ( ( 0.23 * ( M_dump_c / 0.2 ) - 0.2027 * ( M_duri_c / 0.1 ) * ( M_pi_c / 0.1 ) ) / ( ( 1 + M_dump_c / 0.2 ) + ( 1 + M_duri_c / 0.1 ) * ( 1 + M_pi_c / 0.1 ) - 1 ) )

        # (NTD8) Reaction NTD8
        cdef double V69 = self.params[137] * self.params[136] * ( ( 1.2 * ( M_dgmp_c / 0.18 ) - 0.0056 * ( M_dgsn_c / 0.1574 ) * ( M_pi_c / 0.1574 ) ) / ( ( 1 + M_dgmp_c / 0.18 ) + ( 1 + M_dgsn_c / 0.1574 ) * ( 1 + M_pi_c / 0.1574 ) - 1 ) )

        # (NTD6) Reaction NTD6
        cdef double V70 = self.params[139] * self.params[138] * ( ( 1 * ( M_damp_c / 0.32 ) - 0.0158 * ( M_dad_2_c / 0.1341 ) * ( M_pi_c / 0.1341 ) ) / ( ( 1 + M_damp_c / 0.32 ) + ( 1 + M_dad_2_c / 0.1341 ) * ( 1 + M_pi_c / 0.1341 ) - 1 ) )

        # (PUNP5) Reaction PUNP5
        cdef double V71 = self.params[141] * self.params[140] * ( ( 26.123325 * ( M_pi_c / 2.1566 ) * ( M_uri_c / 0.115 ) - 109.8933 * ( M_ura_c / 2.8972 ) * ( M_r1p_c / 0.0137 ) ) / ( ( 1 + M_pi_c / 2.1566 ) * ( 1 + M_uri_c / 0.115 ) + ( 1 + M_ura_c / 2.8972 ) * ( 1 + M_r1p_c / 0.0137 ) - 1 ) )

        # (GLYK) Reaction GLYK
        cdef double V72 = self.params[143] * self.params[142] * ( ( 3280.3273 * ( M_atp_c / 0.3337 ) * ( M_glyc_c / 0.1345 ) - 93.0643 * ( M_adp_c / 0.1719 ) * ( M_glyc3p_c / 3.9462 ) ) / ( ( 1 + M_atp_c / 0.3337 ) * ( 1 + M_glyc_c / 0.1345 ) + ( 1 + M_adp_c / 0.1719 ) * ( 1 + M_glyc3p_c / 3.9462 ) - 1 ) )

        # (ACPS) Reaction ACPS
        cdef double V73 = self.params[145] * self.params[144] * ( ( 0.7524 * ( M_apoACP_c / 0.0041 ) * ( M_coa_c / 0.0161 ) - 0.0482 * ( M_ACP_c / 0.6226 ) * ( M_pap_c / 0.6226 ) ) / ( ( 1 + M_apoACP_c / 0.0041 ) * ( 1 + M_coa_c / 0.0161 ) + ( 1 + M_ACP_c / 0.6226 ) * ( 1 + M_pap_c / 0.6226 ) - 1 ) )

        # (BPNT) Reaction BPNT
        cdef double V74 = self.params[147] * self.params[146] * ( ( 7.1703 * ( M_pap_c / 0.025 ) - 0.3301 * ( M_amp_c / 0.2965 ) * ( M_pi_c / 0.2965 ) ) / ( ( 1 + M_pap_c / 0.025 ) + ( 1 + M_amp_c / 0.2965 ) * ( 1 + M_pi_c / 0.2965 ) - 1 ) )

        # (FAKr) Reaction FAKr
        cdef double V75 = self.params[149] * self.params[148] * ( ( 74.8817 * ( M_atp_c / 1.4528 ) * ( M_fa_c / 0.3014 ) - 228.2947 * ( M_adp_c / 0.0332 ) * ( M_ap_c / 0.0332 ) ) / ( ( 1 + M_atp_c / 1.4528 ) * ( 1 + M_fa_c / 0.3014 ) + ( 1 + M_adp_c / 0.0332 ) * ( 1 + M_ap_c / 0.0332 ) - 1 ) )

        # (ACPPAT) Reaction ACPPAT
        cdef double V76 = self.params[151] * self.params[150] * ( ( 20.7 * ( M_ACP_c / 0.0179 ) * ( M_ap_c / 0.0179 ) - 0.0494 * ( M_ACP_R_c / 0.5602 ) * ( M_pi_c / 0.5602 ) ) / ( ( 1 + M_ACP_c / 0.0179 ) * ( 1 + M_ap_c / 0.0179 ) + ( 1 + M_ACP_R_c / 0.5602 ) * ( 1 + M_pi_c / 0.5602 ) - 1 ) )

        # (APG3PAT) Reaction APG3PAT
        cdef double V77 = self.params[153] * self.params[152] * ( ( 8.4989 * ( M_ap_c / 0.0206 ) * ( M_glyc3p_c / 0.0206 ) - 0.0658 * ( M_1ag3p_c / 0.4856 ) * ( M_pi_c / 0.4856 ) ) / ( ( 1 + M_ap_c / 0.0206 ) * ( 1 + M_glyc3p_c / 0.0206 ) + ( 1 + M_1ag3p_c / 0.4856 ) * ( 1 + M_pi_c / 0.4856 ) - 1 ) )

        # (AGPAT) Reaction AGPAT
        cdef double V78 = self.params[155] * self.params[154] * ( ( 144.4 * ( M_1ag3p_c / 0.0064 ) * ( M_ACP_R_c / 0.0247 ) - 0 * ( M_ACP_c / 0.4044 ) * ( M_pa_c / 0.4044 ) ) / ( ( 1 + M_1ag3p_c / 0.0064 ) * ( 1 + M_ACP_R_c / 0.0247 ) + ( 1 + M_ACP_c / 0.4044 ) * ( 1 + M_pa_c / 0.4044 ) - 1 ) )

        # (DASYN) Reaction DASYN
        cdef double V79 = self.params[157] * self.params[156] * ( ( 6.4829 * ( M_ctp_c / 0.164 ) * ( M_pa_c / 0.0758 ) - 3.3541 * ( M_cdpdag_c / 0.1446 ) * ( M_ppi_c / 0.1446 ) ) / ( ( 1 + M_ctp_c / 0.164 ) * ( 1 + M_pa_c / 0.0758 ) + ( 1 + M_cdpdag_c / 0.1446 ) * ( 1 + M_ppi_c / 0.1446 ) - 1 ) )

        # (PGSA) Reaction PGSA
        cdef double V80 = self.params[159] * self.params[158] * ( ( 4.0472 * ( M_cdpdag_c / 0.0324 ) * ( M_glyc3p_c / 0.0509 ) - 0 * ( M_cmp_c / 0.0416 ) * ( M_pg3p_c / 0.0536 ) ) / ( ( 1 + M_cdpdag_c / 0.0324 ) * ( 1 + M_glyc3p_c / 0.0509 ) + ( 1 + M_cmp_c / 0.0416 ) * ( 1 + M_pg3p_c / 0.0536 ) - 1 ) )

        # (PGPP) Reaction PGPP
        cdef double V81 = self.params[161] * self.params[160] * ( ( 26.212 * ( M_pg3p_c / 0.0011 ) - 0 * ( M_pg_c / 46.068 ) * ( M_pi_c / 46.068 ) ) / ( ( 1 + M_pg3p_c / 0.0011 ) + ( 1 + M_pg_c / 46.068 ) * ( 1 + M_pi_c / 46.068 ) - 1 ) )

        # (CLPNS) Reaction CLPNS
        cdef double V82 = self.params[163] * self.params[162] * ( ( 120.0193 * ( M_pg_c / 2.9902 ) * ( M_pg_c / 2.9902 ) - 0.0005 * ( M_clpn_c / 1.8783 ) * ( M_glyc_c / 1.8783 ) ) / ( ( 1 + M_pg_c / 2.9902 ) * ( 1 + M_pg_c / 2.9902 ) + ( 1 + M_clpn_c / 1.8783 ) * ( 1 + M_glyc_c / 1.8783 ) - 1 ) )

        # (PAPA) Reaction PAPA
        cdef double V83 = self.params[165] * self.params[164] * ( ( 44.3467 * ( M_pa_c / 0.052 ) - 0.1399 * ( M_12dgr_c / 0.3459 ) * ( M_pi_c / 0.3459 ) ) / ( ( 1 + M_pa_c / 0.052 ) + ( 1 + M_12dgr_c / 0.3459 ) * ( 1 + M_pi_c / 0.3459 ) - 1 ) )

        # (PGMT) Reaction PGMT
        cdef double V84 = self.params[167] * self.params[166] * ( ( 0.668 * ( M_g6p_c / 2.8402 ) - 4.5294 * ( M_g1p_c / 6.6039 ) ) / ( ( 1 + M_g6p_c / 2.8402 ) + ( 1 + M_g1p_c / 6.6039 ) - 1 ) )

        # (GALU) Reaction GALU
        cdef double V85 = self.params[169] * self.params[168] * ( ( 36.5547 * ( M_g1p_c / 0.1659 ) * ( M_utp_c / 0.0903 ) - 5.8538 * ( M_ppi_c / 0.0377 ) * ( M_udpg_c / 0.0459 ) ) / ( ( 1 + M_g1p_c / 0.1659 ) * ( 1 + M_utp_c / 0.0903 ) + ( 1 + M_ppi_c / 0.0377 ) * ( 1 + M_udpg_c / 0.0459 ) - 1 ) )

        # (UDPG4E) Reaction UDPG4E
        cdef double V86 = self.params[171] * self.params[170] * ( ( 449.458 * ( M_udpg_c / 29.0953 ) - 243.4781 * ( M_udpgal_c / 0.1181 ) ) / ( ( 1 + M_udpg_c / 29.0953 ) + ( 1 + M_udpgal_c / 0.1181 ) - 1 ) )

        # (UDPGALM) Reaction UDPGALM
        cdef double V87 = self.params[173] * self.params[172] * ( ( 27 * ( M_udpgal_c / 0.3042 ) - 1.5 * ( M_udpgalfur_c / 0.45 ) ) / ( ( 1 + M_udpgal_c / 0.3042 ) + ( 1 + M_udpgalfur_c / 0.45 ) - 1 ) )

        # (DAGGALT) Reaction DAGGALT
        cdef double V88 = self.params[175] * self.params[174] * ( ( 0.4457 * ( M_12dgr_c / 0.5305 ) * ( M_udpgalfur_c / 0.0589 ) - 0 * ( M_galfur12dgr_c / 0.1699 ) * ( M_udp_c / 0.1699 ) ) / ( ( 1 + M_12dgr_c / 0.5305 ) * ( 1 + M_udpgalfur_c / 0.0589 ) + ( 1 + M_galfur12dgr_c / 0.1699 ) * ( 1 + M_udp_c / 0.1699 ) - 1 ) )

        # (NCTPPRT) Reaction NCTPPRT
        cdef double V89 = self.params[177] * self.params[176] * ( ( 0.0433 * ( M_nac_c / 0.00806 ) * ( M_prpp_c / 0.0038 ) - 0 * ( M_nicrnt_c / 0.0346 ) * ( M_ppi_c / 0.18 ) ) / ( ( 1 + M_nac_c / 0.00806 ) * ( 1 + M_prpp_c / 0.0038 ) + ( 1 + M_nicrnt_c / 0.0346 ) * ( 1 + M_ppi_c / 0.18 ) - 1 ) )

        # (NNATr) Reaction NNATr
        cdef double V90 = self.params[179] * self.params[178] * ( ( 3.6682 * ( M_atp_c / 0.1837 ) * ( M_nicrnt_c / 0.0554 ) - 0.0022 * ( M_dnad_c / 0.1804 ) * ( M_ppi_c / 1.0287 ) ) / ( ( 1 + M_atp_c / 0.1837 ) * ( 1 + M_nicrnt_c / 0.0554 ) + ( 1 + M_dnad_c / 0.1804 ) * ( 1 + M_ppi_c / 1.0287 ) - 1 ) )

        # (NADS) Reaction NADS
        cdef double V91 = self.params[181] * self.params[180] * ( ( 0.0454 * ( M_atp_c / 0.0278 ) * ( M_dnad_c / 0.0465 ) * ( M_nh3_c / 0.0871 ) - 0.0104 * ( M_amp_c / 0.2152 ) * ( M_nad_c / 0.2152 ) * ( M_ppi_c / 0.2152 ) ) / ( ( 1 + M_atp_c / 0.0278 ) * ( 1 + M_dnad_c / 0.0465 ) * ( 1 + M_nh3_c / 0.0871 ) + ( 1 + M_amp_c / 0.2152 ) * ( 1 + M_nad_c / 0.2152 ) * ( 1 + M_ppi_c / 0.2152 ) - 1 ) )

        # (NADK) Reaction NADK
        cdef double V92 = self.params[183] * self.params[182] * ( ( 1.36 * ( M_atp_c / 0.4312 ) * ( M_nad_c / 2 ) - 104.4259 * ( M_adp_c / 4.5824 ) * ( M_nadp_c / 0.18 ) ) / ( ( 1 + M_atp_c / 0.4312 ) * ( 1 + M_nad_c / 2 ) + ( 1 + M_adp_c / 4.5824 ) * ( 1 + M_nadp_c / 0.18 ) - 1 ) )

        # (NADHK) Reaction NADHK
        cdef double V93 = self.params[185] * self.params[184] * ( ( 31.3405 * ( M_atp_c / 0.4054 ) * ( M_nadh_c / 2 ) - 104.2685 * ( M_adp_c / 4.624 ) * ( M_nadph_c / 0.3 ) ) / ( ( 1 + M_atp_c / 0.4054 ) * ( 1 + M_nadh_c / 2 ) + ( 1 + M_adp_c / 4.624 ) * ( 1 + M_nadph_c / 0.3 ) - 1 ) )

        # (RBFK) Reaction RBFK
        cdef double V94 = self.params[187] * self.params[186] * ( ( 1.2796 * ( M_atp_c / 0.1681 ) * ( M_ribflv_c / 0.0094 ) - 0.4249 * ( M_adp_c / 0.1385 ) * ( M_fmn_c / 0.1385 ) ) / ( ( 1 + M_atp_c / 0.1681 ) * ( 1 + M_ribflv_c / 0.0094 ) + ( 1 + M_adp_c / 0.1385 ) * ( 1 + M_fmn_c / 0.1385 ) - 1 ) )

        # (FMNAT) Reaction FMNAT
        cdef double V95 = self.params[189] * self.params[188] * ( ( 0.0027 * ( M_atp_c / 0.0173 ) * ( M_fmn_c / 0.0294 ) - 0 * ( M_fad_c / 0.0027 ) * ( M_ppi_c / 0.0281 ) ) / ( ( 1 + M_atp_c / 0.0173 ) * ( 1 + M_fmn_c / 0.0294 ) + ( 1 + M_fad_c / 0.0027 ) * ( 1 + M_ppi_c / 0.0281 ) - 1 ) )

        # (5FTHFPGS) Reaction 5FTHFPGS
        cdef double V96 = self.params[191] * self.params[190] * ( ( 2 * ( M_5fthf_c / 0.0999 ) * ( M_atp_c / 0.23 ) * ( M_atp_c / 0.23 ) * ( M_glu__L_c / 0.3 ) * ( M_glu__L_c / 0.3 ) - 0 * ( M_5fthfglu3_c / 0.1001 ) * ( M_adp_c / 0.1001 ) * ( M_adp_c / 0.1001 ) * ( M_pi_c / 0.1001 ) * ( M_pi_c / 0.1001 ) ) / ( ( 1 + M_5fthf_c / 0.0999 ) * ( 1 + M_atp_c / 0.23 ) * ( 1 + M_atp_c / 0.23 ) * ( 1 + M_glu__L_c / 0.3 ) * ( 1 + M_glu__L_c / 0.3 ) + ( 1 + M_5fthfglu3_c / 0.1001 ) * ( 1 + M_adp_c / 0.1001 ) * ( 1 + M_adp_c / 0.1001 ) * ( 1 + M_pi_c / 0.1001 ) * ( 1 + M_pi_c / 0.1001 ) - 1 ) )

        # (MTHFC) Reaction MTHFC
        cdef double V97 = self.params[193] * self.params[192] * ( ( 185 * ( M_methfglu3_c / 0.1207 ) - 15.1817 * ( M_10fthfglu3_c / 0.0829 ) ) / ( ( 1 + M_methfglu3_c / 0.1207 ) + ( 1 + M_10fthfglu3_c / 0.0829 ) - 1 ) )

        # (FMETTRS) Reaction FMETTRS
        cdef double V98 = self.params[195] * self.params[194] * ( ( 0.45 * ( M_10fthfglu3_c / 0.011 ) * ( M_mettrna_c / 0.01 ) - 0 * ( M_fmettrna_c / 0.001 ) * ( M_thfglu3_c / 0.1 ) ) / ( ( 1 + M_10fthfglu3_c / 0.011 ) * ( 1 + M_mettrna_c / 0.01 ) + ( 1 + M_fmettrna_c / 0.001 ) * ( 1 + M_thfglu3_c / 0.1 ) - 1 ) )

        # (GHMT) Reaction GHMT
        cdef double V99 = self.params[197] * self.params[196] * ( ( 640 * ( M_ser__L_c / 2.3312 ) * ( M_thfglu3_c / 0.0684 ) - 23.173 * ( M_gly_c / 0.1463 ) * ( M_mlthfglu3_c / 0.1463 ) ) / ( ( 1 + M_ser__L_c / 2.3312 ) * ( 1 + M_thfglu3_c / 0.0684 ) + ( 1 + M_gly_c / 0.1463 ) * ( 1 + M_mlthfglu3_c / 0.1463 ) - 1 ) )

        # (MTHFD) Reaction MTHFD
        cdef double V100 = self.params[199] * self.params[198] * ( ( 19.4491 * ( M_mlthfglu3_c / 0.1738 ) * ( M_nadp_c / 0.0303 ) - 3.2 * ( M_methfglu3_c / 0.0575 ) * ( M_nadph_c / 0.0575 ) ) / ( ( 1 + M_mlthfglu3_c / 0.1738 ) * ( 1 + M_nadp_c / 0.0303 ) + ( 1 + M_methfglu3_c / 0.0575 ) * ( 1 + M_nadph_c / 0.0575 ) - 1 ) )

        # (FTHFCL) Reaction FTHFCL
        cdef double V101 = self.params[201] * self.params[200] * ( ( 3.0655 * ( M_5fthfglu3_c / 0.0191 ) * ( M_atp_c / 0.0164 ) - 0 * ( M_adp_c / 0.5231 ) * ( M_methfglu3_c / 0.5231 ) * ( M_pi_c / 0.5231 ) ) / ( ( 1 + M_5fthfglu3_c / 0.0191 ) * ( 1 + M_atp_c / 0.0164 ) + ( 1 + M_adp_c / 0.5231 ) * ( 1 + M_methfglu3_c / 0.5231 ) * ( 1 + M_pi_c / 0.5231 ) - 1 ) )

        # (GHMT2) Reaction GHMT2
        cdef double V102 = self.params[203] * self.params[202] * ( ( 0 * ( M_methfglu3_c / 0.0684 ) - 0 * ( M_5fthfglu3_c / 0.1463 ) ) / ( ( 1 + M_methfglu3_c / 0.0684 ) + ( 1 + M_5fthfglu3_c / 0.1463 ) - 1 ) )

        # (ATPase) Reaction ATPase
        cdef double V103 = self.params[205] * self.params[204] * ( ( 20 * ( M_adp_c / 0.1 ) * ( M_pi_c / 4.2 ) - 72.3333333333333 * ( M_atp_c / 0.6 ) ) / ( ( 1 + M_adp_c / 0.1 ) * ( 1 + M_pi_c / 4.2 ) + ( 1 + M_atp_c / 0.6 ) - 1 ) )

        # (RIBFLVabc) Reaction RIBFLVabc
        cdef double V104 = self.params[207] * self.params[206] * ( ( 1.66 * ( M_atp_c / 0.1 ) * ( 0.003013 / 0.0019 ) - 0 * ( M_adp_c / 0.1 ) * ( M_pi_c / 0.1 ) * ( M_ribflv_c / 0.0019 ) ) / ( ( 1 + M_atp_c / 0.1 ) * ( 1 + 0.003013 / 0.0019 ) + ( 1 + M_adp_c / 0.1 ) * ( 1 + M_pi_c / 0.1 ) * ( 1 + M_ribflv_c / 0.0019 ) - 1 ) )

        # (P5Pabc) Reaction P5Pabc
        cdef double V105 = self.params[209] * self.params[208] * ( ( 1.66 * ( M_atp_c / 0.1 ) * ( 0.006 / 0.0019 ) - 0 * ( M_adp_c / 0.1 ) * ( M_pi_c / 0.1 ) * ( M_pydx5p_c / 0.0019 ) ) / ( ( 1 + M_atp_c / 0.1 ) * ( 1 + 0.006 / 0.0019 ) + ( 1 + M_adp_c / 0.1 ) * ( 1 + M_pi_c / 0.1 ) * ( 1 + M_pydx5p_c / 0.0019 ) - 1 ) )

        # (5FTHFabc) Reaction 5FTHFabc
        cdef double V106 = self.params[211] * self.params[210] * ( ( 1.66 * ( 0.05002 / 0.0019 ) * ( M_atp_c / 0.1 ) - 0 * ( M_5fthf_c / 0.0019 ) * ( M_adp_c / 0.1 ) * ( M_pi_c / 0.1 ) ) / ( ( 1 + 0.05002 / 0.0019 ) * ( 1 + M_atp_c / 0.1 ) + ( 1 + M_5fthf_c / 0.0019 ) * ( 1 + M_adp_c / 0.1 ) * ( 1 + M_pi_c / 0.1 ) - 1 ) )

        # (NACabc) Reaction NACabc
        cdef double V107 = self.params[213] * self.params[212] * ( ( 1.66 * ( M_atp_c / 0.1 ) * ( 0.0161 / 0.0019 ) - 0 * ( M_adp_c / 0.1 ) * ( M_nac_c / 0.0019 ) * ( M_pi_c / 0.1 ) ) / ( ( 1 + M_atp_c / 0.1 ) * ( 1 + 0.0161 / 0.0019 ) + ( 1 + M_adp_c / 0.1 ) * ( 1 + M_nac_c / 0.0019 ) * ( 1 + M_pi_c / 0.1 ) - 1 ) )

        # (COAabc) Reaction COAabc
        cdef double V108 = self.params[215] * self.params[214] * ( ( 1.66 * ( M_atp_c / 0.1 ) * ( 0.002928 / 0.0019 ) - 0 * ( M_adp_c / 0.1 ) * ( M_coa_c / 0.0019 ) * ( M_pi_c / 0.1 ) ) / ( ( 1 + M_atp_c / 0.1 ) * ( 1 + 0.002928 / 0.0019 ) + ( 1 + M_adp_c / 0.1 ) * ( 1 + M_coa_c / 0.0019 ) * ( 1 + M_pi_c / 0.1 ) - 1 ) )

        # (THMPPabc) Reaction THMPPabc
        cdef double V109 = self.params[217] * self.params[216] * ( ( 1.66 * ( M_atp_c / 0.1 ) * ( 0.008015 / 0.0019 ) - 0 * ( M_adp_c / 0.1 ) * ( M_pi_c / 0.1 ) * ( M_thmpp_c / 0.0019 ) ) / ( ( 1 + M_atp_c / 0.1 ) * ( 1 + 0.008015 / 0.0019 ) + ( 1 + M_adp_c / 0.1 ) * ( 1 + M_pi_c / 0.1 ) * ( 1 + M_thmpp_c / 0.0019 ) - 1 ) )

        # (SPRMabc) Reaction SPRMabc
        cdef double V110 = self.params[219] * self.params[218] * ( ( 3 * ( M_atp_c / 0.1 ) * ( 0.1 / 2 ) - 0 * ( M_adp_c / 0.1 ) * ( M_pi_c / 0.1 ) * ( M_sprm_c / 2 ) ) / ( ( 1 + M_atp_c / 0.1 ) * ( 1 + 0.1 / 2 ) + ( 1 + M_adp_c / 0.1 ) * ( 1 + M_pi_c / 0.1 ) * ( 1 + M_sprm_c / 2 ) - 1 ) )

        # (DCYTabc) Reaction DCYTabc
        cdef double V111 = self.params[221] * self.params[220] * ( ( 0.5 * ( M_atp_c / 0.023 ) * ( 0.02203 / 0.0011 ) - 0 * ( M_adp_c / 2.8 ) * ( M_dcyt_c / 0.02 ) * ( M_pi_c / 10 ) ) / ( ( 1 + M_atp_c / 0.023 ) * ( 1 + 0.02203 / 0.0011 ) + ( 1 + M_adp_c / 2.8 ) * ( 1 + M_dcyt_c / 0.02 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (URIabc) Reaction URIabc
        cdef double V112 = self.params[223] * self.params[222] * ( ( 1.5 * ( M_atp_c / 0.023 ) * ( 0.18 / 0.0017 ) - 0 * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) * ( M_uri_c / 0.02 ) ) / ( ( 1 + M_atp_c / 0.023 ) * ( 1 + 0.18 / 0.0017 ) + ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) * ( 1 + M_uri_c / 0.02 ) - 1 ) )

        # (THMDabc) Reaction THMDabc
        cdef double V113 = self.params[225] * self.params[224] * ( ( 1 * ( M_atp_c / 0.023 ) * ( 0.1007 / 0.0017 ) - 0 * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) * ( M_thymd_c / 0.02 ) ) / ( ( 1 + M_atp_c / 0.023 ) * ( 1 + 0.1007 / 0.0017 ) + ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) * ( 1 + M_thymd_c / 0.02 ) - 1 ) )

        # (ADNabc) Reaction ADNabc
        cdef double V114 = self.params[227] * self.params[226] * ( ( 1.5 * ( 0.15 / 0.0019 ) * ( M_atp_c / 0.023 ) - 0 * ( M_adn_c / 0.02 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 0.15 / 0.0019 ) * ( 1 + M_atp_c / 0.023 ) + ( 1 + M_adn_c / 0.02 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (DADNabc) Reaction DADNabc
        cdef double V115 = self.params[229] * self.params[228] * ( ( 1 * ( M_atp_c / 0.023 ) * ( 0.01992 / 0.0019 ) - 0 * ( M_adp_c / 2.8 ) * ( M_dad_2_c / 0.02 ) * ( M_pi_c / 10 ) ) / ( ( 1 + M_atp_c / 0.023 ) * ( 1 + 0.01992 / 0.0019 ) + ( 1 + M_adp_c / 2.8 ) * ( 1 + M_dad_2_c / 0.02 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (GSNabc) Reaction GSNabc
        cdef double V116 = self.params[231] * self.params[230] * ( ( 0.5 * ( M_atp_c / 0.023 ) * ( 0.13 / 0.0019 ) - 0 * ( M_adp_c / 2.8 ) * ( M_gsn_c / 0.02 ) * ( M_pi_c / 10 ) ) / ( ( 1 + M_atp_c / 0.023 ) * ( 1 + 0.13 / 0.0019 ) + ( 1 + M_adp_c / 2.8 ) * ( 1 + M_gsn_c / 0.02 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (DGSNabc) Reaction DGSNabc
        cdef double V117 = self.params[233] * self.params[232] * ( ( 0.5 * ( M_atp_c / 0.023 ) * ( 0.01873 / 0.0019 ) - 0 * ( M_adp_c / 2.8 ) * ( M_dgsn_c / 0.02 ) * ( M_pi_c / 10 ) ) / ( ( 1 + M_atp_c / 0.023 ) * ( 1 + 0.01873 / 0.0019 ) + ( 1 + M_adp_c / 2.8 ) * ( 1 + M_dgsn_c / 0.02 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (Kt6) Reaction Kt6
        cdef double V118 = self.params[235] * self.params[234] * ( ( 3 * ( M_atp_c / 0.03 ) * ( 12.67 / 0.46 ) * ( M_na1_c / 7.47 ) - 0 * ( M_adp_c / 1.5 ) * ( M_k_c / 1.9 ) * ( 19.38 / 12.7 ) * ( M_pi_c / 10 ) ) / ( ( 1 + M_atp_c / 0.03 ) * ( 1 + 12.67 / 0.46 ) * ( 1 + M_na1_c / 7.47 ) + ( 1 + M_adp_c / 1.5 ) * ( 1 + M_k_c / 1.9 ) * ( 1 + 19.38 / 12.7 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (MG2abc) Reaction MG2abc
        cdef double V119 = self.params[237] * self.params[236] * ( ( 22 * ( M_atp_c / 0.1 ) * ( 1.407 / 0.05 ) - 0 * ( M_adp_c / 2.8 ) * ( M_mg2_c / 0.05 ) * ( M_pi_c / 10 ) ) / ( ( 1 + M_atp_c / 0.1 ) * ( 1 + 1.407 / 0.05 ) + ( 1 + M_adp_c / 2.8 ) * ( 1 + M_mg2_c / 0.05 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (CA2abc) Reaction CA2abc
        cdef double V120 = self.params[239] * self.params[238] * ( ( 9.5 * ( M_atp_c / 0.075 ) * ( 0.6803 / 0.0075 ) - 0 * ( M_adp_c / 0.1 ) * ( M_ca2_c / 0.0075 ) * ( M_pi_c / 10 ) ) / ( ( 1 + M_atp_c / 0.075 ) * ( 1 + 0.6803 / 0.0075 ) + ( 1 + M_adp_c / 0.1 ) * ( 1 + M_ca2_c / 0.0075 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (PIabc) Reaction PIabc
        cdef double V121 = self.params[241] * self.params[240] * ( ( 30 * ( M_atp_c / 0.023 ) * ( 134.0 / 0.0031 ) - 0 * ( M_adp_c / 0.654 ) * ( M_pi_c / 0.385 ) * ( M_pi_c / 0.385 ) ) / ( ( 1 + M_atp_c / 0.023 ) * ( 1 + 134.0 / 0.0031 ) + ( 1 + M_adp_c / 0.654 ) * ( 1 + M_pi_c / 0.385 ) * ( 1 + M_pi_c / 0.385 ) - 1 ) )

        # (ARGabc) Reaction ARGabc
        cdef double V122 = self.params[243] * self.params[242] * ( ( 0.2 * ( 4.166 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_arg__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.166 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_arg__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (ASPabc) Reaction ASPabc
        cdef double V123 = self.params[245] * self.params[244] * ( ( 0.2 * ( 0.1128 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_asp__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 0.1128 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_asp__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (GLUabc) Reaction GLUabc
        cdef double V124 = self.params[247] * self.params[246] * ( ( 0.2 * ( 4.255 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_glu__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.255 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_glu__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (GLYabc) Reaction GLYabc
        cdef double V125 = self.params[249] * self.params[248] * ( ( 0.2 * ( 4.333 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_gly_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.333 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_gly_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (ILEabc) Reaction ILEabc
        cdef double V126 = self.params[251] * self.params[250] * ( ( 0.2 * ( 4.076 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_ile__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.076 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_ile__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (ALAabc) Reaction ALAabc
        cdef double V127 = self.params[253] * self.params[252] * ( ( 0.2 * ( 4.14 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_ala__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.14 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_ala__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (ASNabc) Reaction ASNabc
        cdef double V128 = self.params[255] * self.params[254] * ( ( 0.2 * ( 4.0 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_asn__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.0 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_asn__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (LEUabc) Reaction LEUabc
        cdef double V129 = self.params[257] * self.params[256] * ( ( 0.2 * ( 4.229 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_leu__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.229 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_leu__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (GLNabc) Reaction GLNabc
        cdef double V130 = self.params[259] * self.params[258] * ( ( 0.2 * ( 2.0 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_gln__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 2.0 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_gln__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (HISabc) Reaction HISabc
        cdef double V131 = self.params[261] * self.params[260] * ( ( 0.2 * ( 4.048 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_his__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.048 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_his__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (LYSabc) Reaction LYSabc
        cdef double V132 = self.params[263] * self.params[262] * ( ( 0.2 * ( 4.191 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_lys__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.191 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_lys__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (PROabc) Reaction PROabc
        cdef double V133 = self.params[265] * self.params[264] * ( ( 0.2 * ( 4.174 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_pro__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.174 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_pro__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (PHEabc) Reaction PHEabc
        cdef double V134 = self.params[267] * self.params[266] * ( ( 0.2 * ( 1.076 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_phe__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 1.076 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_phe__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (THRabc) Reaction THRabc
        cdef double V135 = self.params[269] * self.params[268] * ( ( 0.2 * ( 1.126 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_thr__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 1.126 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_thr__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (TRPabc) Reaction TRPabc
        cdef double V136 = self.params[271] * self.params[270] * ( ( 0.2 * ( 0.5245 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_trp__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 0.5245 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_trp__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (TYRabc) Reaction TYRabc
        cdef double V137 = self.params[273] * self.params[272] * ( ( 0.2 * ( 4.11 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_tyr__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 4.11 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_tyr__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (VALabc) Reaction VALabc
        cdef double V138 = self.params[275] * self.params[274] * ( ( 0.2 * ( 8.107 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_val__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 8.107 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_val__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (SERabc) Reaction SERabc
        cdef double V139 = self.params[277] * self.params[276] * ( ( 0.2 * ( 1.119 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_ser__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 1.119 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_ser__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (METabc) Reaction METabc
        cdef double V140 = self.params[279] * self.params[278] * ( ( 0.2 * ( 8.05 / 0.03 ) * ( M_atp_c / 2.8 ) - 0 * ( M_met__L_c / 0.03 ) * ( M_adp_c / 2.8 ) * ( M_pi_c / 10 ) ) / ( ( 1 + 8.05 / 0.03 ) * ( 1 + M_atp_c / 2.8 ) + ( 1 + M_met__L_c / 0.03 ) * ( 1 + M_adp_c / 2.8 ) * ( 1 + M_pi_c / 10 ) - 1 ) )

        # (ARGt2r) Reaction ARGt2r
        cdef double V141 = self.params[281] * self.params[280] * ( ( 3.4 * ( 4.166 / 0.03 ) - 0 * ( M_arg__L_c / 0.03 ) ) / ( ( 1 + 4.166 / 0.03 ) + ( 1 + M_arg__L_c / 0.03 ) - 1 ) )

        # (ASPt2pr) Reaction ASPt2pr
        cdef double V142 = self.params[283] * self.params[282] * ( ( 9.1 * ( 0.1128 / 0.03 ) - 0 * ( M_asp__L_c / 0.03 ) ) / ( ( 1 + 0.1128 / 0.03 ) + ( 1 + M_asp__L_c / 0.03 ) - 1 ) )

        # (CYSt2r) Reaction CYSt2r
        cdef double V143 = self.params[285] * self.params[284] * ( ( 3.4 * ( 3.168 / 0.03 ) - 0 * ( M_cys__L_c / 0.03 ) ) / ( ( 1 + 3.168 / 0.03 ) + ( 1 + M_cys__L_c / 0.03 ) - 1 ) )

        # (GLUt2pr) Reaction GLUt2pr
        cdef double V144 = self.params[287] * self.params[286] * ( ( 3.4 * ( 4.255 / 0.03 ) - 0 * ( M_glu__L_c / 0.03 ) ) / ( ( 1 + 4.255 / 0.03 ) + ( 1 + M_glu__L_c / 0.03 ) - 1 ) )

        # (GLYt2r) Reaction GLYt2r
        cdef double V145 = self.params[289] * self.params[288] * ( ( 5.1 * ( 4.333 / 0.03 ) - 0 * ( M_gly_c / 0.03 ) ) / ( ( 1 + 4.333 / 0.03 ) + ( 1 + M_gly_c / 0.03 ) - 1 ) )

        # (ISOt2r) Reaction ISOt2r
        cdef double V146 = self.params[291] * self.params[290] * ( ( 9.1 * ( 4.076 / 0.03 ) - 0 * ( M_ile__L_c / 0.03 ) ) / ( ( 1 + 4.076 / 0.03 ) + ( 1 + M_ile__L_c / 0.03 ) - 1 ) )

        # (ALAt2r) Reaction ALAt2r
        cdef double V147 = self.params[293] * self.params[292] * ( ( 5.1 * ( 4.14 / 0.03 ) - 0 * ( M_ala__L_c / 0.03 ) ) / ( ( 1 + 4.14 / 0.03 ) + ( 1 + M_ala__L_c / 0.03 ) - 1 ) )

        # (ASNt2r) Reaction ASNt2r
        cdef double V148 = self.params[295] * self.params[294] * ( ( 6.1 * ( 4.0 / 0.03 ) - 0 * ( M_asn__L_c / 0.03 ) ) / ( ( 1 + 4.0 / 0.03 ) + ( 1 + M_asn__L_c / 0.03 ) - 1 ) )

        # (LEUt2r) Reaction LEUt2r
        cdef double V149 = self.params[297] * self.params[296] * ( ( 9.1 * ( 4.229 / 0.03 ) - 0 * ( M_leu__L_c / 0.03 ) ) / ( ( 1 + 4.229 / 0.03 ) + ( 1 + M_leu__L_c / 0.03 ) - 1 ) )

        # (GLNt2r) Reaction GLNt2r
        cdef double V150 = self.params[299] * self.params[298] * ( ( 7.1 * ( 2.0 / 0.03 ) - 0 * ( M_gln__L_c / 0.03 ) ) / ( ( 1 + 2.0 / 0.03 ) + ( 1 + M_gln__L_c / 0.03 ) - 1 ) )

        # (HISt2r) Reaction HISt2r
        cdef double V151 = self.params[301] * self.params[300] * ( ( 3.4 * ( 4.048 / 0.03 ) - 0 * ( M_his__L_c / 0.03 ) ) / ( ( 1 + 4.048 / 0.03 ) + ( 1 + M_his__L_c / 0.03 ) - 1 ) )

        # (LYSt2r) Reaction LYSt2r
        cdef double V152 = self.params[303] * self.params[302] * ( ( 9.1 * ( 4.191 / 0.03 ) - 0 * ( M_lys__L_c / 0.03 ) ) / ( ( 1 + 4.191 / 0.03 ) + ( 1 + M_lys__L_c / 0.03 ) - 1 ) )

        # (PROt2r) Reaction PROt2r
        cdef double V153 = self.params[305] * self.params[304] * ( ( 3.4 * ( 4.174 / 0.03 ) - 0 * ( M_pro__L_c / 0.03 ) ) / ( ( 1 + 4.174 / 0.03 ) + ( 1 + M_pro__L_c / 0.03 ) - 1 ) )

        # (PHEt2r) Reaction PHEt2r
        cdef double V154 = self.params[307] * self.params[306] * ( ( 3.4 * ( 1.076 / 0.03 ) - 0 * ( M_phe__L_c / 0.03 ) ) / ( ( 1 + 1.076 / 0.03 ) + ( 1 + M_phe__L_c / 0.03 ) - 1 ) )

        # (THRt2r) Reaction THRt2r
        cdef double V155 = self.params[309] * self.params[308] * ( ( 5.1 * ( 1.126 / 0.03 ) - 0 * ( M_thr__L_c / 0.03 ) ) / ( ( 1 + 1.126 / 0.03 ) + ( 1 + M_thr__L_c / 0.03 ) - 1 ) )

        # (TRPt2r) Reaction TRPt2r
        cdef double V156 = self.params[311] * self.params[310] * ( ( 3.4 * ( 0.5245 / 0.03 ) - 0 * ( M_trp__L_c / 0.03 ) ) / ( ( 1 + 0.5245 / 0.03 ) + ( 1 + M_trp__L_c / 0.03 ) - 1 ) )

        # (TYRt2r) Reaction TYRt2r
        cdef double V157 = self.params[313] * self.params[312] * ( ( 3.4 * ( 4.11 / 0.03 ) - 0 * ( M_tyr__L_c / 0.03 ) ) / ( ( 1 + 4.11 / 0.03 ) + ( 1 + M_tyr__L_c / 0.03 ) - 1 ) )

        # (VALt2r) Reaction VALt2r
        cdef double V158 = self.params[315] * self.params[314] * ( ( 5.1 * ( 8.107 / 0.03 ) - 0 * ( M_val__L_c / 0.03 ) ) / ( ( 1 + 8.107 / 0.03 ) + ( 1 + M_val__L_c / 0.03 ) - 1 ) )

        # (SERt2r) Reaction SERt2r
        cdef double V159 = self.params[317] * self.params[316] * ( ( 6.1 * ( 1.119 / 0.03 ) - 0 * ( M_ser__L_c / 0.03 ) ) / ( ( 1 + 1.119 / 0.03 ) + ( 1 + M_ser__L_c / 0.03 ) - 1 ) )

        # (METt2r) Reaction METt2r
        cdef double V160 = self.params[319] * self.params[318] * ( ( 3.4 * ( 8.05 / 0.03 ) - 0 * ( M_met__L_c / 0.03 ) ) / ( ( 1 + 8.05 / 0.03 ) + ( 1 + M_met__L_c / 0.03 ) - 1 ) )

        # (GLCpts0) Reaction GLCpts0
        cdef double V161 = (10000 * ptsi * M_pep_c) - (4000 * ptsi_P * M_pyr_c)

        # (GLCpts1) Reaction GLCpts1
        cdef double V162 = (200000 * ptsh * ptsi_P) - (8000 * ptsh_P * ptsi)

        # (GLCpts2) Reaction GLCpts2
        cdef double V163 = (61000 * crr * ptsh_P) - (47000 * crr_P * ptsh)

        # (GLCpts3) Reaction GLCpts3
        cdef double V164 = (3900 * ptsg * crr_P) - (310 * ptsg_P * crr)

        # (GLCpts4) Reaction GLCpts4
        cdef double V165 = (1.29 * 42.78 * ptsg_P) - (1e-05 * M_g6p_c * ptsg)

        # (L_LACt2r) Reaction L_LACt2r
        cdef double V166 = 5e-09 * (M_lac__L_c - 0.001) * 3 / 2e-07

        # (PYRt2r) Reaction PYRt2r
        cdef double V167 = 5e-09 * (M_pyr_c - 0.001) * 3 / 2e-07

        # (ACt) Reaction ACt
        cdef double V168 = 5e-09 * (M_ac_c - 0.3) * 3 / 2e-07

        # (SMt) Reaction SMt
        cdef double V169 = 0.0011

        # (PCt) Reaction PCt
        cdef double V170 = 0.00014

        # (TAGt) Reaction TAGt
        cdef double V171 = 0.00031

        # (CHOLt) Reaction CHOLt
        cdef double V172 = 0.00291

        # (FAt) Reaction FAt
        cdef double V173 = 0.0051

        # (GLYCt) Reaction GLYCt
        cdef double V174 = 0.00333

        # (NAt) Reaction NAt
        cdef double V175 = 0.02



        return np.asarray([V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V25, V26, V27, V28, V29, V30, V31, V32, V33, V34, V35, V36, V37, V38, V39, V40, V41, V42, V43, V44, V45, V46, V47, V48, V49, V50, V51, V52, V53, V54, V55, V56, V57, V58, V59, V60, V61, V62, V63, V64, V65, V66, V67, V68, V69, V70, V71, V72, V73, V74, V75, V76, V77, V78, V79, V80, V81, V82, V83, V84, V85, V86, V87, V88, V89, V90, V91, V92, V93, V94, V95, V96, V97, V98, V99, V100, V101, V102, V103, V104, V105, V106, V107, V108, V109, V110, V111, V112, V113, V114, V115, V116, V117, V118, V119, V120, V121, V122, V123, V124, V125, V126, V127, V128, V129, V130, V131, V132, V133, V134, V135, V136, V137, V138, V139, V140, V141, V142, V143, V144, V145, V146, V147, V148, V149, V150, V151, V152, V153, V154, V155, V156, V157, V158, V159, V160, V161, V162, V163, V164, V165, V166, V167, V168, V169, V170, V171, V172, V173, V174, V175])


####
#Definition of main function:
####

    def __call__(self, float t, np.ndarray[DTYPEDBL_t, ndim=1] y):
        return self.compiledCall(t, y)

    cdef compiledCall(self, float t, np.ndarray[DTYPEDBL_t, ndim=1] y):

        cdef np.ndarray[DTYPEDBL_t, ndim=1] fluxes = self.calcFlux_c(t,y)

        # (PGI) Reaction PGI
        cdef double V1 = fluxes[0]

        # (PFK) Reaction PFK
        cdef double V2 = fluxes[1]

        # (FBA) Reaction FBA
        cdef double V3 = fluxes[2]

        # (TPI) Reaction TPI
        cdef double V4 = fluxes[3]

        # (GAPD) Reaction GAPD
        cdef double V5 = fluxes[4]

        # (GAPDP) Reaction GAPDP
        cdef double V6 = fluxes[5]

        # (PGK) Reaction PGK
        cdef double V7 = fluxes[6]

        # (PGM) Reaction PGM
        cdef double V8 = fluxes[7]

        # (ENO) Reaction ENO
        cdef double V9 = fluxes[8]

        # (PYK) Reaction PYK
        cdef double V10 = fluxes[9]

        # (LDH_L) Reaction LDH_L
        cdef double V11 = fluxes[10]

        # (PDH_acald) Reaction PDH_acald
        cdef double V12 = fluxes[11]

        # (PDH_E3) Reaction PDH_E3
        cdef double V13 = fluxes[12]

        # (PTAr) Reaction PTAr
        cdef double V14 = fluxes[13]

        # (ACKr) Reaction ACKr
        cdef double V15 = fluxes[14]

        # (NOX) Reaction NOX
        cdef double V16 = fluxes[15]

        # (TALA) Reaction TALA
        cdef double V17 = fluxes[16]

        # (TKT1) Reaction TKT1
        cdef double V18 = fluxes[17]

        # (TKT2) Reaction TKT2
        cdef double V19 = fluxes[18]

        # (RPE) Reaction RPE
        cdef double V20 = fluxes[19]

        # (RPI) Reaction RPI
        cdef double V21 = fluxes[20]

        # (PRPPS) Reaction PRPPS
        cdef double V22 = fluxes[21]

        # (PPM) Reaction PPM
        cdef double V23 = fluxes[22]

        # (PPM2) Reaction PPM2
        cdef double V24 = fluxes[23]

        # (DRPA) Reaction DRPA
        cdef double V25 = fluxes[24]

        # (ADPT) Reaction ADPT
        cdef double V26 = fluxes[25]

        # (PUNP1) Reaction PUNP1
        cdef double V27 = fluxes[26]

        # (ADK1) Reaction ADK1
        cdef double V28 = fluxes[27]

        # (RNDR1) Reaction RNDR1
        cdef double V29 = fluxes[28]

        # (TRDR) Reaction TRDR
        cdef double V30 = fluxes[29]

        # (PPA) Reaction PPA
        cdef double V31 = fluxes[30]

        # (PUNP2) Reaction PUNP2
        cdef double V32 = fluxes[31]

        # (DADNK) Reaction DADNK
        cdef double V33 = fluxes[32]

        # (DADK) Reaction DADK
        cdef double V34 = fluxes[33]

        # (PYK5) Reaction PYK5
        cdef double V35 = fluxes[34]

        # (PGK2) Reaction PGK2
        cdef double V36 = fluxes[35]

        # (GUAPRT) Reaction GUAPRT
        cdef double V37 = fluxes[36]

        # (PUNP3) Reaction PUNP3
        cdef double V38 = fluxes[37]

        # (GK1) Reaction GK1
        cdef double V39 = fluxes[38]

        # (RNDR2) Reaction RNDR2
        cdef double V40 = fluxes[39]

        # (PYK3) Reaction PYK3
        cdef double V41 = fluxes[40]

        # (PGK3) Reaction PGK3
        cdef double V42 = fluxes[41]

        # (PUNP4) Reaction PUNP4
        cdef double V43 = fluxes[42]

        # (DGSNK) Reaction DGSNK
        cdef double V44 = fluxes[43]

        # (DGK1) Reaction DGK1
        cdef double V45 = fluxes[44]

        # (PYK6) Reaction PYK6
        cdef double V46 = fluxes[45]

        # (PGK4) Reaction PGK4
        cdef double V47 = fluxes[46]

        # (UPPRT) Reaction UPPRT
        cdef double V48 = fluxes[47]

        # (UMPK) Reaction UMPK
        cdef double V49 = fluxes[48]

        # (PYK2) Reaction PYK2
        cdef double V50 = fluxes[49]

        # (RNDR4) Reaction RNDR4
        cdef double V51 = fluxes[50]

        # (CYTK1) Reaction CYTK1
        cdef double V52 = fluxes[51]

        # (PYK4) Reaction PYK4
        cdef double V53 = fluxes[52]

        # (RNDR3) Reaction RNDR3
        cdef double V54 = fluxes[53]

        # (DCYTK) Reaction DCYTK
        cdef double V55 = fluxes[54]

        # (CYTK2) Reaction CYTK2
        cdef double V56 = fluxes[55]

        # (PYK7) Reaction PYK7
        cdef double V57 = fluxes[56]

        # (DURIK1) Reaction DURIK1
        cdef double V58 = fluxes[57]

        # (TMDK1) Reaction TMDK1
        cdef double V59 = fluxes[58]

        # (TMPK) Reaction TMPK
        cdef double V60 = fluxes[59]

        # (PYK8) Reaction PYK8
        cdef double V61 = fluxes[60]

        # (PYK9) Reaction PYK9
        cdef double V62 = fluxes[61]

        # (CTPS2) Reaction CTPS2
        cdef double V63 = fluxes[62]

        # (CTPSDUMP) Reaction CTPSDUMP
        cdef double V64 = fluxes[63]

        # (DCMPDA) Reaction DCMPDA
        cdef double V65 = fluxes[64]

        # (DUTPDP) Reaction DUTPDP
        cdef double V66 = fluxes[65]

        # (NTD5) Reaction NTD5
        cdef double V67 = fluxes[66]

        # (NTD1) Reaction NTD1
        cdef double V68 = fluxes[67]

        # (NTD8) Reaction NTD8
        cdef double V69 = fluxes[68]

        # (NTD6) Reaction NTD6
        cdef double V70 = fluxes[69]

        # (PUNP5) Reaction PUNP5
        cdef double V71 = fluxes[70]

        # (GLYK) Reaction GLYK
        cdef double V72 = fluxes[71]

        # (ACPS) Reaction ACPS
        cdef double V73 = fluxes[72]

        # (BPNT) Reaction BPNT
        cdef double V74 = fluxes[73]

        # (FAKr) Reaction FAKr
        cdef double V75 = fluxes[74]

        # (ACPPAT) Reaction ACPPAT
        cdef double V76 = fluxes[75]

        # (APG3PAT) Reaction APG3PAT
        cdef double V77 = fluxes[76]

        # (AGPAT) Reaction AGPAT
        cdef double V78 = fluxes[77]

        # (DASYN) Reaction DASYN
        cdef double V79 = fluxes[78]

        # (PGSA) Reaction PGSA
        cdef double V80 = fluxes[79]

        # (PGPP) Reaction PGPP
        cdef double V81 = fluxes[80]

        # (CLPNS) Reaction CLPNS
        cdef double V82 = fluxes[81]

        # (PAPA) Reaction PAPA
        cdef double V83 = fluxes[82]

        # (PGMT) Reaction PGMT
        cdef double V84 = fluxes[83]

        # (GALU) Reaction GALU
        cdef double V85 = fluxes[84]

        # (UDPG4E) Reaction UDPG4E
        cdef double V86 = fluxes[85]

        # (UDPGALM) Reaction UDPGALM
        cdef double V87 = fluxes[86]

        # (DAGGALT) Reaction DAGGALT
        cdef double V88 = fluxes[87]

        # (NCTPPRT) Reaction NCTPPRT
        cdef double V89 = fluxes[88]

        # (NNATr) Reaction NNATr
        cdef double V90 = fluxes[89]

        # (NADS) Reaction NADS
        cdef double V91 = fluxes[90]

        # (NADK) Reaction NADK
        cdef double V92 = fluxes[91]

        # (NADHK) Reaction NADHK
        cdef double V93 = fluxes[92]

        # (RBFK) Reaction RBFK
        cdef double V94 = fluxes[93]

        # (FMNAT) Reaction FMNAT
        cdef double V95 = fluxes[94]

        # (5FTHFPGS) Reaction 5FTHFPGS
        cdef double V96 = fluxes[95]

        # (MTHFC) Reaction MTHFC
        cdef double V97 = fluxes[96]

        # (FMETTRS) Reaction FMETTRS
        cdef double V98 = fluxes[97]

        # (GHMT) Reaction GHMT
        cdef double V99 = fluxes[98]

        # (MTHFD) Reaction MTHFD
        cdef double V100 = fluxes[99]

        # (FTHFCL) Reaction FTHFCL
        cdef double V101 = fluxes[100]

        # (GHMT2) Reaction GHMT2
        cdef double V102 = fluxes[101]

        # (ATPase) Reaction ATPase
        cdef double V103 = fluxes[102]

        # (RIBFLVabc) Reaction RIBFLVabc
        cdef double V104 = fluxes[103]

        # (P5Pabc) Reaction P5Pabc
        cdef double V105 = fluxes[104]

        # (5FTHFabc) Reaction 5FTHFabc
        cdef double V106 = fluxes[105]

        # (NACabc) Reaction NACabc
        cdef double V107 = fluxes[106]

        # (COAabc) Reaction COAabc
        cdef double V108 = fluxes[107]

        # (THMPPabc) Reaction THMPPabc
        cdef double V109 = fluxes[108]

        # (SPRMabc) Reaction SPRMabc
        cdef double V110 = fluxes[109]

        # (DCYTabc) Reaction DCYTabc
        cdef double V111 = fluxes[110]

        # (URIabc) Reaction URIabc
        cdef double V112 = fluxes[111]

        # (THMDabc) Reaction THMDabc
        cdef double V113 = fluxes[112]

        # (ADNabc) Reaction ADNabc
        cdef double V114 = fluxes[113]

        # (DADNabc) Reaction DADNabc
        cdef double V115 = fluxes[114]

        # (GSNabc) Reaction GSNabc
        cdef double V116 = fluxes[115]

        # (DGSNabc) Reaction DGSNabc
        cdef double V117 = fluxes[116]

        # (Kt6) Reaction Kt6
        cdef double V118 = fluxes[117]

        # (MG2abc) Reaction MG2abc
        cdef double V119 = fluxes[118]

        # (CA2abc) Reaction CA2abc
        cdef double V120 = fluxes[119]

        # (PIabc) Reaction PIabc
        cdef double V121 = fluxes[120]

        # (ARGabc) Reaction ARGabc
        cdef double V122 = fluxes[121]

        # (ASPabc) Reaction ASPabc
        cdef double V123 = fluxes[122]

        # (GLUabc) Reaction GLUabc
        cdef double V124 = fluxes[123]

        # (GLYabc) Reaction GLYabc
        cdef double V125 = fluxes[124]

        # (ILEabc) Reaction ILEabc
        cdef double V126 = fluxes[125]

        # (ALAabc) Reaction ALAabc
        cdef double V127 = fluxes[126]

        # (ASNabc) Reaction ASNabc
        cdef double V128 = fluxes[127]

        # (LEUabc) Reaction LEUabc
        cdef double V129 = fluxes[128]

        # (GLNabc) Reaction GLNabc
        cdef double V130 = fluxes[129]

        # (HISabc) Reaction HISabc
        cdef double V131 = fluxes[130]

        # (LYSabc) Reaction LYSabc
        cdef double V132 = fluxes[131]

        # (PROabc) Reaction PROabc
        cdef double V133 = fluxes[132]

        # (PHEabc) Reaction PHEabc
        cdef double V134 = fluxes[133]

        # (THRabc) Reaction THRabc
        cdef double V135 = fluxes[134]

        # (TRPabc) Reaction TRPabc
        cdef double V136 = fluxes[135]

        # (TYRabc) Reaction TYRabc
        cdef double V137 = fluxes[136]

        # (VALabc) Reaction VALabc
        cdef double V138 = fluxes[137]

        # (SERabc) Reaction SERabc
        cdef double V139 = fluxes[138]

        # (METabc) Reaction METabc
        cdef double V140 = fluxes[139]

        # (ARGt2r) Reaction ARGt2r
        cdef double V141 = fluxes[140]

        # (ASPt2pr) Reaction ASPt2pr
        cdef double V142 = fluxes[141]

        # (CYSt2r) Reaction CYSt2r
        cdef double V143 = fluxes[142]

        # (GLUt2pr) Reaction GLUt2pr
        cdef double V144 = fluxes[143]

        # (GLYt2r) Reaction GLYt2r
        cdef double V145 = fluxes[144]

        # (ISOt2r) Reaction ISOt2r
        cdef double V146 = fluxes[145]

        # (ALAt2r) Reaction ALAt2r
        cdef double V147 = fluxes[146]

        # (ASNt2r) Reaction ASNt2r
        cdef double V148 = fluxes[147]

        # (LEUt2r) Reaction LEUt2r
        cdef double V149 = fluxes[148]

        # (GLNt2r) Reaction GLNt2r
        cdef double V150 = fluxes[149]

        # (HISt2r) Reaction HISt2r
        cdef double V151 = fluxes[150]

        # (LYSt2r) Reaction LYSt2r
        cdef double V152 = fluxes[151]

        # (PROt2r) Reaction PROt2r
        cdef double V153 = fluxes[152]

        # (PHEt2r) Reaction PHEt2r
        cdef double V154 = fluxes[153]

        # (THRt2r) Reaction THRt2r
        cdef double V155 = fluxes[154]

        # (TRPt2r) Reaction TRPt2r
        cdef double V156 = fluxes[155]

        # (TYRt2r) Reaction TYRt2r
        cdef double V157 = fluxes[156]

        # (VALt2r) Reaction VALt2r
        cdef double V158 = fluxes[157]

        # (SERt2r) Reaction SERt2r
        cdef double V159 = fluxes[158]

        # (METt2r) Reaction METt2r
        cdef double V160 = fluxes[159]

        # (GLCpts0) Reaction GLCpts0
        cdef double V161 = fluxes[160]

        # (GLCpts1) Reaction GLCpts1
        cdef double V162 = fluxes[161]

        # (GLCpts2) Reaction GLCpts2
        cdef double V163 = fluxes[162]

        # (GLCpts3) Reaction GLCpts3
        cdef double V164 = fluxes[163]

        # (GLCpts4) Reaction GLCpts4
        cdef double V165 = fluxes[164]

        # (L_LACt2r) Reaction L_LACt2r
        cdef double V166 = fluxes[165]

        # (PYRt2r) Reaction PYRt2r
        cdef double V167 = fluxes[166]

        # (ACt) Reaction ACt
        cdef double V168 = fluxes[167]

        # (SMt) Reaction SMt
        cdef double V169 = fluxes[168]

        # (PCt) Reaction PCt
        cdef double V170 = fluxes[169]

        # (TAGt) Reaction TAGt
        cdef double V171 = fluxes[170]

        # (CHOLt) Reaction CHOLt
        cdef double V172 = fluxes[171]

        # (FAt) Reaction FAt
        cdef double V173 = fluxes[172]

        # (GLYCt) Reaction GLYCt
        cdef double V174 = fluxes[173]

        # (NAt) Reaction NAt
        cdef double V175 = fluxes[174]

        cdef double DeltaM_ACP_c = +1*V73 -1*V76 +1*V78 + 0
        cdef double DeltaM_ACP_R_c = +1*V76 -1*V78 + 0
        cdef double DeltaM_apoACP_c = -1*V73 + 0
        cdef double DeltaM_trdrd_c = -1*V40 -1*V51 -1*V54 -1*V29 +1*V30 + 0
        cdef double DeltaM_trdox_c = +1*V40 +1*V51 +1*V54 +1*V29 -1*V30 + 0
        cdef double DeltaM_dhlpl_PdhC_c = +1*V12 -1*V13 + 0
        cdef double DeltaM_acdhlpl_PdhC_c = 0
        cdef double DeltaM_lpl_PdhC_c = -1*V12 +1*V13 + 0
        cdef double Deltaptsi_P = +1*V161 -1*V162 + 0
        cdef double Deltaptsi = -1*V161 +1*V162 + 0
        cdef double Deltaptsh_P = +1*V162 -1*V163 + 0
        cdef double Deltaptsh = -1*V162 +1*V163 + 0
        cdef double Deltacrr_P = +1*V163 -1*V164 + 0
        cdef double Deltacrr = -1*V163 +1*V164 + 0
        cdef double Deltaptsg_P = +1*V164 -1*V165 + 0
        cdef double Deltaptsg = -1*V164 +1*V165 + 0
        cdef double DeltaM_g6p_c = -1*V1 -1*V84 +1*V165 + 0
        cdef double DeltaM_f6p_c = +1*V1 -1*V2 -1*V19 -1*V17 + 0
        cdef double DeltaM_atp_c = -1*V129 -1*V2 -1*V130 -1*V131 -1*V132 -1*V133 +1*V7 -1*V134 -1*V135  
        DeltaM_atp_c += +1*V10 -1*V136 -1*V137 -1*V138 -1*V139 +1*V15 -1*V140 -1*V22 -1*V28 -1*V33  
        DeltaM_atp_c += -1*V34 -1*V39 -1*V44 -1*V45 -1*V49 -1*V52 -1*V55 -1*V56 -1*V58 -1*V59  
        DeltaM_atp_c += -1*V60 -1*V63 -1*V64 -1*V72 -1*V75 -1*V90 -1*V91 -1*V92 -1*V93 -1*V94  
        DeltaM_atp_c += -1*V95 -2*V96 -1*V101 +1*V103 -1*V104 -1*V105 -1*V106 -1*V107 -1*V108 -1*V109  
        DeltaM_atp_c += -1*V110 -1*V111 -1*V112 -1*V113 -1*V114 -1*V115 -1*V116 -1*V117 -1*V118 -1*V119  
        DeltaM_atp_c += -1*V120 -1*V121 -1*V122 -1*V123 -1*V124 -1*V125 -1*V126 -1*V127 -1*V128 + 0
        cdef double DeltaM_adp_c = +1*V129 +1*V2 +1*V130 +1*V131 +1*V132 +1*V133 -1*V7 +1*V134 +1*V135  
        DeltaM_adp_c += -1*V10 +1*V136 +1*V137 +1*V138 +1*V139 -1*V15 +1*V140 +2*V28 -1*V29 +1*V33  
        DeltaM_adp_c += +1*V34 +1*V39 +1*V44 +1*V45 +1*V49 +1*V52 +1*V55 +1*V56 +1*V58 +1*V59  
        DeltaM_adp_c += +1*V60 +1*V63 +1*V64 +1*V72 +1*V75 +1*V92 +1*V93 +1*V94 +2*V96 +1*V101  
        DeltaM_adp_c += -1*V103 +1*V104 +1*V105 +1*V106 +1*V107 +1*V108 +1*V109 +1*V110 +1*V111 +1*V112  
        DeltaM_adp_c += +1*V113 +1*V114 +1*V115 +1*V116 +1*V117 +1*V118 +1*V119 +1*V120 +1*V121 +1*V122  
        DeltaM_adp_c += +1*V123 +1*V124 +1*V125 +1*V126 +1*V127 +1*V128 + 0
        cdef double DeltaM_fdp_c = +1*V2 -1*V3 + 0
        cdef double DeltaM_dhap_c = +1*V3 -1*V4 + 0
        cdef double DeltaM_g3p_c = +1*V3 +1*V4 -1*V5 -1*V6 +1*V17 -1*V18 -1*V19 +1*V25 + 0
        cdef double DeltaM_nad_c = -1*V5 +1*V11 -1*V13 +2*V16 +1*V91 -1*V92 + 0
        cdef double DeltaM_pi_c = +1*V129 +1*V130 +1*V131 +1*V132 -1*V5 +1*V133 +1*V134 +1*V135 +1*V136  
        DeltaM_pi_c += +1*V137 +1*V138 +1*V139 +1*V140 -1*V14 -1*V27 +2*V31 -1*V32 -1*V38 -1*V43  
        DeltaM_pi_c += +1*V63 +1*V64 +1*V67 +1*V68 +1*V69 +1*V70 -1*V71 +1*V74 +1*V76 +1*V77  
        DeltaM_pi_c += +1*V81 +1*V83 +2*V96 +1*V101 -1*V103 +1*V104 +1*V105 +1*V106 +1*V107 +1*V108  
        DeltaM_pi_c += +1*V109 +1*V110 +1*V111 +1*V112 +1*V113 +1*V114 +1*V115 +1*V116 +1*V117 +1*V118  
        DeltaM_pi_c += +1*V119 +1*V120 +2*V121 +1*V122 +1*V123 +1*V124 +1*V125 +1*V126 +1*V127 +1*V128 + 0
        cdef double DeltaM_13dpg_c = -1*V36 +1*V5 -1*V7 -1*V42 -1*V47 + 0
        cdef double DeltaM_nadh_c = +1*V5 -1*V11 +1*V13 -2*V16 -1*V93 + 0
        cdef double DeltaM_nadp_c = -1*V100 +1*V30 +1*V92 -1*V6 + 0
        cdef double DeltaM_3pg_c = +1*V36 +1*V6 +1*V7 -1*V8 +1*V42 +1*V47 + 0
        cdef double DeltaM_nadph_c = -1*V30 +1*V100 +1*V93 +1*V6 + 0
        cdef double DeltaM_2pg_c = -1*V9 +1*V8 + 0
        cdef double DeltaM_pep_c = -1*V161 -1*V35 +1*V9 -1*V10 -1*V41 -1*V46 -1*V50 -1*V53 -1*V57  
        DeltaM_pep_c += -1*V61 -1*V62 + 0
        cdef double DeltaM_pyr_c = +1*V161 +1*V35 -1*V167 +1*V41 +1*V10 -1*V11 +1*V46 +1*V50 +1*V53  
        DeltaM_pyr_c += +1*V57 +1*V61 +1*V62 + 0
        cdef double DeltaM_lac__L_c = +1*V11 -1*V166 + 0
        cdef double DeltaM_acald_c = +1*V25 -1*V12 + 0
        cdef double DeltaM_coa_c = -1*V73 +1*V108 -1*V12 +1*V14 + 0
        cdef double DeltaM_accoa_c = +1*V12 -1*V14 + 0
        cdef double DeltaM_actp_c = +1*V14 -1*V15 + 0
        cdef double DeltaM_ac_c = +1*V15 -1*V168 + 0
        cdef double DeltaM_o2_c = 0
        cdef double DeltaM_e4p_c = -1*V17 +1*V19 + 0
        cdef double DeltaM_s7p_c = +1*V17 -1*V18 + 0
        cdef double DeltaM_r5p_c = +1*V18 +1*V21 -1*V22 +1*V23 + 0
        cdef double DeltaM_xu5p__D_c = +1*V18 +1*V19 -1*V20 + 0
        cdef double DeltaM_ru5p__D_c = +1*V20 -1*V21 + 0
        cdef double DeltaM_amp_c = +1*V74 +1*V22 +1*V26 +1*V91 -1*V28 + 0
        cdef double DeltaM_prpp_c = -1*V37 -1*V48 +1*V22 -1*V89 -1*V26 + 0
        cdef double DeltaM_r1p_c = +1*V71 +1*V27 +1*V38 -1*V23 + 0
        cdef double DeltaM_2dr1p_c = +1*V43 +1*V32 -1*V24 + 0
        cdef double DeltaM_2dr5p_c = -1*V25 +1*V24 + 0
        cdef double DeltaM_ade_c = -1*V26 +1*V27 +1*V32 + 0
        cdef double DeltaM_ppi_c = +1*V90 +1*V66 +1*V37 +1*V79 +1*V48 +1*V85 +1*V95 +1*V89 +1*V26  
        DeltaM_ppi_c += +1*V91 -1*V31 + 0
        cdef double DeltaM_adn_c = +1*V114 -1*V27 + 0
        cdef double DeltaM_dadp_c = +1*V34 -1*V35 -1*V36 +1*V29 + 0
        cdef double DeltaM_dad_2_c = -1*V33 +1*V115 +1*V70 -1*V32 + 0
        cdef double DeltaM_damp_c = +1*V33 -1*V34 -1*V70 + 0
        cdef double DeltaM_datp_c = +1*V35 +1*V36 + 0
        cdef double DeltaM_gua_c = +1*V43 -1*V37 +1*V38 + 0
        cdef double DeltaM_gmp_c = +1*V37 -1*V39 + 0
        cdef double DeltaM_gsn_c = +1*V116 -1*V38 + 0
        cdef double DeltaM_gdp_c = -1*V41 -1*V42 +1*V39 -1*V40 + 0
        cdef double DeltaM_dgdp_c = +1*V45 -1*V46 -1*V47 +1*V40 + 0
        cdef double DeltaM_gtp_c = +1*V41 +1*V42 + 0
        cdef double DeltaM_dgsn_c = +1*V117 -1*V43 -1*V44 +1*V69 + 0
        cdef double DeltaM_dgmp_c = +1*V44 -1*V45 -1*V69 + 0
        cdef double DeltaM_dgtp_c = +1*V46 +1*V47 + 0
        cdef double DeltaM_ura_c = +1*V71 -1*V48 + 0
        cdef double DeltaM_ump_c = -1*V49 +1*V48 + 0
        cdef double DeltaM_udp_c = +1*V49 -1*V50 -1*V51 +1*V88 + 0
        cdef double DeltaM_utp_c = +1*V50 -1*V85 -1*V63 + 0
        cdef double DeltaM_dudp_c = +1*V51 -1*V62 + 0
        cdef double DeltaM_cmp_c = -1*V52 +1*V80 + 0
        cdef double DeltaM_cdp_c = +1*V52 -1*V53 -1*V54 + 0
        cdef double DeltaM_ctp_c = -1*V79 +1*V53 +1*V63 + 0
        cdef double DeltaM_dcdp_c = -1*V57 +1*V54 +1*V56 + 0
        cdef double DeltaM_dcyt_c = +1*V111 -1*V55 + 0
        cdef double DeltaM_dcmp_c = -1*V65 +1*V64 +1*V55 -1*V56 + 0
        cdef double DeltaM_dctp_c = +1*V57 + 0
        cdef double DeltaM_duri_c = -1*V58 +1*V68 + 0
        cdef double DeltaM_dump_c = +1*V65 +1*V66 -1*V68 +1*V58 -1*V64 + 0
        cdef double DeltaM_thymd_c = +1*V113 -1*V59 +1*V67 + 0
        cdef double DeltaM_dtmp_c = +1*V59 -1*V60 -1*V67 + 0
        cdef double DeltaM_dtdp_c = +1*V60 -1*V61 + 0
        cdef double DeltaM_dttp_c = +1*V61 + 0
        cdef double DeltaM_dutp_c = -1*V66 +1*V62 + 0
        cdef double DeltaM_gln__L_c = +1*V130 +1*V150 -1*V63 -1*V64 + 0
        cdef double DeltaM_glu__L_c = +1*V144 +1*V124 +1*V64 +1*V63 -2*V96 + 0
        cdef double DeltaM_nh3_c = +1*V65 -1*V91 + 0
        cdef double DeltaM_uri_c = -1*V71 +1*V112 + 0
        cdef double DeltaM_glyc_c = +1*V82 +1*V174 -1*V72 + 0
        cdef double DeltaM_glyc3p_c = -1*V77 -1*V80 +1*V72 + 0
        cdef double DeltaM_pap_c = +1*V73 -1*V74 + 0
        cdef double DeltaM_fa_c = -1*V75 +1*V173 + 0
        cdef double DeltaM_ap_c = +1*V75 -1*V76 -1*V77 + 0
        cdef double DeltaM_1ag3p_c = +1*V77 -1*V78 + 0
        cdef double DeltaM_pa_c = -1*V83 +1*V78 -1*V79 + 0
        cdef double DeltaM_cdpdag_c = +1*V79 -1*V80 + 0
        cdef double DeltaM_pg3p_c = -1*V81 +1*V80 + 0
        cdef double DeltaM_pg_c = +1*V81 -2*V82 + 0
        cdef double DeltaM_clpn_c = +1*V82 + 0
        cdef double DeltaM_12dgr_c = +1*V83 -1*V88 + 0
        cdef double DeltaM_g1p_c = +1*V84 -1*V85 + 0
        cdef double DeltaM_udpg_c = +1*V85 -1*V86 + 0
        cdef double DeltaM_udpgal_c = +1*V86 -1*V87 + 0
        cdef double DeltaM_udpgalfur_c = +1*V87 -1*V88 + 0
        cdef double DeltaM_galfur12dgr_c = +1*V88 + 0
        cdef double DeltaM_nac_c = -1*V89 +1*V107 + 0
        cdef double DeltaM_nicrnt_c = +1*V89 -1*V90 + 0
        cdef double DeltaM_dnad_c = +1*V90 -1*V91 + 0
        cdef double DeltaM_ribflv_c = -1*V94 +1*V104 + 0
        cdef double DeltaM_fmn_c = +1*V94 -1*V95 + 0
        cdef double DeltaM_fad_c = +1*V95 + 0
        cdef double DeltaM_5fthf_c = +1*V106 -1*V96 + 0
        cdef double DeltaM_5fthfglu3_c = -1*V101 +1*V102 +1*V96 + 0
        cdef double DeltaM_methfglu3_c = -1*V97 +1*V100 +1*V101 -1*V102 + 0
        cdef double DeltaM_10fthfglu3_c = +1*V97 -1*V98 + 0
        cdef double DeltaM_mettrna_c = -1*V98 + 0
        cdef double DeltaM_fmettrna_c = +1*V98 + 0
        cdef double DeltaM_thfglu3_c = +1*V98 -1*V99 + 0
        cdef double DeltaM_ser__L_c = -1*V99 +1*V159 +1*V139 + 0
        cdef double DeltaM_gly_c = +1*V145 +1*V99 +1*V125 + 0
        cdef double DeltaM_mlthfglu3_c = +1*V99 -1*V100 + 0
        cdef double DeltaM_pydx5p_c = +1*V105 + 0
        cdef double DeltaM_thmpp_c = +1*V109 + 0
        cdef double DeltaM_sprm_c = +1*V110 + 0
        cdef double DeltaM_na1_c = -1*V118 +1*V175 + 0
        cdef double DeltaM_k_c = +1*V118 + 0
        cdef double DeltaM_mg2_c = +1*V119 + 0
        cdef double DeltaM_ca2_c = +1*V120 + 0
        cdef double DeltaM_arg__L_c = +1*V122 +1*V141 + 0
        cdef double DeltaM_asp__L_c = +1*V123 +1*V142 + 0
        cdef double DeltaM_ile__L_c = +1*V146 +1*V126 + 0
        cdef double DeltaM_ala__L_c = +1*V147 +1*V127 + 0
        cdef double DeltaM_asn__L_c = +1*V148 +1*V128 + 0
        cdef double DeltaM_leu__L_c = +1*V129 +1*V149 + 0
        cdef double DeltaM_his__L_c = +1*V131 +1*V151 + 0
        cdef double DeltaM_lys__L_c = +1*V132 +1*V152 + 0
        cdef double DeltaM_pro__L_c = +1*V153 +1*V133 + 0
        cdef double DeltaM_phe__L_c = +1*V154 +1*V134 + 0
        cdef double DeltaM_thr__L_c = +1*V155 +1*V135 + 0
        cdef double DeltaM_trp__L_c = +1*V156 +1*V136 + 0
        cdef double DeltaM_tyr__L_c = +1*V137 +1*V157 + 0
        cdef double DeltaM_val__L_c = +1*V138 +1*V158 + 0
        cdef double DeltaM_met__L_c = +1*V140 +1*V160 + 0
        cdef double DeltaM_cys__L_c = +1*V143 + 0
        cdef double DeltaM_sm_c = +1*V169 + 0
        cdef double DeltaM_pc_c = +1*V170 + 0
        cdef double DeltaM_tag_c = +1*V171 + 0
        cdef double DeltaM_chsterol_c = +1*V172 + 0


        return np.asarray([DeltaM_ACP_c, DeltaM_ACP_R_c, DeltaM_apoACP_c, DeltaM_trdrd_c, DeltaM_trdox_c, DeltaM_dhlpl_PdhC_c, DeltaM_acdhlpl_PdhC_c, DeltaM_lpl_PdhC_c, Deltaptsi_P, Deltaptsi, Deltaptsh_P, Deltaptsh, Deltacrr_P, Deltacrr, Deltaptsg_P, Deltaptsg, DeltaM_g6p_c, DeltaM_f6p_c, DeltaM_atp_c, DeltaM_adp_c, DeltaM_fdp_c, DeltaM_dhap_c, DeltaM_g3p_c, DeltaM_nad_c, DeltaM_pi_c, DeltaM_13dpg_c, DeltaM_nadh_c, DeltaM_nadp_c, DeltaM_3pg_c, DeltaM_nadph_c, DeltaM_2pg_c, DeltaM_pep_c, DeltaM_pyr_c, DeltaM_lac__L_c, DeltaM_acald_c, DeltaM_coa_c, DeltaM_accoa_c, DeltaM_actp_c, DeltaM_ac_c, DeltaM_o2_c, DeltaM_e4p_c, DeltaM_s7p_c, DeltaM_r5p_c, DeltaM_xu5p__D_c, DeltaM_ru5p__D_c, DeltaM_amp_c, DeltaM_prpp_c, DeltaM_r1p_c, DeltaM_2dr1p_c, DeltaM_2dr5p_c, DeltaM_ade_c, DeltaM_ppi_c, DeltaM_adn_c, DeltaM_dadp_c, DeltaM_dad_2_c, DeltaM_damp_c, DeltaM_datp_c, DeltaM_gua_c, DeltaM_gmp_c, DeltaM_gsn_c, DeltaM_gdp_c, DeltaM_dgdp_c, DeltaM_gtp_c, DeltaM_dgsn_c, DeltaM_dgmp_c, DeltaM_dgtp_c, DeltaM_ura_c, DeltaM_ump_c, DeltaM_udp_c, DeltaM_utp_c, DeltaM_dudp_c, DeltaM_cmp_c, DeltaM_cdp_c, DeltaM_ctp_c, DeltaM_dcdp_c, DeltaM_dcyt_c, DeltaM_dcmp_c, DeltaM_dctp_c, DeltaM_duri_c, DeltaM_dump_c, DeltaM_thymd_c, DeltaM_dtmp_c, DeltaM_dtdp_c, DeltaM_dttp_c, DeltaM_dutp_c, DeltaM_gln__L_c, DeltaM_glu__L_c, DeltaM_nh3_c, DeltaM_uri_c, DeltaM_glyc_c, DeltaM_glyc3p_c, DeltaM_pap_c, DeltaM_fa_c, DeltaM_ap_c, DeltaM_1ag3p_c, DeltaM_pa_c, DeltaM_cdpdag_c, DeltaM_pg3p_c, DeltaM_pg_c, DeltaM_clpn_c, DeltaM_12dgr_c, DeltaM_g1p_c, DeltaM_udpg_c, DeltaM_udpgal_c, DeltaM_udpgalfur_c, DeltaM_galfur12dgr_c, DeltaM_nac_c, DeltaM_nicrnt_c, DeltaM_dnad_c, DeltaM_ribflv_c, DeltaM_fmn_c, DeltaM_fad_c, DeltaM_5fthf_c, DeltaM_5fthfglu3_c, DeltaM_methfglu3_c, DeltaM_10fthfglu3_c, DeltaM_mettrna_c, DeltaM_fmettrna_c, DeltaM_thfglu3_c, DeltaM_ser__L_c, DeltaM_gly_c, DeltaM_mlthfglu3_c, DeltaM_pydx5p_c, DeltaM_thmpp_c, DeltaM_sprm_c, DeltaM_na1_c, DeltaM_k_c, DeltaM_mg2_c, DeltaM_ca2_c, DeltaM_arg__L_c, DeltaM_asp__L_c, DeltaM_ile__L_c, DeltaM_ala__L_c, DeltaM_asn__L_c, DeltaM_leu__L_c, DeltaM_his__L_c, DeltaM_lys__L_c, DeltaM_pro__L_c, DeltaM_phe__L_c, DeltaM_thr__L_c, DeltaM_trp__L_c, DeltaM_tyr__L_c, DeltaM_val__L_c, DeltaM_met__L_c, DeltaM_cys__L_c, DeltaM_sm_c, DeltaM_pc_c, DeltaM_tag_c, DeltaM_chsterol_c])

