# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 11:43:05 2021

@author: 81806
"""


import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
    

def Evap(Q_apply, md, T_wickout, h_evap, n_g, l_g, t_g, w_g, h_out, T_amb, 
         h_cont_echb, h_cont_hb, t_insu_hb, k_insu_hb, h_cont_ec, t_insu_ec, k_insu_ec,
         wick_l, HB_l, HB_w, EC_l, EC_w):
    
    # 蒸発・熱コンダクタンス [W/K]
    G_evap = h_evap * (wick_l - w_g * n_g) * l_g 
    
    # 伝熱面温度:Tec [K]　ただし，Tw.mのときTeとしているので注意
    T_ec = T_wickout + md * HV(T_wickout) / G_evap
    
    # 蒸気生成に用いられる熱量:Q_ev [W]
    Q_ev = md * HV(T_wickout)
    
    ###### ヒータブロック・蒸発器ケース・外界に関する各種熱リーク  ################
    G_hbamb = (1 / 
               ((1 / HB_l * HB_w * h_cont_hb) 
               + t_insu_hb / (k_insu_hb * HB_l * HB_w)
               + 1 / (HB_l * HB_w * h_out) ) )
    
    G_ecamb = (1 / 
               ((1 / ((EC_l * EC_w - HB_l * HB_w) * h_cont_ec)
                + t_insu_ec / (k_insu_ec * (EC_l * EC_w - HB_l * HB_w))
                + 1 / ((EC_l + EC_w - HB_l * HB_w) * h_out)) ) )

    G_hbec = HB_l * HB_w * h_cont_echb
    
    Q_ecamb = G_ecamb * (T_ec - T_amb)
    
    T_hb = (Q_apply + G_hbec * T_ec + G_hbamb * T_amb) / (G_hbec + G_hbamb)
    
    Q_hbamb = G_hbamb * (T_hb - T_amb)
    
    Q_hbec = G_hbec * (T_hb - T_ec)    
    
    ###########################################################################
    # 溝に関するパラメータ
    d_g = 4 * (w_g * t_g) / (2 * (w_g + t_g))     # 溝の水力等価直径 [m]
    
    Re_gr = 4 * (md / n_g) / (np.pi * MuV(T_wickout) * d_g)   # 一つ溝あたりのレイノルズ数
    
    if Re_gr < 2300:    # 層流
        h_gr = 3.66 * kV(T_wickout) / d_g                     # 溝の熱伝達率[W/m^2/K]
        
    else:               # 乱流
        Pr = MuV(T_wickout) * VCp(T_wickout) / kV(T_wickout)
        h_gr = 0.023 * Re_gr**(0.8) * Pr**(0.4)* kV(T_wickout) / d_g
    
    
    # 溝出口の蒸気温度:T_grout [K]
    A = w_g * t_g
    B = w_g 
    C = (h_gr * w_g / (md/n_g) / VCp(T_wickout) * (w_g * t_g + w_g * l_g))  + w_g
    D = w_g * T_wickout + h_gr * w_g * T_ec/ (md/n_g) / VCp(T_wickout) * (w_g * t_g + w_g * l_g)
    
    T_grout = (T_wickout - D/C) * A**(C/B) * (A+B*l_g)**(-C/B) + D/C
    
    # 伝熱面から溝へ与えられる熱量:Qecgr [W]
    Q_ecgr = md * VCp(T_wickout) * (T_grout - T_wickout)
    
    # 溝の圧力損失を計算する場合の等価長さ:L_g [m]
    L_g = l_g
    
    # 溝の圧力損失:dP_gr [Pa]
    if Re_gr < 2300:     # 層流
        dP_gr = 128 * md / n_g * MuV(T_wickout) / VD(T_wickout) / d_g**4 / np.pi * L_g
    
    else:                # 乱流
        f = 0.0791 * Re_gr**(-0.25)
        dP_gr = 32 * f / np.pi**2 * (md/n_g) / (VD(T_wickout) *d_g**5) * L_g
    
    # ウィック出口の圧力: Pwickout [Pa]
    P_wickout = P_sat(T_wickout)
    # 溝出口の圧力 : P_grout [Pa]
    P_grout   = P_wickout - dP_gr
    
    
    return [P_wickout, P_grout, dP_gr, T_ec, T_hb, T_grout, Q_ev, Q_ecgr, Q_hbamb, Q_ecamb, Q_hbec, Re_gr]

    

def line_1D(md, P_in, T_in, T_out, X_in, component, L_x, num_n, ID, OD, theta, h_sink, 
            k_v, k_c, k_l, t_insu, k_insu, h_cont, h_out, wi_cl, P_cri):
    
    ###########################################################################
    # インプットパラメータ 
    ###########################################################################
    
    g = 9.806               # 重力加速度 [m/s^2]
    
    ###########################################################################
    #           forループ前の各種アウトプットパラメータの初期化 
    ###########################################################################
    if component == 1:
        k_line = k_v        #蒸気管の熱伝導率: k_line [W/m/K]を指定
    elif component==2:
        k_line = k_c        #凝縮器の熱伝導率: k_line [W/m/K]を指定
    elif component==3:
        k_line = k_l        #液管の熱伝導率: k_line [W/m/K]を指定
    
    dL_n = L_x/num_n        # 計算ノード間隔: dL_n [-] / 対象機器の長さ: L_x / 対象機器の分割ノード数: num_n
                        
    Q_out = 0                    # 各種管の全体の放熱量: Qout [W]
    new_T = np.zeros(num_n)      # 各種管の各ノードにおける温度: newT [K]
    new_P = np.zeros(num_n)      # 各種管の各ノードにおける圧力: newP [Pa]
    new_X = np.zeros(num_n)      # 各種管の各ノードにおけるクオリティー: newX [-]
    new_V = np.zeros(num_n)      # 各種管の各ノードにおけるボイド率: newV [-]
    
    new_Re    = np.zeros(num_n) 
    new_Re_V  = np.zeros(num_n) 
    new_Re_L  = np.zeros(num_n) 
    
    h_in      = np.zeros(num_n) 
    
    reg       = np.zeros(num_n)  # 凝縮のレジーム記録用
    Gcon      = np.zeros(num_n)  # 各輸送管の各セルの熱コンダクタンスの記録用
    dQout     = np.zeros(num_n)  # 各輸送管の各セルの排熱記録用
    fai_v_r   = np.zeros(num_n)  # LMパラメータ(φv）の各セルの排熱記録用
    dPv2f_r   = np.zeros(num_n)  # dPv2fの各セルの圧損記録用
    
    hL_r      = np.zeros(num_n)  # hLの各セルの圧損記録用
    hNu_r     = np.zeros(num_n)  # hNuの各セルの圧損記録用
    WeGT_r    = np.zeros(num_n)  # WeGTの各セルの圧損記録用
    Z_r       = np.zeros(num_n)  # Zの各セルの圧損記録用
    JG_r      = np.zeros(num_n)  # JGの各セルの圧損記録用
    
    dPall_r   = np.zeros(num_n)  # dPallの各セルの圧損記録用 
    dPf_r     = np.zeros(num_n)  # dPfの各セルの圧損記録用 
    dP2f_r    = np.zeros(num_n)  # dP2fの各セルの圧損記録用
    dPh_r     = np.zeros(num_n)  # dPhの各セルの圧損記録用
    dPa_r     = np.zeros(num_n)  # dPaの各セルの圧損記録用
    dPrat_r   = np.zeros(num_n)  # dPa/dPall割合の各セルの圧損記録用
    T_sur_r   = np.zeros(num_n)  # T_surの各セルの圧損記録用
    
    hin            = 0
    T_sur          = 0
    
    """
    # 仮変数を定義
    hin            = 0
    regime         = 0
    G_line_out     = 0
    G_line_out_sur = 0
    q              = 0
    dP2f           = 0
    dPv2f          = 0
    dPl2f          = 0
    dP_h           = 0
    dPa            = 0
    
    h_L  = 0
    h_LT = 0
    h_Nu = 0
    WeGT = 0
    pr   = 0
    Z    = 0
    G    = 0
    JG   = 0
    PrL  = 0
    ReLT = 0
    ReLO = 0
    judge_1 = 0
    judge_2 = 0
    judge_3 = 0
    judge_4 = 0
    
    htp0        = 0
    htp_minus90 = 0
    
    Re   = 0
    Re_V = 0
    Re_L = 0
    Pr   = 0
    f    = 0
    
    T = 0
    X = 0
    P = 0

    newX           = 0
    newV           = 0
    
    T_sur          = 0
    #T_sur     = np.zeros(num_n)
    """
    
    ###########################################################################
    # 　　　　各種管の各ノードの熱力学パラメータ計算のためのforループ
    ###########################################################################
    for i in range(0, num_n):
        
        if i == 0:
            T = T_in
            X = X_in
            P = P_in
            
            if X == 0:
                V = 0
            elif X == 1:
                V = 1
            else:
                V = (1 + VD(T)/LD(T) * (1-X)/X * 
                     (0.4 + 0.6 * ((VD(T)/LD(T) + 0.4 * (1-X)/X) / (1 + 0.4*(1-X)/X))**0.5))**(-1)
            
        else:
            T = new_T[i-1]
            X = new_X[i-1]
            P = new_P[i-1]
            V = new_V[i-1]


        #######################################################################
        #                     蒸気単相の場合における計算
        #######################################################################
        
        if X == 1:    # 蒸気単相の場合における計算
            Re = 4 * md /(np.pi * MuV(T) * ID)
            Re_V = 0            # 計算とは関係なし（記録用）
            Re_L = 0            # 計算とは関係なし（記録用）
            Pr = MuV(T) * VCp(T) / kV(T)
            regime = 0          #　凝縮レジーム判定用（単相なので意味はない）
            
            #print(T)
            
            if Re < 2300:       # 層流の場合
                f = 16 / Re
                hin = 4.36 * kV(T) / ID
            
            else:               # 乱流の場合
                f = 0.0791 * Re**(-0.25)
                if T > T_out:   # 管内流体が外気より高温な場合(流体自身が冷却される場合）
                                # ただし，Tout: 管外温度(=Tamb)
                    hin = 0.023* Re**(0.8) * Pr**(0.3) * kV(T) / ID
                    
                else:           # 管内流体が外気より低温な場合(流体自身が加熱される場合）
                    hin = 0.023* Re**(0.8) * Pr**(0.4) * kV(T) / ID
                    
            if component==1 or component==3:    # 蒸発管と液管の場合
                
                G_line_out = dL_n / (
                              (1 / (hin * np.pi * ID)) 
                             + (np.log(OD/ID) / (2 * np.pi * k_line))
                             + (1 / (h_cont * np.pi * OD))
                             + (np.log((OD + 2 * t_insu) / OD) / (2 * np.pi * k_insu))
                             + (1 / (h_out * np.pi * (OD + 2 * t_insu)))
                             )
            
            elif component==2:                  # 凝縮器の場合
                
                G_line_out = dL_n / (
                               ((1 / (hin * np.pi * ID)))
                             + ((np.log(4 * OD / ID) / (2 * np.pi * k_line)))      # 4*OD = 4*hi_clなので注意
                             + ((1 / (h_sink * 2 * wi_cl)))
                             )
            
            
            if component==2:
                G_line_out_sur = dL_n * h_sink * 2 * wi_cl
                T_sur = G_line_out * (T - T_out) / G_line_out_sur + T_out
            
            #print("G_lineout : ", G_line_out)
            
            dP_f = 32 * md**2 * f / (np.pi**2 * VD(T) * ID**5) * dL_n
            dP_h = VD(T) * g * dL_n * np.sin(theta * np.pi / 180)
            
            if i == 0:
                new_T[i] = T_in - G_line_out * (T - T_out) / (md * VCp(T))
                q        = G_line_out * (T - T_out)        # 長さdL_n部分の外気への放熱量:q [W]     
                new_P[i] = P_in - dP_f - dP_h               
                           # theta：各種管の傾斜角度[deg]
            
            else:
                new_T[i] = new_T[i-1] - G_line_out * (T - T_out) / (md * VCp(T))
                q        =  G_line_out * (T - T_out)
                new_P[i] = new_P[i-1] - dP_f - dP_h          
               
        
            # 凝縮器の場合，蒸気単相から気液二相への条件判定を行う． 
            if new_P[i] > P_sat(new_T[i]):
                X = 0.999999999
        
    
            if X == 0:
                V = 0
            elif X == 1:
                V = 1
            else:
                V = (1 + VD(T)/LD(T) * (1-X)/X  
                        * (0.4 + 0.6 * ((VD(T)/LD(T) + 0.4 * (1-X)/X) / (1 + 0.4*(1-X)/X))**0.5))**(-1)
            
            new_V[i]   = V
            new_X[i]   = X
            h_in[i]    = hin
            reg[i]     = regime
            Gcon[i]    = G_line_out
            dQout[i]   = q
            fai_v_r[i] = 0
            dPv2f_r[i] = 0
            
            hL_r[i]    = 0
            hNu_r[i]   = 0
            WeGT_r[i]  = 0
            Z_r[i]     = 0
            JG_r[i]    = 0
            
            dPf_r[i]   = dP_f
            dP2f_r[i]  = 0
            dPh_r[i]   = dP_h
            dPa_r[i]   = 0
            dPrat_r[i] = 0
            
            if i == 0:
                dPall_r[i] = P - new_P[i]
            else:
                dPall_r[i] = new_P[i-1] - new_P[i]
               
            if component==2:
                T_sur_r[i] = T_sur
            else:
                T_sur_r[i] = 0
              
        #######################################################################
        #                     気液二相の場合における計算
        #######################################################################
        
        elif X > 0 and X < 1:    # 気液二相の場合における計算　
            
            if (T - T_out) == 0:
                G_line_out = 0
                
            else:
                pr   = P / P_cri                   # 作動流体の臨界圧との比 
                G    = md / (np.pi * ID**2 / 4)    # 質量流量
                JG   = X * G / (g * ID * VD(T) * (LD(T) - VD(T)))**0.5
                Z    = (1/X - 1)**0.8 * pr**0.4
                PrL  = MuL(T) * LCp(T) / kL(T)
                ReLT = G * ID / MuL(T)             # Reynolds number for all mass flowing as liquid  
                ReLO = G * (1 - X) * ID / MuL(T)   # Reynolds number assuming liquid phase flowing alone
                
                h_LT = 0.023 * ReLT**(0.8) * PrL**(0.4) * kL(T) / ID
                
                h_L  = h_LT * (1 + 1.128 * X**0.817 * (LD(T) / VD(T))**0.3685 
                               * (MuL(T) / MuV(T))**0.2363 * (1 - MuV(T) / MuL(T))**2.144 * PrL**(-0.1))
                
                h_Nu = 1.32 * ReLO**(-1/3) * (LD(T) * (LD(T) - VD(T)) * g * kL(T)**3 / MuL(T)**2)**(1/3)
                
                judge_1 = 0.98 * (Z + 0.263)**(-0.62)
                judge_2 = 0.95 * (1.254 + 2.27 * Z**1.249)**(-1)
                judge_3 = 1 / (2.4 * Z + 0.73)
                judge_4 = 0.89 - 0.93 * np.exp(-0.087 * Z**(-1.17))
                
                WeGT    = G**2 * ID / (VD(T) * surf_T(T))
                
                # 以下，-30<=θ<=+90°(水平と鉛直上昇流も範囲内)のときの hin 計算
                if (theta >= -30) and (theta <= 90):
                    if (JG >= judge_1) and (WeGT > 100):    # Regime I
                        
                        hin = h_L
                        regime = 1
                
                    elif (JG <= judge_2) and (WeGT > 20):   # Regime III
                        
                        hin = h_Nu
                        regime = 3
                        
                    else:                                   # Regime II
                        
                        hin = h_L + h_Nu
                        regime = 2
                        
                        
                # 以下，θ=-90°鉛直下降のときの hin 計算    
                elif theta == -90:
                    if (JG >= judge_1) and (WeGT > 100):    # Regime I
                        
                        hin = h_L
                        regime = 1
                
                    elif (JG <= judge_2) and (WeGT > 20):   # Regime III
                        
                        hin = h_Nu
                        regime = 3
                        
                    else:                                   # Regime II
                        
                        hin = h_L + h_Nu
                        regime = 2
                
                    
                # 以下，θが -90<θ<-30°のときの hin の計算
                elif (theta > -90) and (theta < -30):
                    # 以下，水平と仮定した場合のhtp0を計算
                    if (JG >= judge_1) and (WeGT > 100):    # Regime I
                        htp0 = h_L
                        regime = 1
                    elif (JG <= judge_2) and (WeGT > 20):   # Regime III
                        htp0 = h_Nu
                        regime = 3
                    else:                                   # Regime II
                        htp0 = h_L + h_Nu
                        regime = 2
                    
                    # 以下，鉛直下降と仮定した場合のhtp_minus90を計算
                    if (JG >= judge_1) and (WeGT > 100):    # Regime I
                        htp_minus90 = h_L
                        regime = 1
                    elif (JG <= judge_2) and (WeGT > 20):   # Regime III
                        htp_minus90 = h_Nu
                        regime = 3
                    else:                                   # Regime II
                        htp_minus90 = h_L + h_Nu
                        regime = 2
                        
                    hin = htp0 + (htp0 - htp_minus90) * (theta + 30) / 60
                    
                    
                ###################################################################
                # ↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
                ###################################################################
            
                #####   #####   #####   #####   #####   #####   #####   #####

                if component==1 or component==3:        # 蒸発管と液管の場合
                    
                    G_line_out = dL_n / (
                                  (1 / (hin * np.pi * ID)) 
                                 + (np.log(OD/ID) / (2 * np.pi * k_line))
                                 + (1 / (h_cont * np.pi * OD))
                                 + (np.log((OD + 2 * t_insu) / OD) / (2 * np.pi * k_insu))
                                 + (1 / (h_out * np.pi * (OD + 2 * t_insu)))
                                 ) 
                                 # 各種管の長さdL_nにおける熱コンダクタンス:G_lineout [W/K]
                
                
                elif component==2:                      # 凝縮器の場合
                    
                    G_line_out = dL_n / (
                                   (1 / (hin * np.pi * ID))
                                 + (np.log(4 * OD / ID) / (2 * np.pi * k_line))      # 4*OD = 4*hi_clなので注意
                                 + (1 / (h_sink * 2 * wi_cl))
                                 )
                                 # 凝縮管の長さdL_nにおける熱コンダクタンス:G_lineout [W/K]
                    
            
            if component==2:
                G_line_out_sur = dL_n * h_sink * 2 * wi_cl
                T_sur = G_line_out * (T - T_sur) / G_line_out_sur + T_out
                
            
            newX = X - G_line_out * (T - T_out) / (md * HV(T))
            
            
            if newX < 0.0001:
                newX = 0
            elif newX > 1:
                newX = 1
                
            if newX == 0:
                newV = 0
            elif newX == 1:
                newV = 1
            else:
                newV = (1 + VD(T) / LD(T) * (1 - newX) / newX
                        * (0.4 + 0.6 * ((VD(T) / LD(T) + 0.4 * (1 - newX) / newX) 
                        / (1 + 0.4 * (1 - newX) / newX))**0.5))**(-1)
            
            
            
            ###################################################################
            
            Re = 0     # 計算とは関係なし（記録用）
            Re_V = 4 * md * X / (np.pi * MuV(T) * ID)
            Re_L = 4 * md * (1 - X) / (np.pi * MuL(T) * ID)
            
            if (Re_V < 2000) and (Re_L < 2000):
                C = 5
                dPv2f = 128 * md * X * MuV(T) * dL_n / (np.pi * VD(T) * ID**4)          # 見掛け圧力損失
                dPl2f = 128 * md * (1 - X) * MuL(T) * dL_n / (np.pi * LD(T) * ID**4)    # 見掛け圧力損失
            
            elif (Re_V < 2000) and (Re_L >= 2000):
                C = 10
                dPv2f = 128 * md * X * MuV(T) * dL_n / (np.pi * VD(T) * ID**4)  
                dPl2f = 32 * (md * (1 - X))**2 * 0.0791 * Re_L**(-0.25) * dL_n / (np.pi * LD(T) * ID**5)
                
            elif (Re_V >= 2000) and (Re_L < 2000):
                C = 12
                dPv2f = 32 * (md * X)**2 * 0.0791 * Re_V**(-0.25) * dL_n / (np.pi * VD(T) * ID**5)
                dPl2f = 128 * md * (1 - X) * MuL(T) * dL_n / (np.pi * LD(T) * ID**4)
                
            elif (Re_V >= 2000) and (Re_L >= 2000):
                C = 20
                dPv2f = 32 * (md * X)**2 * 0.0791 * Re_V**(-0.25) * dL_n / (np.pi * VD(T) * ID**5)
                dPl2f = 32 * (md * (1 - X))**2 * 0.0791 * Re_L**(-0.25) * dL_n / (np.pi * LD(T) * ID**5)
                
            
            # 二相摩擦損失: dP2f[Pa]の計算（※ロックハート＆マルチネリ相関式）
            if md == 0:
                dP2f = 0
            else:
                dP2f = (1 + C * (dPl2f / dPv2f)**0.5 + dPl2f / dPv2f) * dPv2f
                
            dP_h = (VD(T) * V + LD(T) * (1 - V)) * g * dL_n * np.sin(theta / 180 * np.pi)
            
            
            if newV == 0:
                dPa = 0
            else:
                dPa = G**2 * (  (X**2 / VD(T) / V + (1 - X)**2 / (LD(T) * (1 - V)))
                              - (newX**2 / VD(T) / newV + (1 - newX)**2 / (LD(T) * (1 - newV)))  )
                              
            
            if i == 0:
                new_P[i] = P_in - dP2f - dP_h - dPa
            else:
                new_P[i] = new_P[i-1] - dP2f - dP_h - dPa
            
            q = G_line_out * (T - T_out)
            
            new_V[i] = newV
            new_X[i] = newX
            new_T[i] = T_sat(new_P[i])      # 二相流の場合の管内流体温度:new_T(i)は常に飽和温度と仮定している．
            h_in[i]  = hin
            reg[i]   = regime
            Gcon[i]  = G_line_out
            dQout[i] = q
            fai_v_r[i] = 1 + C * (dPl2f / dPv2f)**0.5 + dPl2f / dPv2f
            dPv2f_r[i] = dPv2f
            
            hL_r[i]   = h_L
            hNu_r[i]  = h_Nu
            WeGT_r[i] = WeGT
            Z_r[i]    = Z
            JG_r[i]   = JG
            
            dPf_r[i]  = 0
            dP2f_r[i] = dP2f
            dPh_r[i]  = dP_h
            dPa_r[i]  = dPa
            
            if i == 0:
                dPall = P - new_P[i]
            else:
                dPall = new_P[i-1] - new_P[i]
                
            dPall_r[i] = dPall
            dPrat_r[i] = dPa / dPall
            
            if component==2:
                T_sur_r[i] = T_sur
            else:
                T_sur_r[i] = 0
            
        #######################################################################
        #                      液単相の場合における計算
        #######################################################################    
                
        elif X == 0:    # 液単相の場合における計算
            
            Re = 4 * md / (np.pi * MuL(T) * ID)
            Re_V = 0
            Re_L = 0
            Pr   = MuL(T) * LCp(T) / kL(T)
            regime = 0
            
            
            if Re < 2300:       # 層流の場合
                f = 16 / Re
                hin = 4.36 * kL(T) / ID
            
            else:               # 乱流の場合
                f = 0.0791 * Re**(-0.25)
                if T > T_out:   # 管内流体が外気より高温な場合(流体自身が冷却される場合）
                                # ただし，Tout: 管外温度(=Tamb)
                    hin = 0.023* Re**0.8 * Pr**0.3 * kL(T) / ID
                    
                else:           # 管内流体が外気より低温な場合(流体自身が加熱される場合）
                    hin = 0.023* Re**0.8 * Pr**0.4 * kL(T) / ID
                    
                
            
            if component==1 or component==3:    # 蒸発管と液管の場合
                
                G_line_out = dL_n / (
                               (1 / (hin * np.pi * ID)) 
                             + (np.log(OD/ID) / (2 * np.pi * k_line))
                             + (1 / (h_cont * np.pi * OD))
                             + (np.log((OD + 2 * t_insu) / OD) / (2 * np.pi * k_insu))
                             + (1 / (h_out * np.pi * (OD + 2 * t_insu)))
                             )
            
            elif component==2:                  # 凝縮器の場合
                
                G_line_out = dL_n / (
                               (1 / (hin * np.pi * ID))
                             + (np.log(4 * OD / ID) / (2 * np.pi * k_line))      # 4*OD = 4*hi_clなので注意
                             + (1 / (h_sink * 2 * wi_cl))
                             )
            
            if component==2:
                G_line_out_sur = dL_n * h_sink * 2 * wi_cl
                T_sur = G_line_out * (T - T_out) / G_line_out_sur + T_out
                
            
            dP_f = 32 * md**2 * f / (np.pi**2 * LD(T) * ID**5) * dL_n
            dP_h = LD(T) * g * dL_n * np.sin(theta * np.pi / 180)     # theta：各種管の傾斜角度[deg]    
            
            if i == 0:
                new_T[i] = T_in - G_line_out * (T - T_out) / (md * LCp(T))
                q        = G_line_out * (T - T_out)
                new_P[i] = P_in - dP_f - dP_h
            
            else:
                new_T[i] = new_T[i-1] - G_line_out * (T - T_out) / (md * LCp(T))
                q        = G_line_out * (T - T_out)
                new_P[i] = new_P[i-1] - dP_f - dP_h
             
            new_V[i]   = 0
            new_X[i]   = 0
            h_in[i]    = hin
            reg[i]     = regime
            Gcon[i]    = G_line_out
            fai_v_r[i] = 0
            dPv2f_r[i] = 0
            
            hL_r[i]    = 0
            hNu_r[i]   = 0
            WeGT_r[i]  = 0
            Z_r[i]     = 0
            JG_r[i]    = 0
            
            dPf_r[i]   = dP_f
            dP2f_r[i]  = 0
            dPh_r[i]   = dP_h
            dPa_r[i]   = 0
            dPrat_r[i] = 0
            
            if component==2:
                T_sur_r[i] = T_sur
            else:
                T_sur_r[i] = 0
            
            if i == 0:
                dPall_r[i] = P - new_P[i]
            else:
                dPall_r[i] = new_P[i-1] - new_P[i]
                 
        #######################################################################
        #######################################################################
        
        new_Re[i]   = Re
        new_Re_V[i] = Re_V
        new_Re_L[i] = Re_L
        dQout[i]    = q
        
        Q_out = Q_out + q    # 漏れ熱量: Qout [W]
        
        
    return [new_P, new_T, new_X, new_V, Q_out, new_Re, new_Re_V, new_Re_L, h_in, reg, Gcon, 
            dQout, fai_v_r, dPv2f_r, hL_r, hNu_r, WeGT_r, Z_r, JG_r,
            dP2f_r, dP2f_r, dPh_r, dPa_r, dPall_r, dPrat_r, T_sur_r]
                
    
    
def cc_void(P_ll, V_vl, V_cl, V_ll, T_wickout, T_vl, T_cl, T_ll, ID_vl, ID_cl, ID_ll,
            N_vl, N_cl, N_ll, L_vl, L_cl, L_ll, p_hi,
            EC_l, EC_w, EC_h, EC_t_x, EC_t_y, wick_l, wick_w, wick_h,
            w_g, l_g, t_g, n_g, V_fluid, N_evap):
    
    ###########################################################################
    #                       インプットパラメータ 
    ###########################################################################
    
    # なし
    
    ###########################################################################
    
    N_vl = len(T_vl)
    N_cl = len(T_cl)
    N_ll = len(T_ll)
    
    m_lines = 0

    # 蒸発管内の流体の（液と蒸気合わせた）総質量の計算
    for i in range(0, N_vl):
        m_lines = m_lines + ((VD(T_vl[i]) * V_vl[i] + LD(T_vl[i]) * (1 - V_vl[i])) 
                             * np.pi * ID_vl**2 * 0.25 * L_vl / N_vl)
    
    # 凝縮管内の流体の（液と蒸気合わせた）総質量の計算
    for i in range(0, N_cl):
        m_lines = m_lines + ((VD(T_cl[i]) * V_cl[i] + LD(T_cl[i]) * (1 - V_cl[i])) 
                             * np.pi * ID_cl**2 * 0.25 * L_cl / N_cl)
    
    # 液管内の流体の（液と蒸気合わせた）総質量の計算
    for i in range(0, N_ll):
        m_lines = m_lines + ((VD(T_ll[i]) * V_ll[i] + LD(T_ll[i]) * (1 - V_ll[i])) 
                             * np.pi * ID_ll**2 * 0.25 * L_ll / N_ll)
        
    # 蒸発器内の（液と蒸気合わせた）総質量: m_evap [kg] の計算
    # 平板型蒸発器用
    m_evap = (LD(T_wickout) * (wick_l * wick_w * wick_h) * p_hi
              * VD(T_wickout) * (w_g * l_g * n_g)) * N_evap
    

    # 封入量から計算した作動流体の総質量: m_fluid [kg]
    # 作動流体の封入量: V_fluid [m^3] ただし，-35 [degC]の場合
    m_fluid = V_fluid * LD(-35 + 273.15)
    
    # CCの総質量: m_cc [kg]
    m_cc = m_fluid - m_lines - m_evap
    
    if m_cc < 0:
        m_cc = 0
    
    # CCの容積: Vol_cc [m^3]
    Vol_cc = (EC_l -2 * EC_t_x) * (EC_w - EC_t_x) * (EC_h - 2 * EC_t_y - wick_h - t_g) * N_evap
    
    
    # CCのボイド率: V_cc [-]
    P_ccin = P_ll[-1]
    T_ccin = T_sat(P_ccin)    
    V_cc   = (m_cc - LD(T_ccin) * Vol_cc) / ((VD(T_ccin) - LD(T_ccin)) * Vol_cc)
    
    if V_cc < 0:
        V_cc = 0
        
    
    return [V_cc, P_ccin, T_ccin, m_fluid, m_lines, m_evap, m_cc]
    
    

def Q_w(n, L_n, Tt_s, T_ec, w_e, l_g, n_g, p_hi, k_wick):
    
    dh = L_n / (n - 1)          # 格子間距離( =L_m /(m-1) )
    a1 = (w_e / 2 / dh) + 1     # 溝と伝熱面の境界（単位：格子番号） 
    a1 = int(a1)
    
    dTc_dy = (Tt_s[1, 1:(a1-1)] - Tt_s[2, 1:(a1-1)]) / dh
    sum_dTc_dy = np.sum(dTc_dy)
    q_w_leak = lam_w(T_ec, p_hi, k_wick) * sum_dTc_dy / (a1 - 1)    # 伝熱面からウィックに伝わる熱流束[W/m2]
    Q_w_leak = q_w_leak * (w_e * l_g * n_g)                 # 伝熱面からウィックに伝わる熱量[W]
    
    return [Q_w_leak, q_w_leak]



def Q_ccinll(md, T_ll, T_ccin):
    
    Q_sub = md * LCp(T_ccin) * (T_ccin - T_ll[-1])
    
    return Q_sub



def lam_w(T, p_hi, k_wick):
    
    k_max = p_hi * kL(T) + (1 - p_hi)*k_wick
    k_min = kL(T) * k_wick /(p_hi * k_wick + (1 - p_hi) * kL(T))
    lam_w = (k_max**0.42) * (k_min**0.58)
    
    return lam_w




###############################################################################
# 計算結果のプロット
###############################################################################
def graphing(wb_filename):
    
    plt.rcParams["font.family"] = "Times New Roman"      #全体のフォントを設定
    plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
    plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
    plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
    plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
    plt.rcParams["xtick.major.width"] = 1.5              #x軸主目盛り線の線幅
    plt.rcParams["ytick.major.width"] = 1.5              #y軸主目盛り線の線幅
    plt.rcParams["xtick.minor.width"] = 1.0              #x軸補助目盛り線の線幅
    plt.rcParams["ytick.minor.width"] = 1.0              #y軸補助目盛り線の線幅
    plt.rcParams["xtick.major.size"] = 7.5               #x軸主目盛り線の長さ
    plt.rcParams["ytick.major.size"] = 7.5               #y軸主目盛り線の長さ
    plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
    plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ
    plt.rcParams["font.size"] = 12                       #フォントの大きさ
    plt.rcParams["axes.linewidth"] = 1.5                 #囲みの太さ
    plt.rcParams['axes.grid'] = True
    
    ###############################################################################


    fig = plt.figure(figsize=(6.5, 5.0), dpi=1000)
    ax1 = fig.add_subplot(1, 1, 1)
    
    plt.rcParams["mathtext.fontset"] = "stix"

    ###############################################################################

    #sheetname = "LHP_result_space"    
    #bookname = sheetname + ".xlsx"

    df = pd.read_excel(wb_filename, skiprows=0)  
    df_array = np.array(df)

    ###############################################################################
    
    
    Q_apply = df_array[:,0]
    T_HB = df_array[:,20]
    T_EC = df_array[:,21]
    T_CON_ave  = df_array[:,31]
    
    P_cap   = df_array[:,18]
    dP_wick = df_array[:,44]
    dP_gr   = df_array[:,45]
    dP_vl   = df_array[:,46]
    dP_cl   = df_array[:,47]
    dP_ll   = df_array[:,48]

    R_LHP = np.zeros_like(T_HB)

    R_LHP = (T_EC - T_CON_ave) / Q_apply

    ###############################################################################
    ## 定常温度の推移グラフ    
        
    fig = plt.figure(figsize=(6.5, 5.0), dpi=1000)
    ax1 = fig.add_subplot(1, 1, 1)

    plt.rcParams["mathtext.fontset"] = "stix"


    ax1.plot(Q_apply, T_HB, "-", color="m")
    ax1.plot(Q_apply, T_HB, ".", label="Heater", color="m")

    ax1.plot(Q_apply, T_EC, "-", color="r")
    ax1.plot(Q_apply, T_EC, ".", label="EVA", color="r")

    ax1.plot(Q_apply, dP_wick, "o-", label="dP Wick", color="green")

    ax1.plot(Q_apply, T_CON_ave, "-", color="b")
    ax1.plot(Q_apply, T_CON_ave, ".", label="CON", color="b")

    ax1.set_xlabel("Heat load, W", fontsize=14)
    ax1.set_ylabel("Temperature, °C", fontsize=14)
    ax1.set_xlim([0, 16])
    ax1.set_ylim([-40, -20])
    ax1.grid()

    """
    ax2 = ax1.twinx()

    ax2.plot(Q_apply, T_EC-T_CON_ave, ":", label="R", color="navy", zorder=1, linewidth=2.0)
    ax2.plot(Q_apply, T_EC-T_CON_ave, "^", label="R", color="navy", zorder=1, linewidth=2.0)
    ax2.set_ylabel("Temperature difference, °C", fontsize=14)
    ax2.set_ylim([0, 20])
    ax2.grid()
    """

    plt.xticks(fontsize=11)

    #handler1, label1 = ax1.get_legend_handles_labels()
    #handler2, label2 = ax2.get_legend_handles_labels()
    #ax1.legend(ncol=2, bbox_to_anchor=(0.01, 1), loc='upper left', frameon=False)

    plt.savefig("LHP_result_Steady.png", bbox_inches="tight", dpi=1000)
    plt.show()
    
    ###############################################################################
    ## 毛細管圧力と圧力損失の推移グラフ

    fig = plt.figure(figsize=(6.5, 5.0), dpi=1000)
    ax1 = fig.add_subplot(1, 1, 1)

    plt.rcParams["mathtext.fontset"] = "stix"

    ax1.plot(Q_apply, P_cap, ".", label=r'$P_{cap}$', color="black", zorder=1, linewidth=2.0)
    ax1.plot(Q_apply, P_cap, "-", color="black", zorder=1, linewidth=2.0)

    ax1.bar(Q_apply, dP_wick, bottom=0, label=r'$\Delta P_{wick}$', color="green", width=0.4)
    ax1.bar(Q_apply, dP_ll, bottom=dP_wick, label=r'$\Delta P_{ll}$', color="blue", width=0.4)
    ax1.bar(Q_apply, dP_cl, bottom=dP_wick+dP_ll , label=r'$\Delta P_{cl}$', color="yellow", width=0.4)
    ax1.bar(Q_apply, dP_vl, bottom=dP_wick+dP_ll+dP_cl , label=r'$\Delta P_{vl}$', color="orange", width=0.4)
    ax1.bar(Q_apply, dP_gr, bottom=dP_wick+dP_ll+dP_cl+dP_vl , label=r'$\Delta P_{gr}$', color="red", width=0.4)

    ax1.set_xlabel("Heat load, W", fontsize=14)
    ax1.set_ylabel("Pressure, Pa", fontsize=14)
    ax1.set_xlim([0, 16])
    ax1.set_ylim([0, 10000])
    ax1.grid()

    """
    ax2 = ax1.twinx()

    ax2.plot(Q_apply, P_cap, "^", label=r'$P_{cap}$', color="black", zorder=1, linewidth=2.0)
    ax2.plot(Q_apply, P_cap, "-", color="black", zorder=1, linewidth=2.0)
    ax2.set_ylabel("Pressure, Pa", fontsize=14)
    ax2.set_ylim([0, 35000])
    ax2.grid()
    """

    plt.xticks(fontsize=14)

    #handler1, label1 = ax1.get_legend_handles_labels()
    #handler2, label2 = ax2.get_legend_handles_labels()
    ax1.legend(ncol=2, bbox_to_anchor=(0.01, 0.65), fontsize=14, loc='upper left', frameon=False)

    plt.savefig("LHP_result_Pcap.png", bbox_inches="tight", dpi=1000)
    plt.show()

    ###############################################################################
    
    
    

    
###############################################################################
# 作動流体の物性値(アンモニア)
###############################################################################
def HV(T):
    # Latent heat for ammonia in the temperture range from 200 to 400K.
    # T [K], HV [J/kg]
    if T < 350:
        HV = (-9.87238213401298000000 *1e-9* T**6 
              +1.47058915143106000000 *1e-5* T**5 
              -9.18993946714153000000 *1e-3* T**4 
              +3.05640971759590000000 *1e0 * T**3 
              -5.74426150475301000000 *1e2 * T**2 
              +5.60395267583865000000 *1e4 * T**1 
              -5.78465373705890000000 *1e5)
    
    else:
        HV = (-5.76355614967383000000 *1e-5* T**6 
              +1.27867330755187000000 *1e-1* T**5 
              -1.18187176683415000000 *1e2 * T**4 
              +5.82537962176215000000 *1e4 * T**3 
              -1.61487357222268000000 *1e7 * T**2 
              +2.38717036084949000000 *1e9 * T**1 
              -1.47007522866971000000 *1e11)
    
    return HV

  
def kL(T):
    # Liquid thermal conductivity for ammonia in the temperture range from 200 to 400K.
    # T [K], kL[W/m/K]
    kL=( 3.59969357043530000000 * 1e-16 * T**6 
        -4.15633795258152000000 * 1e-13 * T**5 
        +4.71813471135146000000 * 1e-11 * T**4 
        +1.08900791944389000000 * 1e-7  * T**3 
        -5.13724520074396000000 * 1e-5  * T**2 
        +5.07474187790046000000 * 1e-3  * T**1 
        +1.00641284828102000000 * 1e0 )
    
    return kL
    

def kV(T):
    # Vapor thermal conductivity for ammonia in the temperture range from 200 to 400K.
    # T [K], kV[W/m/K]
    if T < 350:    
        kV=( 2.29968681030093000000 * 1e-15 * T**6 
            -3.53876775076550000000 * 1e-12 * T**5 
            +2.26877986954517000000 * 1e-9  * T**4
            -7.74064969209584000000 * 1e-7  * T**3
            +1.48520324606629000000 * 1e-4  * T**2
            -1.52162379277920000000 * 1e-2  * T**1 
            +6.69823510176668000000 * 1e-1 )
        
    else:          
        kV=( 2.30770078966377000000 * 1e-11 * T**6 
            -5.09596846670998000000 * 1e-8  * T**5
            +4.68848476341596000000 * 1e-5  * T**4
            -2.30038370864699000000 * 1e-2  * T**3 
            +6.34816575224020000000 * 1e0   * T**2 
            -9.34216805827399000000 * 1e2   * T**1 
            +5.72775330958816000000 * 1e4 )
    
    return kV


def LCp(T):      
    # Liquid Cp for ammonia in the temperture range from 200 to 400K.
    # T [K], LCp[J/kg/K]
    if T < 350:   
        LCp=( 4.00047340631642000000 * 1e-10 * T**6 
             -6.31493785827400000000 * 1e-7  * T**5
             +4.16234019511637000000 * 1e-4  * T**4
             -1.46149838816783000000 * 1e-1  * T**3 
             +2.87542697467256000000 * 1e1   * T**2 
             -2.99394597415278000000 * 1e3   * T**1
             +1.32545347836710000000 * 1e5)
        
    else:         
        LCp=( 2.90332845096497000000 * 1e-5  * T**6 
             -6.46348153723758000000 * 1e-2  * T**5 
             +5.99392121770294000000 * 1e1   * T**4 
             -2.96372668403237000000 * 1e4   * T**3 
             +8.24081781141665000000 * 1e6   * T**2 
             -1.22175046310063000000 * 1e9   * T**1 
             +7.54507795299835000000 * 1e10 )
        
    return LCp

    
def VCp(T):
    # Vapor Cp for ammonia in the temperture range from 200 to 400K.
    # T [K], VCp[J/kg/K]
    if T < 350:       
        VCp=( 6.38451007448256000000 * 1e-10 * T**6 
             -9.75933125332493000000 * 1e-7  * T**5 
             +6.20867976304444000000 * 1e-4  * T**4 
             -2.10030952299074000000 * 1e-1  * T**3
             +3.98696860547729000000 * 1e1   * T**2 
             -4.02789074249468000000 * 1e3   * T**1 
             +1.71164579579488000000 * 1e5 )
        
    else:
        VCp=( 1.71347800179689000000 * 1e-5 * T**6  
             -3.79091705671829000000 * 1e-2 * T**5
             +3.49408323780249000000 * 1e1  * T**4 
             -1.71732280643692000000 * 1e4  * T**3 
             +4.74702053883247000000 * 1e6  * T**2 
             -6.99703820757237000000 * 1e8  * T**1 
             +4.29654085122100000000 * 1e10)
        
    return VCp
    
    
def LD(T):
    # Liquid densityt for ammonia in the temperture range from 200 to 400K.
    # T [K], LD[kg/m3]
    if T < 350:
        LD=( -2.66868904091108000000 * 1e-12 * T**6
             +4.02344287509133000000 * 1e-9  * T**5
             -2.58317259468985000000 * 1e-6  * T**4 
             +9.00353871930187000000 * 1e-4  * T**3 
             -1.80768911395490000000 * 1e-1  * T**2 
             +1.87941960014401000000 * 1e1   * T**1 
             +1.35656461569591000000 * 1e1)
        
    else:
        LD= (-2.08994748035485000000 * 1e-8  * T**6 
             +4.63881726827847000000 * 1e-5  * T**5 
             -4.28951936810615000000 * 1e-2  * T**4 
             +2.11515535938156000000 * 1e1   * T**3 
             -5.86578876742799000000 * 1e3   * T**2 
             +8.67425098148504000000 * 1e5   * T**1 
             -5.34366028151538000000 * 1e7)
        
    return LD


def VD(T):
    # Vapor densityt for ammonia in the temperture range from 200 to 400K.
    # T [K], LD[kg/m3]
    if T < 350:
        VD=  (2.77699409174847000000 * 1e-12 * T**6 
             -4.21142076299227000000 * 1e-9  * T**5 
             +2.68614077033097000000 * 1e-6  * T**4 
             -9.13007872375855000000 * 1e-4  * T**3  
             +1.73801285279068000000 * 1e-1  * T**2  
             -1.75527635130298000000 * 1e1   * T**1  
             +7.34762347834233000000 * 1e2)
        
    else:
        VD=  (1.73015967767920000000 * 1e-8  * T**6 
             -3.83704884464729000000 * 1e-5  * T**5 
             +3.54532538928763000000 * 1e-2  * T**4  
             -1.74688261237741000000 * 1e1   * T**3  
             +4.84102237613608000000 * 1e3   * T**2  
             -7.15396082779866000000 * 1e5   * T**1  
             +4.40430173102304000000 * 1e7)
        
    return VD
    

def MuL(T):
    # Liquid viscosity for ammonia in the temperture range from 200 to 400K.
    # T [K], MuL[Pa s]
    if T < 350:
        MuL=( 7.75690552859348000000 * 1e-17 * T**6 
             -1.39766522433775000000 * 1e-13 * T**5 
             +1.05177529474923000000 * 1e-10 * T**4 
             -4.23802644367013000000 * 1e-8  * T**3 
             +9.66662252625950000000 * 1e-6  * T**2 
             -1.18813765591077000000 * 1e-3  * T**1 
             +6.19885370282222000000 * 1e-2)
        
    else:
        MuL=(-1.32560992659653000000 * 1e-15 * T**6 
             +2.92423299476449000000 * 1e-12 * T**5 
             -2.68767697870695000000 * 1e-9  * T**4 
             +1.31736996143760000000 * 1e-6  * T**3 
             -3.63172961924272000000 * 1e-4  * T**2 
             +5.33886876774134000000 * 1e-2  * T**1 
             -3.26935368216519000000 * 1e0)
        
    return MuL


def MuV(T):
    # Vapor viscosity for ammonia in the temperture range from 200 to 400K.
    # T [K], MuV[Pa s]
    if T < 350:
        MuV=( 2.14262475025136000000 * 1e-19 * T**6 
             -3.22072506175791000000 * 1e-16 * T**5 
             +2.03802614095042000000 * 1e-13 * T**4 
             -6.96357202101342000000 * 1e-11 * T**3 
             +1.35908003907145000000 * 1e-8  * T**2 
             -1.41147977811128000000 * 1e-6  * T**1 
             +6.59679830684079000000 * 1e-5 )
        
    else:
        MuV=( 8.71409333316924000000 * 1e-16 * T**6 
             -1.92084609316123000000 * 1e-12 * T**5 
             +1.76426949993459000000 * 1e-9  * T**4 
             -8.64246059127046000000 * 1e-7  * T**3 
             +2.38135446225062000000 * 1e-4  * T**2 
             -3.49939053630324000000 * 1e-2  * T**1 
             +2.14252912743510000000 * 1e0 )
        
    return MuV
    

    
def P_sat(T):
    # Saturation Pressure for ammonia in the temperture range from 200 to 400K.
    # T [K], Psat[Pa]
    P_sat=( 9.0717585748967300000 * 1e-8 * T**6 
          -1.57624767606723000000 * 1e-4 * T**5 
          +1.16515614969311000000 * 1e-1 * T**4 
          -4.50913119640308000000 * 1e1  * T**3 
          +9.57485000669399000000 * 1e3  * T**2 
          -1.06046740760258000000 * 1e6  * T**1 
          +4.80493435057333000000 * 1e7)
            
    return P_sat

    
def surf_T(T):
    # Surface tension for ammonia in the temperture range from 200 to 400K.
    # T [K], surfT[N/m]
    surf_T =( 2.22680989169017000000 * 1e-16 * T**6 
             -3.57009512955957000000 * 1e-13 * T**5 
             +2.32186811115964000000 * 1e-10 * T**4 
             -7.64228182540059000000 * 1e-8  * T**3 
             +1.29109542970105000000 * 1e-5  * T**2 
             -1.20466411539778000000 * 1e-3  * T**1 
             +1.07377689046167000000 * 1e-1)
    
    return surf_T
    
    
    
def T_sat(P_sat):
    # Saturation temperature for ethanol in the temperture range from 250 to 500K.
    # P [Pa], T_sat[K]
    #R2=0.999999930075995
    T_sat = (-6.64639787828492000000 * 1e-4 * (np.log(P_sat))**6 
             +5.01927020849770000000 * 1e-2 * (np.log(P_sat))**5 
             -1.54712285954372000000 * 1e0  * (np.log(P_sat))**4 
             +2.51629313826281000000 * 1e1  * (np.log(P_sat))**3 
             -2.27612328645078000000 * 1e2  * (np.log(P_sat))**2 
             +1.09619751317666000000 * 1e3  * (np.log(P_sat))**1 
             -2.03382600846871000000 * 1e3)
    
    return T_sat


    