import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import defopt
import re
from scipy import interpolate
import sys
sys.path.insert(1, '/nesi/project/uoo03104/.conda/envs/xesmf_stable_env/lib/python3.7/site-packages/cmcrameri/')
import cm
#add .plot(cmap=cm.hawaii) for cb friendly
import preprocess_xsect as prep
#warnings.filterwarnings('ignore')
def dataframe_to_datetime(d):
    d['Datetime'] = pd.to_datetime(d['date'] + ' ' + d['hour'])
    d = d.set_index('Datetime')
    d = d.drop(['date','hour'], axis=1)
    d['date']=d.index
    return d



def plot_timeseries(*, save_dir: str='/nesi/project/uoo03104/snakemake_output/Taylor200_glac_update1_2yrloop/DEC18', station_name: str='cwg', plot_name: str='precip', date: str='2018-12-01 04:00:00'):
    """
    Plot timeseries of modelled crocus energy balance
    
    @param  save_dir directory to save timeseries png
    @param station_name LTER network name of stream gauge

    """

    cohm = pd.read_table('/nesi/nobackup/uoo03104/validation_data/COHM_MB.txt', delim_whitespace=True, header=0)
    cohm = dataframe_to_datetime(cohm)
    cohm = cohm.loc[cohm.index >= '2018-12-01 13:00:00', :]
    cohm = cohm.loc[cohm.index <= '2019-01-01 12:00:00', :]
    cohm.index = cohm.index.tz_localize('Antarctica/Mcmurdo').tz_convert('UTC')
    cohm["hsnow(obs,m)_scaled"] = cohm["hsnow(obs,m)"] - cohm["hsnow(obs,m)"][0]

    cohm_aws = pd.read_table('/nesi/project/uoo03104/COHM_AWS.txt', delim_whitespace=True, header=0)
    cohm_aws = dataframe_to_datetime(cohm_aws)
    cohm_aws = cohm_aws.loc[cohm_aws.index >= '2018-12-01 13:00:00', :]
    cohm_aws = cohm_aws.loc[cohm_aws.index <= '2019-01-01 12:00:00', :]
    cohm_aws.index = cohm_aws.index.tz_localize('Antarctica/Mcmurdo').tz_convert('UTC')

    c = pd.read_csv('/nesi/nobackup/uoo03104/validation_data/long_aws_cwg.csv',delimiter=',',sep='\t', header=0, skiprows=[0,2,3])
    c = c.set_index('TIMESTAMP')
    c.index = pd.to_datetime(c.index)
    c.index = c.index.tz_localize('Antarctica/Mcmurdo').tz_convert('UTC')
    c = c.astype(float)
    hsnow = pd.DataFrame()
    hsnow["SR50T"] = -1*c["SR50T_Avg"].loc[c.index >= '2021-12-01 13:00:00']
    hsnow["hsnow"] = hsnow["SR50T"] - hsnow["SR50T"][0]

    #precip = hsnow["hsnow"].diff()
    time_mask = (hsnow.index.hour == 00) & (hsnow.index.minute == 00) #filter each day
    precip = hsnow[time_mask].diff()["hsnow"] #daily diff in sfc height
    precip[precip < 0.] = 0.
    precip = precip *100. #convert to swe using 250kg/m3 density
    precip = precip.resample('H').ffill()

    daily_sfcheight = hsnow[time_mask].diff()["hsnow"].cumsum()
    daily_sfcheight.loc["2021-12-02 00:00:00+00:00"] = 0.0
    
    SWin = c["incommingSW_Avg"].resample('D').sum()
    SWout = c["outgoingSW_Avg"].resample('D').sum()
    Alb = SWout/SWin
    Alb = Alb.resample('H').ffill()

    lwd = c["incomingLW_Avg"].loc["2021-12-01 13:00:00+00:00":"2022-01-01 00:00:00+00:00"]
    lwu = c["outgoingLW_Avg"].loc["2021-12-01 13:00:00+00:00":"2022-01-01 00:00:00+00:00"]
    airT = c["AirTC_Avg"] + 273.15
    airT = airT.loc["2021-12-01 13:00:00+00:00":"2022-01-01 00:00:00+00:00"]
    sigma = 5.67e-8
    eps = 0.98
    tsfc_1 = ((1/sigma)*lwu)**(0.25)
    tsfc_98 = ((1/(eps*sigma))*(lwu - lwd + (eps*lwd)))**(0.25)

    df = pd.read_csv(f'{save_dir}/timeseries_ldasout_{station_name}.csv', index_col=0)
    df.index = pd.to_datetime(df.index)
    df.index = df.index.tz_localize('UTC')
    hsnow_m = pd.DataFrame()
    hsnow_m["SNOWH"] = df["SNOWH"].loc[df.index>='2021-12-01 13:00:00']
    hsnow_m["SNOWH_scaled"] = hsnow_m["SNOWH"] - hsnow_m["SNOWH"][0]

    df["ACCPRCP_scaled"] = df["ACCPRCP"] - df["ACCPRCP"][0]
    df["ACSNOM_scaled"] = df["ACSNOM"] - df["ACSNOM"][0]
    df["PSNOWTOTSWE_scaled"] = df["PSNOWTOTSWE"] - df["PSNOWTOTSWE"][0]

    #df["subl"] = (df["LH"]*(3600/2838200))  #old subl 
    df["subl"] = (df["LH"]*(3600/(df['PSNOWRHO0']*2.8345e6))) #same as croc

    df['CUMSUM_ACCPRCP'] = df['ACCPRCP'].cumsum()
    df['CUMSUM_subl'] = df['subl'].cumsum()
    df['CUMSUM_subldrift'] = df['PSNDRIFT'].cumsum()
    df['CUMSUM_ACSNOM'] = df['ACSNOM_scaled'].cumsum()
    
    #en = pd.read_csv('../energybal/middle_energybal.csv', index_col=0)
    #en.index = en.index.tz_localize('UTC')
    #df["subl"] = (en["LH"]*(3600/2838200))
    if plot_name=='albdiognosis':
        df2 = df["2021-12-01 13:00:00+00:00":]
        ftsize=14
        fig, axs = plt.subplots(3, 1, sharex=True, figsize=(18,14))
        fig.suptitle('Albedo Comparison', fontsize=ftsize) 
        axs[0].plot(df2.index, df2["ALBEDO"], color="grey", label="albedo")
        axs[0].legend(loc="upper right")
        axs[1].plot(df2.index, df2["PSNOWRHO0"], color='blue', label='density')
        axs[1].legend(loc="upper right")
        axs[2].plot(df2.index, hsnow_m["SNOWH_scaled"], color='orange', label="snow height")
        axs[2].legend(loc="upper right")
        plt.legend()
        plt.savefig(f'{save_dir}/timeseries_albdiognosis_{station_name}.png')

    if plot_name=='tsfc':
        plt.figure(figsize=[12,7])
        (tsfc_98.resample('H').mean()).plot(label="obs")
        df["TG"].loc[df["TG"] >= "2021-12-01 13:00:00:00"].plot(label='model')
        plt.title("Surface Temperature")
        plt.legend(loc="upper right")
        plt.savefig(f'{save_dir}/timeseries_tsfc_{station_name}.png')

    if plot_name=='precip':
        plt.figure(figsize=[12, 7])
        precip.plot(label='obs')
        #cohm["precip(obs,mmwe)"].plot(label='obs')
        df["ACCPRCP"].loc[df.index >= '2021-12-01 13:00:00'].plot(label='model')
        plt.title('precip')
        plt.ylabel('precip (mmwe)')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_precip_{station_name}.png')
        #plt.show()

    if plot_name=='albedo':
        plt.figure(figsize=[12, 7])
        Alb.loc[Alb.index >= '2021-12-02 00:00:00'].plot(label='obs')
        #cohm_aws["Albedo(obs,-)"].plot(label='obs')
        df["ALBEDO"].loc[df.index >= '2021-12-01 13:00:00'].plot(label='model')
        plt.title('Albedo')
        plt.ylabel('Albedo')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_albedo_{station_name}.png')
        #plt.show()

    if plot_name=='snowheight':
        plt.figure(figsize=[12, 7])
        daily_sfcheight.plot(label='obs')
        #cohm["hsnow(obs,m)_scaled"].plot(label='obs')
        hsnow_m["SNOWH_scaled"].plot(label='model')
        plt.title('Surface height')
        plt.ylabel('surface height (m)')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_snowheight_{station_name}.png')
        #plt.show()

    if plot_name=='icetemp':
        plt.figure(figsize=[12, 7])
        df["PSNOWTEMP0"].plot(label='icetemp0')
        df["PSNOWTEMP1"].plot(label='icetemp1')
        plt.title('Glacier Temperature')
        plt.ylabel('T(K)')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_icetemp_{station_name}.png')
        #plt.show()

    ### check units
    # plt.figure()
    # df["ACSNOM_scaled"].plot(label="runoff")
    # df["subl"].plot(label="sublimation")
    # df["ACCPRCP_scaled"].plot(label="precipitation")
    # df["PSNOWTOTSWE_scaled"].plot(label="total snow swe")
    # plt.ylabel('mmwe')
    # plt.title('Mass Balance Components')
    # plt.legend(loc='upper right')
    # plt.savefig(f'{save_dir}/timeseries_massbal_{station_name}.png')
    # #plt.show()

    if plot_name=='massbal':
        plt.figure(figsize=[12, 7])
        df["CUMSUM_ACSNOM"].plot(label="runoff")
        df["CUMSUM_subl"].plot(label="sublimation")
        df["CUMSUM_subldrift"].plot(label="sublimation_drift")
        df["CUMSUM_ACCPRCP"].plot(label="precipitation")
        df["PSNOWTOTSWE_scaled"].plot(label="total snow swe")
        plt.ylabel('mmwe')
        plt.title('Mass Balance Components')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_massbal_{station_name}.png')
        #plt.show()


    df_heat = df.filter(regex=r"PSNOWHEAT")
    df_liq = df.filter(regex=r"PSNOWLIQ")
    df_dz = df.filter(regex=r"PSNOWDZ")

    df_rph = (df.filter(regex=r"PSNOWRPH").sum(axis=1))/3600 #J/m2 to W/m2
    df_mph = (df.filter(regex=r"PSNOWMPH").sum(axis=1))/3600 #J/m2 to W/m2
    df_swp = (df.filter(regex=r"PSNOWSWP").sum(axis=1)) #W/m2

    
    #calculate heat content in W/m2
    df_sum = pd.DataFrame()
    #df_sum.index = df_heat.index
    for i in range(40):
        df_sum["PSNOWHEATDZ"+str(i)] = (df_heat["PSNOWHEAT"+str(i)])/3600
        #df_sum["PSNOWHEATDZ"+str(i)] = (df_heat["PSNOWHEAT"+str(i)]*df_dz["PSNOWDZ"+str(i)])/3600
    df_sum["SUM_hc"] = df_sum.filter(regex=r"PSNOWHEATDZ").sum(axis=1)
    df_sum["DIFF_hc"] = df_sum["SUM_hc"].diff()
    df_sum["HEATCONTENT"] = df_sum["DIFF_hc"]


    #calculate phase change in W/m2: 334000 J/kg lh of fusion of melt and 3600 to get J/s to J/hr
    for i in range(40):
        df_sum["PSNOWLIQDZ"+str(i)] = (df_liq["PSNOWLIQ"+str(i)]*df_dz["PSNOWDZ"+str(i)]*334000)/3600
    df_sum["SUM_pc"] = df_sum.filter(regex=r"PSNOWLIQDZ").sum(axis=1)
    df_sum["PHASECHANGE"] = df_sum["SUM_pc"].diff()

    #calculate residual
    df['qm'] = df['FSA'] - df['FIRA'] - df['LH'] - df['HFX'] - df['GRDFLX']

    df['FIRA'] = -1*df['FIRA']
    df['LH'] = -1*df['LH']
    df['HFX'] = -1*df['HFX']
    df['GRDFLX'] = -1*df['GRDFLX']


    if plot_name=='energybal':
        plt.figure(figsize=[12, 8])
        df["qm"].plot(label="residual", linewidth=0.5)
        df["FSA"].plot(label="net shortwave radiation", linewidth=0.5)
        df["FIRA"].plot(label="net longwave radiation", linewidth=0.5)
        df["LH"].plot(label="latent heat flux", linewidth=0.7)
        df["HFX"].plot(label="sensible heat flux", linewidth=0.7)
        df["GRDFLX"].plot(label="ground flux", linewidth=0.7)
        df_rph.plot(label="refreezing", linewidth=0.7)
        df_mph.plot(label="melt", linewidth=0.7)
        df_swp.plot(label="SW penetrative radiation", linewidth=0.7)
        df_sum["HEATCONTENT"].plot(label="heat content", linewidth=0.7)
        #df_sum["PHASECHANGE"].plot(label="phase change")

        #df.plot(linewidth=0.7)
        plt.title(f'Energy balance components for {station_name}',fontsize=18)
        plt.ylabel('Energy (W/m2)', fontsize=14)
        plt.xlabel(f'Datetime (UTC)')
        plt.legend(loc='upper right')
        plt.ylim([-300., 450.])
        plt.savefig(f'{save_dir}/timeseries_energybal_{station_name}.png')

    if plot_name=='melt':
        plt.figure(figsize=[12, 8])
        df_rph.plot(label="refreezing", linewidth=0.7)
        df_mph.plot(label="melt", linewidth=0.7)
        df_swp.plot(label="SW penetrative radiation", linewidth=0.7)
        plt.title(f'Energy balance components for {station_name}',fontsize=18)
        plt.ylabel('Energy (W/m2)', fontsize=14)
        plt.xlabel(f'Datetime (UTC)')
        plt.legend(loc='upper right')
        #plt.ylim([-300., 450.])
        plt.savefig(f'{save_dir}/timeseries_melt_{station_name}.png')
    
    if plot_name=='heatcontent':
        plt.figure(figsize=[12, 7])
        df_sum["HEATCONTENT"].plot(label="heat content")
        plt.title('Heat Content')
        plt.ylabel('Heat content (W/m2)')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_heatcontent_{station_name}.png')

    if plot_name=='phasechange':
        plt.figure(figsize=[12, 7])
        df_sum["PHASECHANGE"].plot(label="phase change")
        plt.title('Phase change energy')
        plt.ylabel('Phase change (W/m2)')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_phasechange_{station_name}.png')

#------------------- scripts to generate xsections, vert profiles etc

    if plot_name=='icetD':

        df_snowh, df_dz, df_var = prep.proc_xsection(save_dir)

        
        #calculate each of the heights for each timestep
        z_005 = df_snowh -0.05
        z_010 = df_snowh -0.1
        z_020 = df_snowh -0.2
        z_050 = df_snowh -0.5
        z_100 = df_snowh -1.0
        z_200 = df_snowh -2.0

        dt = pd.DataFrame(columns=["0.05", "0.1", "0.2", "0.5", "1.0", "2.0"], index=z_005.index)

        for i in range(len(z_005)):
            # print(z_005.iloc[i]) #target to interp to
            f = interpolate.interp1d(df_dz.iloc[i].values, df_var.iloc[i].values, bounds_error=False)
            t_005 = f(z_005.iloc[i])
            t_010 = f(z_010.iloc[i])
            t_020 = f(z_020.iloc[i])
            t_050 = f(z_050.iloc[i])
            t_100 = f(z_100.iloc[i])
            t_200 = f(z_200.iloc[i])
            dt.iloc[i] = pd.Series({'0.05':t_005, '0.1':t_010, '0.2':t_020, '0.5':t_050, '1.0':t_100, '2.0':t_200}, dtype=np.float64)

        dt = dt.astype(np.float64)
        dt = dt.resample('d').mean()
        plt.figure(figsize=[12, 7])
        dt.plot()
        plt.savefig(f'{save_dir}/timeseries_icetD_{station_name}.png')

    if plot_name=='icetH':

        df_snowh, df_dz, df_var = prep.proc_xsection(save_dir)


        #calculate each of the heights for each timestep
        #z_005 = df_snowh -0.05
        #z_010 = df_snowh -0.1
        #z_020 = df_snowh -0.2
        #z_050 = df_snowh -0.5
        #z_100 = df_snowh -1.0
        #z_200 = df_snowh -2.0

        #dt = pd.DataFrame(columns=["0.05", "0.1", "0.2", "0.5", "1.0", "2.0"], index=z_005.index)

        #for i in range(len(z_005)):
            # print(z_005.iloc[i]) #target to interp to
            #f = interpolate.interp1d(df_dz.iloc[i].values, df_var.iloc[i].values, bounds_error=False)
            #t_005 = f(z_005.iloc[i])
            #t_010 = f(z_010.iloc[i])
            #t_020 = f(z_020.iloc[i])
            #t_050 = f(z_050.iloc[i])
            #t_100 = f(z_100.iloc[i])
            #t_200 = f(z_200.iloc[i])
            #dt.iloc[i] = pd.Series({'0.05':t_005, '0.1':t_010, '0.2':t_020, '0.5':t_050, '1.0':t_100, '2.0':t_200}, dtype=np.float64)

        #dt = dt.astype(np.float64)
        #dt = dt-273.15
        
        hsnow = hsnow.resample('H').mean()

        hsnow["h_TC1"] = 0.05 + hsnow["hsnow"] #calculate height of sensor as snowpack melts
        hsnow["h_TC2"] = 0.1 + hsnow["hsnow"]
        hsnow["h_TC3"] = 0.2 + hsnow["hsnow"]
        hsnow["h_TC4"] = 0.5 + hsnow["hsnow"]
        hsnow["h_TC5"] = 1.0 + hsnow["hsnow"]
        hsnow["h_TC6"] = 2.0 + hsnow["hsnow"]

        hsnow["h_TC1"] = hsnow["h_TC1"].where(hsnow["h_TC1"]>=0.0, np.nan) #filter if sensor melts out
        hsnow["h_TC2"] = hsnow["h_TC2"].where(hsnow["h_TC2"]>=0.0, np.nan)
        hsnow["h_TC3"] = hsnow["h_TC3"].where(hsnow["h_TC3"]>=0.0, np.nan)
        hsnow["h_TC4"] = hsnow["h_TC4"].where(hsnow["h_TC4"]>=0.0, np.nan)
        hsnow["h_TC5"] = hsnow["h_TC5"].where(hsnow["h_TC5"]>=0.0, np.nan)
        hsnow["h_TC6"] = hsnow["h_TC6"].where(hsnow["h_TC6"]>=0.0, np.nan)

        df_dz = df_dz.loc[hsnow.index[0]:]
        df_var = df_var.loc[hsnow.index[0]:]
       
        df_snowh = df_snowh.loc[hsnow.index[0]:]
        df_snowh.index = df_snowh.index.tz_localize('UTC')
        z_005 = df_snowh - hsnow["h_TC1"] #snow height of each sensor mapped to model
        z_010 = df_snowh - hsnow["h_TC2"]
        z_020 = df_snowh - hsnow["h_TC3"]
        z_050 = df_snowh - hsnow["h_TC4"]
        z_100 = df_snowh - hsnow["h_TC5"]
        z_200 = df_snowh - hsnow["h_TC6"]

        #interpolate for each sensor height in model
        dt = pd.DataFrame(columns=["0.05", "0.1", "0.2", "0.5", "1.0", "2.0"], index=z_005.index)
        z_005 = z_005.loc[:df_snowh.index[-1]]
        for i in range(len(z_005)):
            # print(z_005.iloc[i]) #target to interp to
            f = interpolate.interp1d(df_dz.iloc[i].values, df_var.iloc[i].values, bounds_error=False)
            t_005 = f(z_005.iloc[i])
            t_010 = f(z_010.iloc[i])
            t_020 = f(z_020.iloc[i])
            t_050 = f(z_050.iloc[i])
            t_100 = f(z_100.iloc[i])
            t_200 = f(z_200.iloc[i])
            dt.iloc[i] = pd.Series({'0.05':t_005, '0.1':t_010, '0.2':t_020, '0.5':t_050, '1.0':t_100, '2.0':t_200}, dtype=np.float64)

        dt = dt.astype(np.float64)
        dt = dt-273.15

        #c["TC5_Avg"].loc["2021-12-15 18:00:00":"2021-12-16 19:30:00"] = c["TC5_Avg"].loc["2021-12-15 18:00:00":"2021-12-16 19:30:00"].where(c["TC5_Avg"]<-8.7, np.nan) #filter out weird meltwater spike
        #c["TC4_Avg"].loc["2021-12-15 18:00:00":"2021-12-16 19:30:00"] = c["TC4_Avg"].loc["2021-12-15 18:00:00":"2021-12-16 19:30:00"].where(c["TC4_Avg"]<-5.2, np.nan)
        melt_index = np.argwhere(c["TC1_Avg"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].values>0.0)[0][0] #index where TC1 melts out
        c["TC1_Avg"].loc["2021-12-01 00:00:00":][melt_index:] = np.nan       
        
        melt_index = np.argwhere(c["TC2_Avg"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].values>0.0)[0][0] #index where TC2 melts out
        c["TC2_Avg"].loc["2021-12-01 00:00:00":][melt_index:] = np.nan

        melt_index = np.argwhere(c["TC3_Avg"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].values>0.0)[0][0] #index where TC3 melts out
        c["TC3_Avg"].loc["2021-12-01 00:00:00":][melt_index:] = np.nan

        plt.figure(figsize=(12,8))
        c["TC1_Avg"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='0.05_ob', color='red', linestyle='dotted')
        c["TC2_Avg"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='0.1_ob', color='orange', linestyle='dotted')
        c["TC3_Avg"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='0.2_ob', color='green', linestyle='dotted')
        c["TC4_Avg"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='0.5_ob', color='blue', linestyle='dotted')
        c["TC5_Avg"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='1.0_ob', color='purple', linestyle='dotted')
        c["TC6_Avg"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='2.0_ob', color='brown', linestyle='dotted')
        dt["0.05"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='0.05_Croc', color='red')
        dt["0.1"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='0.1_Croc', color='orange')
        dt["0.2"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='0.2_Croc', color='green')
        dt["0.5"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='0.5_Croc', color='blue')
        dt["1.0"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='1.0_Croc', color='purple')
        dt["2.0"].loc["2021-12-01 00:00:00":"2021-12-31 23:00:00"].plot(label='2.0_Croc', color='brown')
        plt.axhline(y=0.0, color='k', linestyle='dotted', label='XTT')
        plt.title('Glacier Temperature')
        plt.ylabel('T(C)')
        plt.legend(loc='upper left')
        plt.savefig(f'{save_dir}/timeseries_icetH_{station_name}.png')
        #plt.show()


        #dt = dt.astype(np.float64)
        #plt.figure()
        #dt.plot()
        #plt.savefig(f'{save_dir}/timeseries_icetH_{station_name}.png')

    if plot_name=='flow':
        df = prep.preprocess_flow(save_dir)

        ftsize=14
        fig, axs = plt.subplots(4, 1, sharex=True, figsize=(18,14))
        fig.suptitle('Runoff at CWG', fontsize=ftsize)

        axs[0].plot(df.index, df["PSNOWTHRUFAL"], color='blue', label="PTHRUFAL")
        axs[0].legend(loc='upper right')
        axs[1].plot(df.index, df["PSNOWTHRUFAL"].cumsum(), color='blue', label="PTHRUFAL_acc")
        axs[1].legend(loc='upper right')
        axs[2].plot(df.index, df["FLOW_ICE"], color='green', label="FLOW_ICE")
        axs[2].legend(loc='upper right')
        axs[3].plot(df.index, df["FLOW_SNOW"], color='orange', label="FLOW_SNOW")
        axs[3].legend(loc='upper right')

        my_xticks = []
        for i in df.index.values:
            my_xticks.append(i)
        my_xticks2 = [re.sub(r'\:00\:00\.0+$', '', str(d)) for d in my_xticks]

        plt.xticks(rotation=45)
        plt.subplots_adjust(top=0.925, bottom=0.12, left=0.085, right=0.9)
        fig.supylabel('Runoff (mm)')
        plt.savefig(f'{save_dir}/timeseries_flow_{station_name}.png')

    if plot_name=='4panel':
        height,temp,heat,rho,liq, thruf, fsno, fice,melt,refrz = prep.proc_4panel(date, save_dir, station_name)
        fig, axs = plt.subplots(2, 2)
        fig.suptitle(date)
        fig.text(0.04, 0.5, 'Height (m)', va='center', rotation='vertical')
        axs[0, 0].plot(temp, height)
        axs[0, 0].invert_yaxis() 
        axs[0, 0].axvline(x=273.15, color='k', linestyle='dotted', label='XTT')
        axs[0,0].legend(loc='upper right')
        axs[0, 0].set_xlabel('PSNOWTEMP (K)')
        axs[0, 1].plot(rho, height, 'tab:orange')
        axs[0, 1].set_xlabel('PSNOWRHO (kg/m3)')
        axs[0, 1].axvline(x=850., color='k', linestyle='dotted', label='XRHOTHRESHOLD')
        axs[0,1].legend(loc='upper right')
        axs[0, 1].invert_yaxis()
        axs[1, 0].plot(liq, height, 'tab:green',label='liq')
        axs[1, 0].plot(melt, height, 'tab:red',label='melt')
        axs[1, 0].plot(refrz, height, 'tab:blue',label='refrz')
        axs[1,0].legend(loc='upper right')
        axs[1, 0].set_xlabel('PSNOWLIQ (kg/m3)')
        axs[1, 0].invert_yaxis()
        axs[1, 1].plot(heat, height, 'tab:red')
        axs[1, 1].set_xlabel('PSNOWHEAT (J/m2)')
        axs[1, 1].invert_yaxis()
        plt.savefig(f'{save_dir}/timeseries_4panel_{station_name}.png')

    if plot_name=='xsect_top2':
        var_names = ["PSNOWTEMP", "PSNOWRHO", "PSNOWLIQ", "PSNOWHEAT"]
        Z_list = []
        cp_list = []

        for var_name in var_names[2:]:
            df_snowh, df_dz, df_var = prep.proc_xsection(save_dir, station_name, var_name)
            z = np.arange(df_snowh.values.max() - 0.1, df_snowh.values.max(), 0.001)
            #z = np.arange(df_snowh.values.max() - 0.5, df_snowh.values.max(), 0.01)
            z = np.append(np.arange(df_snowh.values.max() - 1.0, df_snowh.values.max() - 0.1, 0.5), z)
            #z = np.append(np.arange(df_snowh.values.max() - 3.5, df_snowh.values.max() - 0.5, 0.5), z)
            z.sort()
            z_rev = z

            dt = pd.DataFrame(columns=["depths", "var"], index=df_dz.index)
            dt2 = pd.DataFrame(columns=z_rev, index=df_dz.index)

            for index, row in df_dz.iterrows(): #iterating over all of the timesteps
                z_real = df_dz.loc[index].to_list() #extract the depth
                #find index of the first element isnan
                try:
                    nan_index = np.argwhere(np.isnan(z_real))[0][0]
                except IndexError:
                    nan_index = -1 #this only happens when init run
                #interp wants monotonically increasing z_real and no NaNs
                z_real_rev = z_real[0:nan_index]
                z_real_rev.reverse()
                t_real_rev = df_var.loc[index][0:nan_index].to_list()
                t_real_rev.reverse() #reverse is an in place operation

                f = interpolate.interp1d(z_real_rev, t_real_rev, bounds_error=False)
                t = f(z_rev)
                dt.loc[index] = pd.Series({'depths':z_rev, 'var':t})
                for ind in range(len(t)):
                    dt2.loc[index][z_rev[ind]] = t[ind]
        
            data=dt2
            x_vals = np.linspace(0, len(data.index), len(data.index), dtype=int)
            y_vals = z_rev
            X, Y = np.meshgrid(x_vals, y_vals, indexing='ij')
            Z = data.values
            Z_list.append(Z)
        
        ftsize=14
        fig, axs = plt.subplots(len(var_names[2:]), 1, sharex=True, figsize=(16,20))
        fig.suptitle("Cross-sections at CWG AWS", fontsize=ftsize)
    
        for r in range(len(var_names[2:])):
            r+=2 #2,3
            l_ax = r-2 # 0,1 for axes
            var_dict = {
            'PSNOWTEMP': [1, cm.vik, 100, "Temperature", "K", [273.15]],  #[levels, cmap, nbins, label, unit]
            'PSNOWRHO': [1, cm.hawaii, 100, "Density", "kg/m3", [850]],
            'PSNOWLIQ': [1, cm.vik, 100, "Liquid content", "mmwe", [0]],
            'PSNOWHEAT': [1, cm.vik, 100, "Heat content", "J/m2", [0]],
            }
            #print(r)
            #breakpoint()
            x = Z_list[l_ax][~pd.isnull(Z_list[l_ax])]
            min_z = x.min()
            max_z = x.max()
            step_z = (max_z - min_z)/100.
            lev = np.arange(min_z, max_z+(2*step_z), step_z)
            cp = axs[l_ax].contourf(X, Y, Z_list[l_ax], cmap=var_dict[var_names[r]][1], levels=lev) #levels=var_dict[var_names[r]][0])
            fig.colorbar(cp, ax=axs[l_ax], label=f'{var_dict[var_names[r]][3]} ({var_dict[var_names[r]][4]})')
            axs[l_ax].plot(X, df_snowh.values, '-k', linewidth=0.1)

        my_xticks = []
        for i in data.index.values:
            my_xticks.append(i)
        my_xticks2 = [re.sub(r'\:00\:00\.0+$', '', str(d)) for d in my_xticks]
        
        plt.xticks(list(range(0,len(data.index),1)), my_xticks2, rotation=45)
        n=var_dict[var_name][2]
        plt.locator_params(axis='x', nbins=n)
        plt.subplots_adjust(top=0.934, bottom=0.145, left=0.063, right=0.99, hspace=0.19, wspace=0.2)
        fig.supylabel('Height (m)')

        plt.savefig(f'{save_dir}/timeseries_xsect_top2_{station_name}.png', bbox_inches='tight')
    plt.close(plt.figure())

    if plot_name=='xsect_top':
        #var_name="PSNOWTEMP"
        var_names = ["PSNOWTEMP", "PSNOWRHO", "PSNOWLIQ", "PSNOWHEAT"]
        Z_list = []
        cp_list = []

        for var_name in var_names[:2]:
            df_snowh, df_dz, df_var = prep.proc_xsection(save_dir, station_name, var_name)
            z = np.arange(df_snowh.values.max() - 0.1, df_snowh.values.max(), 0.001)
            #z = np.arange(df_snowh.values.max() - 0.5, df_snowh.values.max(), 0.01)
            z = np.append(np.arange(df_snowh.values.max() - 1.0, df_snowh.values.max() - 0.1, 0.5), z)
            #z = np.append(np.arange(df_snowh.values.max() - 3.5, df_snowh.values.max() - 0.5, 0.5), z)
            z.sort()
            z_rev = z

            dt = pd.DataFrame(columns=["depths", "var"], index=df_dz.index)
            dt2 = pd.DataFrame(columns=z_rev, index=df_dz.index)

            for index, row in df_dz.iterrows(): #iterating over all of the timesteps
                z_real = df_dz.loc[index].to_list() #extract the depth
                #find index of the first element isnan
                try:
                    nan_index = np.argwhere(np.isnan(z_real))[0][0]
                except IndexError:
                    nan_index = -1 #this only happens when init run
                #interp wants monotonically increasing z_real and no NaNs
                z_real_rev = z_real[0:nan_index]
                z_real_rev.reverse()
                t_real_rev = df_var.loc[index][0:nan_index].to_list()
                t_real_rev.reverse() #reverse is an in place operation

                f = interpolate.interp1d(z_real_rev, t_real_rev, bounds_error=False)
                t = f(z_rev)
                dt.loc[index] = pd.Series({'depths':z_rev, 'var':t})
                for ind in range(len(t)):
                    dt2.loc[index][z_rev[ind]] = t[ind]
        
            data=dt2
            #var_dict = {
            #'PSNOWTEMP': [np.arange(data.min().min(), 274.0, 0.5), cm.vik, len(data.index)/240, "Temperature", "K", [273.15]],  #[levels, cmap, nbins, label, unit]
            #'PSNOWRHO': [np.arange(data.min().min(), data.max().max(), 0.5), cm.hawaii, len(data.index)/240, "Density", "kg/m3", [850]],
            #'PSNOWLIQ': [np.arange(data.min().min(), data.max().max(), 0.5), cm.vik, len(data.index)/240, "Liquid content", "mmwe", [0]],
            #'PSNOWHEAT': [np.arange(data.min().min(), data.max().max(), 0.5), cm.vik, len(data.index)/240, "Heat content", "J/m2", [0]],
            #}
            x_vals = np.linspace(0, len(data.index), len(data.index), dtype=int)
            y_vals = z_rev
            # y_vals = np.linspace(0, len(z), len(z), dtype=int)
            X, Y = np.meshgrid(x_vals, y_vals, indexing='ij')
            Z = data.values
            Z_list.append(Z)
        
        #breakpoint()
        #var_dict = {
        #    'PSNOWTEMP': [np.arange(Z_list[r].min().min(), 274.0, 0.5), cm.vik, len(data.index)/240, "Temperature", "K", [273.15]],  #[levels, cmap, nbins, label, unit]
        #    'PSNOWRHO': [np.arange(Z_list[r].min().min(), Z_list[r].max().max(), 0.5), cm.hawaii, len(data.index)/240, "Density", "kg/m3", [850]],
        #    'PSNOWLIQ': [np.arange(Z_list[r].min().min(), Z_list[r].max().max(), 0.5), cm.vik, len(data.index)/240, "Liquid content", "mmwe", [0]],
        #    'PSNOWHEAT': [np.arange(Z_list[r].min().min(), Z_list[r].max().max(), 0.5), cm.vik, len(data.index)/240, "Heat content", "J/m2", [0]],
        #}
        ftsize=14
        fig, axs = plt.subplots(len(var_names[:2]), 1, sharex=True, figsize=(16,20))
        fig.suptitle("Cross-sections at CWG AWS", fontsize=ftsize)
    
        for r in range(len(var_names[0:2])):
            var_dict = {
            'PSNOWTEMP': [1, cm.vik, 100, "Temperature", "K", [273.15]],  #[levels, cmap, nbins, label, unit]
            'PSNOWRHO': [1, cm.hawaii, 100, "Density", "kg/m3", [850]],
            'PSNOWLIQ': [1, cm.vik, 100, "Liquid content", "mmwe", [0]],
            'PSNOWHEAT': [1, cm.vik, 100, "Heat content", "J/m2", [0]],
            }
            x = Z_list[r][~pd.isnull(Z_list[r])]
            min_z = x.min()
            max_z = x.max()
            step_z = (max_z - min_z)/100.
            #new_max = max_z + (max_z - min_z)/50
            lev = np.arange(min_z, max_z + (2*step_z), step_z)

            cp = axs[r].contourf(X, Y, Z_list[r], cmap=var_dict[var_names[r]][1], levels=lev)   #, levels=var_dict[var_names[r]][0])
            fig.colorbar(cp, ax=axs[r], label=f'{var_dict[var_names[r]][3]} ({var_dict[var_names[r]][4]})')
            axs[r].plot(X, df_snowh.values, '-k', linewidth=0.1)
        my_xticks = []
        for i in data.index.values:
            my_xticks.append(i)
        my_xticks2 = [re.sub(r'\:00\:00\.0+$', '', str(d)) for d in my_xticks]
        
        #plt.subplots_adjust(bottom=1.3)
        plt.xticks(list(range(0,len(data.index),1)), my_xticks2, rotation=45)
        n=var_dict[var_name][2]
        plt.locator_params(axis='x', nbins=n)
        plt.subplots_adjust(top=0.934, bottom=0.145, left=0.063, right=0.99, hspace=0.19, wspace=0.2)
        fig.supylabel('Height (m)')
        #fig.supxlabel('Datetime (UTC)')
        #plt.contour(cp, levels=var_dict[var_name][-1], colors='white')        

        #plt.title(f'Cross Section of the changes in {var_name} for a pixel')
        #plt.xlabel('Datetime')
        #plt.ylabel('Height (m)')
        plt.savefig(f'{save_dir}/timeseries_xsect_top_{station_name}.png', bbox_inches='tight')

if __name__ == "__main__":
    defopt.run(plot_timeseries)
