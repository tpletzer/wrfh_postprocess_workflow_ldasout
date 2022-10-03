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
import preprocess_xsect.py as prep

def dataframe_to_datetime(d):
    d['Datetime'] = pd.to_datetime(d['date'] + ' ' + d['hour'])
    d = d.set_index('Datetime')
    d = d.drop(['date','hour'], axis=1)
    d['date']=d.index
    return d



def plot_timeseries(*, save_dir: str='/nesi/project/uoo03104/snakemake_output/Taylor200_glac_update1_2yrloop/DEC18', station_name: str='cwg', plot_name: str='precip'):
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


    df = pd.read_csv(f'{save_dir}/timeseries_ldasout_{station_name}.csv', index_col=0)
    df.index = pd.to_datetime(df.index)
    df.index = df.index.tz_localize('UTC')
    df["SNOWH_scaled"] = df["SNOWH"] - df["SNOWH"][0]

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

    if plot_name=='precip':
        plt.figure()
        cohm["precip(obs,mmwe)"].plot(label='obs')
        df["ACCPRCP"].plot(label='model')
        plt.title('precip')
        plt.ylabel('precip (mmwe)')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_precip_{station_name}.png')
        #plt.show()

    if plot_name=='albedo':
        plt.figure()
        cohm_aws["Albedo(obs,-)"].plot(label='obs')
        df["ALBEDO"].plot(label='model')
        plt.title('Albedo')
        plt.ylabel('Albedo')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_albedo_{station_name}.png')
        #plt.show()

    if plot_name=='snowheight':
        plt.figure()
        cohm["hsnow(obs,m)_scaled"].plot(label='obs')
        df["SNOWH_scaled"].plot(label='model')
        plt.title('Surface height')
        plt.ylabel('snowh (mm)')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_snowheight_{station_name}.png')
        #plt.show()

    if plot_name=='icetemp':
        plt.figure()
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
        plt.figure()
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
    df_swp = (df.filter(regex=r"PSNOWSWP").sum(axis=1))/3600 #J/m2 to W/m2

    
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
        plt.figure(figsize=[9.0, 5.0])
        df["qm"].plot(label="residual")
        df["FSA"].plot(label="net shortwave radiation")
        df["FIRA"].plot(label="net longwave radiation")
        df["LH"].plot(label="latent heat flux")
        df["HFX"].plot(label="sensible heat flux")
        df["GRDFLX"].plot(label="ground flux")
        df_rph.plot(label="refreezing")
        df_mph.plot(label="melt")
        df_swp.plot(label="SW penetrative radiation")
        df_sum["HEATCONTENT"].plot(label="heat content")
        #df_sum["PHASECHANGE"].plot(label="phase change")

        #df.plot(linewidth=0.7)
        plt.title(f'Energy balance components for {station_name}',fontsize=18)
        plt.ylabel('Energy (W/m2)', fontsize=14)
        plt.xlabel(f'Datetime (UTC)')
        plt.legend(loc='upper right')
        plt.ylim([-300., 450.])
        plt.savefig(f'{save_dir}/timeseries_energybal_{station_name}.png')
    
    if plot_name=='heatcontent':
        plt.figure()
        df_sum["HEATCONTENT"].plot(label="heat content")
        plt.title('Heat Content')
        plt.ylabel('Heat content (W/m2)')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_heatcontent_{station_name}.png')

    if plot_name=='phasechange':
        plt.figure()
        df_sum["PHASECHANGE"].plot(label="phase change")
        plt.title('Phase change energy')
        plt.ylabel('Phase change (W/m2)')
        plt.legend(loc='upper right')
        plt.savefig(f'{save_dir}/timeseries_phasechange_{station_name}.png')

#------------------- scripts to generate xsections, vert profiles etc

    if plot_name=='icetD'

        df_snowh, df_dz, df_var = prep.proc_xsection()

        
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
        dt.resample('d').mean()
        plt.figure()
        dt.plot()
        plt.savefig(f'{save_dir}/timeseries_icetD_{station_name}.png')

    if plot_name=='icetH'

        df_snowh, df_dz, df_var = prep.proc_xsection()

        
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
        plt.figure()
        dt.plot()
        plt.savefig(f'{save_dir}/timeseries_icetH_{station_name}.png')

    plt.close(plt.figure())

if __name__ == "__main__":
    defopt.run(plot_timeseries)
