import xarray as xr
import pandas as pd
import glob
from datetime import datetime
import defopt
import numpy as np


fid_dict = {
    'terminus': [227, 198],  #[i, j]
    'middle': [219, 207],   # [217, 211],
    'highelv': [209, 220],
    'soil': [208, 195],
    }


def LDASOUT_energybal_todf(*, file_dir: str='/nesi/project/uoo03104/code/wrf_hydroCrocus_mahuika/Taylor200_glac_update1_2yrloop/NWM', 
    save_dir: str='/nesi/project/uoo03104/snakemake_output/Taylor200_glac_update1_2yrloop/DEC18', station_name: str='soil'): 
    
    files = glob.glob(f'{file_dir}/*LDASOUT*')
    files = sorted(files)

    end_index = len(files)

    DATETIME = []
    albedo = []
    swdown = []
    lwdown = []
    fira = []
    fsa = []
    sag = []
    lh = []
    grdflx = []
    hfx = []
    rainrate = []
    ugdrnoff = []
    accprecip = []
    snowh = []
    sneqv = []
    qsnow = []
    acsnow = []
    acsnom = []
    qrain = []
    flow_ice = []
    flow_snow = []
    glacierthickness = []
    psnowthrufal = []
    psnowheight = []
    psnowtotswe = []
    psndrift = []

    snowliq = {}
    snowswe = {}
    snowheat = {}
    snowtemp = {}
    snowrho = {}
    snowdz = {}
    snowrefrz = {}
    snowrph = {}
    snowmph = {}
    snowmelt = {}
    snowswp = {}     

    test = xr.open_dataset(files[0],decode_times=False)
    lev_size = test['glacier_levels'].values.size

    for lev in range(lev_size):
        snowliq[lev] = []
        snowswe[lev] = []
        snowheat[lev] = []
        snowtemp[lev] = []
        snowrho[lev] = []
        snowdz[lev] = []
        snowrefrz[lev] = []
        snowrph[lev] = []
        snowmph[lev] = []
        snowmelt[lev] = []
        snowswp[lev] = []

    snliq = {}
    snice = {}
    for lev2 in range(3):
        snliq[lev2] = []
        snice[lev2] = []

    pix_i = fid_dict[station_name][0]
    pix_j = fid_dict[station_name][1]

    for file in files[0:end_index]:
        ds = xr.open_dataset(file,decode_times=False)

        albedo.append(ds['ALBEDO'][:,pix_j,pix_i].values)
        swdown.append(ds['SWFORC'][:,pix_j,pix_i].values)
        lwdown.append(ds['LWFORC'][:,pix_j,pix_i].values)
        fira.append(ds['FIRA'][:,pix_j,pix_i].values)
        fsa.append(ds['FSA'][:,pix_j,pix_i].values)
        sag.append(ds['SAG'][:,pix_j,pix_i].values)
        lh.append(ds['LH'][:,pix_j,pix_i].values)
        grdflx.append(ds['GRDFLX'][:,pix_j,pix_i].values)
        hfx.append(ds['HFX'][:,pix_j,pix_i].values)
        rainrate.append(ds['RAINRATE'][:,pix_j,pix_i].values)
        ugdrnoff.append(ds['UGDRNOFF'][:,pix_j,pix_i].values)
        accprecip.append(ds['ACCPRCP'][:,pix_j,pix_i].values)
        snowh.append(ds['SNOWH'][:,pix_j,pix_i].values)
        sneqv.append(ds['SNEQV'][:,pix_j,pix_i].values)
        qsnow.append(ds['QSNOW'][:,pix_j,pix_i].values)
        acsnow.append(ds['ACSNOW'][:,pix_j,pix_i].values)
        acsnom.append(ds['ACSNOM'][:,pix_j,pix_i].values)
        qrain.append(ds['QRAIN'][:,pix_j,pix_i].values)
        flow_ice.append(ds['FLOW_ICE'][:,pix_j,pix_i].values)
        flow_snow.append(ds['FLOW_SNOW'][:,pix_j,pix_i].values)
        glacierthickness.append(ds['glacier_thickness'][:,pix_j,pix_i].values)
        psnowthrufal.append(ds['PSNOWTHRUFAL'][:,pix_j,pix_i].values)
        psnowheight.append(ds['PSNOWHEIGHT'][:,pix_j,pix_i].values)
        psnowtotswe.append(ds['PSNOWTOTSWE'][:,pix_j,pix_i].values)
        psndrift.append(ds['PSNOWSUBL'][:,pix_j,pix_i].values)

        for l in range(lev_size):
            snowliq[l].append(ds['PSNOWLIQ'][:,pix_j,l,pix_i].values)
            snowswe[l].append(ds['PSNOWSWE'][:,pix_j,l,pix_i].values)
            snowheat[l].append(ds['PSNOWHEAT'][:,pix_j,l,pix_i].values)
            snowtemp[l].append(ds['PSNOWTEMP'][:,pix_j,l,pix_i].values)
            snowrho[l].append(ds['PSNOWRHO'][:,pix_j,l,pix_i].values)
            snowdz[l].append(ds['PSNOWDZ'][:,pix_j,l,pix_i].values) 
            snowrefrz[l].append(ds['PSNOWREFRZ'][:,pix_j,l,pix_i].values)
            snowrph[l].append(ds['PSNOWRPH'][:,pix_j,l,pix_i].values)
            snowmph[l].append(ds['PSNOWMPH'][:,pix_j,l,pix_i].values)
            snowmelt[l].append(ds['PSNOWMELT'][:,pix_j,l,pix_i].values)
            snowswp[l].append(ds['PSNOWSWP'][:,pix_j,l,pix_i].values)    

        for l1 in range(3):
            snliq[l1].append(ds['SNLIQ'][:,pix_j,l1,pix_i].values)
            snice[l1].append(ds['SNICE'][:,pix_j,l1,pix_i].values)

        DATETIME.append(file.split('/')[-1].split('.')[0])

    
    lst = []
    comp = [swdown, albedo, lwdown, fira, fsa, sag, lh, grdflx, hfx, rainrate, ugdrnoff, accprecip, 
            snowh, sneqv, qsnow, acsnow, acsnom, qrain, flow_ice, flow_snow, glacierthickness, psnowthrufal, 
            psnowheight, psnowtotswe, psndrift]
    col_names = ["SWFORC", "ALBEDO", "LWFORC", "FIRA", "FSA", "SAG", "LH", "GRDFLX", "HFX", "RAINRATE", "UGDRNOFF", 
                "ACCPRCP", "SNOWH", "SNEQV", "QSNOW", "ACSNOW", "ACSNOM", "QRAIN", "FLOW_ICE", "FLOW_SNOW", "glacier_thickness", 
                "PSNOWTHRUFAL", "PSNOWHEIGHT", "PSNOWTOTSWE", "PSNDRIFT"]

    for l in range(lev_size):
        col_names.append(f'PSNOWLIQ{l}')
        col_names.append(f'PSNOWSWE{l}')
        col_names.append(f'PSNOWHEAT{l}')
        col_names.append(f'PSNOWTEMP{l}')
        col_names.append(f'PSNOWRHO{l}')
        col_names.append(f'PSNOWDZ{l}')
        col_names.append(f'PSNOWREFRZ{l}')
        col_names.append(f'PSNOWRPH{l}')
        col_names.append(f'PSNOWMPH{l}')
        col_names.append(f'PSNOWMELT{l}')
        col_names.append(f'PSNOWSWP{l}')

        comp.append(snowliq[l])
        comp.append(snowswe[l])
        comp.append(snowheat[l])
        comp.append(snowtemp[l])
        comp.append(snowrho[l])
        comp.append(snowdz[l])
        comp.append(snowrefrz[l])
        comp.append(snowrph[l])
        comp.append(snowmph[l])
        comp.append(snowmelt[l])
        comp.append(snowswp[l])

    for l1 in range(3):
        col_names.append(f'SNLIQ{l1}')
        col_names.append(f'SNICE{l1}')

        comp.append(snliq[l1])
        comp.append(snice[l1])

    for ind in range(len(swdown)):
        lst2 = []
        for val in comp:
            lst2.append(val[ind][0])
        lst.append(lst2)

    df = pd.DataFrame(lst, index=pd.to_datetime(DATETIME), columns=col_names)

    #df_dailycyc = df.groupby([df.index.hour]).mean()

    df.to_csv(f'{save_dir}/timeseries_ldasout_{station_name}.csv')

    #df_dailycyc.to_csv(f'{save_dir}/{station_name}_dailycyc.csv')


if __name__=='__main__':
    defopt.run(LDASOUT_energybal_todf)
