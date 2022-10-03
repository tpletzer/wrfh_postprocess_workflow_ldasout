import glob
import re

configfile: "config.yaml"

# directory containing the netcdf and csv files
FILE_DIR = config['file_dir']
SAVE_DIR = config['save_dir']
#STATION_NAMES = config['station_names'] 
PLOT_NAMES = config['plot_names']

NC_FILES = glob.glob(f"{FILE_DIR}/*LDASOUT*")


rule all:
    input:
        expand(f'{SAVE_DIR}/timeseries_{pname}_middle.png' for pname in PLOT_NAMES),

#rule clean:
#    shell:
#        "rm {SAVE_DIR}/*.csv {SAVE_DIR}/*.png *.html"

rule produceStationData:
    input:
        NC_FILES,
    output:
        "{SAVE_DIR}/timeseries_ldasout_cwg.csv"
    shell:
        "python generate_timeseries_ldasout.py -f {FILE_DIR} --save-dir={SAVE_DIR} --station-name=cwg"

rule createTimeseriesPlot:
    input:
        "{SAVE_DIR}/timeseries_ldasout_cwg.csv"
    output:
        report("{SAVE_DIR}/timeseries_{pname}_cwg.png", category="timeseries plot")
    shell:
        "python plot_timeseries_energybal.py --save-dir={SAVE_DIR} --station-name=cwg --plot-name={wildcards.pname}"

