# wrfh_postprocess_workflow

## Requirements

You will need to have "Snakemake" installed. The easiest is to create a conda environment
```
conda create -p /path/to/snakemake_env
conda activate /path/to/snakemake_env
conda install -c conda-forge -c bioconda snakemake
```

## How to check the syntax of your workflow

```
snakemake --dry-run
```

## How to run the workflow

Using 4 workers (for instance)
```
snakemake -j 4
```

The `output` directory will coontain the output files, including the intermediate files:
```
$ ls output/
1_Scatter.png.nc	2_TimeSeries.png.nc	4_Scatter.png.nc	stat.csv
1_TimeSeries.png.nc	3_Scatter.png.nc	4_TimeSeries.png.nc	wrf-h_report.pdf
2_Scatter.png.nc	3_TimeSeries.png.nc	avg.csv
```

## How to clean up the output directory

```
snakemake -j 1 clean
```

## How to create a report

Once you have run your workflow,
```
snakemake --report report.html
```

## How to submit yyour workflow to mahuika

Edit the `mahuika.json` file. Then type
```
snakemake -j 999 --cluster-config mahuika.json --cluster "sbatch --account={cluster.account} --ntasks={cluster.ntasks} --time={cluster.time}"
```

