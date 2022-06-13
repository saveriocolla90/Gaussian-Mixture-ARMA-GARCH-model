import numpy as np
import pandas as pd
from tabulate import tabulate

fdati = "./RUNEUSDT_4H.csv"
fout = "./RUNEUSDT_4H_conv.txt"
n_dati = 301

#-Read and Convert the index from TimeStamp to DateTime
CSV_read = pd.read_csv(fdati).set_index('Open Time')
CSV_read.index = pd.to_datetime(CSV_read.index/1000, unit='s')	
	
#-Get the timeseries of price return
time_series = CSV_read.iloc[-n_dati:]
time_series = 100.*time_series.Close.pct_change().dropna()
print(time_series)
quit()
#-Create a table
table = [[]]*len(time_series.index)
for i in range(len(time_series.index)):
	table[i] = [time_series.index[i].strftime("%Y%m%d %H%M%S"), round(time_series[i],4)]

with open(fout,"w") as f:
	f.write(tabulate(table, tablefmt="plain"))
print("Input file: %s" % fdati)
print("Time series converted to file: %s" % fout)

