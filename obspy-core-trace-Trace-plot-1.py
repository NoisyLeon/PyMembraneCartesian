from obspy import read
st = read()
tr = st[0]
tr.plot()