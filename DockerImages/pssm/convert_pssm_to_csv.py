import pandas as pd
import sys
pssm=sys.argv[1]
df = pd.read_csv(pssm,skiprows=2,skip_blank_lines=True,  delim_whitespace=True,header=0,
                  names=["WT",u'A', u'R', u'N', u'D', u'C', u'Q', u'E', u'G', u'H', u'I', u'L', u'K',
       u'M', u'F', u'P', u'S', u'T', u'W', u'Y', u'V', u'A.1', u'R.1', u'N.1',
       u'D.1', u'C.1', u'Q.1', u'E.1', u'G.1', u'H.1', u'I.1', u'L.1', u'K.1',
       u'M.1', u'F.1', u'P.1', u'S.1', u'T.1', u'W.1', u'Y.1', u'V.1',"S1","S2"])

df.dropna(axis=0,inplace=True)
col_to_keep = [u'A', u'R', u'N', u'D', u'C', u'Q', u'E', u'G', u'H', u'I', u'L', u'K',u'M', u'F', u'P', u'S', u'T', u'W', u'Y', u'V']
df = df[col_to_keep]
df = df.T

df.to_csv("/home/pssm.csv")
