#!/usr/bin/env python3
import os
import time
import pandas as pd
from sqlalchemy import create_engine

engine = create_engine('mysql+mysqldb://anonymous@ensembldb.ensembl.org/homo_sapiens_core_92_38')

outfilename = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'genes.json')

q = """
select
et.exon_id,
et.transcript_id,
g.stable_id as id,
g.description,
r.name as chr,
g.seq_region_start as start,
g.seq_region_end as end,
e.seq_region_start as exon_start,
e.seq_region_end as exon_end,
t.seq_region_strand as fwdstrand
from exon_transcript et, exon e, gene g, transcript t, seq_region r
where
g.canonical_transcript_id = et.transcript_id and
g.seq_region_id = r.seq_region_id and
r.coord_system_id = 4 and
r.name NOT RLIKE 'CHR' and
et.transcript_id = t.transcript_id and
e.exon_id =et.exon_id
"""

start_time = time.time()

df = pd.read_sql_query(q, engine,index_col='exon_id')
df['exons'] = list(zip(df.exon_start, df.exon_end))
df['fwdstrand'] = df['fwdstrand'].map({1:True,-1:False})
df['tss'] = df.apply(lambda row: row['start'] if row['fwdstrand'] else row['end'], axis=1)
flatdf = pd.DataFrame(df.groupby(['id','description','tss','chr','start','end','fwdstrand'])['exons'].apply(list)).reset_index()
flatdf.set_index('id', inplace=True)
print(flatdf['chr'].value_counts())
flatdf.to_json(outfilename,orient='index')

print("--- %s seconds ---" % (time.time() - start_time))