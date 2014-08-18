select oh, 1DPM_clean_A and resi 909 and name O
select po, 1DPM_clean_A and resi 900 and name O3
select oe, 1DPM_clean_A and resi 900 and name O2

select lpo, tvx and name O1
select loh, tvx and name O3
select loe, tvx and name O2

pair_fit loh,oh, lpo,po, loe,oe
