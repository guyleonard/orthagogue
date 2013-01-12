#!/bin/bash


# -rwxr-xr-x 1 mironov metacenter 31K May  3 18:32 /norstore/user/mironov/workspace/omcl/omcl2.0.3/sample/goodProteins.blast
# ==> /norstore/user/mironov/workspace/omcl/omcl2.0.3/sample/goodProteins.blast <==
# ath|At1g01190_atath|At1g01190_at100.0053500153515350.01107
# ath|At1g01190_atath|At1g01280_at30.444502871295524534964e-55 204
# ath|At1g01190_atath|At1g11600_at28.295093091659531215091e-46 176
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/sample/goodProteins.blast";
./orthAgogue -i $PATH  -nss -o 50 -c 1 -t 0 -p 1 -s '|' 2>err.txt 1>out.txt

# -rw-r--r-- 1 mironov metacenter 11M May  7 12:16 /norstore/user/mironov/workspace/omcl/omcl2.0.3/test/goodProteins.blast
# ==> /norstore/user/mironov/workspace/omcl/omcl2.0.3/test/goodProteins.blast <==
# 559292|P38903|2A5D_YEAST559292|P38903|2A5D_YEAST100.0075700175717570.01557
# 559292|P38903|2A5D_YEAST284812|Q10428|2AD1_SCHPO52.544932043223713535177e-142 500
# 559292|P38903|2A5D_YEAST284812|P78759|2AD2_SCHPO46.455492348194709936143e-135 478
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/test/goodProteins.blast";
./orthAgogue -i $PATH  -nss -o 50 -c 1 -t 1 -p 0 -s '_' 2>>err.txt 1>>out.txt

# -rw-r--r-- 1 mironov metacenter 6.8M May  5 15:10 /norstore/user/mironov/workspace/omcl/omcl2.0.3/2-set/goodProteins.blast
# ==>> /norstore/user/mironov/workspace/omcl/omcl2.0.3/2-set/goodProteins.blast <==
# sce|P38903sce|P38903100.00610001487571487570.01128
# sce|P38903spo|Q1042852.544932043223713535173e-135 478
# sce|P38903spo|P7875948.0849522242237091476143e-127 451
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/2-set/goodProteins.blast";
 ./orthAgogue -i $PATH  -nss -o 50 -c 1 -t 0 -p 1 -s '|' 2>>err.txt 1>>out.txt

# -rw-r--r-- 1 mironov metacenter 54M May  5 15:49 /norstore/user/mironov/workspace/omcl/omcl2.0.3/3-set/goodProteins.blast
# ==>> /norstore/user/mironov/workspace/omcl/omcl2.0.3/3-set/goodProteins.blast <==
# cel|P41932cel|P41932100.0024800124812486e-116 413
# cel|P41932cel|Q2065585.54249342124812481e-99 359
# cel|P41932cel|Q95ZT187.76196231119511965e-80 294
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/3-set/goodProteins.blast";
./orthAgogue -i $PATH  -nss -o 50 -c 1 -t 0 -p 1 -s '|' 2>>err.txt 1>>out.txt

# -rw-r--r-- 1 mironov metacenter 169M May  5 16:31 /norstore/user/mironov/workspace/omcl/omcl2.0.3/4-set/goodProteins.blast
# ==>> /norstore/user/mironov/workspace/omcl/omcl2.0.3/4-set/goodProteins.blast <==
# cel|P41932cel|P41932100.0024800124812481e-115 413
# cel|P41932cel|Q2065585.54249342124812482e-99 359
# cel|P41932dme|P2931081.17239432624272459e-91 330
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/4-set/goodProteins.blast";
./orthAgogue -i $PATH  -nss -o 50 -c 1 -t 0 -p 1 -s '|' 2>>err.txt 1>>out.txt

# -rw-r--r-- 1 mironov metacenter 398M May  5 18:01 /norstore/user/mironov/workspace/omcl/omcl2.0.3/5-set/goodProteins.blast
# ==>> /norstore/user/mironov/workspace/omcl/omcl2.0.3/5-set/goodProteins.blast <==
# ath|P48347ath|P48347100.0025400125412548e-118 421
# ath|P48347ath|B9DFR1100.0025300125312532e-117 420
# ath|P48347ath|F4I1C1100.0024000124012401e-114 410
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/5-set/goodProteins.blast";
./orthAgogue -i $PATH  -nss -o 50 -c 1 -t 0 -p 1 -s '|' 2>>err.txt 1>>out.txt

# -rw-r--r-- 1 mironov metacenter 994M May  6 12:02 /norstore/user/mironov/workspace/omcl/omcl2.0.3/6-set/goodProteins.blast
# ==>> /norstore/user/mironov/workspace/omcl/omcl2.0.3/6-set/goodProteins.blast <==
# ath|P48347ath|P48347100.0025400125412541e-117 421
# ath|P48347ath|B9DFR1100.0025300125312533e-117 420
# ath|P48347ath|F4I1C1100.0024000124012402e-114 410
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/6-set/goodProteins.blast";
./orthAgogue -i $PATH  -nss -o 50 -c 1 -t 0 -p 1 -s '|' 2>>err.txt 1>>out.txt

# -rw-r--r-- 1 mironov metacenter 1.9G May  6 18:33 /norstore/user/mironov/workspace/omcl/omcl2.0.3/7-set/goodProteins.blast 
# ==>> /norstore/user/mironov/workspace/omcl/omcl2.0.3/7-set/goodProteins.blast <==
# ath|P48347ath|P48347100.0025400125412542e-117 421
# ath|P48347ath|B9DFR1100.0025300125312534e-117 420
# ath|P48347ath|F4I1C1100.0024000124012402e-114 410
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/7-set/goodProteins.blast";
./orthAgogue -i $PATH  -nss -o 50 -c 1 -t 0 -p 1 -s '|' 2>>err.txt 1>>out.txt

# -rw-r--r-- 1 mironov metacenter 6.1G May  8 11:06 /norstore/user/mironov/workspace/omcl/omcl2.0.3/sprot/goodProteins.blast
# ==>> /norstore/user/mironov/workspace/omcl/omcl2.0.3/sprot/goodProteins.blast <==
# 11179|HN_NDVC11179|HN_NDVC100.005540018571185710.01089
# 11179|HN_NDVC11183|HN_NDVJ94.9555428018571185710.01034
# 11179|HN_NDVC11182|HN_NDVI92.6055441018571185710.01024
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/sprot/goodProteins.blast";
# ./orthAgogue -i $PATH  -nss -o 50 -c 1 -t 0 -p 1 -s '|' 2>>err.txt 1>>out.txt
./orthAgogue -i "/norstore/user/mironov/workspace/omcl/omcl2.0.3/sprot/goodProteins.blast" -nss -o 50 -c 1 -t 0 -p 1 -s '|' 2>>err.txt 1>>out.txt
# -rw-r--r-- 1 mironov metacenter 1.7G May 11 14:55 /norstore/user/mironov/workspace/omcl/omcl2.0.3/gexo/goodProteins.blast
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/gexo/goodProteins.blast";
# ==>> /norstore/user/mironov/workspace/omcl/omcl2.0.3/gexo/goodProteins.blast <==
# human|P31946human|P31946100.0024600124612465e-118 422
# human|P31946human|B5BU2499.5924610124612469e-118 421
# human|P31946rat|P3521398.3724640124612464e-115 412
./orthAgogue -i $PATH  -nss -o 50 -c 1 -t 0 -p 1 -s '|' 2>>err.txt 1>>out.txt

# /norstore/user/mironov/workspace/omcl/omcl2.0.3/sprot/goodProteins.blast