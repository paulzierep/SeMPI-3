#!/usr/bin/env bash
# Download the data SeMPI-3 is built and benchmarked on.
#
#   MIBiG 4.0 JSON      ~30 MB   curated clusters, compounds and modules
#   MIBiG antiSMASH 8   ~64 MB   antiSMASH output for every MIBiG entry
#
# Both expand to about 400 MB. antiSMASH itself is NOT needed: MIBiG
# publishes its output.
set -euo pipefail

BASE_URL="https://dl.secondarymetabolites.org/mibig"
RAW_DIR="$(dirname "$0")/../data/raw"
mkdir -p "$RAW_DIR"
cd "$RAW_DIR"

echo "==> MIBiG 4.0 JSON"
[ -f mibig_json_4.0.tar.gz ] || curl -fL -o mibig_json_4.0.tar.gz "$BASE_URL/mibig_json_4.0.tar.gz"
mkdir -p json && tar xzf mibig_json_4.0.tar.gz -C json

echo "==> MIBiG 4.0 antiSMASH 8 output"
[ -f mibig_antismash_4.0.tar.bz2 ] || curl -fL -o mibig_antismash_4.0.tar.bz2 \
    "$BASE_URL/mibig_antismash_4.0_gbk_as8b1.tar.bz2"
mkdir -p antismash && tar xjf mibig_antismash_4.0.tar.bz2 -C antismash

echo "==> PARAS A-domain training set"
# 2,873 adenylation domains with curated substrates and 34-residue binding
# pocket signatures, from Terlouw et al. (PARAS/PARASECT). CC BY 4.0.
mkdir -p ../external/paras
[ -f ../external/paras/dataset_paras.tsv ] || curl -fL -o ../external/paras/dataset_paras.tsv \
    "https://raw.githubusercontent.com/BTheDragonMaster/parasect/main/condensed/data/dataset_paras.tsv"

echo "==> done"
du -sh json antismash ../external
