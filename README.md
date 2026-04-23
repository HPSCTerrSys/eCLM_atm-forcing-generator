# eCLM atmospheric forcing generator

[![docs](https://github.com/HPSCTerrSys/eCLM_atm-forcing-generator/actions/workflows/doc.yml/badge.svg)](https://github.com/HPSCTerrSys/eCLM_atm-forcing-generator/actions/workflows/doc.yml)

## Introduction

This repository shows how to generate atmospheric forcings for eCLM
simulations.

## Usage / Documentation

Please check the documentation at https://hpscterrsys.github.io/eCLM_atm-forcing-generator/INDEX.html

## Docker usage

To run an opinionated script to download and process either ERA5 or SEAS5 run the following:

Ensure CDSAPI_KEY & CDSAPI_URL are both set to appropriate values and domain.lnd.DE-RuS.240717.nc is your domainfile.

Build the image

```bash
docker build -t atm-forcing .
```

```bash
mkdir -p  data/2026-01
docker run -e CDSAPI_KEY=$CDSAPI_KEY -e CDSAPI_URL=$CDSAPI_URL  -it  -v $(pwd)domain.lnd.DE-RuS.240717.nc:/home/nonroot/domain.nc -v $(pwd)/data/2026-01:/home/nonroot/2026-01   atm-forcing ERA5 2026 01
```

or for SEAS5

```bash
mkdir -p  data/2026-01
docker run -e CDSAPI_KEY=$CDSAPI_KEY -e CDSAPI_URL=$CDSAPI_URL  -it  -v $(pwd)domain.lnd.DE-RuS.240717.nc:/home/nonroot/domain.nc -v $(pwd)/data/2026-01:/home/nonroot/2026-01   atm-forcing SEAS5 2026 01
```

The resulting data can be found at data/2026-01.

## License
eCLM atmospheric forcing generator is open source software and is licensed under the [MIT-License](https://github.com/HPSCTerrSys/eCLM_atm-forcing-generator/blob/master/LICENSE).
