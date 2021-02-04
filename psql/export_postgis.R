#!/usr/bin/R

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :GPLv3
# ==============================================================================

## Libraries
library(readr)
library(sf)
library(here)
library(glue)

## Load database
adan_db <- read_delim(here("data", "baobabs", "data_Adansonia.csv"), delim=",")

## As simple features
adan_sf <- st_as_sf(adan_db, coords=c("Long", "Lat"), crs=4326)

## Export sf as gpk
write_sf(adan_sf, here("data", "baobabs", "adansonia_occ.gpkg"))

## Load a single layer GeoPackage into a Postgis database
## Postgis extension must be activated for the database
## SSH tunnel: ssh fdb -L 5432:localhost:5432
psql_db <- "ecology"
psql_host <- "localhost"
psql_port <- "5432"
psql_user <- Sys.getenv("PSQL_USER")
psql_pwd <- Sys.getenv("PSQL_PASSWORD")
psql_tab <- "adansonia_occ"
f <- here("data", "baobabs", "adansonia_occ.gpkg")
cmd <- glue("ogr2ogr -f \"PostgreSQL\" \\
            PG:\"dbname='{psql_db}' host='{psql_host}' port='{psql_port}' \\
            user='{psql_user}' password='{psql_pwd}'\" {f} \\
            -nln \"{psql_tab}\"")
system(cmd)

# End of file