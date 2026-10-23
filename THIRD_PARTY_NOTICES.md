# Third-party inputs and software

The MIT licence in [LICENSE](LICENSE) covers this repository's original code. It does not grant rights to
external datasets, model matrices, or helper code. These inputs are not redistributed here. Obtain the exact
study versions through their providers, record their terms and checksums, and cite the underlying sources.
Provider pages below were checked on 2026-10-23; their current downloads do not establish the versions used
in the study. Required paths and remaining provenance gaps are in [INPUTS.md](docs/INPUTS.md).

| Input | Source and acquisition | Scope of this notice |
|---|---|---|
| OCIM2-48L circulation model | DeVries (2022), [doi:10.1029/2021GL096018](https://doi.org/10.1029/2021GL096018). Obtain the exact model matrix and MATLAB dependencies from the model provider/author. | `OCIM2_48L_CTL.mat` and external helper functions are absent. Their redistribution permission is not inferred from the paper. `ocim/co2model_batch.m` is the study integration code. |
| World Ocean Atlas 2009 | NOAA, [WOA09](https://www.nodc.noaa.gov/OC5/WOA09/pr_woa09.html). | Obtain the model-prepared `woa09si.mat` and `woa09po4.mat` with their preparation provenance. The NOAA download page alone does not identify those MAT files. |
| Seabed lithology | [dbSEABED, University of Colorado](https://instaar.colorado.edu/~jenkinsc/dbseabed/kml/); the study scripts consume prepared raster/cache inputs. | The public KML product is not the raster input used here. Exact dbSEABED/HubOcean assets, processing history and asset-specific terms still need to accompany the input archive. |
| Sediment organic-carbon stock | Atwood, Witt, Mayorga, Hammill and Sala (2020), [Global Patterns in Marine Sediment Carbon Stocks](https://doi.org/10.3389/fmars.2020.00165). | Obtain the actual raster edition used by the study, including the applicable corrections and uncertainty layers. A citation to the article does not establish the identity of a local raster. |
| Longhurst provinces | Flanders Marine Institute (2009), Longhurst Provinces, via [Marine Regions sources](https://www.marineregions.org/sources.php). | Obtain the boundary dataset, retain its supplied attribution and terms, and record the downloaded edition. |
| Land mask | The external land-mask file named by the relevant scripts. | Its exact source and edition have not been established in this reconstruction. No attribution or redistribution permission is invented here; these must be recovered with the study input archive. |
| GEBCO bathymetry | [GEBCO gridded bathymetry](https://www.gebco.net/data-products/gridded-bathymetry-data). | Step 7 names a GEBCO 2026 NetCDF input. Obtain that edition and follow its accompanying terms and attribution. |
| AIS positions and vessel metadata | Global Fishing Watch / the study's data-access arrangement; [data access](https://globalfishingwatch.org/data-download/) and [public API licence documentation](https://api-doc.globalfishingwatch.org/our-apis/documentation/docs/license-rate-limits). | No AIS positions are included. Public API terms are not treated as permission to redistribute the study's private AIS extract. Recover and follow the agreement applicable to that extract. |

R, Python and MATLAB dependencies retain their own licences. The dependency lists in
[ENVIRONMENT.md](docs/ENVIRONMENT.md) describe imports found in the source; they do not relicense dependencies
or supply the original execution environment. No OCIM or other third-party licence is asserted to be MIT.
