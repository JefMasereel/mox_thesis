current RSP pipeline uses MATLAB scripts in this order
- find_breaths	read RSP signals, apply peak annotation, estimate breaths, extract features
--> logdata_rspc (-c suffix used to distinguish results from alternative annotation methods)

additional scripts to save results for use in python (visualization tools etc)
- export_rsp_raw	save raw RSP segments to .mat file
- export_rsp_annotated	save processed RSP segments to .mat file