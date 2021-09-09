current BVP pipeline uses MATLAB scripts in this order
- init_process 		same iterables used across all 3 modalities
- make_rdeco_in		read and format data for easier use in R-DECO
- make_rdeco_out	copy script to data directory and run R-DECO
- make_validfids	derive all fiducial points from clean peak annotations
- review_validfids	optional, can be used to review and/or refine annotations
- make_features		extract basic pulse rate features from fiducial points directly
			build PRV signals from fiducial points, and extract PRV power features