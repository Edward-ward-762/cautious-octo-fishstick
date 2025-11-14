##Pull Pipeline
nextflow pull Edward-ward-762/cautious-octo-fishstick -r main

##Run Pipeline
nextflow run Edward-ward-762/cautious-octo-fishstick \
-r main \
--qc_input "default" \
--mask_input "default"
