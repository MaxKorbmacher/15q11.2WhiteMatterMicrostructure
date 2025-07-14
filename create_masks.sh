# FSL cmds to create masks

fslmaths /Users/max/fsl/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz -thr 2 -uthr 2 corpus_callosum_body.nii.gz

fslmaths /Users/max/fsl/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz -thr 37 -uthr 38 cing_hippo.nii.gz

fslmaths /Users/max/fsl/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz -thr 47 -uthr 48 UnciFasc.nii.gz
