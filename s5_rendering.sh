#!/bin/bash



sub=sub-S006
base_dir=/bcbl/home/home_g-m/lmengxing/TESTDATA/
fs_dir=$base_dir/ThaTractExampleFigure/$sub/ses-T01/output/flywheel/v0/output/RTP/fs/
SR_tck=$base_dir/rtp-pipeline_4.3.5d/analysis-52/$sub/ses-T01/output/LKN76_clean.tck
MR_tck=$base_dir/ThaTractExampleFigure/$sub/ses-T01/output/LKN44_clean.tck
AR_tck=$base_dir/ThaTractExampleFigure/$sub/ses-T01/output/LKN38_clean.tck
OR_tck=$base_dir/ThaTractExampleFigure/$sub/ses-T01/output/LKN28_clean.tck
DT_tck=$base_dir/rtp-pipeline_4.3.5d/analysis-32/$sub/ses-T01/output/LKN67_clean.tck





3dcalc -a $fs_dir/ROIs/Right-VPL.nii.gz -expr "step(a)" \
            -prefix  $fs_dir/ROIs/Right-VPL_bin.nii.gz -overwrite
3dcalc -a $fs_dir/ROIs/Left-Dentate.nii.gz -expr "step(a)" \
            -prefix  $fs_dir/ROIs/Left-Dentate_bin.nii.gz -overwrite


vglrun mrview -load $fs_dir/brainmask.nii.gz \
				-overlay.load $fs_dir/ROIs/R_S1_dil-1.nii.gz \
				-overlay.colour 0,255,0 \
				-overlay.intensity 0.2,1 \
				-overlay.threshold_min 0.2 \
				-overlay.threshold_max 1 \
				-overlay.interpolation 1 \
                -overlay.load $fs_dir/ROIs/Right-VPL_bin.nii.gz \
                -overlay.colour 255,0,0 \
                -overlay.intensity 0.2,1 \
                -overlay.threshold_min 0.2 \
                -overlay.threshold_max 1 \
                -overlay.interpolation 1 \
                -tractography.load $SR_tck \
                -tractography.lighting 1 \
                -overlay.load $fs_dir/ROIs/R_Primary_Motor_Cortex_dil-1.nii.gz \
                -overlay.colour 0,255,0 \
                -overlay.intensity 0.2,1 \
                -overlay.threshold_min 0.2 \
                -overlay.threshold_max 1 \
                -overlay.interpolation 1 \
               -overlay.load $fs_dir/ROIs/Right-VLa_AND_Right-VLp_dil-1.nii.gz \
                -overlay.colour 255,0,0 \
                -overlay.intensity 0.2,1 \
                -overlay.threshold_min 0.2 \
                -overlay.threshold_max 1 \
                -overlay.interpolation 1 \
                -tractography.load $MR_tck \
                -tractography.lighting 1 \
                -overlay.load $fs_dir/ROIs/Left-Dentate_bin.nii.gz \
                -overlay.colour 128,0,128 \
                -overlay.intensity 0.2,1 \
                -overlay.threshold_min 0.2 \
                -overlay.threshold_max 1 \
                -overlay.interpolation 1 \
                -tractography.load $DT_tck \
                -tractography.lighting 1 \
                -overlay.load $fs_dir/ROIs/Right-V1_AND_Right-V2_dil-1.nii.gz \
                -overlay.colour 0,255,0 \
                -overlay.intensity 0.2,1 \
                -overlay.threshold_min 0.2 \
                -overlay.threshold_max 1 \
                -overlay.interpolation 1 \
                -overlay.load $fs_dir/ROIs/Right-LGN_dil-1.nii.gz \
                -overlay.colour 255,0,0 \
                -overlay.intensity 0.2,1 \
                -overlay.threshold_min 0.2 \
                -overlay.threshold_max 1 \
                -overlay.interpolation 1 \
                -tractography.load $OR_tck \
                -tractography.lighting 1 \
                -overlay.load $fs_dir/ROIs/R_Primary_Auditory_Cortex_dil-1.nii.gz \
                -overlay.colour 0,255,0 \
                -overlay.intensity 0.2,1 \
                -overlay.threshold_min 0.2 \
                -overlay.threshold_max 1 \
                -overlay.interpolation 1 \
                -overlay.load $fs_dir/ROIs/Right-MGN_dil-1.nii.gz \
                -overlay.colour 255,0,0 \
                -overlay.intensity 0.2,1 \
                -overlay.threshold_min 0.2 \
                -overlay.threshold_max 1 \
                -overlay.interpolation 1 \
                -tractography.load $AR_tck \
	            -tractography.lighting 1 \
 			-mode 1 \
				-voxel 127,77,127 \
				-noannotations -fullscreen &

