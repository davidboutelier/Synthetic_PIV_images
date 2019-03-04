# Synthetic_PIV_images
Scripts to make synthetic images for PIV

The MATLAB scripts in this repo are based on Thielicke and Stamhuis PIVlab. I have modified thescript to make shear zones instead of vortex.

## Vortex
Simple test to check vorticity

## Pure shear
derived from the vortex. Diffrent sign for one component. It produces pure shear at the center (no vorticity). 

## Static shear zone
v-component of velocity changes across a zone parallel to the y-axis. v is set for both sides of the shear zone, and an erf function is used to smooth the velocity from one side to the other across the shear zone.

## Moving shear zone
Similar to static shear zone, but all particles alos have a u homogeneous component of velocity and the shear zone is advected along x with u.

Thielicke, W and Stamhuis, E J 2014 PIVlab – Towards User-friendly, Affordable and Accurate Digital Particle Image Velocimetry in MATLAB. Journal of Open Research Software, 2: e30, DOI: http://dx.doi.org/10.5334/jors.bl
