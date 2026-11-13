# Soma-free 3D neurite tracing

Work in `/app` with the supplied `data/` and write your submission to `output/`.

Reconstruct distal neuronal processes from a light-sheet fluorescence microscopy block that contains no cell body. Choose seeds and directions from the image, apply suitable denoising or enhancement, and trace the three-dimensional neurite skeleton. Preserve branches and avoid merging nearby but unconnected fibers.

The input, `data/extractedBlock-27.tif`, is a zlib-compressed 16-bit grayscale TIFF with 200 × 200 × 200 voxels. The physical spacing is 0.525 × 0.525 × 1.5 µm in x, y and z. Write `output/reconstruction.swc` in standard space-delimited SWC format, with seven fields per row: node ID, type, x, y, z, radius and parent ID. Use type 2 for axon, 3 for dendrite, 4 for apical dendrite and 0 for unknown type. Use parent ID −1 at each root. Report coordinates and radii to four decimals. Physical coordinates in µm are preferred. If reporting voxel coordinates, include the header `# units: voxel`.
