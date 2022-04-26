# Layered_surface_nerve_segmentation
Segmentation of nerve axons and myelin structures from CT scan

## Metrics
Myelin density: 68.13%  
Average nerve area: 1587.5 px²  
Average axon area: 505.96 px²  
Average myelin area: 1081.53 px²  
Average nerve radius: 22.23 px  
Average axon radius: 12.25 px  
Average myelin thickness: 9.98 px  

## Development
For each of the selected nerves, the image is unwrapped (polar coordinate mapping represented as a rectangular image) around its center. The unwrapped image presents a dark line spaning the whole width of the image, corresponding to the myelin layer around the axon. Using the layered-surface detection principile to detect two surfaces in the unwrapped images (from bright to dark and from dark to bright), the inner and outer surfaces of the myelin layered can be detected. Solving the max-flow - min-cut problem of the corresponding graph of the image surfaces with delta = 0.3 and wrapping active, the pixel coordinates of each surfaces are obtained. The surfaces are then transformed back to cartesian coordinates and the shape is saved. This process is repeated for every center and over the first 350 slices of the volume. 

## Visualization
![image](https://user-images.githubusercontent.com/44910949/165318853-7fd0d295-84e6-4a79-b93c-3a21e25b8745.png)
