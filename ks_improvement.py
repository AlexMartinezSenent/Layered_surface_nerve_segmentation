#%%
import skimage.io
import numpy as np
import scipy.interpolate
import scipy.ndimage
import matplotlib.pyplot as plt
import slgbuilder
import cv2 as cv
from mpl_toolkits import mplot3d



def unwrapp_image(image, center,angles):
    r = 50
    radii = np.arange(r) + 1 #radial coordinate for unwrapping
    Y = center[0] + np.outer(radii,np.cos(angles))
    X = center[1] + np.outer(radii,np.sin(angles))
    F = scipy.interpolate.interp2d(np.arange(image.shape[0]), np.arange(image.shape[1]), image)
    p_coords = np.c_[Y.ravel(),X.ravel()]
    U = np.array([F(p[0],p[1]) for p in p_coords])
    U = U.reshape((r,a)).astype(float)

    # fig, ax = plt.subplots(1,2)
    # ax[0].imshow(image, cmap='gray')
    # ax[1].imshow(U, cmap='gray')
    # ax[0].plot(center[0], center[1],"r.")
    # ax[0].plot(p_coords[:,0], p_coords[:,1],"b.")
    return U

def get_boundary_layers(unwrapped_image):
    layers = [slgbuilder.GraphObject(0*unwrapped_image), slgbuilder.GraphObject(0*unwrapped_image)] # no on-surface cost
    helper = slgbuilder.MaxflowBuilder()
    helper.add_objects(layers)
    helper.add_layered_boundary_cost()
    helper.add_layered_region_cost(layers[0], 255-unwrapped_image, unwrapped_image)
    helper.add_layered_region_cost(layers[1], unwrapped_image, 255-unwrapped_image)
    helper.add_layered_smoothness(delta=0.3, wrap=True)  
    helper.add_layered_containment(layers[0], layers[1], min_margin=1, max_margin=10)

    helper.solve()
    segmentations = [helper.what_segments(l).astype(int) for l in layers]
    segmentation_lines = [s.shape[0] - np.argmax(s[::-1,:], axis=0) - 1 for s in segmentations]
    return segmentation_lines


#%%%%%%%%%%%%%%%%%%%%
image = skimage.io.imread('nerves/nerves_part.tiff')
a = 180 # number of angles for unfolding
angles = np.arange(a)*2*np.pi/a
#%%%
# centers = {
#     "1": [167, 218],
#     "2": [130, 260],
#     "3": [100, 140],
#     "4": [170, 120],
#     "5": [280, 120],
#     "6": [230, 280]
# }
centers = {         ##### Nina's centers
    "1": [100, 220 ],
    "2": [132, 260 ],
    "3": [165, 215 ],
    "4": [218, 207 ],
    "5": [95, 138 ],
    "6": [240, 50 ]
}
colors = {
    "1": "red",
    "2": "blue",
    "3": "green",
    "4": "yellow",
    "5": "orange",
    "6": "purple"
}
contours = {}
radii = {}
areas = {}
fig, ax = plt.subplots(figsize=(10,10))
ax.imshow(image[0],cmap="gray")
for im in image[0:350]:
    for key,center in zip(centers.keys(),centers.values()):
        U = unwrapp_image(im,center,angles)
        seg_lines = get_boundary_layers(U)
        coord_seg_lines = [None,None]
        for ind,line in enumerate(seg_lines):
            Y = center[0] + np.multiply(line,np.cos(angles))
            X = center[1] + np.multiply(line,np.sin(angles))
            coord_seg_lines[ind] = [Y,X]
        if key not in contours.keys():
            radii[key] = [seg_lines]
            contours[key] = [coord_seg_lines]
        else:
            contours[key].append(coord_seg_lines)
            radii[key].append(seg_lines)
        centers[key] = [np.mean(coord_seg_lines[0][0]), np.mean(coord_seg_lines[0][1])]

# %%
h = np.ones(len(angles))
fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection="3d")
for key, cont in zip(contours.keys(),contours.values()):
    for z,shapes in enumerate(cont):
        for py,px in shapes:
            ax.plot3D(py,px,z*h,color=colors[key],alpha=0.3)
ax.view_init(elev=30, azim=-40)
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.grid(True)
plt.savefig("layered_surface_6_nerves.png", bbox_inches="tight")
# ax.plot3D(Y,X,z*np.ones((len(Y))), color=colors[key], alpha=0.1

#%%
inner_radius = np.empty(350)
outer_radius = np.empty(350)
radii_metrics = {}
for key, rad in zip(radii.keys(),radii.values()):
    for z, radius in enumerate(rad):
        outer_radius[z] = np.mean(radius[0])
        inner_radius[z] = np.mean(radius[1])
    radii_metrics[key] = [np.mean(inner_radius), np.mean(outer_radius)]
avg_inner_radius = np.mean([radii_metrics[key][0] for key in radii_metrics.keys()])
avg_outer_radius = np.mean([radii_metrics[key][1] for key in radii_metrics.keys()])

#%% 
outer_area = np.empty(350)
inner_area = np.empty(350)
area_metrics = {}
for key, cont in zip(contours.keys(),contours.values()):
    for z,shapes in enumerate(cont):
        outer_px = shapes[0][1]
        outer_py = shapes[0][0]
        inner_px = shapes[1][1]
        inner_py = shapes[1][0]
        outer_area[z] = 0.5*np.abs(np.dot(outer_px,np.roll(outer_py,1))-np.dot(outer_py,np.roll(outer_px,1)))
        inner_area[z] = 0.5*np.abs(np.dot(inner_px,np.roll(inner_py,1))-np.dot(inner_py,np.roll(inner_px,1)))
    area_metrics[key] = [np.mean(outer_area), np.mean(inner_area)]
avg_inner_area = np.mean([area_metrics[key][1] for key in area_metrics.keys()])
avg_outer_area = np.mean([area_metrics[key][0] for key in area_metrics.keys()])


#%%
print(f"Myelin density: {round((avg_outer_area-avg_inner_area)*100/avg_outer_area,2)}%")
print(f"Average nerve area: {round(avg_outer_area,2)} px²")
print(f"Average axon area: {round(avg_inner_area,2)} px², Average myelin area: {round(avg_outer_area-avg_inner_area,2)} px²")
print(f"Average nerve radius: {round(avg_outer_radius,2)} px")
print(f"Average axon radius: {round(avg_inner_radius,2)} px, Average myelin thickness: {round(avg_outer_radius-avg_inner_radius,2)} px")


# Myelin density: 68.13%
# Average nerve area: 1587.5 px²
# Average axon area: 505.96 px², Average myelin area: 1081.53 px²
# Average nerve radius: 22.23 px
# Average axon radius: 12.25 px, Average myelin thickness: 9.98 px

# 220 100
# 215 135
# 220 170
# 200 220
# 140 100
# 160 130


# GARBAGE
# %%
# unwrapped_image = unwrapp_image(image[1],(100,140),angles)
# # %%
# layer = slgbuilder.GraphObject(unwrapped_image)
# helper = slgbuilder.MaxflowBuilder()
# helper.add_object(layer)
# helper.add_layered_boundary_cost()
# helper.add_layered_smoothness(delta=1, wrap=False)

# helper.solve()
# segmentation = helper.what_segments(layer)
# segmentation_line = segmentation.shape[0] - np.argmax(segmentation[::-1,:], axis=0) - 1

# plt.imshow(unwrapped_image)
# plt.plot(segmentation_line, 'r')
# #%%%%%%%%
# Y = 230 + np.multiply(segmentation_line,np.cos(angles))
# X = 280 + np.multiply(segmentation_line,np.sin(angles))
# plt.imshow(image[1],cmap='gray')
# plt.plot(Y,X,"r-")


# # %%
# layers = [slgbuilder.GraphObject(0*unwrapped_image), slgbuilder.GraphObject(0*unwrapped_image)] # no on-surface cost
# helper = slgbuilder.MaxflowBuilder()
# helper.add_objects(layers)
# helper.add_layered_boundary_cost()
# helper.add_layered_region_cost(layers[0], 255-unwrapped_image, unwrapped_image)
# helper.add_layered_region_cost(layers[1], unwrapped_image, 255-unwrapped_image)
# helper.add_layered_smoothness(delta=0.2, wrap=True)  
# helper.add_layered_containment(layers[0], layers[1], min_margin=1, max_margin=10)

# helper.solve()
# segmentations = [helper.what_segments(l).astype(float) for l in layers]
# segmentation_lines = [s.shape[0] - np.argmax(s[::-1,:], axis=0) - 1 for s in segmentations]

# plt.imshow(unwrapped_image)
# for line in segmentation_lines:
#     print(line)
#     plt.plot(line, 'r')

# #%%
# for line in segmentation_lines:
#     Y = 230 + np.multiply(line,np.cos(angles))
#     X = 280 + np.multiply(line,np.sin(angles))
#     plt.imshow(image[1],cmap='gray')
#     plt.plot(Y,X,"r-")

