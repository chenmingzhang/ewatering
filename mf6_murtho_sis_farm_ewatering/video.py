# this code needs to install opencv

import cv2
import os

# Path to the folder containing your time-lapse images
image_folder = '//puffball.labs.eait.uq.edu.au/uqczhan2/Desktop/modflow6/ewatering_flopy'

# Output video file name
output_video = 'output.mp4'

# Frame rate (fps)
fps = 2.5

# Get a list of image files in the folder
images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
images.sort()  # Ensure images are in the correct order

# Get the dimensions of the first image (assuming all images have the same dimensions)
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

# Define the codec and create a VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
video = cv2.VideoWriter(output_video, fourcc, fps, (width, height))

for image in images:
    img_path = os.path.join(image_folder, image)
    frame = cv2.imread(img_path)
    video.write(frame)

cv2.destroyAllWindows()
video.release()

print(f"Time-lapse video '{output_video}' created successfully.")