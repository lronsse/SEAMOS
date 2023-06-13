import holoocean
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image as im
import time
import scipy
import cv2

pixel_brightness_threshold = 6.5e-02

# Json approach seaweed system from the side
#"location": [0, 15, -1],
#"rotation": [0.0, 0.0, -90]

# Json passing across under the system
#"location": [0, -1, -3],
#"rotation": [0.0, 0.0, 30.0]

#### GET SONAR CONFIG
scenario = "ExampleLevel-HoveringSonarCamera"
config = holoocean.packagemanager.get_scenario(scenario)
config = config['agents'][0]['sensors'][-1]["configuration"]
azi = config['Azimuth']
minR = config['RangeMin']
maxR = config['RangeMax']
binsR = config['RangeBins']
binsA = config['AzimuthBins']



#### RUN SIMULATION
def array_to_image(array, name):
# Converts and array of data from the sonar into an image
    array[array < 0 ] = 0


    array = array * 255
    array = np.rot90(array, 2)
    # creating image object of
    # above array
    data = im.fromarray(array)
    data = data.convert("L")

    # saving the final output
    # as a PNG file
    data.save('images/' + name + '.png')

def seaweed_recognition_scaling(array):
# Image analysis by multiplying the brightness of a pixel with the distance from the vehicle
    new = np.zeros(np.shape(array))
    for i in range(np.shape(array)[0]):
        for j in range(np.shape(array)[1]):
             new[j, i] = (np.shape(array)[1]-j) * array[j, i]
    return new

def seaweed_recognition_first_point(array):
# Creates an array where each columns has only one bright pixel and that is the first pixel above some threshold
    new_image = np.zeros(np.shape(array))
    point_array = np.zeros(np.shape(array)[1])
    for i in range(np.shape(array)[1]):
        a = 0
        j = 0
        while j < np.shape(array)[1] and a == 0:
            if array[j,i] > pixel_brightness_threshold:
                new_image[j,i] = 255
                point_array[i] = j
                a = 1
            j += 1

    return new_image, point_array

def square_to_fan(array):
# Converts the rays of the sonar into a fan like thing
    x = np.zeros(len(array))
    y = np.zeros(len(array))
    i = 0
    while i < len(array):
        alpha = (180 - azi/2 + (i * azi/binsA) )* np.pi/180
        y[i] = -array[i]*np.cos(alpha)
        x[i] = array[i] * np.sin(alpha)
        i += 1
    # Sorts the x values descending and y together with them
    arr1inds = x.argsort()
    x_out = x[arr1inds[::-1]]
    y_out = y[arr1inds[::-1]]

    return x_out, y_out

def polyfit(x, y, degree):
    # Get parameters for a linear fit of the sonar data
    results = {}

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results['determination']

def av_corrected(x, y):
    # Remove all values above the average
    coeff_y = 4
    y = y/ coeff_y
    average = np.average(y)
    y_corrected = []
    x_corrected = []
    for j in range(len(y)):
        if y[j] < average:
            y_corrected.append(y[j])
            x_corrected.append(x[j])
    x_corrected = np.array(x_corrected)
    y_corrected = np.array(y_corrected)
    n = len(x_corrected)

    return x_corrected, y_corrected, n, average

def split_domain(m, x, y):
    # Split the domain into smaller one a create values to create separate plots for every of them.
    # Also average of every small domain to get the distance
    # And whether the mean squared error (mse) and the number of points passing the 'average test' is good enough to plot it
    # -> reason to believe the system is there

    line_size = 20
    mse_cutoff = 0.5
    n_cutoff = 10
    should_plot = np.zeros(m)
    midpoints = np.zeros(m)
    output_line = np.zeros((m, 2, line_size))
    x_corrected0, y_corrected0, n, average = av_corrected(x, y)
    for i in range(m):
        left_bound = int(len(x)*i/m)
        right_bound = int(len(x)*(i+1)/m-1)
        x_use = x[left_bound: right_bound]
        y_use = y[left_bound: right_bound]

        #print(i, x[0], x[-1], left_bound, right_bound, x_use[0], x_use[-1])
        x_corrected, y_corrected, n, average = av_corrected(x_use, y_use)
        a, b, r, p, std = scipy.stats.linregress(x_corrected, y_corrected)
        y_fit = a * x_corrected + b
        mse = np.square(y_fit - y_corrected).mean()
        midpoints[i] = y_corrected.mean()
        #print(i, average, mse, n)
        if mse < mse_cutoff and n > n_cutoff:
            should_plot[i] = 1

        xseq = np.linspace(x_use[0], x_use[-1], num=line_size)
        yseq =  b + a * xseq

        for k in range(line_size):
            output_line[i, 0, k] = xseq[k]
            output_line[i, 1, k] = yseq[k]
    return output_line, x_corrected0, y_corrected0, should_plot, midpoints

def rot2eul(rotation_matrix):
    # Matrix rotation
    beta = -np.arcsin(rotation_matrix[2, 0])
    alpha = np.arctan2(rotation_matrix[2, 1] / np.cos(beta), rotation_matrix[2, 2] / np.cos(beta))
    gamma = np.arctan2(rotation_matrix[1, 0] / np.cos(beta), rotation_matrix[0, 0] / np.cos(beta))
    return np.array((alpha, beta, gamma)) / np.pi * 180

def position():
    # Getting position of the drone
    truth_state = state['PoseSensor']
    truth_xy_states = [[], []]
    truth_titles = ["pos z", "vel", "roll", "pitch", "yaw"]
    truth_states = [[] for _ in range(len(truth_titles))]
    truth_xy_states[0].append(truth_state[0, 3])
    truth_xy_states[1].append(truth_state[1, 3])
    truth_states[0].append(truth_state[2, 3])
    x = truth_xy_states[0]
    y = truth_xy_states[1]
    z = truth_states[0]
    velocities = state['VelocitySensor']
    truth_states[1].append((velocities[0] ** 2 + velocities[1] ** 2 + velocities[2] ** 2) ** 0.5)

    truth_euler = rot2eul(truth_state[0:3, 0:3])
    truth_states[2].append(truth_euler[0])
    truth_states[3].append(truth_euler[1])
    truth_states[4].append(truth_euler[2])

    alpha = truth_states[2]
    beta = truth_states[3]
    gamma = truth_states[4]

    return x, y, z, alpha, beta, gamma

def run_sonar():
    # Get a fan shaped array from sonar
    sonar_image = state['ImagingSonar']
    adjusted_image, point_array = seaweed_recognition_first_point(sonar_image)
    array_to_image(sonar_image, "original" + str(i))
    array_to_image(adjusted_image, "adjusted" + str(i))

    x, y = square_to_fan(point_array)

    # Split domain and plot
    output, x_corrected, y_corrected, should_plot, midpoints = split_domain(5, x, y)
    fig, ax = plt.subplots(figsize=(9, 6))
    plt.xlim([-100, 100])
    plt.ylim([-100, 100])
    print(should_plot)
    for k in range(len(output)):
        if should_plot[k] == 1:
            plt.plot(output[k, 0], output[k, 1], color='orange')
            print(should_plot * midpoints)

    ax.scatter(x_corrected, y_corrected, s=60, alpha=0.7, edgecolors="k")
    # x, y, z, roll, pitch, yaw = position()

    plt.savefig('images/plot' + str(i))
    time.sleep(0.2)
    plt.close()
    return should_plot * midpoints

def camera_noise(img):
    mean = 0
    stddev = 180
    noise = np.zeros(img.shape, np.uint8)
    cv2.randn(noise, mean, stddev)

    # Add noise to image
    img = cv2.add(img, noise)
    return img
def run_camera():
    img = state["LeftCamera"]

    # Generate random Gaussian noise
    img = camera_noise(img)

    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    ret, thresh = cv2.threshold(gray, 110, 255, 0)
    contours, hierarchy = cv2.findContours(thresh, 1, 2)
    # print("Number of contours detected:", len(contours))

    for cnt in contours:
        print(cnt)
        img = cv2.drawContours(img, [cnt], -1, (0, 255, 255), 3)
        x1, y1 = cnt[0][0]
        approx = cv2.approxPolyDP(cnt, 0.1 * cv2.arcLength(cnt, True), True)
        if 4 < len(approx) < 10:
            x, y, w, h = cv2.boundingRect(cnt)
            ratio = float(w) / h
            if w * h > 10:
                pass
                # cv2.putText(img, 'Rectangle', (x1, y1), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 0), 2)
    #cv2.imshow("Shapes", img)
    #cv2.waitKey(0)
    #cv2.destroyAllWindows()
    mask = np.ones(img.shape[:2], dtype="uint8") * 255

    # Draw the contours on the mask
    cv2.drawContours(mask, contours, -1, 0, -1)

    # remove the contours from the image and show the resulting images
    img = cv2.bitwise_and(img, img, mask=mask)
    #cv2.imshow("Mask", mask)
    #cv2.imshow("After", img)
    #cv2.waitKey(0)
    name = 'images/image' + str(i) + '.png'
    cv2.imwrite(name, img)
    #cv2.waitKey(0)
    cv2.destroyAllWindows()

# Cutoffs for the mean squared error and number of points that have to pass the 'average test' to decide that there is a system there
mse_cutoff = 10
n_cutoff = 28

# Direction of motion
command = np.array([0,0,0,0,10,10,10,10])

# Start simulation
with holoocean.make(scenario) as env:
    for i in range(3000):
        env.act("auv0", command)
        state = env.tick()

        if 'ImagingSonar' in state:
            where_system = run_sonar()
            if "LeftCamera" in state:
                run_camera()
            #plt.ioff()
            #plt.show()


print("Finished Simulation!")