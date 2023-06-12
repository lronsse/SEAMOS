import holoocean
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image as im
import time
import scipy

pixel_brightness_threshold = 6.5e-02

#### GET SONAR CONFIG
scenario = "PierHarbor-HoveringImagingSonar"
config = holoocean.packagemanager.get_scenario(scenario)
config = config['agents'][0]['sensors'][-1]["configuration"]
azi = config['Azimuth']
minR = config['RangeMin']
maxR = config['RangeMax']
binsR = config['RangeBins']
binsA = config['AzimuthBins']

#### GET PLOT READY
#plt.ion()
#fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(8,5))
#ax.set_theta_zero_location("N")
#ax.set_thetamin(-azi/2)
#ax.set_thetamax(azi/2)

#theta = np.linspace(-azi/2, azi/2, binsA)*np.pi/180
#r = np.linspace(minR, maxR, binsR)
#T, R = np.meshgrid(theta, r)
#z = np.zeros_like(T)

#plt.grid(False)
#plot = ax.pcolormesh(T, R, z, cmap='gray', shading='auto', vmin=0, vmax=1)
#plt.tight_layout()
#fig.canvas.draw()
#fig.canvas.flush_events()

#### RUN SIMULATION
def array_to_image(array, name):

    array[array < 0 ] = 0


    array = array * 255
    array = np.rot90(array, 2)
    # creating image object of
    # above array
    data = im.fromarray(array)
    data = data.convert("L")

    # saving the final output
    # as a PNG file
    data.save(name + '.png')

def seaweed_recognition_scaling(array):
    new = np.zeros(np.shape(array))
    for i in range(np.shape(array)[0]):
        for j in range(np.shape(array)[1]):
             new[j, i] = (np.shape(array)[1]-j) * array[j, i]
    return new

def seaweed_recognition_first_point(array):
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
    x = np.zeros(len(array))
    y = np.zeros(len(array))
    i = 0
    while i < len(array):
        alpha = (180 - azi/2 + (i * azi/binsA) )* np.pi/180
        y[i] = -array[i]*np.cos(alpha)
        x[i] = array[i] * np.sin(alpha)
        i += 1
    return x, y

def polyfit(x, y, degree):
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

mse_cutoff = 10
n_cutoff = 28
command = np.array([0,0,0,0,20,20,20,20])
with holoocean.make(scenario) as env:
    for i in range(3000):
        env.act("auv0", command)
        state = env.tick()

        if 'ImagingSonar' in state:
            sonar_image = state['ImagingSonar']
            adjusted_image, point_array = seaweed_recognition_first_point(sonar_image)
            #print(adjusted_image)
            array_to_image(sonar_image, "original" + str(i))
            array_to_image(adjusted_image, "adjusted" + str(i))
            x, y = square_to_fan(point_array)

            # Initialize layout
            fig, ax = plt.subplots(figsize=(9, 6))

            # Fit linear regression via least squares with numpy.polyfit
            # It returns an slope (b) and intercept (a)
            # deg=1 means linear fit (i.e. polynomial of degree 1)

            x_corrected, y_corrected, n, average = av_corrected(x, y)

            # Add scatterplot
            ax.scatter(x_corrected, y_corrected, s=60, alpha=0.7, edgecolors="k")

            a, b, r, p, std = scipy.stats.linregress(x_corrected, y_corrected)
            y_fit = a * x_corrected + b

            mse = np.square(y_fit - y_corrected).mean()
            rsq = polyfit(x, y, 1)
            print(i, average, mse, n)
            # Create sequence of 100 numbers from 0 to 100
            xseq = np.linspace(min(x), max(x), num=100)



            # Plot regression line
            ax.plot(xseq, b + a * xseq, color="k", lw=2.5);
            plt.xlim([-40, 40])
            plt.ylim([0, 60])
            if mse < mse_cutoff and n >= n_cutoff:
                plt.title('Sth can be seen')
            else:
                plt.title('Not really')

            plt.savefig('plot'+ str(i))
            time.sleep(0.5)
            plt.close()


print("Finished Simulation!")
plt.ioff()
plt.show()