import holoocean
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image as im
import time
import scipy

pixel_brightness_threshold = 6.5e-02

# Json approach seaweed system from the side
#"location": [0, 15, -1],
#"rotation": [0.0, 0.0, -90]

#### GET SONAR CONFIG
scenario = "ExampleLevel-HoveringSonar"
config = holoocean.packagemanager.get_scenario(scenario)
config = config['agents'][0]['sensors'][-1]["configuration"]
azi = config['Azimuth']
minR = config['RangeMin']
maxR = config['RangeMax']
binsR = config['RangeBins']
binsA = config['AzimuthBins']



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
    data.save('images/' + name + '.png')

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

    arr1inds = x.argsort()
    x_out = x[arr1inds[::-1]]
    y_out = y[arr1inds[::-1]]

    return x_out, y_out

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

def split_domain(m, x, y):
    line_size = 20
    mse_cutoff = 0.5
    n_cutoff = 10
    should_plot = np.zeros(m)
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
        print(i, average, mse, n)
        if mse < mse_cutoff and n > n_cutoff:
            should_plot[i] = 1

        xseq = np.linspace(x_use[0], x_use[-1], num=line_size)
        yseq =  b + a * xseq
        for k in range(line_size):
            output_line[i, 0, k] = xseq[k]
            output_line[i, 1, k] = yseq[k]
    return output_line, x_corrected0, y_corrected0, should_plot


mse_cutoff = 10
n_cutoff = 28
command = np.array([0,0,0,0,10,10,10,10])
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
            output, x_corrected, y_corrected, should_plot = split_domain(5, x, y)
            fig, ax = plt.subplots(figsize=(9, 6))
            plt.xlim([-100, 100])
            plt.ylim([0, 100])
            print(should_plot)
            for k in range(len(output)):
                if should_plot[k] == 1:
                    plt.plot(output[k,0], output[k,1], color='orange')
                #print(i, output[i,0], output[i,1])
            ax.scatter(x_corrected, y_corrected, s=60, alpha=0.7, edgecolors="k")


            plt.savefig('images/plot'+ str(i))
            time.sleep(0.2)
            plt.close()

            #x_corrected, y_corrected, n, average = av_corrected(x, y)

            # Add scatterplot
            #ax.scatter(x_corrected, y_corrected, s=60, alpha=0.7, edgecolors="k")

            #a, b, r, p, std = scipy.stats.linregress(x_corrected, y_corrected)
            #y_fit = a * x_corrected + b

            #mse = np.square(y_fit - y_corrected).mean()
            #rsq = polyfit(x, y, 1)
            #print(i, average, mse, n)


            # Create sequence of 100 numbers from 0 to 100
            #xseq = np.linspace(min(x), max(x), num=100)



            # Plot regression line


            #if mse < mse_cutoff and n >= n_cutoff:
            #    plt.title('Sth can be seen')
            #else:
            #    plt.title('Not really')


            #time.sleep(0.5)
            #plt.close()


print("Finished Simulation!")
plt.ioff()
plt.show()