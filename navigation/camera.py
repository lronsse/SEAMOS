import holoocean, cv2
import numpy as np
scenario = "ExampleLevel-HoveringSonarCamera"
env = holoocean.make(scenario)
env.act('auv0', [0,0,0,0,10,10,10,10])

for _ in range(1000):
    state = env.tick()

    if "LeftCamera" in state:
        img = state["LeftCamera"]


        # Generate random Gaussian noise
        mean = 0
        stddev = 250
        noise = np.zeros(img.shape, np.uint8)
        cv2.randn(noise, mean, stddev)

        # Add noise to image
        img = cv2.add(img, noise)



        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

        ret, thresh = cv2.threshold(gray, 110, 255, 0)
        contours, hierarchy = cv2.findContours(thresh, 1, 2)
        #print("Number of contours detected:", len(contours))

        for cnt in contours:
            img = cv2.drawContours(img, [cnt], -1, (0, 255, 255), 3)
            x1, y1 = cnt[0][0]
            approx = cv2.approxPolyDP(cnt, 0.1 * cv2.arcLength(cnt, True), True)
            if 4 < len(approx) < 10:
                x, y, w, h = cv2.boundingRect(cnt)
                ratio = float(w) / h
                if w*h > 10 :
                    pass
                    #cv2.putText(img, 'Rectangle', (x1, y1), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 0), 2)

            cv2.imshow("Shapes", img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()