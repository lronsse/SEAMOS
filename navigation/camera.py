import holoocean, cv2

env = holoocean.make("Dam-HoveringCamera")
env.act('auv0', [10,10,10,10,0,0,0,0])

for _ in range(200):
    state = env.tick()

    if "LeftCamera" in state:
        img = state["LeftCamera"]
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

        ret, thresh = cv2.threshold(gray, 110, 255, 0)
        contours, hierarchy = cv2.findContours(thresh, 1, 2)
        #print("Number of contours detected:", len(contours))

        for cnt in contours:
            img = cv2.drawContours(img, [cnt], -1, (0, 255, 255), 3)
            x1, y1 = cnt[0][0]
            approx = cv2.approxPolyDP(cnt, 0.1 * cv2.arcLength(cnt, True), True)
            if 3 < len(approx) < 10:
                x, y, w, h = cv2.boundingRect(cnt)
                ratio = float(w) / h
                if w*h > 10 :
                    pass
                    #cv2.putText(img, 'Rectangle', (x1, y1), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 0), 2)

                cv2.imshow("Shapes", img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()