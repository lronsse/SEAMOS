import holoocean
import numpy as np
from pynput import keyboard

pressed_keys = list()
speed = 25
force = 25

truth_titles = ["pos z", "vel", "roll", "pitch", "yaw"]

truth_xy_states = [[], []]
truth_states = [[] for _ in range(len(truth_titles))]

def on_press(key):
    global pressed_keys
    if hasattr(key, 'char'):
        pressed_keys.append(key.char)
        pressed_keys = list(set(pressed_keys))

def on_release(key):
    global pressed_keys
    if hasattr(key, 'char'):
        pressed_keys.remove(key.char)

listener = keyboard.Listener(
    on_press=on_press,
    on_release=on_release)
listener.start()

def parse_keys(keys, val, speed):
    command = np.zeros(8)
    if 'i' in keys:
        command[0:4] += val
    if 'k' in keys:
        command[0:4] -= val
    if 'j' in keys:
        command[[4, 7]] += val
        command[[5, 6]] -= val
    if 'l' in keys:
        command[[4, 7]] -= val
        command[[5, 6]] += val

    if 'w' in keys:
        command[4:8] += val
    if 's' in keys:
        command[4:8] -= val
    if 'a' in keys:
        command[[4, 6]] += val
        command[[5, 7]] -= val
    if 'd' in keys:
        command[[4, 6]] -= val
        command[[5, 7]] += val


    return command
scenario = "ExampleLevel-HoveringSonar"

with holoocean.make(scenario) as env:
#with holoocean.make("Dam-Hovering") as env:
    while True:
        if 'q' in pressed_keys:
            break
        command = parse_keys(pressed_keys, force, speed)

        #send to holoocean
        env.act("auv0", command)
        state = env.tick()

        truth_state = state['PoseSensor']

        print(truth_state[0, 3])
        # truth_xy_states[1].append(truth_state[1, 3])
        # truth_states[0].append(truth_state[2, 3])
        #
        # print(truth_xy_states[0])