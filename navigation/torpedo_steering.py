import holoocean
import numpy as np
from pynput import keyboard

pressed_keys = list()
speed = 100
force = 1000

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
    command = np.zeros(5)
    if 'i' in keys:
        command[[0]] -= val
        command[[2]] += val
    if 'k' in keys:
        command[[0]] += val
        command[[2]] -= val
    if 'j' in keys:
        command[[1]] -= val
        command[[3]] += val
    if 'l' in keys:
        command[[1]] += val
        command[[3]] -= val

    if 'y' in keys:
        command[[4]] += speed
    if 'h' in keys:
        command[[4]] -= speed


    return command

with holoocean.make("SimpleUnderwater-Torpedo") as env:
#with holoocean.make("Dam-Hovering") as env:
    while True:
        if 'q' in pressed_keys:
            break
        command = parse_keys(pressed_keys, force, speed)

        #send to holoocean
        env.act("auv0", command)
        state = env.tick()