import ModelSolutions as MS
from RayTracing import RayTracing
import json
import importlib
import codecs
import sys
sys.stdout.reconfigure(encoding='utf-8')
open("defn.py", "w").close()
import defn

PROCEDURE_ID = 0 # 0: main(), 1: test(TEST_I, TEST_J, TEST_SHOW), 2: initialize_all()
TEST_I, TEST_J, TEST_SHOW = 0, 0, True
FILENAME = "dump_o4-mini_screenshots.txt"
NUMBER_OF_RUNS = 10

SHAPES = [
    "Sphere",
    "Cube",
    "Square",
    "Cylinder",
    "Two Spheres",
    "Disc",
    "Quad with Circle Hole",
    "Plane with Triangle Hole",
    "Ellipsoid",
    "Cut Sphere",
    "Sphere with Sphere Hole",
    "Sphere with Cylinder Hole",
]

SOLUTIONS = [
    MS.shape0,
    MS.shape1,
    MS.shape2,
    MS.shape3,
    MS.shape4,
    MS.shape5,
    MS.shape6,
    MS.shape7,
    MS.shape8,
    MS.shape9,
    MS.shape10,
    MS.shape11,
]
ORIGINS = [
    [[0, 0, -3], [-1, 3, -2], [2, 2, 3], [0, 0, 10]],
    [[0, 0, -3], [-1, 3, -2], [2, 2, 3], [0, 0, 10]],
    [[0, 0, -3], [-1, 3, -2], [2, 2, 3], [0, 0, 10]],
    [[0, 0, -3], [-1, 3, -2], [0, 5, -1], [0, 0, 10]],
    [[0, 0, -3], [-1, 3, -2], [4, 0, -3], [0, 0, 10]],
    [[0, 0, -3], [-1, 2.5, -2], [-3, -3, -2], [0, 0, 10]],
    [[0, 0, -3], [-1, 3, -2], [2, 2, 3], [0, 0, 10]],
    [[0, 0, -3], [-1, 3, -2], [0, -1, 5], [0, 0, 10]],
    [[0, 0, -3], [-1, 3, -2], [10, 10, 11], [0, 0, 10]],
    [[0, 0, -3], [-1, 3, -2], [3, 3, 0], [0, 0, 10]],
    [[0, 0, -3], [-1, 3, -2], [3, 3, 0], [0, 0, 10]],
    [[0, 0, -3], [-1, 3, -2], [3, 1, 0], [0, 0, 10]],
]
DIRECTIONS = [
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [-0.57735026919, -0.57735026919, -0.57735026919], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [-0.57735026919, -0.57735026919, -0.57735026919], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [-0.57735026919, -0.57735026919, -0.57735026919], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [0, -0.98058067569, 0.19611613513], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [-0.8, 0, 0.6], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [0.57735026919, 0.57735026919, 0.57735026919], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [-0.57735026919, -0.57735026919, -0.57735026919], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [0, 0, -1], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [-0.57735026919, -0.57735026919, -0.57735026919], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [-0.70710678118, -0.70710678118, 0], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [-0.70710678118, -0.70710678118, 0], [0, 0, 1]],
    [[0, 0, 1], [0.33333333333, -0.66666666667, 0.66666666667], [-0.94868329805, -0.31622776601, 0], [0, 0, 1]],
]

def test(i, j, show=False):
    tester = RayTracing(SOLUTIONS[i])
    return tester.test(ORIGINS[i][j], DIRECTIONS[i][j], show=show)

def initialize_all():
    for i in range(len(SOLUTIONS)):
        tester = RayTracing(SOLUTIONS[i])
        for j in range(len(ORIGINS[i])):
            tester.test(ORIGINS[i][j], DIRECTIONS[i][j], show=False, save=True, name=f"Results/model_{i}_test{j}")

def extract_python(item):
    a = item.find("```python") + 10
    if a == 9:
        a = 0
        b = len(item)
    else:
        b = item.find("```", a)
    c = item.find("def ") + 4
    if c == 3:
        return item[a:b], "???"
    d = item.find("(", c)
    if d == -1:
        return item[a:b], "???"
    return item[a:b], item[c:d]

def same(i, j, k):
    count = 0
    try:
        with open(f"Results/model_{i}_test{j}.json", "r") as file1:
            with open(f"Results/python_{i}.{k}_test{j}.json", "r") as file2:
                mat1 = json.load(file1)
                mat2 = json.load(file2)
                if len(mat1) != len(mat2):
                    return False
                if len(mat1[-1]) != len(mat2[-1]):
                    return False
                for y in range(len(mat1)):
                    for x in range(len(mat1[0])):
                        if abs(mat1[y][x] - mat2[y][x]) > 0.001:
                            count += 1
                            if count > 100:
                                return False
    except:
        return False
    return True

if PROCEDURE_ID == 0:
    for i in range(len(SOLUTIONS)):
        for j in range(len(ORIGINS[i])):
            with open(f"Results/model_{i}_test{j}.json", "w") as file_model:
                tester = RayTracing(SOLUTIONS[i])
                json.dump(tester.test(ORIGINS[i][j], DIRECTIONS[i][j], show=False, save=True, name=f"Results/model_{i}_test{j}"), file_model)
    tester = RayTracing(SOLUTIONS[0])
    error = [0] * len(SOLUTIONS)
    with open(FILENAME, "r") as file_read: # CHANGE FILE NAME FOR DIFFERENT MODEL
        stuff = json.load(file_read)
        for i in range(len(SOLUTIONS)):
            for k in range(NUMBER_OF_RUNS):
                item = stuff[i][k]
                func, name = extract_python(item)
                print(f"{i}.{k}")
                try:
                    file_defn = codecs.open("defn.py", "w", "utf-8")
                    file_defn.write(func)
                    file_defn.close()
                    importlib.reload(defn)
                    exec(f"tester = RayTracing(defn.{name})")
                    for j in range(len(ORIGINS[i])):
                        with open(f"Results/python_{i}.{k}_test{j}.json", "w") as file_test:
                            json.dump(tester.test(ORIGINS[i][j], DIRECTIONS[i][j], show=False, save=True, name=f"Results/python_{i}.{k}_test{j}"), file_test)
                except Exception as e:
                    print(e)
                    error[i] += 1
                    continue

    correct = [0] * len(SOLUTIONS)
    for i in range(len(SOLUTIONS)):
        incorrect = 0
        for k in range(NUMBER_OF_RUNS):
            for j in range(len(ORIGINS[i])):
                if not same(i, j, k):
                    incorrect += 1
                    break
        correct[i] = NUMBER_OF_RUNS - incorrect

    for i in range(len(SOLUTIONS)):
        print(f"{SHAPES[i]}: Correct Results = {correct[i]}, Incorrect Results = {NUMBER_OF_RUNS - correct[i] - error[i]}, Error = {error[i]}.")

    print(f"error: {error}")
    print(f"correct: {correct}")
elif PROCEDURE_ID == 1:
    test(TEST_I, TEST_J, TEST_SHOW)
elif PROCEDURE_ID == 2:
    initialize_all()