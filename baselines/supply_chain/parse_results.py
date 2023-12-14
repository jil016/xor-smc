import numpy as np


if __name__ == "__main__":
    res = np.array([[85.0, 90.0, 85.0, 90.0, 85.0],
                    [90.0, 85.0, 85.0, 90.0, 85.0],
                    [85.0, 90.0, 85.0, 90.0, 90.0],
                    [85.0, 85.0, 85.0, 90.0, 85.0],
                    [90.0, 85.0, 85.0, 90.0, 90.0],
                    [85.0, 90.0, 85.0, 85.0, 85.0],
                    [90.0, 90.0, 85.0, 90.0, 90.0],
                    [85.0, 85.0, 30.0, 85.0, 85.0],
                    [85.0, 90.0, 90.0, 85.0, 85.0],
                    [85.0, 90.0, 90.0, 85.0, 85.0]])

    print(np.sum(res, axis=0) / 10.0)
