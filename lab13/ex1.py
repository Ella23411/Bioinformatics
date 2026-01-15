import numpy as np

def predict_discrete_steps(A, x0, steps=5):
   
    A = np.asarray(A, dtype=float)
    x = np.asarray(x0, dtype=float).reshape(-1)

    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"A must be square. Got shape {A.shape}.")
    n = A.shape[0]
    if x.shape[0] != n:
        raise ValueError(f"x0 length must match A size {n}. Got {x.shape[0]}.")

    states = np.zeros((steps + 1, n), dtype=float)
    states[0] = x

    for k in range(steps):
        x = np.dot(A, x)
        states[k + 1] = x

    return states

if __name__ == "__main__":
    
    A = [
        [0.9, 0.1],
        [0.2, 0.8]
    ]
    x0 = [1, 0]

    print("Matrix A:")
    print(np.asarray(A, dtype=float))
    print("\nInitial vector x0:")
    print(np.asarray(x0, dtype=float))

    states = predict_discrete_steps(A, x0, steps=5)

    print("\nPredicted states:")
    for k, xk in enumerate(states):
        print(f"Step {k}: {xk}")