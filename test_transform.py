import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
from numpy import exp, sin, cos, sqrt, conj

theta = 1 / 2 * np.arcsin(2 / np.sqrt(5))

mat_for_d = np.array([[1/2, 3**0.5/2], [3**0.5/2, -1/2]])
print(scipy.linalg.sqrtm(mat_for_d))
U = np.array([[sin(theta), -1j * cos(theta)], [cos(theta), 1j * sin(theta)]])
V = np.array([[cos(theta), -1j * sin(theta)], [sin(theta), 1j * cos(theta)]])
S = np.array([[sqrt(5) / 2, 0], [0, sqrt(5) / 2]])
A = (U @ S @ conj(V).T)
S = np.array([[1 / 2, 0], [0, 1 / 2]])

B = (U @ S @ V.T)
print(A)
print(B)
print(A@B.T - (A@B.T).T)
print(A@A.T - B@B.T)

# A = np.array([[1, -1 / 2], [1 / 2, 1]])
# B = np.array([[0, 1 / 2], [1 / 2, 0]])

print(np.linalg.svd(A))
print(np.linalg.svd(B))

T = 10
NT = 100
gamma = 1
v = np.load('v_func.npy')
times = np.linspace(0, 30, 400)
dt = times[2] - times[1]
sigma_u = 1
u = exp(-times ** 2 / (2 * sigma_u ** 2))
u = u / (np.linalg.norm(u) * dt ** 0.5)
print(np.linalg.norm(u) ** 2 * dt)
second_term = np.zeros_like(times)

for i, tprime in enumerate(times):
    # second_term[i] = sum([v[j].conjugate() * gamma * np.exp(gamma / 2 * (tprime - t))
    #                     * cosh(epsilon_int(tprime, t)) for j, t in enumerate(times[i:])])*dt
    second_term[i] = sum(
        np.array([v[j] * gamma * exp(gamma / 2 * (tprime - t)) for j, t in enumerate(times)])[times - tprime > 0]) * dt
plt.plot(times, v - second_term, linestyle='dashed')
plt.plot(u)
plt.show()
